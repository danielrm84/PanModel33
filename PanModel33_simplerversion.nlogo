;**************************************************************************************
;          Project Name: Evolutionary adaptive responses to rapid environmental change
;                                  UNIVERSITY OF POTSDAM
;
; Author(s):       Kahl S., Paraskevopoulou, S., Folkertsma R., and Romero-Mujalli D.
;
; Wrote the code:  Daniel Romero Mujalli
;
; Written:         27. 6.2016
; Last update:     29. 5.2019
;
; Type of model:   Individual-based model (IBM)
;
; Summary:         The purpose of the model is:
;                  to study the local adaptation ability (genetic changes and
;              phenotypic plasticity) of different kinds of organisms (life history
;              strategies) under scenarios of environmental change
;
;
;              NOTES / COMMENTS / QUESTIONS:
;
;
;**************************************************************************************
;**************************************************************************************
;                            MIT License
;
;                   Copyright (c) 2019 Daniel Romero Mujalli
;
;   Permission is hereby granted, free of charge, to any person obtaining a copy of
;   this software and associated documentation files (the "Software"), to deal in the
;   Software without restriction, including without limitation the rights to use, copy,
;   modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
;   and to permit persons to whom the Software is furnished to do so, subject to the
;   following conditions:
;
;   The above copyright notice and this permission notice shall be included in all
;   copies or substantial portions of the Software.
;
;   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
;   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
;   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
;   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
;   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
;
;**************************************************************************************
; load r-extension
;extensions [r]

;global parameters of the model:
globals
[
  genetic-mean          ; initial genetic mean of the population
  env-effect-mean       ; mean random environmental effect on phenotype (standard model and
                        ; implicit genetics only)
  initial-env-optimum   ; the initial phenotypic optimum as given by the environment
  size-of-environment   ; the preferred size of the environment (distance from the center)
  female-male-ratio     ; controls the sex ratio in the population. Example: a value of 0.7
                        ; means 70% females and 30% males in the population
  env-effect-variance   ; variance of environmental effect (standard model and implicit
                        ; genetics only)
  extinction?           ; boolean: true if the population is extinct
  pop-gv                ; population level additive genetic variance
  density-compensation  ; governs the population dynamics defining the life strategy
  strength-selection    ; strength of selection used in fitness function
  mu-dist-mut-effects   ; mean distribution of effect size from % of beneficial-mutations
  ;time-limit            ; desired time limit (in generations)
]

; properties of the individuals
turtles-own
[
  phenotype            ; phenotype 'z'
  genetic-component    ; genetic component 'a' of the phenotype
  environmental-effect ; environmental effect 'e'
  fitness              ; fitness wi of individual i
  fecundity            ; fecundity lambda of individual i
  stage                ; whether "adult" or "juvenile", as string
  sex                  ; sex: "female" or "male"
  reproduced?          ; boolean
  dna-strain1          ; first strain of the chromosome
  dna-strain2          ; second strain of the chromosome
]

; properties of the environment (patches)
patches-own
[
  optimum              ; phenotypic optimum tita as given by the environment
  noise                ; variance of the environmental optimum
  degree-maladaptation ; degree of maladaptation
  time-to-extinction   ; store the time when the population went extinct
  mean-env-optimum     ; the mean environmental optimum
]






;**************************************************************************************
;        TO SETUP
;**************************************************************************************
; This procedure set the values for the global parameters of the model, for the optimum
; phenotype as given by the environment (here the environment consists of a group of
; patches with green color). Then each individual is created with own values for the
; genetic and environmental components. These two components are used to calculate
; the phenotye. Finally, the distribution of phenotypes of the population is plotted.
to setup

  ; clear interface, reset ticks
  clear-all
  ;r:clear
  reset-ticks

  ; set values for global parameters
  set-global-parameters

  ; set the environment:
  ; patch at the center (0,0) and those in radius 10 from the center are use as
  ; the environment to test for local adaptation
  ; more than one patch were selected for better visualization of the world
  ;***************************************************************************
  ; THIS SHOULD BE CHANGED TO ACCOUNT FOR ENVIRONMENTAL HETEROGENEITY OR
  ; FRAGMENTATION
  ;***************************************************************************
  ask patch 0 0
  [
    set pcolor green
    set optimum initial-env-optimum
    set noise 0
    set degree-maladaptation compute-degree-maladaptation
    set time-to-extinction 0
    set mean-env-optimum initial-env-optimum

    ask patches with [distance myself <= size-of-environment]
    [
      set pcolor green
      set optimum initial-env-optimum
      set noise 0
      set degree-maladaptation compute-degree-maladaptation
      set time-to-extinction 0
      set mean-env-optimum initial-env-optimum
    ]
  ]
  ;***************************************************************************

  ; create a population of N individuals,
  ; initialize individuals (the turtles) and
  ; calculate phenotype (update phenotype)
  create-turtles population-size
  [
    set color blue
    move-to one-of patches with [pcolor = green]
  ]

  ; set initial conditions according to the method for simulating genetics
  if-else (how-genetics? = "implicit")
  [ ask turtles [ set-initial-conditions-implicit-genetics ] ]
  [
    if-else (how-genetics? = "explicit")
    [ ask turtles [ set-initial-conditions-explicit-genetics ] ]
    [ error "undefined method for modelling genetics" stop ] ;catch exception
  ]

  ; update phenotype
  if-else (how-plasticity? = "standard-model")
  [ ask turtles [ update-phenotype ] ] ;standard model
  [ ask turtles [ update-phenotypic-response] ] ; considers plasticity

  ; update mean of the distribution of effect size
  ;if-else (beneficial-mutations != 0.5)
  ;[set mu-dist-mut-effects mean-dist (beneficial-mutations) (mut-effect-size)]
  ;[set mu-dist-mut-effects 0]
  set mu-dist-mut-effects 0

  ; calculate fitness of each individual
  if-else( fitness-function = "Bjoerklund2009")
  [
    ask turtles [ check-fitness-Bjoerklund2009 ]
    ; DEBUG:
    ;print "fitness according to Bjoerklund"
  ]
  [ ; else, negative-exponential
    if-else (fitness-function = "negative-exponential")
    [ ask turtles [ check-fitness-negative-exponential ] ]
    [ error "fitness function not identified" stop ] ; else, error message
  ]

  ; update fitness-plot
  plot-mean-fitness

  ;update-output:
  update-output

end






;**************************************************************************************
;        TO GO
;**************************************************************************************
to go

  ; iterations are counted at the beginning of the loop.
  ; Each iteration means one generation. The model
  ; uses non-overlapping generations:
  tick

  ; Environmental scenario:
  ; update-environment: updates optimum phenotype as given by the-environment according
  ; to the scenario of environmental change.
  ; If specified, the environmental change starts after a given number of
  ; iterations/generations (mutation selection balance/mutation-selection drift balance)
  ;**********************************************************************
  ; NEEDS TO BE CHANGED IF EXTENDING THE MODEL TO SIMULATE HETEROGENEOUS /
  ; FRAGMENTED LANDSCAPES
  ;**********************************************************************
    ask patch 0 0 [ update-environment ]
  ;**********************************************************************

  ; update-stage, reproductive status and sex of individuals:
  ask turtles [   update-stage-sex  ]

  ; calculate fitness of each individual
  if-else( fitness-function = "Bjoerklund2009")
  [
    ask turtles [ check-fitness-Bjoerklund2009 ]
    ; DEBUG:
    ;print "fitness according to Bjoerklund"
  ]
  [ ; else, negative-exponential
    if-else (fitness-function = "negative-exponential")
    [ ask turtles [ check-fitness-negative-exponential ] ]
    [ error "fitness function not identified" stop ] ; else, error message
  ]

  ; update fitness-plot
  plot-mean-fitness

  ; calcualte fecundity of each individual
  ask turtles [ set-fecundity-Bjoerklund ]

  ; reproduction
  ; store population-level additive genetic variance
    if (how-genetic-variance = "population-level")
    [
      if-else ( count turtles > 1)
      [
        set pop-gv variance [genetic-component] of turtles with [stage = "adult"]
      ]
      [ set pop-gv 0 ] ; else, if only one turtle
    ]
  ; lottery polygyny: females reproduce only once, but males can repeat
  ask turtles with [stage = "adult" and sex = "female" ]
  [
    if-else (how-genetics? = "implicit")  [ reproduce-implicit-genetics ]
    [ if-else (how-genetics? = "explicit")[ reproduce-explicit-genetics ]
      [ error "undefined method of how-genetics?" stop ] ; catch exception
    ]
  ]

  ; mortality of adults
  ask turtles with [ stage = "adult" ] [ die ]

  ; update phenotype
  if-else (how-plasticity? = "standard-model")
  [ ask turtles [ update-phenotype ] ] ; standard model
  [ ask turtles [ update-phenotypic-response] ] ; account for plasticity

  ; compute the degree of maladaptation according to Björklund et al (2009)
  ;**********************************************************************
  ; NEEDS TO BE CHANGED IF EXTENDING THE MODEL TO SIMULATE HETEROGENEOUS /
  ; FRAGMENTED LANDSCAPES
  ;**********************************************************************
  ask patch 0 0
  [
    ; DEBUG:
    ;write "degree-maladaptation before: " print degree-maladaptation
    let new-value compute-degree-maladaptation
    let old degree-maladaptation

    set degree-maladaptation old + new-value
    ; DEBUG:
    ;write "degree-maladaptation after: " print degree-maladaptation

    ask patches with [distance myself <= size-of-environment]
    [ set degree-maladaptation old + new-value ]
  ]

  ; check time-to-extinction:
  if (count turtles < 1 and extinction? = false)
  [
    ask patches with [ pcolor = green ] [ set time-to-extinction (ticks - time-to-balance) ]
    set extinction? true
  ]
  ;**********************************************************************

  ; update-output:
  update-output

  ; end simulation:
  if (extinction? = true or (ticks - time-to-balance) > time-limit) [ stop ]

end






;**************************************************************************************
;        TO SET-GLOBAL-PARAMETERS
;**************************************************************************************
to set-global-parameters

  set genetic-mean         0.0
  set env-effect-mean      0.0 ; standard model and implicit genetics only
  set initial-env-optimum  0.0
   ;**********************************************************************
  ; NEEDS TO BE CHANGED IF EXTENDING THE MODEL TO SIMULATE HETEROGENEOUS /
  ; FRAGMENTED LANDSCAPES
  ;**********************************************************************
  set size-of-environment 10
  ;**********************************************************************
  set female-male-ratio    0.5 ; 0.5 means random set of sex
  set extinction? false
  ; set the variance of the random environmental effect "ve" on the development of the
  ; trait. If how-plasticity? = "standard-model", then ve is computed according to given
  ; h2 and gv as in Bjoerklund et al (2009).
  if (how-plasticity? = "standard-model")
  [ set env-effect-variance   compute-env-effect-variance (heritability) (genetic-variance)]

  ; set level of density compensation, which is simulated as in Bjoerklund et al (2009),
  ; according to the density dependence effeect selected by the user: this impacts the
  ; population dynamics. For example, monotonic or oscillatory population dynamics.

  set density-compensation density-dependence-effect

  ; DEBUG:
  ;write "density dependence effect: " print density-dependence-effect
  ;write "value of density compensation: " print density-compensation


  ; set the strength of selection according to fitness function and type of organism
  ; whether specialist, moderate, generalist
  if-else (fitness-function = "Bjoerklund2009")
  [
    ; values are set as in Bjoerklund et al (2009)
    if-else (type-organism = "specialist") [ set strength-selection 10 ]
    [ ;else
      if-else (type-organism = "moderate") [ set strength-selection 20 ]
      [;else
        if-else (type-organism = "generalist") [ set strength-selection 40 ]
        [;else, catch exception
          error "Unidentified type of organism ..." stop
        ]
      ]
    ]
  ]
  [;else: exponential function
    ; values are set as in Burger and Lynch (1995): strength of selection
    if-else (fitness-function = "negative-exponential")
    [
      if-else (type-organism = "specialist") [ set strength-selection 1 ]
      [ ;else
        if-else (type-organism = "moderate") [ set strength-selection 2.2 ]
        [;else
          if-else (type-organism = "generalist") [ set strength-selection 3.2 ]
          [;else, catch exception
            error "Unidentified type of organism ..." stop
          ]
        ]
      ]
    ]
    [;else; catch exception
      error "Unidentified fitness function" stop
    ]
  ]

  ; DEBUG
  ;write "strength of selection: " print strength-selection
  ;write "env-effect-variance: " print env-effect-variance

  ; end of parameter values

end







;**************************************************************************************
;        TO SET-INITIAL-CONDITIONS-IMPLICIT-GENETICS
;**************************************************************************************
; This function set the initial conditions for both, the genetic and environmental
; components of the phenotype for each individual (i.e., turtle).
; Here the genetics is implicitly modelled.
; Important!: random-normal function uses the standard deviation = sqrt(variance) of the
; distribution.
; The function works in a turtle context, example: ask turtles [ set-initial-cond... ]
to set-initial-conditions-implicit-genetics

  ; genetic component 'a':
  set genetic-component random-normal genetic-mean sqrt (genetic-variance)

  ; environmental component 'e':
  set environmental-effect random-normal env-effect-mean sqrt (env-effect-variance)

  ; DEBUG:
  ;write "turtle_" write who print ": "
  ;write "genetic component: " print genetic-component
  ;write "environ component: " print environmental-effect

end







;**************************************************************************************
;        TO SET-INITIAL-CONDITIONS-EXPLICIT-GENETICS
;**************************************************************************************
; this function sets the initial conditions for the genetic and environmental component
; of each individual turtle.
; The genetics is modelled explicitly, considering the number of loci.
; Three methods for the initialization of allele values were implemented. Only one must
; be active.
; alleles can be either:
; - real allele values set according to normal distribution as in Vincenzi (2014)
; - real allele values randomly set
; - random set of ones and zeros (on / off)
; (loci effect on phenotype is assumed to be additive).
; This function works in a turtle context, example: ask turtles [set-initial-condi... ]
to set-initial-conditions-explicit-genetics

  ; set allele values for each locus. Notice that the two dna strains are simulated
  ; separately

  ; population is assumed locally adapted. Allele values from a normal distribution
  ; as in Vincenzi (2014)
  ; N(0, va), where mean = 0 = initial environmental optimum, and
  ; va is the additive genetic variance per locus at the start of the simulation, given by
  ; VG (genetic-variance) / 2*L, L is the number-of-loci, and VG is the additive genetic
  ; variance of the quantitative trait at the start of the simulation
  ; (THIS METHOD is based on the continuum-of-alleles model of Kimura 1970: an introduction
  ; to population genetics theory (cited in Vincenzi 2014) )
  ; Vincenzi (2014) simulate additive effect through the addition of allele values
  set dna-strain1 n-values number-of-loci [ random-normal 0
                                           sqrt (genetic-variance / (2 * number-of-loci)) ]
  set dna-strain2 n-values number-of-loci [ random-normal 0
                                           sqrt (genetic-variance / (2 * number-of-loci)) ]

  ; alleles take values from uniform distribution in range (-1, 1)
  ;set dna-strain1 n-values number-of-loci [ (random-float 2) - 1 ]
  ;set dna-strain2 n-values number-of-loci [ (random-float 2) - 1 ]

  ; another method: alleles take values of 1 or 0 (on / off)
  ; this method might not be appropriated for small number of loci
  ; set dna-strain1 n-values number-of-loci [ random 2 ]
  ; set dna-strain2 n-values number-of-loci [ random 2 ]

  ; DEBUG
  ;type "dna-strain1: " print dna-strain1
  ;type "dna-strain2: " print dna-strain2
  ;show item (number-of-loci - 1) dna-strain1

  ; set the value for the genetic component 'a' according to the genome
  ; loci effect is assumed to be additive
  set genetic-component (sum dna-strain1) + (sum dna-strain2)

  ; set environmental component 'e':
  set environmental-effect random-normal env-effect-mean sqrt (env-effect-variance)

  ; DEBUG:
  ;write "turtle_" write who print ": "
  ;write "genetic component: " print genetic-component
  ;write "environ component: " print environmental-effect
end






;**************************************************************************************
;        TO UPDATE-PHENOTYPE
;**************************************************************************************
; This function calculates (updates) the phenotype 'z' of the turtle following
; z = a + e
; where 'a' is the genetic-component and 'e' the environmental-effect on the phenotype
; The function works in a turtle context, example: ask turtles [ update-phenotype ]
to update-phenotype

  set phenotype genetic-component + environmental-effect

  ;DEBUG:
    ;write "phenotype z of turtle " write who print ":"
    ;print phenotype
    ;write "genetics: " print genetic-component
    ;write "env-effect: " print environmental-effect

end






;**************************************************************************************
;        TO UPDATE-PHENOTYPIC-RESPONSE
;**************************************************************************************
; updates the phenotype z according to the plastic response (selected method).
; The methods are similar, differing only in the parameter b, which is a sinusoidal
; function in method 2 and 3. This imposes limits to the plastic response of z depending on
; the amount of change that the organism is currently experimenting.
; Mehod 1 is based on the typical reaction norm approach, with a small modification.
; Methods assume a reference environment Q* where no plasticity is needed, and
; the phenotype z develops according to the genetic component a.
; 1) method 1: linear reaction norm
; z = a + bQt; where a is the genetic component (also breeding value); Qt the optimum
; phenotype as given by the environment; and b governs the plastic response (the slope)
;
; 2) method 2 and 3: sinusoidal reaction norms
; z = a + bQt; b = sin(abs Qt - a) These methods accounts for adaptive phenotypic
; plasticity, and consider that the plastic response has its limits
;
; 3) method 4: random plasticity
; common approach (e.g., Vincenzi 2014), where a value is draw randomly from a normal
; distribution
to update-phenotypic-response

  let a genetic-component
  let Qt [optimum] of patch-here
  let b 0

   if-else (how-plasticity? = "linear-RN" ) ; method 1
   [
     ;set b 1
     set b slope
   ]
   [ ;else method 2
     if-else (how-plasticity? = "adaptive-sinusoidal")
     [
       ; this method uses the amount of change defined as the departure
       ; from the reference environment Q* (Q* = Qt, such that Qt - a = 0)
       set Qt Qt - a
       ; in this method, 1 / slope affects the plastic ability of the organism,
       ; such that small slope less plasticity, large slope, high plasticity,
       ; just as for the linear reaction norm
       let omega 1 / slope
       ;
       ; In Netlogo sin function assumes angle is given in degrees
       if (abs (omega * Qt) < pi)
       [
         let angle (abs (omega * Qt) * 180) / pi
         set b sin (angle)
       ]
     ] ; no negative effect
     [ ; else method 3
      if-else (how-plasticity? = "adaptive-logistic")
      [
        ; this method is similar to method 2, but assumes that the plastic
        ; response saturates after certain amount of change (similar to
        ; logistic function
        set Qt Qt - a
        ; in this method, 1 / slope affects the plastic ability of the organism,
        ; such that small slope less plasticity, large slope, high plasticity,
        ; just as for the linear reaction norm
        let omega 1 / slope
        ; check condition for saturation of the plastic response
        if (abs (omega * Qt) > pi / 2)
        [
          set Qt (pi / (2 * omega)) * (abs Qt) / Qt ; keep the sign
        ]
        ; In Netlogo sin function assumes angle is given in degrees
        let angle (abs (omega * Qt) * 180) / pi
        set b sin (angle)
      ]
      [ ;else method 4
        if-else (how-plasticity? = "random-noise")
        [
          set b random-normal 0 sqrt 1
          ; random plasticity is assumed noise around the breeding value a
          ; thus
          set Qt 1
        ]
        [; else catch exception
          error "Unidentified method of plasticity ..." stop
        ]
      ]
     ]
   ]

  ; update plastic response
  set phenotype a + (b * Qt)

  ; DEBUG:
   ;write "phenotype of" write who write ": " print phenotype
   ;write "a: " print a
   ;write "b: " print b
   ;write "Qt: " print Qt


end






;**************************************************************************************
;        TO UPDATE-STAGE-SEX
;**************************************************************************************
; update the stage of the individuals to adulthood, and set sex randomly.
; The function works in a turtle context, example: ask turtles [ update-stage-sex ]
to update-stage-sex

  ; set stage to adulthood:
  set stage "adult"
  ; set reproduced? to false
  set reproduced? false
  ; set sex randomly:
  if-else (random-float 1 <= female-male-ratio)
  [ set sex "female" ]
  [ set sex  "male"  ] ; else, male

end






;**************************************************************************************
;        TO CHECK-FITNESS-BJOERKLUND2009
;**************************************************************************************
; Fitness wi of individual i is calculated according to Björklund et al. (2009):
; wi = 1 - [(zi - tita(t))^2]/gamma
; where tita(t) is the optimum phenotype as given by the environment at time t, and
; gamma is the strength of selection
; the maximum fitness is 1
; The function works in a turtle context
to check-fitness-Bjoerklund2009

  if-else([pcolor] of patch-here = green)
  [
    set fitness (1 - (((phenotype - [optimum] of patch-here) ^ 2) / strength-selection))
  ]
  [ error "turtle outside green patches" stop ] ; change this for heterogeneous landscapes

  ; DEBUG
  ;write "fitness of turtle " write who write ": " print fitness
  ;write "optimum of patch-here: " print [optimum] of patch-here
  if (fitness > 1)
  [
    write "phenotype of turtle : " write who write ": " print phenotype
    write "optimum: " print [optimum] of patch-here
    write "gamma: " print strength-selection
    write "fitness: " print fitness
    error "warning! fitness value is greater than 1" stop
  ]

end






;**************************************************************************************
;        TO CHECK-FITNESS-NEGATIVE-EXPONENTIAL
;**************************************************************************************
; Fitness wi of individual i is calculated according to a negative exponential function
; as in Burger & Lynch (1995):
; wi = exp [ - (zi - tita(t))^2]/2*gamma^2
; where tita(t) is the optimum phenotype as given by the environment at time t, and
; gamma is the strength of selection.
; wi is in range (0, 1) => maximum fitness of 1
; The function works in a turtle context
to check-fitness-negative-exponential

  let tita 0

  if-else([pcolor] of patch-here = green)
  [ set tita [optimum] of patch-here ]
  [ error "turtle outside green patches!" stop]

  let gamma strength-selection

  set fitness exp ( - ( (phenotype - tita) ^ 2) / (2 * (gamma) ^ 2) )

  ; DEBUG
  ;write "fitness of turtle " write who write ": " print fitness
  ;write "optimum of patch-here: " print [optimum] of patch-here

  ; catch error:
  if (fitness > 1)
  [
    write "phenotype of turtle : " write who write ": " print phenotype
    write "optimum: " print [optimum] of patch-here
    write "gamma: " print strength-selection
    write "fitness: " print fitness
    error "warning! fitness value is greater than 1" stop
  ]

end






;**************************************************************************************
;        TO SET-FECUNDITY-BJOERKLUND
;**************************************************************************************
; Fecundity lambda is calculated according to Björklund et al. (2009) for each
; individual:
; lambda = wi*exp[alpha(1 - N/K)]
; where wi is the fitness of individual i, alpha is the density-compensation, N the
; population size, K the carrying capacity and exp the exponential function
; The function works in a turtle context
to set-fecundity-bjoerklund

  let alpha density-compensation
  let N count turtles with [ stage = "adult" ]
  let K carrying-capacity
  let wi fitness

  set fecundity wi * ( exp ( alpha * (1 - N / K ) ) )

  ; DEBUG
  ;write "fecundity of turtle " write who write ": " print fecundity
  ;write "population size: " print N

end







;**************************************************************************************
;        TO REPRODUCE-IMPLICIT-GENETICS
;**************************************************************************************
; This function implements sexual reproduction according to Björklund et al. (2009).
; Female individuals randomly select a partner of opposite sex to mate. The fitness of
; a pair is the sum of the fitness value of the two parents. Then fecundity is calculated
; considering density dependence. The number of offspring is drawn from a poisson
; distribution centered on the value of fecundity.
;
; Genetics is implicit according to the infinitesimal model of quantitative genetics which
; assumes that traits are affected by a large number of loci of additive effects.
; Therefore trait inheritance can be approximated using a normal distribution with mean
; parents trait value, and variance: half the genetic variance.
; Mutations and recombination are implicit.
;
; The function works in a turtle context, example: ask turtles [ reproduce ]
to reproduce-implicit-genetics

  ; SEXUAL REPRODUCTION:
  ; pick a random partner of opposite sex:
  let my-partner one-of turtles with [ stage = "adult" and
                                       sex = "male"   ]

  ; test whether there is a partner available to reproduce:
  if ( my-partner != nobody )
  [
    ; DEBUG:
    ;write "me "print who
    ;write "partner " print [who] of my-partner
    ;type "my-sex: " print sex
    ;type "partner sex: " print [sex] of my-partner
    ;if(reproduced? = true) [ print "me"]
    ;if([reproduced?] of my-partner = true) [print "partner"]

    ;create a list containing the genetic-component of both parents:
    let list-genetic-parents list (genetic-component) ([genetic-component] of my-partner)

    ;DEBUG:
    ;write "my-genetic-component: " print genetic-component
    ;write "my-partner-genetic-comp: " print [genetic-component] of my-partner
    ;write "list-genetic-components-of-parents: " print list-genetic-parents

    ; calculate the mean parental genetic-component which will be inherited by offspring
    let genetic-mean-of-parents mean ( list-genetic-parents )

    ;DEBUG:
    ;write "my genetic-component: " print genetic-component
    ;write "partner genetic-component: " print [genetic-component] of my-partner
    ;write "genetic-mean: " print genetic-mean-of-parents

    ; calculate genetic variance among parents:
    let genetic-variance-of-parents variance ( list-genetic-parents )
    ; DEBUG:
    ;write "genetic-variance-of-parents: " print genetic-variance-of-parents
    ;type "adults "print count turtles with [stage = "adult"]
    ;type "total " print count turtles
    ; set the genetic-variance for the offspring according to the selected method
    let my-genetic-variance 0

    if-else ( how-genetic-variance = "parental-level")
    [
      set my-genetic-variance (genetic-variance-of-parents);parental-level
                                                           ;(Bjoerklund2009)
    ]
    [; else,
      if-else (how-genetic-variance = "parameter")
      [ set my-genetic-variance genetic-variance ] ; as in Reed et al. (2011):
      ;the population-level additive genetic variance is an input parameter (set by user)
      [; else
        if-else (how-genetic-variance = "population-level")
        [
        ; method as in Vincenzi & Piotti (2014): genetic variance is the total
        ; additive genetic variance for the trait at the population level
          set my-genetic-variance pop-gv
        ]
        ; else, catch exception
        [ error "undefined method for genetic variance" stop ]
      ]
    ]
    ; the genetic variance is simulated as in Vincenzi et al (2012). The infinitesimal
    ; model is modified to account for the decline in additive genetic variance and the
    ; new input of variation through mutation (mutational-variance).
    ; This method also accounts for offspring additive genetic variance equals half the
    ; additive genetic variance
    set my-genetic-variance (1 / 2) * (my-genetic-variance + mutational-variance)

    ; calculate the fecundity of the breeding pair pair-fecundity as in Björklund et al.
    ; Björklund(2009) assumed that the fitness of the pair was the sum of the fitness
    ; values of the two parents (wsum) which is equivalent to the sum of their fecundities.
    ; Given that wsum is used, the population is saturated above the carrying capacity.
    let fecundity-of-pair ( fecundity + [ fecundity ] of my-partner )

    ;DEBUG:
    ;write "my-fecundity: " print fecundity
    ;write "fecundity-of-partner: " print [fecundity] of my-partner
    ;write "fecundity-of-pair: " print fecundity-of-pair

    let number-of-offspring random-poisson fecundity-of-pair
    ; DEBUG
    ;write "nr. offspring: " print number-of-offspring

    if (number-of-offspring >= 1) ; reproduction occurs
    [
      hatch number-of-offspring
      [ ; according to netlogo dictionary, each new turtle inherits all its
        ; variables, including its location, from its parent, except [who]

        ; set stage as juvenile
        set stage "juvenile"
        ; set random orientation of the newborn and move to a green patch
        set heading random 360
        move-to one-of patches with [pcolor = green]

        set genetic-component random-normal genetic-mean-of-parents ; mean
                                sqrt my-genetic-variance

        ; DEBUG
        ;type "genetic-component of offspring: " print genetic-component
        ;type "genetic-variance: "  print my-genetic-variance

         ; environmental component 'e':
        ; Set the variance of environmental effect ve according to method of heritability
        ; if h2 = fixed, compute the variance of environmental effect according to the
        ; genetic variance and level of heritability
        ; else, ve is a constant parameter
        let ve 0 ; initialization
        if (how-plasticity? = "standard-model")
        [ set ve compute-env-effect-variance (heritability) (my-genetic-variance)]

        set environmental-effect random-normal env-effect-mean sqrt ve

        ; DEBUG:
        ;write "genetic-variance: " print my-genetic-variance
        ;write "env-effect-variance: " print ve
      ] ; end of inheritance

   ] ; end of hatching n offspring ( reproduction)

    ; update reproduced? status of the partner to true
    ask my-partner [ set reproduced? true ]
    ;type "genetic-component of partner " print [genetic-component] of my-partner
  ]; end of if statement: whether a partner is available for reproduction

  ; update reproduced? of this turtle or individual to true:
    set reproduced? true
    ;type "my-genetic-component " print genetic-component

end






;**************************************************************************************
;        TO REPRODUCE-EXPLICIT-GENETICS
;**************************************************************************************
; This function implements sexual reproduction and lottery polygyny.
; Females randomly select a partner of opposite sex to mate. The fitness of a pair
; is the sum of the fitness value of the two parents. Then fecundity is calculated
; considering density dependence. The number of offspring is drawn from a poisson
; distribution centered on the value of fecundity.
;
; Genetics is explicit. This means that loci are explicitly modelled. Loci effect on
; trait value is assumed to be additive.
; Recombination is simulated through the random selection of parental allele values.
; Then, mutations take place with probability mut-rate.
;
; The function works in a turtle context, example: ask turtles [ reproduce ]
to reproduce-explicit-genetics

  ; store the loci (dna strains) of this turtle (parent 1)
  let my-loci1 dna-strain1
  let my-loci2 dna-strain2

  ; SEXUAL REPRODUCTION:
  ; pick a random partner of opposite sex:
  let my-partner one-of turtles with [ stage = "adult" and
                                       sex = "male" ]

  ; test whether there is a partner available to reproduce:
  if ( my-partner != nobody )
  [
    ; DEBUG:
    ;write "me "print who
    ;write "partner " print [who] of my-partner
    ;type "my-sex: " print sex
    ;type "partner sex: " print [sex] of my-partner
    ;if(reproduced? = true) [ print "me"]
    ;if([reproduced?] of my-partner = true) [print "partner"]

    ; store loci of my-partner (parent 2)
    let partner-loci1 [dna-strain1] of my-partner
    let partner-loci2 [dna-strain2] of my-partner

    ;DEBUG:
    ;write "Loci parent 1: " write my-loci1 write " and " print my-loci2
    ;write "Loci parent 2: " write partner-loci1 write " and " print partner-loci2

    ; calculate the fecundity of the breeding pair pair-fecundity as in Björklund et al.
    ; Björklund(2009) assumed the fitness of the pair as the sum of the fitness
    ; values of the two parents (wsum) which is equivalent to the sum of their fecundities.
    let fecundity-of-pair ( fecundity + [ fecundity ] of my-partner )

    ;DEBUG:
    ;write "my-fecundity: " print fecundity
    ;write "fecundity-of-partner: " print [fecundity] of my-partner
    ;write "fecundity-of-pair: " print fecundity-of-pair

    ; set the genetic variance. Only for the calculation of environmental variance
    ; according to the heritability
    let genetic-variance-of-parents variance list (genetic-component)
                                                  ([genetic-component] of my-partner)
    let my-genetic-variance 0

    if-else ( how-genetic-variance = "parental-level")
    [
      set my-genetic-variance (genetic-variance-of-parents); parental-level
                                                           ; (Bjoerklund2009)
    ]
    [; else,
      if-else (how-genetic-variance = "parameter")
      [ set my-genetic-variance genetic-variance ] ; as in Reed et al. (2011):
        ;the population-level additive genetic variance is an input parameter (set by user)
      [; else
        if-else (how-genetic-variance = "population-level")
        [
        ; method as in Vincenzi & Piotti (2014): genetic variance is the total
        ; additive genetic variance for the trait at the population level
          set my-genetic-variance pop-gv
        ]
        ; else, catch exception
        [ error "undefined method for genetic variance" stop ]
        ]
    ]
    ; additive genetic variance of offspring equals half the additive genetic variance:
    set my-genetic-variance (1 / 2) * (my-genetic-variance)

    let number-of-offspring random-poisson fecundity-of-pair
    ; DEBUG
    ;write "nr. offspring: " print number-of-offspring

    if (number-of-offspring >= 1) ; reproduction occurs
    [
      hatch number-of-offspring
      [
        ; set stage as juvenile
        set stage "juvenile"
        set color blue
        ; set random orientation of the newborn and move to a green patch
        set heading random 360
        move-to one-of patches with [pcolor = green]


        ;*********************************************************************************
        ; INHERITANCE GENETIC COMPONENT OF PHENOTYPE
        ;*********************************************************************************
        ; inheritance with recombination and mutations, explicit genetics:

        ;1st mutation takes place during the creation of gametes inside each parent
        ; Each locus is checked for mutations according to the given mutation rate.
        ; The effect-size of the mutation is implemeted using a normal distribution
        ; N(0, variance), where the variance controls for the effect-size (Vincenzi 2014).
        ;
        ; A mutation can be explicitly simulated per locus or as a single mutation per
        ; strain with probability number-of-loci * mut-rate-per-locus, as in Vincenzi (2014)
        ; the effect-size of mutations can be modified in the future. For example,
        ; Vincenzi (2014) based the effect-size of mutations on a normal distribution
        ; N(0, effect-size). Through modifications of the variance of the distribution
        ; one can control the effect size of the mutation. Note that using the normal
        ; distribution, the expected value is the mean (in the example, 0 effect).

        ;DEBUG:
        ;write "parent1 strain1: " print my-loci1
        ;write "parent1 strain2: " print my-loci2

        ;parent1
        let l 0 ; counter
        ; copy strains info to avoid modifications of the original variable
        let strain1-parent1 my-loci1
        let strain2-parent1 my-loci2
        while [l < number-of-loci]
        [
          if (random-float 1 < mut-rate-per-locus) ; mutation occurs
          [ ; dna-strain1 = my-loci1 = strain1-parent1

            let dummy random-normal mu-dist-mut-effects sqrt mut-effect-size

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
           ;if([mean-env-optimum] of patch-here < phenotype) ; (optimum push for smaller trait value)
           ;[ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the desired probability of beneficial mutations
           ;] ; else if ([mean-env-optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain1-parent1
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain1-parent1 )
            set strain1-parent1 replace-item l strain1-parent1 dummy

            ; DEBUG
            ;write "value after: " print item l strain1-parent1
          ]
          if (random-float 1 < mut-rate-per-locus) ; mutation occurs
          [ ; dna-strain2 = my-loci2 = strain2-parent1

            let dummy random-normal mu-dist-mut-effects sqrt mut-effect-size

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            ;if([mean-env-optimum] of patch-here < phenotype) ; optimum at the left
                                                             ;(optimum push for smaller trait value)
            ;[ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the desired probability of beneficial mutations
            ;] ; else if ([mean-env-optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain2-parent1
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain2-parent1 )
            set strain2-parent1 replace-item l strain2-parent1 dummy

            ; DEBUG
            ;write "value after: " print item l strain2-parent1
          ]
          set l (l + 1)
        ]; while loop

        ;parent2
        set l 0 ; counter
        let strain1-parent2 partner-loci1
        let strain2-parent2 partner-loci2
        while [l < number-of-loci]
        [
          if (random-float 1 < mut-rate-per-locus) ; mutation occurs
          [ ; dna-strain1 = partner-loci1= strain1-parent2

            let dummy random-normal mu-dist-mut-effects sqrt mut-effect-size

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            ;if([mean-env-optimum] of patch-here < phenotype) ; optimum at the left
                                                             ;(optimum push for smaller trait value)
            ;[ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the probability of beneficial mutations
            ;] ; else if ([mean-env-optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain1-parent2
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain1-parent2 )
            set strain1-parent2 replace-item l strain1-parent2 dummy

            ; DEBUG
            ;write "value after: " print item l strain1-parent2
          ]
          if (random-float 1 < mut-rate-per-locus) ; mutation occurs
          [ ; dna-strain2 = partner-loci2 = strain2-parent2

            let dummy random-normal mu-dist-mut-effects sqrt mut-effect-size

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            ;if([mean-env-optimum] of patch-here < phenotype) ; optimum at the left
                                                             ; (optimum push for smaller trait value)
            ;[ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the probability of beneficial mutations
            ;] ; else if ([mean-env-optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain2-parent2
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain2-parent2 )
            set strain2-parent2 replace-item l strain2-parent2 dummy

            ; DEBUG
            ;write "value after: " print item l strain2-parent2
          ]
          set l (l + 1)
        ]; while loop

        ; 2nd inheritance with recombination:
        set l 0 ; counter
        set dna-strain1 n-values number-of-loci [0] ; initialization
        set dna-strain2 n-values number-of-loci [0]
        while [l < number-of-loci]
        [
          ; offspring strain 1, from parent 1
          set dna-strain1 replace-item l dna-strain1 (one-of (list item l strain1-parent1
                                                                   item l strain2-parent1))
          ; offspring strain 2, from parent 2
          set dna-strain2 replace-item l dna-strain2 (one-of (list item l strain1-parent2
                                                                   item l strain2-parent2))
          set l (l + 1)
        ] ; while loop
        ; DEBUG
        ;write "offspring strain1: " print dna-strain1
        ;write "offspring strain2: " print dna-strain2

        ; end of inheritance explicit genetics
        ;***********************************************************************************

        ;*********************************************************************************

        ; update the genetic component of the offspring
        set genetic-component (sum dna-strain1) + (sum dna-strain2)

        ; DEBUG
        ;type "genetic-component: " print genetic-component

        ; environmental component 'e':
        ; Set the variance of environmental effect ve according to method of heritability
        ; if h2 = fixed, compute the variance of environmental effect according to the
        ; genetic variance and level of heritability
        ; else, ve is a constant parameter
        let ve 0 ; initialization
        if (how-plasticity? = "standard-model")
        [ set ve compute-env-effect-variance (heritability) (my-genetic-variance)]

        set environmental-effect random-normal env-effect-mean sqrt ve

        ; DEBUG:
        ;write "genetic-variance: " print my-genetic-variance
        ;write "env-effect-variance: " print ve
      ] ; end of inheritance

   ] ; end of hatching n offspring ( reproduction)

    ; update reproduced? status of the partner to true
    ask my-partner [ set reproduced? true ]

  ]; end of if statement: whether a partner is available for reproduction

  ; update reproduced? of this turtle or individual to true:
    set reproduced? true

end





;**************************************************************************************
;        TO UPDATE-ENVIRONMENT
;**************************************************************************************
; Currently, there are two main environmental scenarios: climate-change and cyclic,
; respectively.
; In climate change scenario,
; the mean environmental optimum changes at rate r every iteration. For the model, each
; iteration is equivalent to a generation
; The optimum Q updates according to:
; 1) Qt = Q0 + rt (Directional deterministic)
; 2) Qt = Q0 + rt; Qt* = Qt + E, (Directional stochastic)
; where E is the type of noise defined as:
; Et = aE(t-1) + b*Dt, where D = N(0, 1); here, a is not the genetic component!
; Here b is assumed to be b = sqrt VEt*(1 - a^2), as in Schwager et al (2006), where
; VEt is the environmental variance at time t; and the
; parameter a (set by the user), determines the strength of the autocorrelation
; (based on Ripa and Lundberg 1996, and Schwager et al. 2006)
;      a = 0 (white noise)
;  0 < a < 1 (red noise) Björklund et al used a = 0.7
; -1 < a < 0 (blue noise) Björklund et al used a = -0.7
;
; The parameter VEt can change or not in time, depending on the rate of change k of
; the variance. Thus, VEt = VE + kt, where VE is the inital environmental variance.
; This consideration allows to simulate the increased probability of extreme events
; as predicted by climatic IPCC scenarios, as in Vincenzi et al (2012).
; The function works in a patch context, example: ask patches [ update-optimum ]
;
; In the cyclic environment,
; The optimum Qt changes according to a sinusoidal function that considers two
; parameters:
; the amplitude A, and the period T, thus:
; Qt = A.sin(2.pi.t / T), according to Burger and Krall (2004)
to update-environment

  ; initial optimum of the environment Q0 = tita-0 according to Björklund2009
  let tita-0 0
  let new-optimum 0
  let new-noise 0
  let VE env-variance
  let k rate-change-of-env-variance
  let autocorr level-autocorr
  let b sqrt (1 - (autocorr ^ 2)) ; scaling factor of the variance according to
                           ; Schwager et al (2006)

  ; first catch exception:
  if (scenario? != "climate-change" and scenario? != "cyclic")
  [ error "update-environment: undefined scenario of environment" stop ]

  ; update the optimum according to the scenario of environment
  if (scenario? = "climate-change")
  [
    ; set parameter values for the scenario climate-change
    let r rate-change-of-optimum

    ; update moving mean envirnomental optimum
    if (ticks > time-to-balance)
    [ set new-optimum tita-0 + (r * (ticks - time-to-balance) ) ]; 1) directional deterministic
  ]

  if (scenario? = "cyclic") ; according to Burger and Krall 2004
  [
    ; set parameter values for the cyclic scenario
    let A amplitude ; amplitude of the wave
    let T period ; period of the wave T = 1 / fr; fr is the frequency

    ; update the environmental optimum
    set new-optimum ( A * sin ( (2 * pi * (ticks - time-to-balance)) / T) )
  ]

  ; check for increasing environmetal variance and update the environmental variance
  if (ticks > time-to-balance)
  [set VE VE + (k * (ticks - time-to-balance) )]

  ; update the value of the mean-env-optimum of the patch before adding stochastic noise
  set mean-env-optimum new-optimum

  ; apply stochasticity around the environmental optimum
  ; stochastic environment with potentially different colour of noise
  ; depending on the level of outocorrelation a (level-autocorr)
  ;type "noise before " print noise
  set new-noise (autocorr * noise) + (b * (sqrt VE) * (random-normal 0 sqrt 1))
  set new-optimum new-optimum + new-noise
  ;type "noise after " print noise

  ; update optimum and noise of the patch
  set optimum new-optimum
  set noise new-noise

  ;**********************************************************************
  ; NEEDS TO BE CHANGED IF EXTENDING THE MODEL TO SIMULATE HETEROGENEOUS /
  ; FRAGMENTED LANDSCAPES
  ;**********************************************************************
  ask patches with [distance myself <= size-of-environment]
  [
    set optimum new-optimum
    set noise new-noise
  ]
  ;**********************************************************************

  ; DEBUG
  ;type "new optimum: " print optimum
  ;type "new noise: "   print noise
  ;type "patch optimum " print optimum
  ;type "patch noise: " print noise
  ;type "env-variance: " print env-variance
  ;type "new env-variance: " print VE

end






;**************************************************************************************
;        TO-REPORT COMPUTE-ENV-EFFECT-VARIANCE
;**************************************************************************************
; This function reports the corresponding value of environmental effect variance given the
; specified values of heritability and genetic variance
; h2: narrow-sense heritability
; gv: additive genetic variance
to-report compute-env-effect-variance [h2 gv]

  ; catch exception h2 <= 0:
  if-else (h2 <= 0)
  [ error "heritability muss be greater than 0 (h2 > 0)" ]
  [ report (gv / h2 ) - gv ] ; else

end






;**************************************************************************************
;        TO-REPORT COMPUTE-DEGREE-MALADAPTATION
;**************************************************************************************
; this function compute the degree of maladaptation according to Björklund et al (2009).
; defined as:
; the sum(from i = 1 to t = 30) of [(average_z - opt)^2] / gamma
to-report compute-degree-maladaptation

 if-else ( count turtles > 0 )
 [
   let z mean [ phenotype ] of turtles
   let opt [ optimum ] of patch 0 0
   let gamma strength-selection

   ; DEBUG:
   ;write "mean z: " print z
   ;write "optimum: " print opt

   ; compute the degree-maladaptation
   report ( ( z - opt ) ^ 2 ) / gamma
 ]
 [ report 0 ] ; else (extinct population)

end





;**************************************************************************************
;        TO UPDATE-OUTPUT
;**************************************************************************************
; this function plots the distribution of phenotypes in the population
to update-output

  ;type "total " print count turtles
  ;type "adults " print count turtles with [stage = "adult"]
  ;type "juveniles " print count turtles with [stage = "juvenile"]

  ; plot time series:
  set-current-plot "Time series"
  set-current-plot-pen "turtles"
  plot count turtles

  ; plot environmental optimum vs mean phenotypic response
  set-current-plot "environmental optimum vs mean phenotypic response"
  set-current-plot-pen "optimum"
  plot [optimum] of patch 0 0
  set-current-plot-pen "mean z"
  if (count turtles > 0)
  [ plot mean [phenotype] of turtles ]


  ; plot genetic variance
  set-current-plot "genotypic-phenotypic-variance"
  set-current-plot-pen "genotype"
  let dummy 0
  if(count turtles > 1)
  [ set dummy variance [genetic-component] of turtles ]
  plot dummy

  ; plot phenotypic variance in same plot of genetic-variance
  set-current-plot-pen "phenotype"
  set dummy 0
  if (count turtles > 1)
  [ set dummy variance [phenotype] of turtles ]
  plot dummy


  ; store the highest current phenotipic value in the population
  ; this value will be used below to set the range of the x axis
  ;let highest-value [phenotype] of max-one-of turtles [phenotype]
  set-current-plot "frequency distribution"


  ;set-plot-x-range -1 end-point-of-environment
  set-plot-pen-mode 1        ; (0 for lines, 1 for bars)
  set-histogram-num-bars 10
  ; plot the phenotypic frequency in the population:
  set-current-plot-pen "phenotype"
  histogram [phenotype] of turtles
  ; plot the genotypic frequency in the population:
  set-current-plot-pen "genotype"
  histogram [genetic-component] of turtles
  ; plot the environmental optimum
  set-current-plot-pen "env-optimum"
  plot-pen-reset
  set-plot-pen-mode 2
  ; plot a vertical line to show the current optimum phenotype in the environment
  let n 0
  while [n < carrying-capacity]
  [
    plotxy [optimum] of one-of patches with [pcolor = green] n
    set n n + 0.1
  ]

end






;**************************************************************************************
;        TO PLOT-MEAN-FITNESS
;**************************************************************************************
; plot mean fitness of turtles
to plot-mean-fitness

  ;type "total " print count turtles
  ;type "adults " print count turtles with [stage = "adult"]
  ;type "juveniles " print count turtles with [stage = "juvenile"]

  set-current-plot "degree of local adaptation"
  let dummy 0
  if (count turtles > 0 )
  [ set dummy mean [fitness] of turtles ]
  plot dummy

end






;**************************************************************************************
;        TO-REPORT MEAN-DIST
;**************************************************************************************
; this reporter returns the mean for the normal distribution "mean-dist"of mutation
; effect size according to the given proportion "p" of beneficial mutations (input
; parameter) and the mutation-effect-size or mutational variance "var". The mean-dist
; is used as the mean of the normal distribution from which a random variable is drawn.
; This reporter use the r-extension for Netlogo as it takes advantage from the qnorm
; function from the stats library in r
;to-report mean-dist [p var]
;
;  let x 0 ; point in the pdf. beneficial mutations are assumed to be any point => 0
;  r:put "p" p
;  r:put "sdev" sqrt var
;
;  report x - (r:get "sdev * qnorm(1 - p, mean = 0, sd = 1)")
;
;end






;**************************************************************************************
;        CHANGE LOG
;**************************************************************************************
;
;**************************************************************************************
;        MODEL VERSION: _V1.1
;**************************************************************************************
; The implementation of beneficial mutations was modified (see below) so it works for
; cyclic/seasonal environmental change
; a mutation in the direction of the optimum occurs with probability p =
; % beneficial mutations / 100
; if([mean-env-optimum] of patch-here < phenotype) ; optimum at the left (optimum push
;                                                  ; for smaller trait value)
; [ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
;                          ; account for the desired probability of beneficial mutations
; ] ; else if ([mean-env-optimum] of patch-here > phenotype) no need of further modification
@#$#@#$#@
GRAPHICS-WINDOW
1006
18
1186
199
-1
-1
5.212121212121212
1
10
1
1
1
0
1
1
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
7
10
71
43
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
44
398
191
frequency distribution
NIL
freq
-10.0
10.0
0.0
10.0
true
true
"" ""
PENS
"phenotype" 1.0 0 -13791810 true "" ""
"env-optimum" 1.0 2 -2674135 true "" ""
"genotype" 1.0 0 -16777216 true "" ""

BUTTON
73
10
177
43
Go! one itera
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
179
10
295
43
Run the model!
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
191
398
338
Time series
time in generations
N
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"turtles" 1.0 0 -13791810 true "" ""

SLIDER
803
36
959
69
carrying-capacity
carrying-capacity
10
1000
1000.0
10
1
NIL
HORIZONTAL

SLIDER
803
70
959
103
population-size
population-size
10
1000
1000.0
10
1
NIL
HORIZONTAL

CHOOSER
803
104
959
149
type-organism
type-organism
"specialist" "moderate" "generalist"
1

INPUTBOX
802
427
961
487
genetic-variance
0.2
1
0
Number

SLIDER
401
193
602
226
rate-change-of-optimum
rate-change-of-optimum
0
1
0.02
0.01
1
NIL
HORIZONTAL

CHOOSER
400
97
604
142
scenario?
scenario?
"climate-change" "cyclic"
0

SLIDER
678
320
795
353
heritability
heritability
0.1
1
1.0
0.1
1
NIL
HORIZONTAL

PLOT
963
345
1249
487
genotypic-phenotypic-variance
time in generations
genetic variation
0.0
2.0
0.0
1.0
true
true
"" ""
PENS
"genotype" 1.0 0 -16777216 true "" ""
"phenotype" 1.0 0 -13791810 true "" ""

CHOOSER
803
184
959
229
fitness-function
fitness-function
"Bjoerklund2009" "negative-exponential"
1

INPUTBOX
400
255
510
315
time-to-balance
0.0
1
0
Number

CHOOSER
403
442
534
487
how-genetic-variance
how-genetic-variance
"parameter" "parental-level" "population-level"
1

SLIDER
607
63
792
96
env-variance
env-variance
0
2
1.0
0.1
1
NIL
HORIZONTAL

PLOT
963
201
1249
345
degree of local adaptation
time in generations
fitness
0.0
5.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

CHOOSER
400
316
510
361
how-genetics?
how-genetics?
"implicit" "explicit"
1

INPUTBOX
535
427
649
487
mut-rate-per-locus
0.001
1
0
Number

INPUTBOX
649
427
801
487
mut-effect-size
0.2
1
0
Number

SLIDER
535
393
649
426
number-of-loci
number-of-loci
1
50
1.0
1
1
NIL
HORIZONTAL

TEXTBOX
405
12
796
47
....................................... ECOLOGY .......................................
14
54.0
1

TEXTBOX
519
269
793
303
...................... EVOLUTION ........................
14
93.0
1

SLIDER
297
10
398
43
time-limit
time-limit
10
1000
150.0
10
1
NIL
HORIZONTAL

SLIDER
607
97
792
130
level-autocorr
level-autocorr
-0.9
0.9
0.0
0.1
1
NIL
HORIZONTAL

SLIDER
607
131
792
164
rate-change-of-env-variance
rate-change-of-env-variance
0
0.1
0.0
0.01
1
NIL
HORIZONTAL

TEXTBOX
464
75
537
93
Scenario:
14
54.0
1

TEXTBOX
814
10
954
28
Type of Organism
14
14.0
1

TEXTBOX
548
363
792
381
.............. Explicit genetics .................
14
93.0
1

TEXTBOX
411
363
537
397
Implicit genetics
14
93.0
1

INPUTBOX
403
381
534
441
mutational-variance
0.01
1
0
Number

PLOT
7
337
398
487
environmental optimum vs mean phenotypic response
time
optimum
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"optimum" 1.0 0 -10899396 true "" ""
"mean z" 1.0 0 -14454117 true "" ""

CHOOSER
803
294
959
339
how-plasticity?
how-plasticity?
"standard-model" "random-noise" "linear-RN" "adaptive-sinusoidal" "adaptive-logistic"
0

SLIDER
803
340
959
373
slope
slope
0.5
2
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
667
42
766
60
Stochasticity
14
54.0
1

TEXTBOX
406
170
792
188
Directional / Climate change\t        Cyclic environmental change
14
54.0
1

SLIDER
615
192
787
225
amplitude
amplitude
0
4
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
615
227
787
260
period
period
0.1
4
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
405
42
605
63
....... The Environment .............
14
54.0
1

TEXTBOX
681
299
802
317
if standard model,
14
93.0
1

SLIDER
803
150
959
183
density-dependence-effect
density-dependence-effect
0
3
0.5
0.1
1
NIL
HORIZONTAL

TEXTBOX
806
268
953
286
            Plasticity
14
53.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
