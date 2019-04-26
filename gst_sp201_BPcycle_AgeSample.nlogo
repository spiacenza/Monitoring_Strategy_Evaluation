
;;; Green Sea Turtle IBM - The impact of individual variability in demographic rates on population recovery and Monitoring Strategy Evaluation
;extensions [ stats ]

;; Define global parameters
globals [
  year 
  timesHere
  timesHere2
  timesGo
  timesAge
  maxTicks
  numberturtles
  numberHatchlings
  numberPelagicJuveniles
  numberNeriticJuveniles
  numberSubadults
  numberNeophytes
  numberAdults
  rh ; intrinsic rate of increase for Ricker function for hatchling production
  Climate ; Climate index to simulate good and poor environmental conditions annually
  q-norm ; normal random variable for logit-normal distribution for selecting detection stochastic variable
  q ; catachability coefficient for monitoring nesters
  q-std-dev ; Standard deviation used for setting q
  scale-factor ; scaling factor for SuperIndividual of hatchlings and pelagicjuveniles for intialization of SAD
  SampleSize
  countNesters
  countNesters2
  Nesters
  Monitored-Nesters  
  Catch ; number of turtles removed during harvest
  Catch_hatchlings
  BreedProb
  Nests
  Monitored-Nests
  Monitored-Nests-rand
 ]

;; Define life history stages
breed [ Hatchlings Hatchling ] ; Hatchlings
breed [ PelagicJuveniles PelagicJuvenile ] ; Pelagic juveniles
breed [ NeriticJuveniles NeriticJuvenile ] ; neritic juveniles
breed [ Subadults Subadult ] ; subadults
breed [ Neophytes Neophyte ] ; Neophyte Nesters
breed [ Adults Adult ] ; adults
;; Define turtle specific parameters
turtles-own [ age AgeClass reprostatus reprostatus2 times-nested remigs remigs_old ClutchFrequency ClutchSize newHatchlings AgeMaturity
  ] ; ; reprostatus = nester or skip-nester, newHatchlings = reproduce new hatchlings , times-nested = count how many times in lifetime nested, remigs=how many years since last nesting, AgeClass = Hatchling, PelagicJuvenile NeriticJuvenile Subadult Neophyte Adult

;; Initialization procedures
to setup
  clear-all
  ;set Monitored-Nesters no-turtles
  ifelse DetectA = 3.841832 [set q-std-dev 2.417518] [ifelse DetectA = -1.643709E-9 [set q-std-dev 0.41646] 
    [ifelse DetectA = -2.216764 [set q-std-dev 0.221504] [set q-std-dev 0]]] ;ifelse DetectA >= 0.94 [set q-std-dev 0.025] [set q-std-dev (0.2165 * DetectA)]
  set q-norm random-normal DetectA q-std-dev ;while [(q <= 0) or (q > 1.0) ] [ set q random-normal DetectA q-std-dev ]; random-normal 0.79 0.17 while [(q <= 0.3) or (q >= 1.0) ] [ set q random-normal 0.75 0.25 ]   ; Catchability coefficient = p(t.) top ranked model MSORB Mark
  set q 1 / (1 + (exp (- q-norm)))
  set scale-factor 5.4 ; per life table
  set rh rhA ;0.15
  set numberHatchlings round (499311 / scale-factor) ; based on Stable age distribution from age structured Matrix projection model and Hatchlings and Pel. Juveniles scaled for Super-Individuals
  set numberPelagicJuveniles round (508129 / scale-factor); 
  set numberNeriticJuveniles 161297
  set numberSubadults 73388
  set numberNeophytes 1416
  set numberAdults 18909
  set maxTicks 601;601
  set year 0
  set timesHere 0
  set timesHere2 0
  set timesGo 0
  set timesAge 0
  set Climate random-float 1.0
   ifelse BreedProbVariability [set BreedProb (0.27 * (sin (50 * year)) + 0.3)] [set BreedProb BreedProbA];random-gamma BPGammaShape 19.06008193 [set BreedProb random-normal BreedProbA 0.114 while [(BreedProb < 0) or (BreedProb >= 1.01) ] [ set BreedProb random-normal BreedProbA 0.114 ]] [set BreedProb BreedProbA]
  ;set breedprob random-normal BreedProbA 0.114 while [(breedprob < 0) or (breedprob >= 1.01)] [set breedprob random-normal 0.251900417 0.114083441]

    
 ;; Initializing individuals for start of model run
  set-default-shape turtles "turtle"

 create-Hatchlings numberHatchlings [ setxy random-xcor random-ycor 
   set color magenta - 1
   set age 0 
   set times-nested 0
   set remigs 0
   set remigs_old remigs
   ifelse AgeMatVariability [
   set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] ; 30 (median of 40 and 22) - 17 + 1 = 14 Zug 2001 and Van Houtan et al. 2014
  ifelse ClutchFreqVariability [
   set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (MSORD Program Mark, with Sec Seasons >=4)
   ifelse ClutchSizeVariability [
   set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA] 
   if age = AgeMaturity [set breed Neophytes ]
  if age > AgeMaturity  [set breed Adults]
  ifelse (breed = Adults) or (breed = Neophytes) [set reprostatus "nester"]
  [set reprostatus "skip-nester"]
  set reprostatus2 reprostatus
     ]

  create-PelagicJuveniles numberPelagicJuveniles [ setxy random-xcor random-ycor 
   set color magenta - 1
   set age 2 + random 1 
   set times-nested 0
   set remigs 0
   set remigs_old remigs
   ifelse AgeMatVariability [
   set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] ; 30 (median of 40 and 22) - 17 + 1 = 14 Zug 2001 and Van Houtan et al. 2014
   ifelse ClutchFreqVariability [
   set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
   ifelse ClutchSizeVariability [
   set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA] 
  if age = AgeMaturity [set breed Neophytes ]
  if age > AgeMaturity  [set breed Adults]
  ifelse (breed = Adults) or (breed = Neophytes) [set reprostatus "nester"]
  [set reprostatus "skip-nester"]
  set reprostatus2 reprostatus
     ]

create-NeriticJuveniles numberNeriticJuveniles [ setxy random-xcor random-ycor 
   set color magenta - 1
   set age 4 + random 7 
   set times-nested 0
   set remigs 0
   set remigs_old remigs
   ifelse AgeMatVariability [
   set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] ; 30 (median of 40 and 22) - 17 + 1 = 14 Zug 2001 and Van Houtan et al. 2014
   ifelse ClutchFreqVariability [
   set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
   ifelse ClutchSizeVariability [
   set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA]   if age = AgeMaturity [set breed Neophytes ]
  if age > AgeMaturity  [set breed Adults]
  ifelse (breed = Adults) or (breed = Neophytes) [set reprostatus "nester"]
  [set reprostatus "skip-nester"]
  set reprostatus2 reprostatus
     ]

create-Subadults numberSubadults [ setxy random-xcor random-ycor 
   set color magenta - 1
   ;set AgeClass "Adult"
   set age 11 
   set times-nested 0
   set remigs 0
   set remigs_old remigs
   ifelse AgeMatVariability [
   set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] ; 30 (median of 40 and 22) - 17 + 1 = 14 Zug 2001 and Van Houtan et al. 2014
   ifelse ClutchFreqVariability [
   set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
   ifelse ClutchSizeVariability [
   set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA] ; mean 98 - 85 = 13 scaled to surivorship to neritic juvenile stage, avg survivorship = 8.5 see survival excel worksheet
  if age = AgeMaturity [set breed Neophytes ]
  if age > AgeMaturity  [set breed Adults]
  ifelse (breed = Adults) or (breed = Neophytes) [set reprostatus "nester"]
  [set reprostatus "skip-nester"]
  set reprostatus2 reprostatus
     ]

create-Neophytes numberNeophytes [ setxy random-xcor random-ycor 
   set color black
   ;set AgeClass "Adult"
   ifelse AgeMatVariability [
   set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] 
   set age AgeMaturity 
   set times-nested 0
   set remigs 0
   set remigs_old remigs
   ifelse ClutchFreqVariability [
   set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
   ifelse ClutchSizeVariability [
  set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA] ; mean 98 - 85 = 13 scaled to surivorship to neritic juvenile stage, avg survivorship = 8.5 see survival excel worksheet
  if age = AgeMaturity [set breed Neophytes ]
  if age > AgeMaturity  [set breed Adults]
  ifelse (breed = Adults) or (breed = Neophytes) [set reprostatus "nester"]
  [set reprostatus "skip-nester"]
  set reprostatus2 reprostatus
     ]
  
  
create-Adults numberAdults [ setxy random-xcor random-ycor 
   set color magenta - 1
   ifelse AgeMatVariability [
   set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] 
   set age 28 + random 57 
   set times-nested random 10
   set remigs random 4
   set remigs_old remigs
   ifelse ClutchFreqVariability [
   set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
   ifelse ClutchSizeVariability [
   set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA]    if age = AgeMaturity [set breed Neophytes ]
  if age > AgeMaturity  [set breed Adults]
  ifelse (breed = Adults) or (breed = Neophytes) [set reprostatus "nester"]
  [set reprostatus "skip-nester"]
  set reprostatus2 reprostatus
     ]

set countNesters2 count turtles with [ reprostatus = "nester"]
set SampleSize round(q * countNesters2)
set Monitored-Nesters no-turtles ;(n-of SampleSize turtles with [reprostatus = "nester"])
 
ask patches [ set pcolor cyan + 3 ]

reset-ticks
  
;Open an output file for testing and analysis
; First, delete it instead of appending it

;if (file-exists? "Individual_turtle_output.csv")
;[
;  carefully
;  [file-delete "Individual_turtle_output.csv"]
;  [print error-message]
;]

;file-open "Individual_turtle_output.csv"
;file-type "id,"
;file-type "Dead,"
;file-type "Tick,"
;file-type "Age,"
;file-type "Breed,"
;file-type "ReproStatus,"
;file-type "ReproStatus2,"
;file-type "Times_Nested,"
;file-type "New_Hatchlings,"
;file-type "AgeMaturity,"
;file-type "ClutchFrequency,"
;file-type "ClutchSize,"
;file-print "Remigs"
;file-close


; write summary data on turtle states to output file
;Open an output file for testing and analysis
; First, delete it instead of appending it
;if (file-exists? "GST_IBM_output.csv")
;[
;  carefully
;  [file-delete "GST_IBM_output.csv"]
;  [print error-message]
;]
;file-open "GST_IBM_output.csv"
;file-type "tick,"
;file-type "Climate,"
;file-type "SampleSize,"
;file-type "Monitored-Nesters,"
;file-type "Hatchlings_Nest,"
;file-type "Hatchlings_Nest_StD,"
;file-type "Mean_Times_Nested,"
;file-type "Times_Nested_StD,"
;file-type "Mean_Remigs,"
;file-type "Remigs_StD,"
;file-type "Mean_AgeMaturity,"
;file-type "AgeMaturity_StD,"
;file-type "Mean_ClutchFrequency,"
;file-type "ClutchFrequency_StD,"
;file-type "Mean_ClutchSize,"
;file-type "ClutchSize_StD,"
;file-type "Hatchlings_Nest_Nster,"
;file-type "Hatchlings_Nest_StD_Nster,"
;file-type "Mean_Times_Nested_Nster,"
;file-type "Times_Nested_StD_Nster,"
;file-type "Mean_Remigs_Nster,"
;file-type "Remigs_StD_Nster,"
;file-type "Mean_AgeMaturity_Nster,"
;file-type "AgeMaturity_StD_Nster,"
;file-type "Mean_ClutchFrequency_Nster,"
;file-type "ClutchFrequency_StD_Nster,"
;file-type "Mean_ClutchSize_Nster,"
;file-type "ClutchSize_StD_Nster,"
;file-type "Mean_Age_Hatchlings,"
;file-type "Age_Hatchlings_StD,"
;file-type "Mean_Age_Pel_Juvenile,"
;file-type "Age_Pel_Juvenile_StD,"
;file-type "Mean_Age_Ner_Juvenile,"
;file-type "Age_Ner_Juvenile_StD,"
;file-type "Mean_Age_Subadults,"
;file-type "Age_Subadults_StD,"
;file-type "Mean_Age_Neophytes,"
;file-type "Age_Neophytes_StD,"
;file-type "Mean_Age_Adults,"
;file-type "Age_Adults_StD,"
;file-type "Count_Hatchlings,"
;file-type "Count_Pel_Juvenile,"
;file-type "Count_Ner_Juvenile,"
;file-type "Count_Subadults,"
;file-type "Count_Neophytes,"
;file-type "Count_Adults,"
;file-type "Count_Nesters,"
;file-type "Count Nests," 
;file-type "Mean_Hatchlings_Nest_Mon-Nesters,"
;file-type "Hatchlings_Nest_Mon-Nesters_StD,"
;file-type "Mean_Times_Nested_Mon-Nesters,"
;file-type "Times_Nested_Mon-Nesters_StD,"
;file-type "Mean_Remigs_Mon-Nesters,"
;file-type "Remigs_Mon-Nesters_StD,"
;file-type "Mean_Age_Maturity_Mon-Nesters,"
;file-type "Age_Maturity_Mon-Nesters_StD,"
;file-type "Mean_Clutch_Freq_Mon-Nesters,"
;file-type "Clutch_Freq_Mon-Nesters_StD,"
;file-type "Mean_Clutch_Size_Mon-Nesters,"
;file-type "Clutch_Size_Mon-Nesters_StD,"
;file-type "Mean_Age_Mon-Nesters,"
;file-type "Age_Mon-Nesters_StD,"
;file-type "Count_Neophytes_Mon-Nesters,"
;file-type "Count_Adults_Mon-Nesters,"
;file-print "Count_Nests_Monitored" 
;file-close

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; start procedures ;;;;;;;;;;;;;;;;;;;;;;;;;

to go
tick
set year year + 1
if ticks >= maxTicks [stop]
if not any? turtles [ stop ]
 set q-norm random-normal DetectA q-std-dev ;while [(q <= 0) or (q > 1.0) ] [ set q random-normal DetectA q-std-dev ]; random-normal 0.79 0.17 while [(q <= 0.3) or (q >= 1.0) ] [ set q random-normal 0.75 0.25 ]   ; Catchability coefficient = p(t.) top ranked model MSORB Mark
  set q 1 / (1 + (exp (- q-norm)));set q random-normal DetectA q-std-dev while [(q <= 0) or (q > 1.0) ] [ set q random-normal DetectA q-std-dev ]; random-normal 0.79 0.17 while [(q <= 0.3) or (q >= 1.0) ] [ set q random-normal 0.75 0.25 ]   ; Catchability coefficient = p(t.) top ranked model MSORB Mark
; if ticks < 175 [set q random-normal DetectA 0.171 while [(q <= 0) or (q > 1.0) ] [ set q random-normal DetectA 0.25 ]] if (ticks >= 175)  [set q ((0.9428 * ln ticks) - 4.6572) while [(q <= 0) or (q > 1.0) ] [ set q ((0.9428 * ln ticks) - 4.6572)]] ;p increasing over time
;if ticks < 175 [set q random-normal DetectA 0.171 while [(q <= 0) or (q > 1.0) ] [ set q random-normal DetectA 0.25 ]] if (ticks >= 175)  [set q ((0.1309 * ln countNesters2) - 0.158) while [(q <= 0) or (q > 1.0) ] [ set q ((0.9428 * ln ticks) - 4.6572)] ]  ;p function of nester abundance
ifelse BreedProbVariability [set BreedProb (0.27 * (sin (50 * year)) + 0.3)] [set BreedProb BreedProbA];[set BreedProb random-gamma BPGammaShape 19.06008193] [set BreedProb BreedProbA] ;[set BreedProb random-normal BreedProbA 0.114 while [(BreedProb < 0) or (BreedProb >= 1.01) ] [ set BreedProb random-normal BreedProbA 0.114 ]] [set BreedProb BreedProbA]
;set breedprob random-normal BreedProbA 0.114 while [(breedprob <= 0) or (breedprob >= 1.01)] [set breedprob random-normal 0.251900417 0.114083441]
set Climate random-float 1.0
Climate_BP
Harvest
Harvest_Nests
ask turtles [
grow_old
set timesGo timesGo + 1
survive
     ]
ask Adults [Nesters?] 
set Nesters (turtles with [reprostatus = "nester"])
set countNesters count Nesters ; used for Hatchling Production equation
ask Nesters [Hatchling_Production]
;ask Neophytes [Hatchling_Production]
;ask Nesters [Remig_DD_Climate]
;ask Neophytes [Remig_DD_Climate]
;ask Adults [Remig_2]
ask Nesters [reproduce]
;ask Neophytes [reproduceNeos]
Monitor
;update-outputs
end

to Climate_BP
 ifelse ClimateBPVariability [ if Climate > ClimateThresh [
 set BreedProb BreedProb * BP-Climate-influence]] [] 
end

to grow_old
set timesAge timesAge + 1
set age age + 1
 if age > 85 [die] 
 if age <= 1 [set breed Hatchlings]
 if (age > 1) and (age < 4) [set breed PelagicJuveniles]
 if (age >= 4) and (age < 11) [set breed NeriticJuveniles]
 if (age >= 11) and (age < AgeMaturity) [set breed Subadults]
 if (age = AgeMaturity) [set breed Neophytes set reprostatus "nester" set reprostatus2 "nester" set color black] ;and (times-nested < 1)
 if (age > AgeMaturity) [set breed Adults set color magenta - 1] ; and (times-nested >= 1)
end

to Harvest
ifelse Harvest_on
    [
if ticks < 200 [] 
if (ticks >= 200) and (ticks < 251) [ set Catch (round ((count turtles with [age >= 11]) / F)) (ask n-of Catch turtles with [age >= 11] [die] ) 
   ] ;set Catch_hatchlings (round (count hatchlings) / F_hatch) (ask n-of Catch_hatchlings hatchlings [die] )
if ticks >= 250 [] 
] 
 [] 
end


to Harvest_Nests
ifelse Harvest_Nests_on
    [
if ticks < 200 [] 
if (ticks >= 200) and (ticks < 251) [ set Catch_hatchlings ((round (count hatchlings) / F_hatch)) (ask n-of Catch_hatchlings hatchlings [die]) 
   ] ;set Catch_hatchlings (round (count hatchlings) / F_hatch) (ask n-of Catch_hatchlings hatchlings [die] )
if ticks >= 250 [] 
] 
 [] 
end

to survive
;ifelse (ticks > 100) and (Climate < 0.9) [  
; set juvenile2Survival 0.824
; set subadultSurvival 0.876
; set adultSurvival 0.929] 
; [set juvenile2Survival 0.7416 ; 0.9*0.824
; set subadultSurvival 0.7884 ; 0.9*0.876
; set adultSurvival 0.8361 ]; .9*0.929 ;  ]]
;if AgeClass = "Hatchling" [if random-float 1.0 > hatchlingSurvival[ die]]  ;;hatchling survival 0.786 (sd=19.2, Niethammer et al. 1997) (0.71 (Van Buskirk and Crowder 1994 FFS HI)
;if AgeClass = "PelagicJuvenile" [if random-float 1.0 > juvenile1Survival [ die]] ;;oceanic juvenile survival 0.8804 (0.835 - 0.927 95 % CI Chaloupka and Limpus 2005 S GBR)
if breed = NeriticJuveniles [if random-float 1.0 > juvenile2Survival [ die]] ;; neritic juvenile survival 0.8804 (0.835 - 0.927 95 % CI Chaloupka and Limpus 2005 S GBR, 5-18 yrs)
if breed = Subadults [if random-float 1.0 > subadultSurvival [ die]] ;; subadult and adult survival 0.8474, ( 0.79-0.91 95% CI, Chaloupka and Limpus 2005 S GBR, 18-35 years)
if breed =  Neophytes [set timesHere timesHere + 1 if random-float 1.0 > adultSurvival [die]]
if breed = Adults [set timesHere timesHere + 1 if random-float 1.0 > adultSurvival [die]]
stop
end

to Nesters?
ifelse remigs = 0 [
   set reprostatus "skip-nester"  set remigs remigs + 1]  
[ 
ifelse random-float 1.0 < BreedProb [set reprostatus "nester"] [set reprostatus "skip-nester" set remigs remigs + 1]]
stop
end

to Hatchling_Production
set newHatchlings round(ClutchSize * ClutchFrequency) * exp (rh * ( 1 - (countNesters / K_nesters))) ;; Ricker Function N(t+1) = (CS*CF) *exp(rh[1-(N(t)/K)])  Alt form of Logistic: X[t+1] = X[t] + r * X[t] * (1 – X[t] / K) density dependence if number hatchlings produce per capita
   if newHatchlings < 0 [set newHatchlings 0] 
end  
 
to reproduce 
; ifelse reprostatus = "skip-nester" 
;   [ set remigs remigs + 1 stop ]
   ;[
     hatch-Hatchlings newHatchlings [ ; inherited by new
      set age 0
      set color blue - 1
      setxy random-xcor random-ycor 
      set reprostatus "skip-nester"
      set reprostatus2 "skip-nester"
      set times-nested 0
      set remigs 0
      set remigs_old remigs
      ifelse AgeMatVariability [
      set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA] ; 30 (median of 40 and 22) - 17 + 1 = 14 Zug 2001 and Van Houtan et al. 2014
      ifelse ClutchFreqVariability [
      set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
      ifelse ClutchSizeVariability [
     set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA] 
           ]
    set times-nested times-nested + 1
    set remigs_old remigs
    set remigs 0
    stop ;] 
end

to reproduceNeos 
 
;ifelse Climate > 0.9 [ ]
;[
;set reprostatus "nester"
      hatch-Hatchlings newHatchlings [
      set age 0
      set color green - 1
      setxy random-xcor random-ycor 
      set reprostatus "skip-nester"
      set reprostatus2 "skip-nester"
      set times-nested 0
      set remigs 0
      set remigs_old remigs
      ifelse AgeMatVariability [
      set AgeMaturity 19 + random-poisson AgeMatRand] [set AgeMaturity AgeMatA]
      ifelse ClutchFreqVariability [
      set ClutchFrequency random-poisson ClutchFreqA] [set ClutchFrequency ClutchFreqA] ; ; mean 4-1 = 3; mean = 3.37304065, SD=1.254475826 (mark, with Sec Seasons >=4)
      ifelse ClutchSizeVariability [
     set ClutchSize random-poisson ClutchSizeA] [set ClutchSize ClutchSizeA] 
                 ]
    set times-nested times-nested + 1
    set remigs_old remigs
    set remigs 0
    stop 
end


to Monitor  
  set Monitored-Nesters no-turtles
  set countNesters2 count Nesters ;turtles with [ reprostatus = "nester"]
  set SampleSize round(q * countNesters2)
  ;set Monitored-Nesters (max-n-of SampleSize Nesters [ClutchFrequency]) ; Monitored-Nesters sampled based on turtles with greatest clutch frequency (better chance of being observed on nesting beach, turtles with [reprostatus = "nester"]
  ;set Monitored-Nesters (max-n-of SampleSize Nesters [ClutchFrequency]) ;
  set Monitored-Nesters (max-n-of SampleSize Nesters [Age]) ; Monitored-Nesters sampled based on oldest turtles  (older turtles have great site fidelity/larger (easier to spot) better chance of being observed on nesting beach
  set Nests sum [ClutchFrequency] of Nesters
  set Monitored-Nests sum [ClutchFrequency] of Monitored-Nesters;round(q * Nests)
  set Monitored-Nests-rand round(q * Nests)
end 

to update-outputs
  if ticks > 100 [
 file-open "Individual_turtle_output.csv"
 ask Monitored-Nesters ;turtles with [(breed = Adults) or (breed = Neophytes) ] ;or (breed = Adults)
 [
  file-type (word who ",")
  file-type (word (is-turtle? nobody) ",")
  file-type (word ticks ",")
  file-type (word age ",")
  file-type (word breed ",")
  file-type (word reprostatus ",")
  file-type (word reprostatus2 ",")
  file-type (word times-nested ",")
  file-type (word newHatchlings ",")
  file-type (word AgeMaturity ",")
  file-type (word ClutchFrequency ",")
  file-type (word ClutchSize ",")
  file-print remigs
 ]
file-close
  ]
  
;if ticks > 150 [
;file-open "GST_IBM_output.csv"
;file-type (word ticks ",")
;file-type (word Climate ",")
;file-type (word SampleSize ",")
;file-type (word count Monitored-Nesters ",")
;file-type (word mean [newHatchlings] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")  
;file-type (word standard-deviation [newHatchlings] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word mean [times-nested] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word standard-deviation [times-nested] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word mean [remigs] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word standard-deviation [remigs] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word mean [AgeMaturity] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word standard-deviation [AgeMaturity] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word mean [ClutchFrequency] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word standard-deviation [ClutchFrequency] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word mean [ClutchSize] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word standard-deviation [ClutchSize] of turtles with [(breed = Adults) or (breed = Neophytes)] ",")
;file-type (word mean [newHatchlings] of turtles with [reprostatus = "nester"] ",")
;file-type (word standard-deviation [newHatchlings] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [times-nested] of turtles with [reprostatus = "nester"] ",")
;file-type (word standard-deviation [times-nested] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [remigs] of turtles with [reprostatus = "nester"] ",")
;file-type (word standard-deviation [remigs] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [AgeMaturity] of turtles with [reprostatus = "nester"] ",")
;file-type (word standard-deviation [AgeMaturity] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [ClutchFrequency] of turtles with [reprostatus = "nester"] ",")
;file-type (word standard-deviation [ClutchFrequency] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [ClutchSize] of turtles with [reprostatus = "nester"] ",")
;file-type (word standard-deviation [ClutchSize] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [age] of Hatchlings ",")
;file-type (word standard-deviation [age] of Hatchlings ",")
;file-type (word mean [age] of PelagicJuveniles ",")
;file-type (word standard-deviation [age] of PelagicJuveniles ",")
;file-type (word mean [age] of NeriticJuveniles ",")
;file-type (word standard-deviation [age] of NeriticJuveniles ",")
;file-type (word mean [age] of Subadults ",")
;file-type (word standard-deviation [age] of Subadults ",")
;file-type (word mean [age] of Neophytes ",")
;file-type (word standard-deviation [age] of Neophytes ",")
;file-type (word mean [age] of Adults ",")
;file-type (word standard-deviation [age] of Adults ",")
;file-type (word count Hatchlings ",")
;file-type (word count PelagicJuveniles ",")
;file-type (word count NeriticJuveniles ",")
;file-type (word count Subadults ",")
;file-type (word count Neophytes ",")
;file-type (word count Adults ",")
;file-type (word count turtles with [reprostatus = "nester"] ",")
;file-type (word sum [ClutchFrequency] of turtles with [reprostatus = "nester"] ",")
;file-type (word mean [newHatchlings] of Monitored-Nesters ",")
;file-type (word standard-deviation [newHatchlings] of Monitored-Nesters ",")
;file-type (word mean [times-nested] of Monitored-Nesters  ",")
;file-type (word standard-deviation [times-nested] of Monitored-Nesters  ",")
;file-type (word mean [remigs] of Monitored-Nesters  ",")
;file-type (word standard-deviation [remigs] of Monitored-Nesters  ",")
;file-type (word mean [AgeMaturity] of Monitored-Nesters ",")
;file-type (word standard-deviation [AgeMaturity] of Monitored-Nesters ",")
;file-type (word mean [ClutchFrequency] of Monitored-Nesters  ",")
;file-type (word standard-deviation [ClutchFrequency] of Monitored-Nesters  ",")
;file-type (word mean [ClutchSize] of Monitored-Nesters  ",")
;file-type (word standard-deviation [ClutchSize] of Monitored-Nesters  ",")
;file-type (word mean [age] of Monitored-Nesters ",")
;file-type (word standard-deviation [age] of Monitored-Nesters ",")
;file-type (word count Monitored-Nesters with [breed = "Neophyte"] ",")
;file-type (word count Monitored-Nesters with [breed = "Adult"] ",")
;file-print sum [ClutchFrequency] of Monitored-Nesters
;file-close
;]
end
@#$#@#$#@
GRAPHICS-WINDOW
1128
253
1373
504
10
10
10.5
1
10
1
1
1
0
0
0
1
-10
10
-10
10
1
1
1
ticks
30.0

BUTTON
8
10
63
43
NIL
setup
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
70
10
133
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

MONITOR
662
585
734
630
Count Adults
count Adults
0
1
11

SLIDER
3
110
149
143
hatchlingSurvival
hatchlingSurvival
0
1
1
0.001
1
3
HORIZONTAL

SLIDER
3
144
149
177
juvenile1Survival
juvenile1Survival
0
1
1
0.001
1
3
HORIZONTAL

SLIDER
3
213
148
246
subadultSurvival
subadultSurvival
0
1
0.876
0.001
1
3
HORIZONTAL

SLIDER
2
280
148
313
adultSurvival
adultSurvival
0
1
0.929
0.001
1
3
HORIZONTAL

SLIDER
2
314
148
347
K_nesters
K_nesters
0
10000
500
1
1
NIL
HORIZONTAL

MONITOR
514
584
584
629
Subadults
count Subadults
0
1
11

MONITOR
523
536
616
581
Count PelJuve
count PelagicJuveniles
1
1
11

MONITOR
734
588
815
633
Hatchlings/Nest
mean [newhatchlings] of turtles with [reprostatus = \"nester\"]
1
1
11

SLIDER
3
179
148
212
juvenile2Survival
juvenile2Survival
0
1
0.824
0.001
1
3
HORIZONTAL

MONITOR
598
537
694
582
Count NerJuve
count NeriticJuveniles
17
1
11

MONITOR
693
446
750
491
Nesters
count turtles with [ reprostatus = \"nester\"]
1
1
11

PLOT
536
153
921
273
Juveniles
ticks
N(t)
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Ner_Juve" 1.0 0 -2674135 true "" "plot count NeriticJuveniles"
"Pel_Juve" 1.0 0 -12087248 true "" "plot count PelagicJuveniles"
"Subadults" 1.0 0 -13345367 true "" "plot count Subadults"

BUTTON
139
10
202
43
NIL
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

PLOT
537
274
922
394
hatchlings
ticks
hatchlings
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count Hatchlings"

MONITOR
752
446
820
491
Skip-nesters
count Adults with [ reprostatus = \"skip-nester\"] + count Neophytes with [ reprostatus = \"skip-nester\"]
0
1
11

PLOT
536
10
919
153
Adults
ticks
N(t)
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Adults" 1.0 0 -10899396 true "" "plot count Adults"
"Neophytes" 1.0 0 -10141563 true "" "plot count Neophytes"
"Total" 1.0 0 -13345367 true "" "plot count turtles with [(breed = Adults) or (breed = Neophytes)]"

PLOT
932
11
1314
131
Nester Abundance
ticks
Nester
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Nesters" 1.0 0 -13345367 true "" "plot count Adults with [reprostatus = \"nester\"]"
"Monitored Nesters" 1.0 0 -2674135 true "" "plot count Monitored-Nesters"

PLOT
151
111
515
246
Age Histogram
Age
Frequency
3.0
85.0
0.0
250.0
true
true
"set-histogram-num-bars 20" "set-histogram-num-bars 20"
PENS
"Adults" 5.0 1 -10899396 true "" "histogram [age] of Adults"
"Neophytes" 5.0 1 -10141563 true "" "histogram [age] of Neophytes"
"Subadult" 5.0 1 -13345367 true "" "histogram [age] of Subadults"
"Ner_Juve" 5.0 1 -2674135 true "" "histogram [age] of NeriticJuveniles"

MONITOR
607
493
694
538
Nester/Skip-nester
count turtles with [ reprostatus = \"nester\"] / ( count Adults with [ reprostatus = \"skip-nester\"] + count Neophytes with [reprostatus = \"skip-nester\"])
5
1
11

MONITOR
611
445
692
490
Remigs years
mean [remigs] of Adults with [reprostatus = \"nester\"]
2
1
11

MONITOR
717
398
831
443
Mean Times Nested/female
mean [times-nested] of Adults
4
1
11

SLIDER
2
247
148
280
neophyteSurvival
neophyteSurvival
0
1
0.929
0.001
1
NIL
HORIZONTAL

MONITOR
583
584
661
629
Count Neos
count Neophytes
0
1
11

SWITCH
147
44
293
77
RemigVariability
RemigVariability
1
1
-1000

SWITCH
293
44
439
77
ClutchFreqVariability
ClutchFreqVariability
0
1
-1000

MONITOR
533
487
606
532
ClutchFreq
mean [ClutchFrequency] of turtles with [reprostatus = \"nester\"]
1
1
11

SWITCH
3
77
148
110
ClutchSizeVariability
ClutchSizeVariability
0
1
-1000

MONITOR
697
538
800
583
Nesters (Ad & Neo)
count turtles with [ reprostatus = \"nester\" ]
0
1
11

MONITOR
529
444
605
489
Neos Times Nested
mean [times-nested] of Neophytes
4
1
11

SWITCH
204
10
350
43
Harvest_on
Harvest_on
0
1
-1000

MONITOR
649
398
716
443
Adult Age
mean [age] of Adults
0
1
11

MONITOR
557
398
650
443
Neophyte Age
Mean [age] of Neophytes
0
1
11

MONITOR
695
494
768
539
Clutch Size
Mean [ClutchSize] of turtles with [reprostatus = \"nester\"]
3
1
11

MONITOR
454
536
521
581
Count Hatchlings
count Hatchlings
0
1
11

MONITOR
499
398
556
443
Climate
Climate
1
1
11

MONITOR
152
249
226
294
SampleSize
SampleSize
0
1
11

MONITOR
228
248
312
293
q (Catchability)
q
4
1
11

MONITOR
229
294
327
339
Monitored Nesters
count Monitored-Nesters
0
1
11

MONITOR
451
489
532
534
AgeMaturity
mean [AgeMaturity] of Neophytes
4
1
11

MONITOR
456
584
513
629
Turtles
count turtles
0
1
11

MONITOR
329
295
386
340
# Nests
Nests
1
1
11

SWITCH
148
77
294
110
AgeMatVariability
AgeMatVariability
0
1
-1000

MONITOR
449
443
526
488
Neo:Nester
count Neophytes / count turtles with [reprostatus = \"nester\"]
3
1
11

MONITOR
374
249
459
294
NIL
adultSurvival
3
1
11

SLIDER
152
340
296
373
AgeMatA
AgeMatA
0
50
30
1
1
NIL
HORIZONTAL

SLIDER
151
373
296
406
ClutchFreqA
ClutchFreqA
0
25
4
1
1
NIL
HORIZONTAL

SLIDER
152
406
296
439
ClutchSizeA
ClutchSizeA
0
100
8
1
1
NIL
HORIZONTAL

TEXTBOX
1043
262
1138
353
Survival Rates:\nAdult       0.929\nNeophyte    0.929\nSubadult    0.876\nNer Juve    0.824\nPel Juve    0.80\nHatchling   0.35
10
0.0
1

SLIDER
296
373
442
406
AgeMatRand
AgeMatRand
0
50
11
1
1
NIL
HORIZONTAL

SLIDER
2
348
148
381
DetectA
DetectA
-5.0
5.0
0.75
0.05
1
NIL
HORIZONTAL

MONITOR
315
249
372
294
Catch
Catch
0
1
11

MONITOR
388
295
475
340
Monitored Nests
Monitored-Nests
0
1
11

SLIDER
3
380
148
413
F
F
0
100
15
1
1
NIL
HORIZONTAL

MONITOR
461
249
535
294
Breed Prob
BreedProb
3
1
11

MONITOR
153
294
226
339
Nesters
count Nesters
0
1
11

PLOT
933
131
1317
251
Nests
ticks
Nests
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Nests Total" 1.0 0 -13840069 true "" "plot Nests"
"Mon. Nests" 1.0 0 -8630108 true "" "plot Monitored-Nests"

SLIDER
297
340
442
373
BreedProbA
BreedProbA
0
1
0.252
1
1
NIL
HORIZONTAL

SWITCH
294
77
440
110
BreedProbVariability
BreedProbVariability
0
1
-1000

MONITOR
475
295
535
340
Real. BP
count Nesters / (count Neophytes + count Adults)
3
1
11

SWITCH
2
44
147
77
ClimateBPVariability
ClimateBPVariability
0
1
-1000

PLOT
3
481
291
631
Breed Prob
Tick
Breed Prob
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 0.1 0 -16777216 true "" "plot BreedProb"

SLIDER
296
406
443
439
BPGammaShape
BPGammaShape
0
100
4.80124
0.05
1
NIL
HORIZONTAL

SLIDER
4
447
148
480
rhA
rhA
0
5
0.15
0.1
1
NIL
HORIZONTAL

SLIDER
296
439
443
472
ClimateThresh
ClimateThresh
0
1
0.9
0.05
1
NIL
HORIZONTAL

SLIDER
151
438
296
471
BP-Climate-influence
BP-Climate-influence
0
1
0.75
0.05
1
NIL
HORIZONTAL

TEXTBOX
937
261
1079
443
Base Model Settings:\nK_nesters 500 (600)\nF 15 (20)\nDetectA 0.75\nAgeMatA 30\nClutchFreqA 4\nClutchSizeA 8\nBPGammaShape 4.80124\nBreedProbA 0.252\nAgeMatRand 11\nBP-Climate_influence 0.75\nrhA 0.15\nClimateThresh 0.90\n
11
0.0
1

TEXTBOX
940
450
1090
484
CHECK BASE MODEL INPUTS!!
14
15.0
1

SLIDER
4
414
148
447
F_hatch
F_hatch
0
100
1.5
1
1
NIL
HORIZONTAL

SWITCH
350
10
496
43
Harvest_Nests_on
Harvest_Nests_on
1
1
-1000

MONITOR
299
479
379
524
Catch hatch
Catch_hatchlings
0
1
11

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

egg
false
0
Circle -7500403 true true 96 76 108
Circle -7500403 true true 72 104 156
Polygon -7500403 true true 221 149 195 101 106 99 80 148

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

turtle 2
true
0
Polygon -10899396 true false 132 85 134 64 130 57 132 32 151 22 171 33 173 57 169 65 172 87
Polygon -10899396 true false 165 210 195 210 240 240 240 255 210 285 195 255
Polygon -10899396 true false 90 210 60 240 60 255 90 285 105 255 120 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 209 225 155 275 143 275 89 225 74 135 88 99
Polygon -16777216 true false 151 21 142 26 159 26
Line -16777216 false 160 35 164 51
Polygon -16777216 true false 161 38 162 46 167 45
Line -16777216 false 169 63 156 69
Line -16777216 false 134 64 144 69
Line -16777216 false 143 69 156 69
Polygon -16777216 true false 139 38 138 46 133 45
Line -16777216 false 140 35 136 51
Polygon -10899396 true false 195 90 225 75 255 75 270 90 285 120 285 150 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 45 75 30 90 15 120 15 150 60 105 75 105 90 105

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
NetLogo 5.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="GST SA_Detect" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="351"/>
    <exitCondition>turtles = 0</exitCondition>
    <metric>ticks</metric>
    <metric>Climate</metric>
    <metric>BreedProb</metric>
    <metric>count Nesters / (count Neophytes + count Adults)</metric>
    <metric>Catch</metric>
    <metric>Catch_hatchlings</metric>
    <metric>q</metric>
    <metric>count turtles</metric>
    <metric>count Adults</metric>
    <metric>count Adults + count Neophytes</metric>
    <metric>count Neophytes</metric>
    <metric>count Subadults</metric>
    <metric>count NeriticJuveniles</metric>
    <metric>count PelagicJuveniles</metric>
    <metric>count Hatchlings</metric>
    <metric>count Nesters</metric>
    <metric>count Nesters / (0.25)</metric>
    <metric>count Nesters with [breed = Neophytes]</metric>
    <metric>count Nesters with [breed = Adults ]</metric>
    <metric>Nests</metric>
    <metric>Nests / (0.25 * 4)</metric>
    <metric>count Monitored-Nesters</metric>
    <metric>(count Monitored-Nesters) / (0.25 * DetectA)</metric>
    <metric>Monitored-Nests</metric>
    <metric>Monitored-Nests / (0.25 * 4 * DetectA)</metric>
    <metric>count Monitored-Nesters with [breed = Neophytes]</metric>
    <metric>count Monitored-Nesters with [breed = Adults ]</metric>
    <metric>count Neophytes / (count Adults + count Neophytes)</metric>
    <metric>count Neophytes / count Nesters</metric>
    <metric>(count Neophytes / count Nesters) * 0.25</metric>
    <metric>count Monitored-Nesters with [breed = Neophytes] / count Monitored-Nesters</metric>
    <metric>(count Monitored-Nesters with [breed = Neophytes] / count Monitored-Nesters) * 0.25</metric>
    <metric>Mean [ClutchSize] of Adults</metric>
    <metric>standard-deviation [ClutchSize] of Adults</metric>
    <metric>mean [newHatchlings] of Adults</metric>
    <metric>standard-deviation [newhatchlings] of Adults</metric>
    <metric>mean [times-nested] of Adults</metric>
    <metric>standard-deviation [times-nested] of Adults</metric>
    <metric>max [times-nested] of Adults</metric>
    <metric>mean [remigs] of Adults</metric>
    <metric>standard-deviation [remigs] of Adults</metric>
    <metric>mean [remigs_old] of Adults</metric>
    <metric>standard-deviation [remigs_old] of Adults</metric>
    <metric>mean [ClutchFrequency] of Adults</metric>
    <metric>standard-deviation [ClutchFrequency] of Adults</metric>
    <metric>mean [AgeMaturity] of Adults</metric>
    <metric>standard-deviation [AgeMaturity] of Adults</metric>
    <metric>mean [age] of Adults</metric>
    <metric>standard-deviation [age] of Adults</metric>
    <metric>Mean [ClutchSize] of Nesters</metric>
    <metric>standard-deviation [ClutchSize] of Nesters</metric>
    <metric>mean [newHatchlings] of Nesters</metric>
    <metric>standard-deviation [newhatchlings] of Nesters</metric>
    <metric>mean [times-nested] of Nesters</metric>
    <metric>standard-deviation [times-nested] of Nesters</metric>
    <metric>max [times-nested] of Nesters</metric>
    <metric>mean [remigs] of Nesters</metric>
    <metric>standard-deviation [remigs] of Nesters</metric>
    <metric>mean [remigs_old] of Nesters</metric>
    <metric>standard-deviation [remigs_old] of Nesters</metric>
    <metric>mean [ClutchFrequency] of Nesters</metric>
    <metric>standard-deviation [ClutchFrequency] of Nesters</metric>
    <metric>mean [AgeMaturity] of Nesters</metric>
    <metric>standard-deviation [AgeMaturity] of Nesters</metric>
    <metric>mean [age] of Nesters</metric>
    <metric>standard-deviation [age] of Nesters</metric>
    <metric>count Nesters with [breed = Neophytes] / count Nesters</metric>
    <metric>mean [newHatchlings] of Monitored-Nesters</metric>
    <metric>standard-deviation [newHatchlings] of Monitored-Nesters</metric>
    <metric>mean [times-nested] of Monitored-Nesters</metric>
    <metric>standard-deviation [times-nested] of Monitored-Nesters</metric>
    <metric>max [times-nested] of Monitored-Nesters</metric>
    <metric>mean [remigs] of Monitored-Nesters</metric>
    <metric>standard-deviation [remigs] of Monitored-Nesters</metric>
    <metric>mean [remigs_old] of Monitored-Nesters</metric>
    <metric>standard-deviation [remigs_old] of Monitored-Nesters</metric>
    <metric>mean [AgeMaturity] of Monitored-Nesters</metric>
    <metric>standard-deviation [AgeMaturity] of Monitored-Nesters</metric>
    <metric>mean [ClutchFrequency] of Monitored-Nesters</metric>
    <metric>standard-deviation [ClutchFrequency] of Monitored-Nesters</metric>
    <metric>mean [ClutchSize] of Monitored-Nesters</metric>
    <metric>standard-deviation [ClutchSize] of Monitored-Nesters</metric>
    <metric>mean [age] of Monitored-Nesters</metric>
    <metric>standard-deviation [age] of Monitored-Nesters</metric>
    <enumeratedValueSet variable="DetectA">
      <value value="-2.216764"/>
      <value value="-1.643709E-9"/>
      <value value="3.841832"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
