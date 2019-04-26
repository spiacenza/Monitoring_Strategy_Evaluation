# Monitoring_Strategy_Evaluation
Code for Monitoring Strategy Evaluation and associated experimental treatments (Piacenza et al. (2019):
Fathoming sea turtles: Monitoring Strategy Evaluation to improve conservation status assessments

Susan E. Piacenza1,3,4, Paul M. Richards2, and Selina S. Heppell1
1 Department of Fisheries and Wildlife, Oregon State University, Corvallis, OR, 97330, USA
2 NOAA NMFS – Southeast Fisheries Science Center, Miami, FL, 33149, USA
3 Present Address: Department of Biology, University of West Florida, Pensacola, FL 32514
4 Email: Susiepiacenza@gmail.com
Ecological Applications In press)

Each file corresponds to a different impact treatment (BP Cycle, Low and High Neritic Juvenile Impacts) and 
sampling treatment(random, age-biased, clutch frequency-biased). Code is for NetLogo 5.1.0. 
Newer versions of NetLogo may not work as well.
The GUI indicates the based model settings. Always confirm base model settings are correct in input boxes. 
Experiments were run by using BehaviorSpace using the GST SA_Detect for the 
BP Cycle treatment (with ["DetectA" -2.216764 -1.643709E-9 3.841832]; correspond to 10%,50% and 90%, respectively,
when converted from the log normalscale) and GST SA_Detect_FLowHigh (with ["DetectA" -2.216764 -1.643709E-9 3.841832]
["F" 2 10], which corresponds to 50% ad 10% harvest rate of juveniles/year). When running BehaviorSpace it works best to use the Table Output option (which continuously writes output as it is produced, and thus if model run doesn't complete there is still an output file to view). If running on a computer with multiple processer cores, especially dual core computers, the number of parallels runs should be (Nprocessor cores/2) - 1; e.g. for a computer with 8 dual processor cores, 3 parallel runs is optimal (in my experience and after trial and error).
