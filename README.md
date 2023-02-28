# ArielNikas-VEWaning
This code is for reproducing the work found in "Estimating Waning of Vaccine Effectiveness: A Simulation Study" published in Clinical Infectious Diseases. (Clin Infect Dis. 2023 Feb 8;76(3):479-486. doi: 10.1093/cid/ciac725.)

For VE waning project focused on comparing methods using an ABM. This constitutes 3 main parts:

1)Julia code for running the simulation, the code gives a single example simulation for 100-0% waning and 30 day spread vaccination. Please note that this example simulation will creat two files begining with "Example-30Day" in your current working directory. To find your current working directory use the pwd() command in julia.  Included is the basic simulation used through most of the paper. Julia files end in .jl while files that describe packages, versions, and the general environment are found in the .toml files. 

2)R code for analyzing the simulation. Included is a .csv file for an example simulation. Example-Analysis.Rmd will (in order) give the true vaccine efficacy for this simulation, run you through how to do Levels 1-3(SR), run you through how to do Levels 1-3(SR), compare all Level 3 methods, find a heat map for min days and events, and compare non-optimized and optimized partitions. This code is setup to work with Batch1-30Day23.csv (gives information on an individual basis) and Batch1-30Day123.csv (gives day by day exposure information) but to reduce computational time generally uses larger bins than found in the paper. This may alter some qualitative aspects but not any of the comparisons or quantitative interpretations. Please note that if minimum bins are set very small the code may take a while to run.

3)Files for figures from both the manuscript and supplement.
