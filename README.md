# ArielNikas-VEWaning
For VE waning project focused on comparing methods using an ABM. This constitutes 3 main parts:

1)Julia code for running the simulation. Included is the basic simulation used through most of the paper. Julia files end in .jl while files that describe packages, versions, and the general environment are found in the .toml files. 

2)R code for analyzing the simulation. Included is a .csv file for an example simulation. Example-Analysis.Rmd will (in order) give the true vaccine efficacy for this simulation, run you through how to do Levels 1-3(SR), then the next chunk gives an example run through of the different Level 3 methods. The last chunk runs through how to find the optimal bins based off of either weighted error or AIC. This code is setup to work with Batch1-30Day23.csv (gives information on an individual basis) and Batch1-30Day123.csv (gives day by day exposure information) and assumes that they will be placed on your desktop. Please note that if minimum bins are set very small the code may take a while to run.

3)Files for the paper, supplement, and figures.
