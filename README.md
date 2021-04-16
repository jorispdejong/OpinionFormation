# Opinion Formation
This repo contains an R project that models the dynamics of a sociological reaction network of opinion formation using both Krylov Subspace Approximation (KSA) and Monte Carlo simulations with the Gillespie algorithm.
The simulations were part of an Mathematics Master's assignment, so for supplementary information about e.g. the KSA method, the Gillespie algorithm, or the sociological model, see the report.

### KSA
The file KSA.R runs the KSA method and will output the files KSA_liberal.png, KSA_totalitarian.png, KSA_liberal.gif, KSA_totalitarian.gif to the 'results' folder if they don't exist yet. Also it will create the RData objects KSA_probabilities_liberal.RData and KSA_probabilities_totalitarian.RData if they don't exist yet, which contains the simulation output. The R script will check if the RData objects exists and load them, to avoid having to run the full simulation.

### Monte Carlo Gillespie
The file MCGillespie.R runs the Monte Carlo Gillspie method and will output the files MCGillespie_liberal.png, MCGillespie_totalitarian.png to the 'results' folder if they don't exist yet. Also it will create the RData objects MCG_Yt_liberal.RData and MCG_Yt_totalitarian.RData if they don't exist yet, which contains the simulation output. The R script will check if the RData objects exists and load them, to avoid having to run the full simulation.

### Functions.R
The Functions.R file contains two helper functions. The first one computes the propensity vector and the second function computes w=exp(tA)v, where A is a large sparse matrix. The second function is rewritten from the matlab package Expokit.