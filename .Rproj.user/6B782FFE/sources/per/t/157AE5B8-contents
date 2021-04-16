###########################################################################
###########      Solving the Master Equation Numerically       ############
###########           EXample: Opinion formation               ############
###########        Monte Carlo Simulations using Gillespie     ############
###########          Authors: Joris de Jong                    ############
###########           Last Modified: 15/04/2021                ############
###########################################################################
# clear global environment
rm(list = ls(all.names = TRUE))

# libraries
suppressPackageStartupMessages(library('pracma')) # interp1()

# source
source('Functions.R')

# reproducibility
set.seed(6)

# show messages
verbose <- T

########################
### MODEL PARAMETERS ###
########################
L <- 40 # number of individuals is 2L
kappa1 <- 20/L # reaction parameter 1
kappa2 <- 40/L # reaction parameter 2

# liberal system: ELASPED TIME = 38.44 sec
# totalitarian system: ELASPED TIME = 34.56 sec
liberal <- F

if(liberal){
  a1 <- 0/L
  a2 <- 0.5/L
  a3 <- 0.5/L
  
  file_name_Yt_obj <- 'MCG_Yt_liberal.RData'
}else{
  a1 <- 1.5/L
  a2 <- 1/L
  a3 <- -0.25/L
  
  file_name_Yt_obj <- 'MCG_Yt_totalitarian.RData'
}

total_time <- 15 # simulation time in days
dt <- 0.1 # time step
time <- seq(0,total_time,dt) # vector with accumulating time steps

##################
### SIMULATION ###
##################
# check if the R object probs already exists, then we simply load it without running the simulation
if(file.exists(file_name_Yt_obj)){
  load(file_name_Yt_obj)
}else{
  if(verbose) print('Running Monte Carlo Simulations using Gillespie algorithm... Please wait')
  
  # measure elapsed time
  system.time(
    {
      #################
      ### GILLESPIE ###
      #################
      # input parameters
      S <- 1000 # number of simulations
      M <- 4 # number of reactions
      N <- 2 # number of species
      
      Y0 <- c(0,0) # starting state vector for each simulations
      
      V <- matrix(c(1,0,-1,0, 0, 1, 0, -1), N,M) # net effect matriX of the system
      
      # simulate a trajectory of public/private opinions over time using Gillespie algorithm
      Yt_data <- array(NA, dim=c(S,N,length(time))) # empty array for collecting the data of each simulation
      pb <- txtProgressBar(1,S) # progress bar
      for(s in 1:S){
        # initialize variables
        t <- 0
        t_collector <- t
        Yt <- Y0
        Yt_collector <- matrix(Yt, nrow=2)
        while(t < total_time){
          # determine which reactions can still occur
          possible_reactions <- 1:M
          if(Yt[1] == L) possible_reactions <- possible_reactions[-c(1)]
          if(Yt[1] == -L) possible_reactions <- possible_reactions[-c(2)]
          if(Yt[2] == L) possible_reactions <- possible_reactions[-c(3)]
          if(Yt[2] == -L) possible_reactions <- possible_reactions[-c(4)]
          
          # compute propensity vector
          alpha <- propensityVector(Yt)[possible_reactions]
          
          # simulate waiting time
          t_add <- rexp(1, sum(alpha)) 
          t <- t + t_add
          t_collector <- c(t_collector, t)
          
          # sample which reaction will be triggered
          m <- sample(possible_reactions, 1, prob = alpha/sum(alpha))
          
          # update the state vector
          Yt <- Yt + V[,m]
          Yt_collector <- cbind(Yt_collector, Yt)
        }
        
        # get Y at every time step by interpolating, then collect the output
        Yt_data[s,1,] <- interp1(t_collector, Yt_collector[1,], xi = time, method = 'nearest')
        Yt_data[s,2,] <- interp1(t_collector, Yt_collector[2,], xi = time, method = 'nearest')
        
        # update the progress bar
        setTxtProgressBar(pb,s)
      }
      
      # save Yt data as R object 
      save(Yt_data, file = file_name_Yt_obj)
    }
  )
}


####################
### PLOT RESULTS ###
####################
# counts the number of occurrences of each state vector after 15 daXs
vals <- Yt_data[,,length(time)]
bins <- -L:L
counts <- sapply(1:ncol(vals), function(j) sapply(bins, function(i) sum(vals[,j]==i)))

# compute the eXpected mean and standard deviation over all simulations
expected_mean <- apply(Yt_data, c(2,3), mean)
expected_std <- apply(Yt_data, c(2,3), std)

if(liberal){
  png_file_name <- 'Gillespie_liberal.png'
}else{
  png_file_name <- 'Gillespie_totalitarian.png'
}

file_needs_closing <- F
if(!file.exists(png_file_name)){
  # open file destination for png
  png(filename = png_file_name, width = 1000, height = 640) 
  
  file_needs_closing <- T
}

# plot the results
par(mfrow=c(2,2), mar=c(5.1, 4.1, 2.1, 2.1), oma=c(0,1,4,0)) # set plotting window

# histograms showing the net public/private opinion on the final time step
plot(bins, counts[,1]/sum(counts[,1]), type = 'h', lwd = 3,
     xlab = 'net public opinion', ylab = 'probability', cex.lab=1.2)
plot(bins, counts[,2]/sum(counts[,2]), type = 'h', lwd = 3,
     xlab = 'net private opinion', ylab = 'probability', cex.lab=1.2)
title('Probability distribution at t=15 days', line=-1, cex.main=1.17, outer = T)

# plots showing the eXpected net public/private opinion over time
plot(time, expected_mean[1,]+expected_std[1,], type = 'l',
     xlab = 'time (days)', ylab = 'Net public opinion', cex.lab=1.2,
     lwd=2, ylim=c(-L,L))
lines(time, expected_mean[1,]-expected_std[1,], lwd=2, col=2)
plot(time, expected_mean[2,]+expected_std[2,], type = 'l',
     xlab = 'time (days)', ylab = 'Net private opinion', cex.lab=1.2,
     lwd=2, ylim=c(-L,L))
lines(time, expected_mean[2,]-expected_std[2,], lwd=2, col=2)
title('Expected net public/private opinion over time', line=-25.5, cex.main=1.17, outer = T)

# main title
title(paste0('Monte Carlo Gillespie - ', ifelse(liberal,'Liberal ', 'Totalitarian '), 'System'), 
      line=1, outer = T, cex.main=1.8)

# close png file
if(file_needs_closing) dev.off()
