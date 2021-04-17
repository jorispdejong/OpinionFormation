###########################################################################
###########      Solving the Master Equation Numerically       ############
###########           Example: Opinion formation               ############
###########        Krylov Subspace Approximation Method        ############
###########          Authors: Joris de Jong                    ############
###########           Last Modified: 15/04/2021                ############
###########################################################################
# clear global environment
rm(list = ls(all.names = TRUE))

# libraries
suppressPackageStartupMessages(library('Matrix')) # sparseMatrix()

# source
source('./Functions.R')

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

# liberal system: ELASPED TIME = 40.75 sec
# totalitarian system: ELASPED TIME = 38.64 sec
liberal <- T

if(liberal){
  a1 <- 0/L
  a2 <- 0.5/L
  a3 <- 0.5/L
  
  file_name_probabilities_obj <- 'KSA/KSA_probabilities_liberal.RData'
}else{
  a1 <- 1.5/L
  a2 <- 1/L
  a3 <- -0.25/L
  
  file_name_probabilities_obj <- 'KSA/KSA_probabilities_totalitarian.RData' 
}

total_time <- 15 # simulation time in days
dt <- 0.1 # time step in days
time <- seq(0,total_time,dt) # vector with accumulating time steps

##################
### SIMULATION ###
##################
# check if the R object probabilities already exists, then we simply load it without running the simulation
if(file.exists(file_name_probabilities_obj)){
  load(file_name_probabilities_obj)
}else{
  if(verbose) print('Running Krylov subspace approximation method... Please wait')
  
  # measure elapsed time
  system.time(
    {
      ############################################
      ### KRYLOV SUBSPACE APPROXIMATION METHOD ###
      ############################################
      # simulation parameters
      m <- 10 # Krylov subspace dimension
      K <- (2*L+1) # dimension of the Q matrix
      K2 <- K*K # total number of states
      
      # create matrix P by creating row and column indices with the corresponding values
      # later these vectors will be used to create a sparse matrix
      row <- rep(0,5*K)
      col <- rep(0,5*K)
      value <- rep(0,5*K)
      count <- 1
      for(j in 1:K){
        j <- 1
        for(i in 1:K){
          i <- 2
          # find propensities
          propensities <- propensityVector(c(j-(L+1), i-(L+1)))
          
          # add value to diagonal
          diag_location <- (j-1)*K+i
          row[count] <- diag_location
          col[count] <- diag_location
          value[count] <- -sum(propensities)
          count <- count+1
          
          # add propensity for reaction 1
          if(propensities[1]>0){
            row[count] <- j*K+i # location with one more y1
            col[count] <- diag_location
            value[count] <- propensities[1]
            count <- count+1
          }
          
          # add propensity for reaction 2
          if(propensities[2]>0){
            row[count] <- (j-2)*K+i # location with one less y1
            col[count] <- diag_location
            value[count] <- propensities[2]
            count <- count+1
          }
          
          # add propensity for reaction 3
          if(propensities[3]>0){
            row[count] <- (j-1)*K+i+1 # location with one more y2
            col[count] <- diag_location
            value[count] <- propensities[3]
            count <- count+1
          }
          
          # add propensity for reaction 4
          if(propensities[4]>0){
            row[count] <- (j-1)*K+i-1 # location with one less y2
            col[count] <- diag_location
            value[count] <- propensities[4]
            count <- count+1
          }
        }
      }
      
      # Save the required matrix as a sparse matrix (saves a lot of memory ~300Mb)
      P <- sparseMatrix(row[value!=0], col[value!=0], x=value[value!=0], dims=c(K2,K2))
      
      # empty array in which will hold the probabilities over time
      probabilities <- array(0, dim = c(K, K, length(time)))
      
      # Set initial condition
      p_old <- matrix(0,nrow(P),1)
      p_old[L*K+L+1,] <- 1 # x1=0, x2=0, so probability is 1
      probabilities[,,1] <- matrix(p_old,K,K,byrow = T)
      
      # Solve p_new = exp(dt*P)*p_old recursively
      pb <- txtProgressBar(2, length(time)) # progress bar
      for(t in 2:length(time)){
        # The Matlab package EXPOKIT's function expv() was rewritten in R,
        # which solves w = exp(t*A)*v and returns w
        p_new <- expv(t=dt, A=P,v=p_old,tol=1e-6,m=m)
        p_old <- p_new
        
        # store probability
        probabilities[,,t] <- matrix(p_old,K,K, byrow = T)
        
        # update progress bar
        setTxtProgressBar(pb, t) 
      }
      
      # save probabilities as R object 
      save(probabilities, file = file_name_probabilities_obj)
    }
  )
}


####################
### PLOT RESULTS ###
####################
# get the probabilities of net public/private opinion on the final time step
probabilities_public <- rowSums(probabilities[,,length(time)])
probabilities_private <- colSums(probabilities[,,length(time)])
bins <- -L:L

# compute the expected mean and standard deviation
Y1_probabilities <- apply(probabilities, c(1,3), sum)
Y2_probabilities <- apply(probabilities, c(2,3), sum)

Y1_mean <- rowSums(sapply(1:length(bins), function(i) bins[i]*Y1_probabilities[i,]))
Y2_mean <- rowSums(sapply(1:length(bins), function(i) bins[i]*Y2_probabilities[i,]))

Y1_std <- sqrt(rowSums(sapply(1:length(bins), function(i) ((bins[i]-Y1_mean)^2)*Y1_probabilities[i,])))
Y2_std <- sqrt(rowSums(sapply(1:length(bins), function(i) ((bins[i]-Y1_mean)^2)*Y2_probabilities[i,])))

# create folder 'results' that will contain the resulting plots
dir.create('KSA/results', showWarnings = F) # if already exists it only gives a warning (which we surpress)

# save a png file with the results if the file does not exist already
if(liberal){
  png_file_name <- 'KSA/results/KSA_liberal.png'
}else{
  png_file_name <- 'KSA/results/KSA_totalitarian.png'
}

file_needs_closing <- F
if(!file.exists(png_file_name)){
  # open png file to store the results
  png(filename = png_file_name, width = 1000, height = 640)
  file_needs_closing <- T
}

# plot the results
par(mfrow=c(2,2), mar=c(5.1, 4.1, 2.1, 2.1), oma=c(0,1,4,0)) # set plotting window

# histograms showing the net public/private opinion after 15 days
plot(bins, probabilities_public, type = 'h', lwd = 3,
     xlab = 'net public opinion', ylab = 'probability', cex.lab=1.2)
lines(bins, probabilities_public)
plot(bins, probabilities_private, type = 'h', lwd = 3,
     xlab = 'net private opinion', ylab = 'probability', cex.lab=1.2)
lines(bins, probabilities_private)
title('Probability distribution at t=15 days', line=-1, cex.main=1.17, outer = T)

# plots showing the expected net public/private opinion over time
plot(seq(0, total_time, dt), Y1_mean+Y1_std, type = 'l',
     ylab = 'Net public opinion', xlab = 'time (days)', cex.lab=1.2,
     lwd=2, ylim=c(-L,L))
lines(seq(0, total_time, dt), Y1_mean-Y1_std, lwd=2, col=2)
plot(seq(0, total_time, dt), Y2_mean+Y2_std, type = 'l',
     ylab = 'Net private opinion', xlab = 'time (days)', cex.lab=1.2,
     lwd=2, ylim=c(-L,L))
lines(seq(0, total_time, dt), Y2_mean-Y2_std, lwd=2, col=2)
title('Expected net public/private opinion over time', line=-25.5, cex.main=1.17, outer = T)

# main title
title(paste0('Krylov Subspace Approximation - ', ifelse(liberal,'Liberal ', 'Totalitarian '), 'System'), 
      line=1, outer = T, cex.main=1.8)

# close png file
if(file_needs_closing) dev.off()


# save a gif with the probabilities in the x1-x2 plane over time (if the file does not exist already)
if(liberal){
  gif_file_name <- 'KSA_liberal.gif'
}else{
  gif_file_name <- 'KSA_totalitarian.gif'
}
if(!file.exists(paste0('KSA/results/',gif_file_name))){
  # Create 3d animation
  suppressPackageStartupMessages(library('animation')) # make animation of plots
  suppressPackageStartupMessages(library('plot3D')) # 3d surface plot
  
  if(verbose) print('Creating GIF... Please wait')
  saveGIF(
    {
      pb <- txtProgressBar(1, length(time)) # progress bar
      for(i in 1:length(time)){
        # surface plot of probabilities at time i
        persp3D(bins, bins, probabilities[,,i], 
                theta = 20, phi = 25, expand = 0.4, ticktype = 'detailed', 
                cex.lab=1.2, nticks = 4,
                xlab = 'net public opinion', ylab = 'net private opinion', zlab = 'probability', 
                main = paste0('Krylov subspace approximation - ', 
                              ifelse(liberal,'Liberal ', 'Totalitarian '), 'System', 
                              ', t = ',
                              sprintf('%.1f',time[i]),
                              ' (days)'))
        
        # update progress bare
        setTxtProgressBar(pb, i)
      }
    }, 
    movie.name = gif_file_name, 
    interval = 0.1, ani.width = 800, ani.height = 550
  ) 
}
# move file to results folder 
# (required work around, because gif cannot be exported straight to results folder)
file.rename(from = gif_file_name, to = paste0('KSA/results/',gif_file_name))
