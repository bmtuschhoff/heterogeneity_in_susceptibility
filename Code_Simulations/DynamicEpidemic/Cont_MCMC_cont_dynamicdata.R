# Estimate parameters for heterogeneity in suscepitibility given contact (HiS) in the continuous case using simulated contact tracing data from a dynamic epidemic (dynam_epi_model_HiST_contactTracing_cont.cpp) and Metropolis-Hastings MCMC
# Estimating the parameters that determine the distribution of individuals' risks of infection (k, theta)
# Using previously generated data to continue running MCMC


install.packages("poisbinom", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("iterators", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("mvtnorm", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')

library("poisbinom", lib="../R_libs")
library("iterators", lib="../R_libs")
library("mvtnorm", lib="../R_libs")
library("utils", lib="../R_libs")


### read in contact tracing data (num inf in each contact network) and set relevant params

Focal <- 51
k <- 0.591715976331361
theta <- 0.626097060979711

contact_data <- read.csv("contact_data_dynamepi_k0.591715976331361_theta0.626097060979711_F51.csv", row.names=1)
inf2_prev <- contact_data$Focal_infected
inf2_naive <- contact_data$Naive_infected
exp2_prev <- contact_data$Focal_exposed
exp2_naive <- contact_data$Naive_exposed


set.seed(3)

##############################################################################
### fit k and theta to data
sims <- 600

# function to calculate log-likelihood
Log_like <- function(kval, thetaval, sims, inf2_prev, inf2_naive, exp2_prev, exp2_naive)
{
  # calculate new thetas for previously exposed groups after first round of infection
  theta_prev <- thetaval / (1 + thetaval)

  prev <- matrix(NA, nrow=sims, ncol=length(inf2_prev))
  naive <- matrix(NA, nrow=sims, ncol=length(inf2_prev))
  
  # pull risks from gamma distributions and calculate probs for likelihood of data for each contact network
  # multiple simulations to account for randomness of rgamma and probabilities
  for(i in (1:sims))
  {
      for(j in (1:length(inf2_prev)))
      {
        prev[i,j] <- dpoisbinom(inf2_prev[j], 1-exp(-rgamma(exp2_prev[j],shape=kval,scale=theta_prev)))
        naive[i,j] <- dpoisbinom(inf2_naive[j], 1-exp(-rgamma(exp2_naive[j],shape=kval,scale=thetaval)))
      }
  }

  
  # average across simulations, result=vector
  mean_prev <- apply(prev, 2, mean)
  mean_naive <- apply(naive, 2, mean)
  
  # computational errors sometimes give prob -1e-18 ~ 0
  mean_prev[mean_prev < 0] <- 0
  mean_naive[mean_naive < 0] <- 0
  
  # calculate log likelihood of the data given this combination of k and theta
  L <- sum(log(mean_prev)) + sum(log(mean_naive))

  return(L)
}




# set arbitrary value for new prior
NewPrior <- 0

# set length of final MCMC chain
steps <- 6e5

OldL <- -Inf
while(is.na(OldL) || OldL == -Inf)
{
  # set original parameter set
  OldKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=matrix(c(0.01,-0.008,
                                                          -0.008,0.05), nrow=2, byrow=T)))
  Oldk <- OldKTheta[1]
  Oldtheta <- OldKTheta[2]

  # calculate log-likelihood for original parameter set
  OldL <- Log_like(Oldk, Oldtheta, sims, inf2_prev, inf2_naive, exp2_prev, exp2_naive)
}


# calculate old prior, just prior for k, other priors flat/improper
OldPrior <- dexp(Oldk, rate=2, log=T)

# calculate old posterior
OldPost <- OldL + OldPrior

# initialize matrix to hold state data
CurrentState <- matrix(nrow=floor(steps/100), ncol=3)

# set name of file to send output to
#fileName <- paste("ContMC2gpsData_pA",p_a,"_pB",p_b,"_fA",f_a,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,".csv", sep="")
fileName <- paste("ContMCContdataDynamEpi_k",k,"_theta",theta,"_F",Focal,"_len",steps,"_sims",sims,"seed3.csv",sep="")

for (i in 1:steps)
{
  # propose new parameter set -- proposal distribution = Lognormal
  NewKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=matrix(c(0.01,-0.008,
                                                          -0.008,0.05), nrow=2, byrow=T)))
  Newk <- Oldk * NewKTheta[1]
  Newtheta <- Oldtheta * NewKTheta[2]
  
  # calculate log-likelihood of new parameter set
  NewL <- Log_like(Newk, Newtheta, sims, inf2_prev, inf2_naive, exp2_prev, exp2_naive)

  # calculate new prior, just prior for k, other priors flat/improper
  NewPrior <- dexp(Newk, rate=2, log=T)
  
  # calculate new posterior
  NewPost <- NewL + NewPrior
  
  # if new posterior better than old posterior, jump to new parameter set
  if(NewPost - OldPost > 0)
  {
    Oldk <- Newk
    Oldtheta <- Newtheta
    OldL <- NewL
    OldPost <- NewPost
  }
  # if new posterior is worse than old posterior, jump to new parameter set with probability
  else if(exp(NewPost - OldPost) > runif(1))
  {
    Oldk <- Newk
    Oldtheta <- Newtheta
    OldL <- NewL
    OldPost <- NewPost
  }
  
  # record current state and log-likelihood for every 100th observation
  if(i %% 100 == 0)
    CurrentState[i/100,] <- c(Oldk, Oldtheta, OldL)
  
  # write current state (since past time output) to file at every 1500th observation
  if(i %% 1500 == 0)
    write.table(CurrentState[(((i-1500)/100 +1):(i/100)),], file=fileName, append=T, sep=",", col.names=NA, row.names=T)
}


final_i <- steps %% 1500
if(final_i == 0)
{
    write.table(CurrentState[(((steps-1500)/100 +1):(steps/100)),], file=fileName, append=T, sep=",", col.names=NA, row.names=T)
} else
{
    write.table(CurrentState[(((steps-final_i)/100 +1):(steps/100)),], file=fileName, append=T, sep=",", col.names=NA, row.names=T)
}
