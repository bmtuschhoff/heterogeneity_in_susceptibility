# Estimate parameters for heterogeneity in suscepitibility given contact (HiS) in the presence of false negatives in the continuous case using simulated contact tracing data and Metropolis-Hastings MCMC
# Estimating the parameters that determine the distribution of individuals' risks of infection (k, theta)
# CORRECTING FOR FALSE NEGATIVES - switch to ABC method with C++ code

install.packages("iterators", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("mvtnorm", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("Rcpp", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("RcppGSL", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')

library("iterators", lib="../R_libs")
library("mvtnorm", lib="../R_libs")
library("Rcpp", lib="../R_libs")
library("RcppGSL", lib="../R_libs")

.libPaths("../R_libs")


# set C++ function as source code
sourceCpp("Cont_MCMCL_FNeg.cpp")


### read in contact tracing data (num inf in each contact network) and set relevant params
Focal <- 1000
N <- 5
N2 <- 5
k <- 1
theta <- 3
false_neg <- 0.1

inf2_prev <- unlist(read.csv("inf2_prevFNeg12_k1_theta3_F1000_N5_N25_len6e+05_sims600_fneg0.1seed.csv", row.names=1), use.names=F)
inf2_naive <- unlist(read.csv("inf2_naiveFNeg12_k1_theta3_F1000_N5_N25_len6e+05_sims600_fneg0.1seed.csv", row.names=1), use.names=F)

prev_inf <- sum(inf2_prev)
naive_inf <- sum(inf2_naive)

set.seed(3)

##############################################################################
### fit k and theta to data

# set arbitrary value for new prior
NewPrior <- 0

# set length of final MCMC chain and number of sims during likelihood calculation
steps <- 6e6
sims <- 100

OldL <- -Inf
while(is.na(OldL) || OldL == -Inf)
{
  # set original parameter set
  OldKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=matrix(c(0.01,-0.008,
                                                          -0.008,0.05), nrow=2, byrow=T)))
  Oldk <- OldKTheta[1]
  Oldtheta <- OldKTheta[2]
  
  # calculate log-likelihood for original parameter set
  OldL <- likelihood(Oldk, Oldtheta, Focal, N2, prev_inf[1], naive_inf[1], false_neg)
}

# calculate old prior, just prior for k, other priors flat/improper
OldPrior <- dexp(Oldk, rate=2, log=T)

# calculate old posterior
OldPost <- OldL + OldPrior

# initialize matrix to hold state data
CurrentState <- matrix(nrow=floor(steps/100), ncol=3)

# set name of file to send output to
fileName <- paste("ContMCFNeg12adj_k",k,"_theta",theta,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"_fneg",false_neg,"Cseed.csv", sep="")

for (i in 1:steps)
{
  # propose new parameter set -- proposal distribution = Lognormal
  NewKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=matrix(c(0.01,-0.008,
                                                          -0.008,0.05), nrow=2, byrow=T)))
  Newk <- Oldk * NewKTheta[1]
  Newtheta <- Oldtheta * NewKTheta[2]
  
  # calculate log-likelihood of new parameter set
  NewL <- likelihood(Newk, Newtheta, Focal, N2, prev_inf[1], naive_inf[1], false_neg)
  
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
