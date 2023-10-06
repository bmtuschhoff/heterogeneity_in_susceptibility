# Estimate parameters for heterogeneity in suscepitibility given contact (HiS) in the discrete case (with 2 types) using simulated contact tracing data and Metropolis-Hastings MCMC with ABC
# Estimating the parameters that determine the distribution of individuals' risks of infection (p_A, p_B, f_A)
# Using previously generated data to continue running MCMC


install.packages("coda", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("mcmcplots", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')


library("coda", lib="../R_libs/")
library("mcmcplots", lib="../R_libs/")


### read in contact tracing data (num inf in each contact network) and set relevant params

## for testing effect of wrong underlying model
#Focal <- 50
#N <- 5
#N2 <- 5
#k <- 0.591715976331361
#theta <- 0.626097060979711

#inf2_prev <- unlist(read.csv("inf2_prev_k0.591715976331361_theta0.626097060979711_F50_N5_N25_len6e+05_sims600seed.csv", row.names=1), use.names=F)
#inf2_naive <- unlist(read.csv("inf2_naive_k0.591715976331361_theta0.626097060979711_F50_N5_N25_len6e+05_sims600seed.csv", row.names=1), use.names=F)

#prev_inf <- sum(inf2_prev)
#naive_inf <- sum(inf2_naive)


## for discrete generated data
Focal <- 5000
N <- 100
N2 <- 100
p_a <- 0.74809046024617
p_b <- 0.125477384938458
f_a <- 0.2

inf2_prev <- unlist(read.csv("inf2_prev_2gps_pA0.74809046024617_pB0.125477384938458_fA0.2_F5000_N100_N2100seed.csv", row.names=1), use.names=F)
inf2_naive <- unlist(read.csv("inf2_naive_2gps_pA0.74809046024617_pB0.125477384938458_fA0.2_F5000_N100_N2100seed.csv", row.names=1), use.names=F)

prev_inf <- sum(inf2_prev)
naive_inf <- sum(inf2_naive)

set.seed(3)

##############################################################################
### fit p_a, p_b, and f_a_old to data using MCMC with Metropolis-Hastings algorithm and ABC

sims <- 100

# function to calculate log-likelihood
Log_Like <- function(pa_val, pb_val, fa_val, Focal, N2, sims, prev_inf, naive_inf)
{
  # avg/expected fraction that are type a in the 2nd round
  fa2 <- (fa_val*(1-pa_val)) / (fa_val*(1-pa_val) + (1-fa_val)*(1-pb_val))
  
  # number of focal individuals that are type A at end of 1st round
  FocalA <- rbinom(sims, Focal, fa2)
  
  # number of type A and B focal individuals infected in 2nd round
  FocalAinf <- rbinom(sims, FocalA, pa_val)
  FocalBinf <- rbinom(sims, Focal-FocalA, pb_val)
  
  # number of non-focal individuals that are type A at beginning of 2nd round
  NonFocalA <- rbinom(sims, (N2-1)*Focal, fa_val)
  
  # number of type A and B non-focal individuals infected in 2nd round
  NonFocalAinf <- rbinom(sims, NonFocalA, pa_val)
  NonFocalBinf <- rbinom(sims, (N2-1)*Focal - NonFocalA, pb_val)
  
  # total number of individuals infected in 2nd round, focal and non-focal
  TotFocalInf <- FocalAinf + FocalBinf
  TotNonFocalInf <- NonFocalAinf + NonFocalBinf
  
  # calc fraction of simulations where TotFocalInf=prev_inf and TotNonFocalInf=naive_inf
  LFocal <- sum(abs(TotFocalInf - prev_inf) / prev_inf <= 0.0) / sims
  LNonFocal <- sum(abs(TotNonFocalInf - naive_inf) / naive_inf <= 0.0) / sims
  
  # log likelihood of the data given this combination of parameters
  L <- log(LFocal) + log(LNonFocal)
  
  return(L)
}


# calculate old prior
OldPrior <- 0
# set arbitrary value for new prior
NewPrior <- 0

# set length of final MCMC chain
steps <- 30e6

OldL <- -Inf
while(OldL == -Inf)
{
  # set original parameter set
  Oldpa <- runif(1)
  Oldpb <- runif(1, 0, Oldpa)
  Oldfa <- runif(1)
  
  # calculate log-likelihood for original parameter set
  OldL <- Log_Like(Oldpa, Oldpb, Oldfa, Focal, N2, sims, prev_inf, naive_inf)
}

# calculate old posterior
OldPost <- OldL + OldPrior

# initialize matrix to hold state data
CurrentState <- matrix(nrow=floor(steps/1500), ncol=4)

# set name of file to send output to
filename <- paste("2gpsMC2gpsData_pA",p_a,"_pB",p_b,"_fA",f_a,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"tol0seed.csv", sep="")

for(i in 1:steps)
{
  # propose new parameter set -- proposal distribution = Unif(0,1) for all parameters
  Newpa <- runif(1, 0, 1)
  Newpb <- runif(1, 0, Newpa)
  Newfa <- runif(1, 0, 1)
  
  # calculate log-likelihood of new parameter set
  NewL <- Log_Like(Newpa, Newpb, Newfa, Focal, N2, sims, prev_inf, naive_inf)
  
  # calculate new prior
  NewPrior <- 0
  
  # calculate new posterior
  NewPost <- NewL + NewPrior
  
  # if new posterior better than old posterior, jump to new parameter set
  if(NewPost - OldPost > 0)
  {
    Oldpa <- Newpa
    Oldpb <- Newpb
    Oldfa <- Newfa
    OldL <- NewL
    OldPost <- NewPost
  }
  # if new posterior is worse than old posterior, jump to new parameter set with probability
  else if(exp(NewPost - OldPost) > runif(1))
  {
    Oldpa <- Newpa
    Oldpb <- Newpb
    Oldfa <- Newfa
    OldL <- NewL
    OldPost <- NewPost
  }
  
  # record current state and log-likelihood for every 1500th observation
  if(i %% 1500 == 0)
    CurrentState[i/1500,] <- c(Oldpa, Oldpb, Oldfa, OldL)

  # write current state (since past time output) to file at every 22,500th observation
  if(i %% 22500 == 0)
    write.table(CurrentState[(((i-22500)/1500 +1):(i/1500)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}

final_i <- steps %% 22500
if(final_i == 0)
{
    write.table(CurrentState[(((steps-22500)/1500 +1):(steps/1500)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
} else
{
    write.table(CurrentState[(((steps-final_i)/1500 +1):(steps/1500)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}


