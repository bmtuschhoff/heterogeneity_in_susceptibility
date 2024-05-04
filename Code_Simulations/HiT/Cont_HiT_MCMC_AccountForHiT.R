# Estimate parameters for heterogeneity in suscepitibility given contact (HiS) and heterogeneity in transmission given contact in the continuous case using simulated contact tracing data and Metropolis-Hastings MCMC
# Estimating the parameters that determine the distributions of individuals' risks of infection (k, theta) and forces of infection (m)
# Using previously generated data to run MCMC

install.packages("mvtnorm", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
library("mvtnorm", lib="../R_libs")

Focal <- 1000
N <- 5
N2 <- 5
k <- 0.591715976331361
theta <- 0.626097060979711
m <- 0.5
phi <- 1/m

inf1 <- unlist(read.csv("inf1_HiT_k0.591715976331361_theta0.626097060979711_F1000_N5_N25_len6e+05_sims600_m0.5phi2seed3.csv", row.names=1), use.names=F)
inf2_naive <- unlist(read.csv("inf2_naiveHiT_k0.591715976331361_theta0.626097060979711_F1000_N5_N25_len6e+05_sims600_m0.5phi2seed3.csv", row.names=1), use.names=F)
inf2_prev <- unlist(read.csv("inf2_prevHiT_k0.591715976331361_theta0.626097060979711_F1000_N5_N25_len6e+05_sims600_m0.5phi2seed3.csv", row.names=1), use.names=F)

# calculate probability and risk of naive indivs being infected
p_naive <- sum(inf2_naive) / (Focal * (N2 - 1))


set.seed(3)

##############################################################################
### fit m and phi to data
sims <- 6000  # sims is number of networks to run for comparing with real data

# function to calculate log-likelihood for estimating HiS and HiT params
Log_like <- function(kval, thetaval, mval, Focal, N2, sims, inf1, inf2_naive, inf2_prev)
{
  inf1_mcmc <- rep(NA, sims)
  inf2_naive_mcmc <- rep(NA, sims)
  inf2_prev_mcmc <- rep(NA, sims)

  # track number focal indivs/networks
  foc <- 0

  # pull risks from gamma distributions and calculate probs for likelihood of data
  # multiple simulations to account for randomness of rgamma and probabilities
  # sort so that highest prob of inf goes with group with most inf
  risk1 <- c()	# initialize vector for risks of focal indivs
  j <- 1	# initialize counter var to save risks
  f <- 1	# initalize counter var to save number inf in each contact network

      # collect focal indivs - simulate contact networks until get enough focal indivs (not inf)
      while(foc < sims)
      {
        ## first exposure events
        # set up relative risk distributions for first infection rounds
        risk <- rgamma(N2, shape=kval, scale=thetaval)

	# draw transmission/force of infection parameter
	lamb <- rgamma(1, shape=mval, scale=1/mval)

        # probability of being inf
        prob_I <- 1 - exp(-lamb*risk)

        # infect individuals
	num_notI1 <- 0		# number of indivs not inf in this contact network
	# use bernoulli for each indiv to check if infected, store risk for those not infected
        for(i in (1:(N2)))
        {
          I <- rbinom(1, 1, prob_I[i])
          if(I == 0)  # indiv not infected
          {
            risk1[j] <- risk[i]
            j <- j + 1
	    foc <- foc + 1
	    num_notI1 <- num_notI1 + 1
          }
        }
	# save num inf in each network and repeat network if multiple focal indivs from this network...
	for(x in 1:num_notI1)
	{
	  if(f <= sims)
	  {
	  	inf1_mcmc[f] <- N - num_notI1
	  	f <- f + 1
	  }
	}

      }


      for(f in (1:sims))
      {      
        # create risk and prob of not infected distribution for naive individuals
        risk0 <- rgamma(N2-1, shape=kval, scale=thetaval)

	# draw transmission/force of infection parameter
	lamb2 <- rgamma(1, shape=mval, scale=1/mval)
        
        # probabilities of getting infected
        prob_I0 <- 1 - exp(-lamb2*risk0)
        prob_I1 <- 1 - exp(-lamb2*risk1[f]) # risk1[f] is focal indiv's risk
        
        # simulate second round of infection
	# focal indiv, 1 if inf
        I1 <- rbinom(1,1,prob_I1)
        
	# naive indivs
        risk01 <- c()
        j <- 1
        for(i in (1:(N2-1)))
        {
          I <- rbinom(1,1,prob_I0[i])
          if(I == 0)  # not infected
          {
            risk01[j] <- risk0[i]
            j <- j + 1
          }
        }
        
        
        # number infected in 2nd round
	inf2_naive_mcmc[f] <- N2 - 1 - length(risk01)
	inf2_prev_mcmc[f] <- I1
        ##########################################################################################################
      }    

  
   num_in_net <- seq(0, N2-1, 1)	# sequence of number of potential individuals inf in a network

   # set up dataframe with number networks with x inf 1st network, y inf 2nd network naive, z inf 2nd network focal
   numI_bin <- data.frame(matrix(0, nrow=length(num_in_net)^2*2, ncol=5))
   colnames(numI_bin) <- c("num_networks_sim", "num_networks_data", "focal", "round1", "round2")
   # set up x inf round1, y inf round2, z inf focal
   numI_bin$focal <- c(rep(0, nrow(numI_bin)/2), rep(1, nrow(numI_bin)/2))
   numI_bin$round1 <- rep(rep(num_in_net, each=length(num_in_net)), 2)
   numI_bin$round2 <- rep(num_in_net, length(num_in_net)*2)

   # for each simulated data point, add that network into numI_bin where it matches x, y, and z
   for(f in 1:sims)
   {
	for(row in 1:nrow(numI_bin))
	{
	  if(inf2_prev_mcmc[f] == numI_bin$focal[row] && inf1_mcmc[f] == numI_bin$round1[row] && inf2_naive_mcmc[f] == numI_bin$round2[row])
	  {
		numI_bin$num_networks_sim[row] <- numI_bin$num_networks_sim[row] + 1
		break
	  }
	}
   }

   # divide number of networks by num sims to get probability of each scenario
   numI_bin$num_networks_sim <- (numI_bin$num_networks_sim + 1) / (sims + 2)

   # for each real data point, add that network into numI_bin where it matches x, y, and z
   for(f in 1:Focal)
   {
	for(row in 1:nrow(numI_bin))
	{
	  if(inf2_prev[f] == numI_bin$focal[row] && inf1[f] == numI_bin$round1[row] && inf2_naive[f] == numI_bin$round2[row])
	  {
		numI_bin$num_networks_data[row] <- numI_bin$num_networks_data[row] + 1
		break
	  }
	}
   }


  # calculate log likelihood of the data
  L <- dmultinom(numI_bin$num_networks_data, prob=numI_bin$num_networks_sim, log=T)
  
  return(L)
}






# set arbitrary value for new prior
NewPrior <- 0

# set length of final MCMC chain
steps <- 6e5

# set up covariance matrix
cov_matrix <- matrix(c(0.01,-0.008,
                      -0.008,0.05), nrow=2, byrow=T)


OldL <- -Inf
while(is.na(OldL) || OldL == -Inf)
{
  # set original parameter set
  OldKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=cov_matrix))
  Oldk <- OldKTheta[1]
  Oldtheta <- OldKTheta[2]
  Oldm <- exp(rnorm(1, mean=0, sd=0.2))

  # calculate log-likelihood for original parameter set
  OldL <- Log_like(Oldk, Oldtheta, Oldm, Focal, N2, sims, inf1, inf2_naive, inf2_prev)
}



# calculate old prior, just prior for k, other priors flat/improper
OldPrior <- dexp(Oldm, rate=2, log=T) + dexp(Oldk, rate=1, log=T) + dexp(Oldtheta, rate=0.2, log=T)

# calculate old posterior
OldPost <- OldL + OldPrior

# initialize matrix to hold state data
CurrentState <- matrix(nrow=floor(steps/1), ncol=5)

# set name of file to send output to
fileName <- paste("ContMCHiTEstHiT_k",k,"_theta",theta,"_m",m,"_phi",phi,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"msep3seed.csv", sep="")

for (i in 1:steps)
{
  ################### update m ################################
  # propose new parameter m -- proposal distribution = Lognormal
  NewM_prop <- exp(rnorm(1, mean=0, sd=0.2))
  Newm <- Oldm * NewM_prop

  # calculate log-likelihood of new parameter set
  NewL <- Log_like(Oldk, Oldtheta, Newm, Focal, N2, sims, inf1, inf2_naive, inf2_prev)

  # calculate new prior
  NewPrior <- dexp(Newm, rate=2, log=T) + dexp(Oldk, rate=1, log=T) + dexp(Oldtheta, rate=0.2, log=T)

  # calculate new posterior
  NewPost <- NewL + NewPrior
  
  # if new posterior better than old posterior, jump to new parameter set
  if(NewPost - OldPost > 0)
  {
    Oldm <- Newm
    OldL <- NewL
    OldPost <- NewPost
  }
  # if new posterior is worse than old posterior, jump to new parameter set with probability
  else if(exp(NewPost - OldPost) > runif(1))
  {
    Oldm <- Newm
    OldL <- NewL
    OldPost <- NewPost
  }


  ################### update k, theta ################################
  # propose new parameter set for k and theta -- proposal distribution = Lognormal
  NewKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=cov_matrix))
  Newk <- Oldk * NewKTheta[1]
  Newtheta <- Oldtheta * NewKTheta[2]

  # calculate log-likelihood of new parameter set
  NewL <- Log_like(Newk, Newtheta, Oldm, Focal, N2, sims, inf1, inf2_naive, inf2_prev)

  # calculate new prior
  NewPrior <- dexp(Oldm, rate=2, log=T) + dexp(Newk, rate=1, log=T) + dexp(Newtheta, rate=0.2, log=T)
  
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

  
  # record current state and log-likelihood
  CurrentState[i,] <- c(Oldk, Oldtheta, Oldm, 1/Oldm, OldL)

  # write current state (since past time output) to file at every 100th observation
  if(i %% 100 == 0)
    write.table(CurrentState[(((i-100)/1 +1):(i)),], file=fileName, append=T, sep=",", col.names=NA, row.names=T)
}

final_i <- steps %% 100
if(final_i == 0)
{
    write.table(CurrentState[(((steps-100)/1 +1):(steps/1)),], file=fileName, append=T, sep=",", col.names=NA, row.names=T)
} else
{
    write.table(CurrentState[(((steps-final_i)/1 +1):(steps/1)),], file=fileName, append=T, sep=",", col.names=NA, row.names=T)
}
