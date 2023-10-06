# Estimate parameters for heterogeneity in suscepitibility given contact (HiS) in the presence of false negatives in the continuous case using simulated contact tracing data and Metropolis-Hastings MCMC
# Estimating the parameters that determine the distribution of individuals' risks of infection (k, theta)

install.packages("poisbinom", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("iterators", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("mvtnorm", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')

library("poisbinom", lib="../R_libs")
library("iterators", lib="../R_libs")
library("mvtnorm", lib="../R_libs")
library("utils", lib="../R_libs")

arg <- commandArgs(trailingOnly=T)

# number of focal individuals
Focal <- as.numeric(arg[1])

# number of sims to run
num_sims <- 1

# total initial population size
N <- 5
# population size for second exposure event
N2 <- 5

# rate of false negatives
false_neg <- as.numeric(arg[2])

# cv of risk and expected fraction of naive individuals infected
cvs <- 1
exp_frac_inf <- 0.75

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0, length(cvs)*length(exp_frac_inf), 5)

set.seed(3)

combo <- 1
for(cv in cvs)
{
  for(f_inf in exp_frac_inf)
  {
    # determine k and theta to use for gamma dist
    k <- 1 / (cv)^2
    theta <- ((1 - f_inf)^(-1/k) - 1)
        
    # initialize vectors to store likelihood tests later on
    ratios <- rep(0, num_sims)

    # vector to track focal indivs
    foc <- rep(0, num_sims)

    # number of infected individuals 1st and 2nd rounds
    inf2_prev <- matrix(NA, Focal, num_sims)
    inf2_naive <- matrix(NA, Focal, num_sims)

    
    for(sim in (1:num_sims))
    {
      risk1 <- c()	# initialize vector for risks of focal indivs
      j <- 1	# initialize counter var to save risks

      # collect focal indivs - simulate contact networks until get enough focal indivs (not inf)
      while(foc[sim] < Focal)
      {
        ## first exposure events
        # set up relative risk distributions for first infection rounds
	if(cv == 0)
		risk <- rep(-log(1-f_inf), N)
	else
        	risk <- rgamma(N, shape=k, scale=theta)

        # probability of being inf
        prob_I <- 1 - exp(-risk)

        # infect individuals
	# use bernoulli for each indiv to check if infected, store risk for those not infected
        for(i in (1:(N)))
        {
          I <- rbinom(1, 1, prob_I[i])
          if(I == 1)  # indiv infected, check if false neg
          {
            check_fn <- rbinom(1, 1, false_neg)
            if(check_fn == 1)  # indiv is false neg -- could be focal with prob 0 of inf
            {
              risk1[j] <- 0
              j <- j + 1
	      foc[sim] <- foc[sim] + 1
            }
          }

          if(I == 0)  # indiv not infected
          {
            risk1[j] <- risk[i]
            j <- j + 1
	    foc[sim] <- foc[sim] + 1
          }
        }
      }



      for(f in (1:Focal))
      {     
        # create risk and prob of not infected distribution for naive individuals
	if(cv == 0)
		risk0 <- rep(-log(1-f_inf), N2-1)
	else
        	risk0 <- rgamma(N2-1, shape=k, scale=theta)

        # probabilities of getting infected
        prob_I0 <- 1 - exp(-risk0)
        prob_I1 <- 1 - exp(-risk1[f]) # risk1[f] is focal indiv's risk
        
        # simulate second round of infection
        # focal indiv, 1 if inf
        I1 <- rbinom(1, 1, prob_I1)
        
        # naive indivs
        risk01 <- c()
        j <- 1
        for(i in (1:(N2-1)))
        {
          I <- rbinom(1, 1, prob_I0[i])
          if(I == 0)  # indiv not infected
          {
            risk01[j] <- risk0[i]
            j <- j + 1
          }
        }
        
        
        # number indivs infected in 2nd round, remove some indivs b/c false neg
        inf2_naive[f,sim] <- rbinom(1, N2-1-length(risk01), 1-false_neg)
        inf2_prev[f,sim] <- rbinom(1, I1, 1-false_neg)

        
        ##########################################################################################################
      }
    }
    
    # sum across all Focal groups for number individuals infected in 2nd round for naive vs focal indivs, vector len=sims
    naive_inf <- apply(inf2_naive, 2, sum)
    prev_inf <- apply(inf2_prev, 2, sum)
    
    # calculate observed probabilities of inf across all focal indivs as best estimates
    p_naive_obs <- naive_inf / (Focal*(N2 - 1))
    p_prev_obs <- prev_inf / Focal
    p_avg_obs <- (naive_inf + prev_inf) / (Focal*N2)
    
    
    # calculate log-likelihoods for homogeneity vs heterogeneity
    Lhom <- dbinom(naive_inf, Focal*(N2-1), p_avg_obs, log=T) + dbinom(prev_inf, Focal, p_avg_obs, log=T)
    Lhet <- dbinom(naive_inf, Focal*(N2-1), p_naive_obs, log=T) + dbinom(prev_inf, Focal, p_prev_obs, log=T)
    
    # calculate likelihood ratio test statistics
    ratios <- -2*(Lhom - Lhet)
    
    # find power = percent of times the ratio test is greater than critical value 3.84 (df=1)
    power <- sum(ratios > 3.84, na.rm=T)/length(ratios) * 100

    # record power, cv, f_inf in table
    powers_table[combo,] <- c(cv, f_inf, power, mean(p_naive_obs), mean(p_prev_obs))
    combo <- combo + 1
  }
}

powers_table <- as.data.frame(powers_table)
colnames(powers_table) <- c("cv", "f_inf", "power", "p_naive", "p_prev")

### write contact networks (num inf) to files to keep running MCMC chain with last accepted params
sims <- 600
steps <- 6e5
naiveFile <- paste("inf2_naiveFNeg12_k",k,"_theta",theta,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"_fneg",false_neg,"seed.csv", sep="")
prevFile <- paste("inf2_prevFNeg12_k",k,"_theta",theta,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"_fneg",false_neg,"seed.csv", sep="")
write.table(inf2_naive, file=naiveFile, sep=",")
write.table(inf2_prev, file=prevFile, sep=",")


##############################################################################
### fit k and theta to data
sims <- 600

# function to calculate log-likelihood
Log_like <- function(kval, thetaval, Focal, N2, sims, inf2_prev, inf2_naive)
{
  # calculate new thetas for previously exposed groups after first round of infection
  theta_prev <- thetaval / (1 + thetaval)

  prev <- matrix(NA, nrow=sims, ncol=Focal)
  naive <- matrix(NA, nrow=sims, ncol=Focal)
  
  # pull risks from gamma distributions and calculate probs for likelihood of data
  # multiple simulations to account for randomness of rgamma and probabilities
  for(i in (1:sims))
  {
      prev[i,] <- dpoisbinom(inf2_prev, 1-exp(-rgamma(1,shape=kval,scale=theta_prev)))
      naive[i,] <- dpoisbinom(inf2_naive, 1-exp(-rgamma(N2-1,shape=kval,scale=thetaval)))
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
  OldL <- Log_like(Oldk, Oldtheta, Focal, N2, sims, inf2_prev, inf2_naive)
}

# calculate old prior, just prior for k, other priors flat/improper
OldPrior <- dexp(Oldk, rate=2, log=T)

# calculate old posterior
OldPost <- OldL + OldPrior

# initialize matrix to hold state data
CurrentState <- matrix(nrow=floor(steps/100), ncol=3)

# set name of file to send output to
fileName <- paste("ContMCFNeg12_k",k,"_theta",theta,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"_fneg",false_neg,"seed.csv", sep="")


for (i in 1:steps)
{
  # propose new parameter set -- proposal distribution = Lognormal
  NewKTheta <- exp(rmvnorm(1, mean=c(0,0), sigma=matrix(c(0.01,-0.008,
                                                          -0.008,0.05), nrow=2, byrow=T)))
  Newk <- Oldk * NewKTheta[1]
  Newtheta <- Oldtheta * NewKTheta[2]
  
  # calculate log-likelihood of new parameter set
  NewL <- Log_like(Newk, Newtheta, Focal, N2, sims, inf2_prev, inf2_naive)

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