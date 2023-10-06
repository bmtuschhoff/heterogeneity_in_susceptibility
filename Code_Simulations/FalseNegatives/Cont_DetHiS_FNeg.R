# Calculate power to detect heterogeneity in suscepitibility given contact (HiS) in the presence of false negatives in the continuous case across parameter combinations using simulated contact tracing data and a likelihood ratio test
# Test if can tell difference between naive individuals and previously exposed (focal) individuals
# If there is a difference, there is HiS in the population.
# If there is not a difference, there is not HiS.


# number of focal individuals
Focal <- 200

# number of sims to run to determine power
num_sims <- 1000

# total initial population size
N <- 5
# population size for second exposure event
N2 <- 5

# rate of false negatives
false_neg <- 0.1

# graph parameters
cvs <- seq(0,3,0.02)
exp_frac_inf <- seq(0.02,0.98,0.02)

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0, length(cvs)*length(exp_frac_inf), 5)

# file to save data to
filename <- paste("ContPowersFNeg12_F",Focal,"_N",N,"_N2",N2,"_sims",num_sims,"_fneg",false_neg,"_c3e0.98seed.csv",sep="")

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

    # number of infected individuals in 2nd exposure event
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
	# risk is gamma distributed, but if no HiS, all individuals have same risk based on f_inf
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
        
        # simulate second exposure events
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

  # record powers in file
  write.table(powers_table[((combo-49):(combo-1)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}

