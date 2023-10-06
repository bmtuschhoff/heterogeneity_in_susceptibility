# Calculate power to detect heterogeneity in suscepitibility given contact (HiS) in the presence of heterogeneity in transmission in the continuous case across parameter combinations using simulated contact tracing data and a likelihood ratio test
# Test if can tell difference between naive individuals and previously exposed (focal) individuals
# If there is a difference, there is HiS in the population.
# If there is not a difference, there is not HiS.


arg <- commandArgs(trailingOnly=T)

# number of focal individuals
Focal <- as.numeric(arg[1])

# number of sims to run to determine power
num_sims <- 1000

# total initial population size
N <- 5
# population size for second round of infection
N2 <- 5

# graph parameters
cvs <- seq(0,3,0.02)
exp_frac_inf <- seq(0.02,0.98,0.02)

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0, length(cvs)*length(exp_frac_inf), 5)

# file to save data to
filename <- paste("ContPowersHiT_F",Focal,"_N",N,"_N2",N2,"_sims",num_sims,"_c03_e0.98_HiTm0.5phi2Lambseed.csv",sep="")

set.seed(3)

# function to numerically integrate to find Lamb
integ_over_lamb <- function(lamb, m, phi, k, theta, Lamb)
{
  (1 + lamb*Lamb*theta)^(-k) * (1 / (gamma(m) * phi^m)) * lamb^(m-1) * exp(-lamb/phi)
}

# function to find avg force of inf Lamb with root finding
find_Lamb_cont <- function(Lamb, m, phi, k, theta, f_inf)
{
  1 - f_inf - integrate(integ_over_lamb, lower=0, upper=Inf, m=m, phi=phi, k=k, theta=theta, Lamb=Lamb)$value
}


combo <- 1
for(cv in cvs)
{
  for(f_inf in exp_frac_inf)
  {
    # determine k and theta to use for gamma dist
    k <- 1 / (cv)^2
    theta <- ((1 - f_inf)^(-1/k) - 1)
        
    # initialize vectors to store likelihood tests later on
    ratios <- rep(0,num_sims)
    
    # vector to track focal indivs
    foc <- rep(0,num_sims)

    # number of infected individuals in 2nd exposure event
    inf2_prev <- matrix(NA, Focal, num_sims)
    inf2_naive <- matrix(NA, Focal, num_sims)

    # calculate avg force of infection Lamb
    if(cv == 0)
      Lamb <- ((1 - f_inf)^(-1/0.5) - 1) / (-log(1-f_inf)*2)
    else
    {
      if(f_inf == 0.98)
        Lamb <- uniroot(find_Lamb_cont, lower=0, upper=200, m=0.5, phi=2, k=k, theta=theta, f_inf=f_inf)$root
      else
        Lamb <- uniroot(find_Lamb_cont, lower=0, upper=200, m=0.5, phi=2, k=k, theta=theta, f_inf=f_inf)$root
    }

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

	# draw transmission/force of infection parameter
	lamb <- rgamma(1, shape=1/2, scale=2)

        # probability of being inf
        prob_I <- 1 - exp(-lamb*risk*Lamb)

        # infect individuals
	# use bernoulli for each indiv to check if infected, store risk for those not infected
        for(i in (1:(N)))
        {
          I <- rbinom(1, 1, prob_I[i])
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


	# draw transmission/force of infection parameter
	lamb2 <- rgamma(1, shape=1/2, scale=2)

        # probabilities of getting infected
        prob_I0 <- 1 - exp(-lamb2*risk0*Lamb)
        prob_I1 <- 1 - exp(-lamb2*risk1[f]*Lamb) # risk1[f] is focal indiv's risk
        
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
               

        # number indivs infected in 2nd round
        inf2_naive[f,sim] <- N2 - 1 - length(risk01)
        inf2_prev[f,sim] <- I1
        
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

