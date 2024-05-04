# Generate contact tracing data to estimate parameters for heterogeneity in suscepitibility given contact (HiS) in the presence of heterogeneity in transmission in the continuous case using Metropolis-Hastings MCMC
# Estimating the parameters that determine the distribution of individuals' risks of infection (k, theta) and forces of infection (m)

arg <- commandArgs(trailingOnly=T)

# number of focal individuals
#Focal <- as.numeric(arg[1])
Focal <- 1000

# number of sims to run
num_sims <- 1

# total initial population size
N <- 5
# population size for second exposure event
N2 <- 5

# cv of risk and expected fraction of naive individuals infected
cvs <- 1.3
exp_frac_inf <- 0.25

# level of heterogeneity in transmission
m <- 0.5

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0, length(cvs)*length(exp_frac_inf), 5)

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
    ratios <- rep(0, num_sims)

    # vector to track focal indivs
    foc <- rep(0, num_sims)
    
    # number of infected individuals 1st and 2nd rounds
    inf1 <- matrix(NA, Focal, num_sims)
    inf2_prev <- matrix(NA, Focal, num_sims)
    inf2_naive <- matrix(NA, Focal, num_sims)

    # calculate avg force of infection Lamb
    Lamb <- uniroot(find_Lamb_cont, lower=0, upper=300, m=m, phi=1/m, k=k, theta=theta, f_inf=f_inf)$root
    
    for(sim in (1:num_sims))
    {  
      risk1 <- c()	# initialize vector for risks of focal indivs
      j <- 1	# initialize counter var to save risks
      f <- 1	# initialize counter var to save number inf in each contact network

      # collect focal indivs - simulate contact networks until get enough focal indivs (not inf)
      while(foc[sim] < Focal)
      {
        ## first exposure events
        # set up relative risk distributions for first infection rounds
        risk <- rgamma(N, shape=k, scale=theta)

	# draw transmission/force of infection parameter
	lamb <- rgamma(1, shape=m, scale=1/m)

        # probability of being inf
        prob_I <- 1 - exp(-lamb*risk*Lamb)

        # infect individuals
	num_notI1 <- 0		# number of indivs not inf in this contact network
	# use bernoulli for each indiv to check if infected, store risk for those not infected
        for(i in (1:(N)))
        {
          I <- rbinom(1, 1, prob_I[i])
          if(I == 0)  # indiv not infected
          {
            risk1[j] <- risk[i]
            j <- j + 1
	    foc[sim] <- foc[sim] + 1
	    num_notI1 <- num_notI1 + 1
          }
        }
        # have repeat of contact network if multiple indivs not inf i.e. multiple indivs will become focal indivs in next round so can keep track of which network goes with which focal indiv...
	for(x in 1:num_notI1)
	{
	  if(f <= Focal)
	  {
	  	inf1[f,sim] <- N - num_notI1
	  	f <- f + 1
	  }
	}
      }


      for(f in (1:Focal))
      {      
        # create risk and prob of not infected distribution for naive individuals
        risk0 <- rgamma(N2-1, shape=k, scale=theta)

	# draw transmission/force of infection parameter
	lamb2 <- rgamma(1, shape=m, scale=1/m)
        
        # probabilities of getting infected
        prob_I0 <- 1 - exp(-lamb2*risk0*Lamb)
        prob_I1 <- 1 - exp(-lamb2*risk1[f]*Lamb) # risk1[f] is focal indiv's risk
        
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
	inf2_naive[f,sim] <- N2 - 1 - length(risk01)
	inf2_prev[f,sim] <- I1
        ##########################################################################################################
      }      
    }

    # sum across all Focal groups for number individuals infected in 2nd round for naive vs focal indivs, vector len=sims
    naive_inf <- apply(inf2_naive, 2, sum)
    prev_inf <- apply(inf2_prev, 2, sum)

    # calculate observed probabilities of inf across all focal indivs as best estimates
    p_naive <- naive_inf / (Focal*(N2 - 1))
    p_prev <- prev_inf / Focal
    p_avg <- (naive_inf + prev_inf) / (Focal*N2)

    # calculate log-likelihoods for homogeneity vs heterogeneity
    Lhom <- dbinom(naive_inf, Focal*(N2-1), p_avg, log=T) + dbinom(prev_inf, Focal, p_avg, log=T)
    Lhet <- dbinom(naive_inf, Focal*(N2-1), p_naive, log=T) + dbinom(prev_inf, Focal, p_prev, log=T)

    # calculate likelihood ratio test statistics
    ratios <- -2*(Lhom - Lhet)
    
    # find power = percent of times the ratio test is greater than critical value 3.84 (df=1)
    power <- sum(ratios > 3.84, na.rm=T)/length(ratios) * 100

    # record power, cv, f_inf in table
    powers_table[combo,] <- c(cv, f_inf, power, mean(p_naive), mean(p_prev))
    combo <- combo + 1
  }
}

powers_table <- as.data.frame(powers_table)
colnames(powers_table) <- c("cv", "f_inf", "power", "p_naive", "p_prev")

### write contact networks (num inf) to files to run MCMC chain
sims <- 600 # param in MCMC
steps <- 6e5 # param in MCMC
inf1File <- paste("inf1_HiT_k",k,"_theta",theta,"_m",m,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"seed3.csv", sep="")
naiveFile <- paste("inf2_naiveHiT_k",k,"_theta",theta,"_m",m,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"seed3.csv", sep="")
prevFile <- paste("inf2_prevHiT_k",k,"_theta",theta,"_m",m,"_F",Focal,"_N",N,"_N2",N2,"_len",steps,"_sims",sims,"seed3.csv", sep="")
write.table(inf1, file=inf1File, sep=",")
write.table(inf2_naive, file=naiveFile, sep=",")
write.table(inf2_prev, file=prevFile, sep=",")


