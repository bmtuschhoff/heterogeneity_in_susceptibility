# Calculate power to detect heterogeneity in transmission given contact (HiT) in the continuous case across different levels of HiT in the presence of heterogeneity in susceptibility using simulated contact tracing data
# Test if number of naive individuals infected in each contact network follows a binomial distribution
# If so, there is not evidence for HiT.
# If not, there is evidence for HiT.

arg <- commandArgs(trailingOnly=T)

# number of focal individuals
Focal <- as.numeric(arg[1])

# number of sims to run
num_sims <- 1000

# total initial population size
N <- 10
# population size for second exposure event
N2 <- 10

# set cv, exp_frac_inf for HiS
cv <- 1.3
f_inf <- 0.25

# values for shape of HiT gamma distribution to check
ms <- seq(0.1, 6, 0.1)

# initialize table to record results
binom_dist_table <- matrix(0, length(ms), 5)

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


set.seed(3)
combo <- 1
for(m in ms)
{
    # determine k and theta to use for gamma dist
    k <- 1 / (cv)^2
    theta <- ((1 - f_inf)^(-1/k) - 1)
      
    # initialize vectors to store likelihood tests later on
    ratios <- rep(0,num_sims)
    
    # vector to track focal indivs
    foc <- rep(0,num_sims)

    # number of infected individuals 1st and 2nd rounds
    inf2_prev <- matrix(NA, Focal, num_sims)
    inf2_naive <- matrix(NA, Focal, num_sims)

    # calculate avg force of infection Lamb
    Lamb <- uniroot(find_Lamb_cont, lower=0, upper=200, m=m, phi=1/m, k=k, theta=theta, f_inf=f_inf)$root
    
    for(sim in (1:num_sims))
    {
      risk1 <- c()	# initialize vector for risks of focal indivs
      j <- 1	# initialize counter var to save risks

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
        risk0 <- rgamma(N2-1, shape=k, scale=theta)

	# draw transmission/force of infection parameter
	lamb2 <- rgamma(1, shape=m, scale=1/m)

        # probabilities of getting infected
        prob_I0 <- 1 - exp(-lamb2*risk0*Lamb)
        prob_I1 <- 1 - exp(-lamb2*risk1[f]*Lamb) # focal indiv
        
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
               

        # number indivs infected in 2nd round
        inf2_naive[f,sim] <- N2 - 1 - length(risk01)
        inf2_prev[f,sim] <- I1
        
        ##########################################################################################################
      }
    }

  # use naive indivs to decide if there is het in transmission
  # simulate for log-likelihoods assuming no HiT and check if data likelihood is within potential likelihoods
  
  ## calculate prob of being infected for naive indivs
  # number of naive infected individuals in each focal group is inf2_naive, matrix Focal x num_sims
  # sum across Focal indiv groups for number individuals obs inf in 2nd round for naive indivs, vector len=num_sims
  naive_inf <- apply(inf2_naive, 2, sum)
   
  # calculate observed probability of inf across all focal indivs as best estimate, vector length=num_sims
  p_naive <- naive_inf / (Focal*(N2 - 1))

  # number naive indivs in each group is N2-1


  ### method to detect HiT
  not_binom <- rep(TRUE, num_sims) # initalize vector detailing whether data shows HiT (T) or not (F) and therefore is binom dist
  for(sim in 1:num_sims)
  {
    potential_Ls <- rep(0,10000) # initialize vector of likelihoods calculated from resimulating data with p_naive
    
    for(lhood in 1:10000)
    {
      # simulate one time number of naive indivs infected in each of Focal networks with prob of inf given by p_naive
      resim_data <- rbinom(Focal, N2-1, p_naive[sim]) # vector of # inf in each network, length=Focal
      
      # calculate log-likelihood of re-simulated data, float
      potential_Ls[lhood] <- sum(dbinom(resim_data, N2-1, p_naive[sim], log=T))
    }
    
    # likelihood of actual data seen with p_naive
    data_L <- sum(dbinom(inf2_naive[,sim], N2-1, p_naive[sim], log=T))
    
    # 95% CI for potential likelihoods
    quantLs <- quantile(potential_Ls, probs=c(0.025,0.975))
    
    # check if likelihood of data fits in 95% CI
    # if data likelihood is in CI, no HiT so is binom dist
    # if data likelihood is not in CI, HiT so not binom dist, leave as True
    if(data_L >= quantLs[1] && data_L <= quantLs[2])
      not_binom[sim] <- FALSE
  }
  
   
  # percent of simulations in which data shows HiT so is not binom dist
  # power to detect HiT
  not_binom <- 100 * sum(not_binom) / num_sims

  binom_dist_table[combo,] <- c(m, f_inf, cv, not_binom, mean(p_naive))
  combo <- combo + 1
}

binom_dist_table <- as.data.frame(binom_dist_table)
colnames(binom_dist_table) <- c("m", "f_inf", "cv", "not_binom_dist", "p_naive")

# record in file
filename <- paste("ContPowersDetHiT_F",Focal,"_N",N,"_N2",N2,"_sims",num_sims,"_cv1.3_e0.25m6seed.csv",sep="")
write.csv(binom_dist_table, file=filename)