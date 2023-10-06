# Generate contact tracing data to estimate parameters for heterogeneity in suscepitibility given contact (HiS) in the discrete case (with 2 types) using Metropolis-Hastings MCMC with ABC
# Estimating the parameters that determine the distribution of individuals' risks of infection (p_A, p_B, f_A)

arg <- commandArgs(trailingOnly=T)

# number of focal individuals
Focal <- as.numeric(arg[1])

# fraction of individuals that are type A in original, naive individual
f_a_old <- 0.2

# cv of risk and expected fraction of naive individuals infected
cvs <- 1.3
exp_frac_inf <- 0.25

# total initial population size
N <- as.numeric(arg[2])
# population size for second exposure event
N2 <- as.numeric(arg[2])

# number of simulations to run
num_sims <- 1

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0,length(cvs)*length(exp_frac_inf),10)

set.seed(3)

# function to find r_b with root finding
find_r_b <- function(r_b, f_a, cv, f_inf)
{
  abs(((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))) - r_b)*sqrt(f_a*(1 - f_a))) / ((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))))*f_a + r_b*(1 - f_a)) - cv)
}

# function to find power to detect het in susc
find_power <- function(p_a, p_b)
{
  # initialize vectors and matrices to hold data
  f <- rep(0, num_sims)		# track num focal indivs
  a_focal <- rep(0, num_sims)	# track num type A focal indivs
  b_focal <- rep(0, num_sims)	# track num type B focal indivs
  
  f_a_new <- rep(NA, num_sims)	# fraction of uninfected pop that is type A after 1st exposure event
  
  Na_new <- rep(0, Focal)	# num type A indivs per network in 2nd exposure event
  Nb_new <- rep(0, Focal)	# num type B indivs per network in 2nd exposure event
  
  Ia_new <- matrix(0, Focal, num_sims)	# num type A naive indivs infected in 2nd exposure event
  Ib_new <- matrix(0, Focal, num_sims)	# num type B naive indivs infected in 2nd exposure event
  I <- matrix(0, Focal, num_sims)	# num focal indivs infected in 2nd exposure event
 

  # collect focal indivs - simulate contact networks until get enough focal indivs (not inf)
  for(sim in 1:num_sims)
  {
    while(f[sim] < Focal)
    {
      ## first exposure events
      # number individuals of each type, add up to total pop size N, based on fraction that are type a
      Na <- rbinom(1, N, f_a_old)
      Nb <- N - Na

      # infect individuals
      Ia <- rbinom(n=1, size=Na, prob=p_a)
      Ib <- rbinom(n=1, size=Nb, prob=p_b)

      # check if any indivs not inf -- if not, save as focal indivs
      if(Na - Ia > 0)
      {
        a_focal[sim] <- a_focal[sim] + Na-Ia
        f[sim] <- f[sim] + Na-Ia
      }
      if(Nb - Ib > 0)
      {
        b_focal[sim] <- b_focal[sim] + Nb-Ib
        f[sim] <- f[sim] + Nb-Ib
      }
    }

    # calculate new fraction of type a after 1st infection round
    f_a_new[sim] <- a_focal[sim] / (a_focal[sim] + b_focal[sim])
  }

 
  ## second exposure events
  for(f in (1:Focal))
  {
    # number naive individuals of each type, add up to N2-1
    Na_new[f] <- rbinom(1, N2-1, f_a_old)
    Nb_new[f] <- N2 - 1 - Na_new[f]
    
    # infection of naive individuals
    Ia_new[f,] <- rbinom(n=num_sims, size=Na_new[f], prob=p_a)
    Ib_new[f,] <- rbinom(n=num_sims, size=Nb_new[f], prob=p_b)
    
    # infection of prev exposed indiv
    for(j in (1:num_sims))
    {
      if(!is.nan(f_a_new[j]))
      {
        # use unif RV to decide if choose type a or b, check for infection, record type
        if(runif(1,0,1) < f_a_new[j])
        {
          I[f,j] <- rbinom(1,1,p_a)
        }
        else
        {
          I[f,j] <- rbinom(1,1,p_b)
        }
      }
    }
    
  }
  
  # sum across all Focal groups for pop sizes and number infected individuals in 2nd round for naive vs prev exposed indiv
  naive_pop <- sum(Na_new) + sum(Nb_new) # integer
  naive_inf <- apply(Ia_new, 2, sum) + apply(Ib_new, 2, sum) # vector, length=num_sims
  prev_inf <- apply(I, 2, sum) # vector, length=num_sims


  ##############################################
  ### write data to files for running MCMC chain
  naiveFile <- paste("inf2_naive_2gps_pA",p_a,"_pB",p_b,"_fA",f_a_old,"_F",Focal,"_N",N,"_N2",N2,"seed.csv", sep="")
  prevFile <- paste("inf2_prev_2gps_pA",p_a,"_pB",p_b,"_fA",f_a_old,"_F",Focal,"_N",N,"_N2",N2,"seed.csv", sep="")
  write.table(Ia_new+Ib_new, file=naiveFile, sep=",")
  write.table(I, file=prevFile, sep=",")
  ##############################################

  
  # calculate probabilities across all focal individuals as best estimates
  p_naive <- naive_inf / naive_pop
  p_prev <- prev_inf / (Focal - rem_focals)
  p_avg <- (naive_inf + prev_inf) / (naive_pop + (Focal - rem_focals))
  
  # calculate log-likelihoods for homogeneity vs heterogeneity
  Lhom <- dbinom(naive_inf, naive_pop, p_avg, log=T) + dbinom(prev_inf, Focal-rem_focals, p_avg, log=T)
  Lhet <- dbinom(naive_inf, naive_pop, p_naive, log=T) + dbinom(prev_inf, Focal-rem_focals, p_prev, log=T)
  
  # calculate likelihood ratio test statistics
  ratio_tests <- -2*(Lhom - Lhet)
  
  # find power = percent of times ratio test is greater than critical value 3.84 (df=1)
  power <- sum(ratio_tests > 3.84, na.rm=T)/length(ratio_tests)*100
  
  # return power, cv, f_inf, sds, means
  return(c(f_inf, cv, power, sd(ratio_tests)/mean(ratio_tests), mean(p_avg), sd(p_avg), mean(p_naive), sd(p_naive), mean(p_prev), sd(p_prev), naive_inf, prev_inf))
}





i <- 1

for(cv in cvs)
{
  for(f_inf in exp_frac_inf)
  {
    if(cv == 0)		# no HiS, all indivs have same risk based on f_inf
    {
      p_a <- f_inf
      p_b <- f_inf
      ra <- -log(1 - p_a)
      rb <- -log(1 - p_b)
    } else
    {
      # find risk for type B indivs numerically solving for root
      if(f_inf < f_a_old)
      {
        if(f_inf < 1-f_a_old)
        {
          rb <- optimize(find_r_b, interval=c(0,-log(1 - (f_inf / (1 - f_a_old)))), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        } else
        {
          rb <- optimize(find_r_b, interval=c(0,50), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        }
      } else
      {
        if(f_inf < 1-f_a_old)
        {
          rb <- optimize(find_r_b, interval=c(-log(1 - ((f_inf - f_a_old) / (1 - f_a_old))),-log(1 - (f_inf / (1 - f_a_old)))), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        } else
        {
          rb <- optimize(find_r_b, interval=c(-log(1 - ((f_inf - f_a_old) / (1 - f_a_old))),50), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        }
      }
      # calculate risk for type A indivs
      ra <- -log(1 - (1/f_a_old)*(f_inf - (1 - exp(-rb))*(1 - f_a_old)))

      # back calculate for p_a and p_b
      p_a <- 1 - exp(-ra)
      p_b <- 1 - exp(-rb)
    }
    
    # record NAs if error in finding ra, rb
    if(p_a < 0 || p_a > 1 || p_b < 0 || p_b > 1)
    {
      powers_table[i,] <- c(f_inf, cv, NA, NA, NA, NA, NA, NA, NA, NA)
      i <- i + 1
    }
    else
    {
      results <- find_power(p_a, p_b)
      naive_inf <- results[11]
      prev_inf <- results[12]
      
      # record cv, f_inf, power, sds, means in table
      powers_table[i,] <- results[1:10]
      i <- i + 1
    }
  }
}

