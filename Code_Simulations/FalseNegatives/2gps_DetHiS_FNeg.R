# Calculate power to detect heterogeneity in suscepitibility given contact (HiS) in presence of false negatives in the discrete case (with 2 types) across parameter combinations using simulated contact tracing data and a likelihood ratio test
# Test if can tell difference between naive individuals and previously exposed (focal) individuals
# If there is a difference, there is HiS in the population.
# If there is not a difference, there is not HiS.

arg <- commandArgs(trailingOnly=T)

# number of focal individuals
Focal <- as.numeric(arg[1])

# fraction of individuals that are type A in original, naive individual
f_a_old <- 0.5

# rate of false negatives
false_neg <- as.numeric(arg[2])

# graph parameters
cvs <- seq(0,3,0.02)
exp_frac_inf <- seq(0.02,0.98,0.02)

# total initial population size
N <- 5
# population size for second exposure event
N2 <- 5

# number of simulations to run
num_sims <- 1000

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0,length(cvs)*length(exp_frac_inf),10)

# file to save data to
filename <- paste("2gpPowersFalseNeg12_fA",f_a_old,"_F",Focal,"_N",N,"_N2",N2,"_fneg",false_neg,"_sims",num_sims,"seed.csv",sep="")

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
  num_fneg <- rep(0, num_sims)	# track num false neg focal indivs
  
  f_a_new <- rep(NA, num_sims)	# fraction of uninfected pop that is type A after 1st exposure event
  f_b_new <- rep(NA, num_sims)	# fraction of uninfected pop that is type B after 1st exposure event
    
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

      # number of indivs that appear as false negs, assume 100% immunity
      fneg_sim <- rbinom(n=1, size=Ia+Ib, prob=false_neg)

      num_fneg[sim] <- num_fneg[sim] + fneg_sim
      f[sim] <- f[sim] + fneg_sim
    }

    # calculate new fractions of type a and b after 1st infection round
    f_a_new[sim] <- a_focal[sim] / (a_focal[sim] + b_focal[sim] + num_fneg[sim])
    f_b_new[sim] <- b_focal[sim] / (a_focal[sim] + b_focal[sim] + num_fneg[sim])
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
        # use unif RV to decide if choose type a or b or c, check for infection, record type
        rand_unif <- runif(1,0,1)
        
        if(rand_unif < f_a_new[j]) # type a
        {
          I[f,j] <- rbinom(1,1,p_a)
        } else if(rand_unif < f_a_new[j] + f_b_new[j]) # type b
        {
          I[f,j] <- rbinom(1,1,p_b)
        } else # type c, pC=0
        {
          I[f,j] <- 0
        }
      }
    }
    
  }
  

  # sum across all Focal groups for pop sizes and number infected individuals in 2nd round for naive vs prev exposed indiv
  naive_pop <- sum(Na_new) + sum(Nb_new) # integer
  naive_inf <- apply(Ia_new, 2, sum) + apply(Ib_new, 2, sum) # vector, length=num_sims
  prev_inf <- apply(I, 2, sum) # vector, length=num_sims
  
  
  # remove inf indivs b/c false neg tests
  naive_inf <- rbinom(n=num_sims, size=naive_inf, prob=1-false_neg)
  prev_inf <- rbinom(n=num_sims, size=prev_inf, prob=1-false_neg)
  
  
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

  # record powers in file
  write.table(powers_table[((i-49):(i-1)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}

