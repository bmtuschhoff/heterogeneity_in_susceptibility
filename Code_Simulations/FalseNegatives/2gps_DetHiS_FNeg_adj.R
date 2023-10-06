# Calculate power to detect heterogeneity in suscepitibility given contact (HiS) in presence of false negatives in the discrete case (with 2 types) across parameter combinations using simulated contact tracing data and a likelihood ratio test
# CORRECTING FOR FALSE NEGATIVES
# Test if can tell difference between naive individuals and previously exposed (focal) individuals
# If there is a difference, there is HiS in the population.
# If there is not a difference, there is not HiS.


install.packages("Rcpp", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
install.packages("RcppGSL", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')

library("Rcpp", lib="../R_libs")
library("RcppGSL", lib="../R_libs")

.libPaths("../R_libs")


# number of focal individuals
Focal <- 200

# fraction of individuals that are type A in original, naive individual
f_a_old <- 0.5

# rate of false negatives
false_neg <- 0.1

# graph parameters
cvs <- seq(0,3,0.02)
exp_frac_inf <- seq(0.02,0.98,0.02)

# total initial population size
N <- 5
# population size for second exposure event
N2 <- 5

# number of simulations to run
num_sims <- 1000

# set C++ function as source code
sourceCpp("2gps_DetHiSL_FNeg.cpp")

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0,length(cvs)*length(exp_frac_inf),3)

# file to save data to
filename <- paste("2gpPowersFalseNeg12adj_fA",f_a_old,"_F",Focal,"_N",N,"_N2",N2,"_fneg",false_neg,"_sims",num_sims,"seed.csv",sep="")

set.seed(3)

# function to find r_b with root finding
find_r_b <- function(r_b, f_a, cv, f_inf)
{
  abs(((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))) - r_b)*sqrt(f_a*(1 - f_a))) / ((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))))*f_a + r_b*(1 - f_a)) - cv)
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
      powers_table[i,] <- c(f_inf, cv, NA)
      i <- i + 1
    }
    else
    {
      power <- find_power(p_a, p_b, f_a_old, false_neg)
      
      # record cv, f_inf, power in table
      powers_table[i,] <- c(f_inf, cv, power)
      i <- i + 1
    }
  }

  # record powers in file
  write.table(powers_table[((i-49):(i-1)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}
