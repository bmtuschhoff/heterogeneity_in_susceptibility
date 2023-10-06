# Calculate power to detect heterogeneity in suscepitibility given contact (HiS) in the presence of false negatives in the continuous case across parameter combinations using simulated contact tracing data and a likelihood ratio test
### CORRECTING FOR FALSE NEGATIVES
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

# number of sims to run to determine power
num_sims <- 1000

# total initial population size
N <- 5
# population size for second exposure event
N2 <- 5

# rate of false negatives
false_neg <- 0.1

# set C++ function as source code
sourceCpp("Cont_DetHiSL_FNeg.cpp")

# graph parameters
cvs <- seq(0,3,0.02)
exp_frac_inf <- seq(0.02,0.98,0.02)

# initialize table for saving power, cv, f_inf,...
powers_table <- matrix(0, length(cvs)*length(exp_frac_inf), 3)

# file to save data to
filename <- paste("ContPowersFNeg12adj_F",Focal,"_N",N,"_N2",N2,"_sims",num_sims,"_fneg",false_neg,"c3e0.98seed.csv",sep="")

set.seed(3)

combo <- 1
for(cv in cvs)
{
  for(f_inf in exp_frac_inf)
  {
    # determine k and theta to use for gamma dist
    k <- 1 / (cv)^2
    theta <- ((1 - f_inf)^(-1/k) - 1)

    # cv=0, flag k so know to change in C++ code
    if(k == Inf)
      k <- -1.0
    
    # determine power to detect HiS in C++ fntn
    power <- find_power(k, theta, f_inf, false_neg)
    
    # record power, cv, f_inf in table
    powers_table[combo,] <- c(cv, f_inf, power)
    combo <- combo + 1
  }

  # record powers in file
  write.table(powers_table[((combo-49):(combo-1)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}

