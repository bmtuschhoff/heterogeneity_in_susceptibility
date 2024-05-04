# Calculate power to detect heterogeneity in suscepitibility given contact (HiS) in the continuous case across parameter combinations using simulated contact tracing data and a likelihood ratio test
# Test if can tell difference between naive individuals and previously exposed (focal) individuals
# If there is a difference, there is HiS in the population.
# If there is not a difference, there is not HiS.
# Looking at ability to detect HiS with "naive" indivs that have been previously exposed so not actually naive
## generate number of exposures for naive indivs and determine their risk from gamma dist, for focal indivs add 1 exposure to each of generated exposures
## number of previous exposures is Poisson i.e., exposed at a constant rate

arg <- commandArgs(trailingOnly=T)

# number of focal individuals
num_focal <- as.numeric(arg[1])

# number of individuals per contact network
N <- 5

# number of naive individuals
num_naive <- num_focal * (N - 1)

# number of sims to run to determine power
num_sims <- 1000

# avg number of previous exposures for "naive" indivs
avg_num_exps <- 1

# graph parameters for level of HiS and expected fraction of truly naive individuals infecte
cvs <- seq(0, 3, 0.02)
exp_frac_inf <- seq(0.02, 0.98, 0.02)

# initialize table to save cv, f_inf, power
powers_table <- matrix(NA, nrow=length(cvs)*length(exp_frac_inf), ncol=5)

# file to save data to
filename <- paste("ContPowersMultiExps_F",num_focal,"_N",N,"_AvgPrevExp",avg_num_exps,"_sims",num_sims,"s3.csv",sep="")

set.seed(3)

## determine power for all param combos
combo <- 1  # counter variable for saving data in table
for(cv in cvs)
{
  for(f_inf in exp_frac_inf)
  {
    # determine k and theta to use for gamma dist, theta0=0 previous exposures
    k <- 1 / (cv)^2
    theta0 <- ((1 - f_inf)^(-1/k) - 1)

    # initialize vector to store likelihood tests later on
    ratios <- rep(0, num_sims)

    # number of total infected individuals in each simulation
    inf_naive <- rep(NA, num_sims)
    inf_focal <- rep(NA, num_sims)

    # determine numbers of individuals infected in each simulation
    for(sim in 1:num_sims)
    {
      ## naive individuals
      # determine number of previous exposures and risk of infection
      num_exps_naive <- rpois(num_naive, lambda=avg_num_exps)  # number of previous exposures
      risks_naive <- rgamma(num_naive, shape=k, scale=theta0/(1 + num_exps_naive*theta0))  # risk of infection for that indiv

      if(cv == 0)
	risks_naive <- rep(-log(1 - f_inf), num_naive)

      # number of naive indivs infected
      inf_naive[sim] <- sum(rbinom(num_naive, 1, 1-exp(-risks_naive)))


      ## focal individuals
      # determine number of previous exposures and risk of infection
      num_exps_focal <- rpois(num_focal, lambda=avg_num_exps) + 1  # number of previous exposures
      risks_focal <- rgamma(num_focal, shape=k, scale=theta0/(1 + num_exps_focal*theta0))  # risk of infection

      if(cv == 0)
	risks_focal <- rep(-log(1 - f_inf), num_focal)

      # number of focal indivs infected
      inf_focal[sim] <- sum(rbinom(num_focal, 1, 1-exp(-risks_focal)))
    }

    # calculate probabilities of infection, vector len=num_sims
    p_naive <- inf_naive / num_naive
    p_focal <- inf_focal / num_focal
    p_avg <- (inf_naive + inf_focal) / (num_naive + num_focal)

    # calculate log-likelihoods for homogeneity vs heterogeneity
    Lhom <- dbinom(inf_naive, num_naive, p_avg, log=T) + dbinom(inf_focal, num_focal, p_avg, log=T)
    Lhet <- dbinom(inf_naive, num_naive, p_naive, log=T) + dbinom(inf_focal, num_focal, p_focal, log=T)

    # calculate likelihood ratio test statistics
    ratios <- -2 * (Lhom - Lhet)

    # find power = percent of times the ratio test is greater than critical value 3.84 (df=1)
    power <- sum(ratios > 3.84, na.rm=T) / num_sims * 100

    # record cv, f_inf, power in table
    powers_table[combo,] <- c(cv, f_inf, power, mean(p_naive), mean(p_focal))
    combo <- combo + 1
  }

  # record powers in file
  write.table(powers_table[((combo-49):(combo-1)),], file=filename, append=T, sep=",", col.names=NA, row.names=T)
}

