# Analyze results of MCMC for estimating params of continuous case (k, theta)
# Run SIR dynamics with estimated params and calculate 95% CIs

library(coda)
library(mcmcplots)

# read in MCMC results
fileName <- "ContMCFNeg12adj_k1_theta3_F1000_N5_N25_len6e+05_sims100_fneg0.1Cseed.csv"
fileName <- "ContMCHiT_k0.591715976331361_theta0.626097060979711_F5000_N5_N25_len6e+05_sims1000_m0.5phi2seed.csv"
fileName <- "ContMCFNeg12_k0.591715976331361_theta0.626097060979711_FI1_FI21_F1000_N5_N25_len6e+05_sims600_fneg0.2.csv"
fileName <- "ContMCContdata_k0.591715976331361_theta0.626097060979711_F5000_N5_N25_len6e+05_sims600N100subseed.csv"
fileName <- "ContMC2gpsData_pA0.4999258614314_pB7.41385686000573e-05_fA0.5_F5000_N5_N25_len6e+05_sims600seed.csv"
results <- read.csv(fileName, header=T)
results <- results[,-1]
results <- results[-seq(16,21328,16),]
colnames(results) <- c("k","theta","L")
results$k <- as.numeric(results$k)
results$theta <- as.numeric(results$theta)
results$L <- as.numeric(results$L)
results <- na.omit(results)

# find burn-in period
rmeanplot(results)
BurnIn <- 2000

# prepare data
results <- results[(BurnIn+1):nrow(results),]
resultsMCMC <- mcmc(results)

# plot and look at results to check for convergence
plot(resultsMCMC)


# save results data frame to more easily access later
resultsc5000N5c1.3e0.25fA0.5wd <- results
resultsc1000N5c1e0.75fN0.1adjrd <- results


############################ SIR Dynamics ###############################
# run SIR model with different levels of heterogeneity in susceptibility

# function to run SIR model and return data frame of number of individuals in each class over time
SIR <- function(N,I_init,beta,gamma,mu,C,EndTime,deltat)
{
  S <- numeric(EndTime/deltat +1)
  I <- numeric(EndTime/deltat +1)
  R <- numeric(EndTime/deltat +1)
  
  S[1] <- N - I_init
  I[1] <- I_init
  R[1] <- 0
  
  t <- seq(0, EndTime, deltat)
  
  for(i in 1:(EndTime/deltat))
  {
    S[i+1] <- S[i] + (-beta*S[i]*I[i] *(S[i]/S[1])^(C^2) - mu*S[i])*deltat
    I[i+1] <- I[i] + (beta*S[i]*I[i] *(S[i]/S[1])^(C^2) - gamma*I[i] - mu*I[i])*deltat
    R[i+1] <- R[i] + (gamma*I[i] - mu*R[i])*deltat
  }
  
  return(data.frame(S,I,R))
}

# define parameters for SIR dynamics
pop_size <- 20010	# initial population size (S+I)
I_init <- 10		# initial number infected
R0 <- 3.0
gam <- 0.1		# recovery rate
mu <- 0			# death rate
EndTime <- 150		# total time
deltat <- 0.001		# time step size

# true parameters dictating HiS distribution
k <- 0.591715976331361
theta <- 0.626097060979711
k <- 1
theta <- 3

# average probability of infection
p <- 1 - (1 + theta)^(-k)

# contact rate
c <- (gam * R0) / (p * (pop_size-I_init))


# run SIR dynamics for true parameters with HiS
actual_results <- SIR(pop_size, I_init, bet=c*p, gam, mu, C=1/sqrt(k), EndTime, deltat)

# run SIR dynamics for true parameters with no HiS (homogeneity)
nohet_results <- SIR(pop_size, I_init, bet=c*p, gam, mu, C=0, EndTime, deltat)


## run SIR dynamics with HiS for estimated parameters from MCMC
# sample 1000 parameter sets from MCMC
set.seed(3)
params <- sample((1:nrow(resultsMCMC)), 1000, replace=F)

# run SIR model for each param set
# initialize matrix to hold all samples
S_samples <- matrix(nrow=EndTime/deltat + 1, ncol=1000)
I_samples <- matrix(nrow=EndTime/deltat + 1, ncol=1000)
n <- 1
for(param in params)
{
  k <- results$k[param]
  theta <- results$theta[param]
  p <- 1 - (1 +theta)^(-k)	# avg probability of infection
  c <- (gam * R0) / (p * (pop_size-I_init))  # contact rate
  
  # save results for S and I separately
  res <- SIR(pop_size, I_init, bet=c*p, gam, mu, C=1/sqrt(k), EndTime, deltat)
  S_samples[,n] <- (res$S) / (res$S[1])
  I_samples[,n] <- (res$I) / pop_size
  n <- n + 1
}

## find 95% CIs for SIR dynamics at each time point
# initialize vectors to hold quantiles and medians for S and I indivs
qs_2.5I <- c()
qs_97.5I <- c()
medMCMCI <- c()
qs_2.5S <- c()
qs_97.5S <- c()
medMCMCS <- c()

# at each time point, find 2.5% and 97.5% quantiles of S and I over all 1000 param sets
for(i in 1:(EndTime/deltat + 1))
{
  qs_2.5I[i] <- quantile(I_samples[i,], probs=0.025)
  qs_97.5I[i] <- quantile(I_samples[i,], probs=0.975)
  medMCMCI[i] <- median(I_samples[i,])
  
  qs_2.5S[i] <- quantile(S_samples[i,], probs=0.025)
  qs_97.5S[i] <- quantile(S_samples[i,], probs=0.975)
  medMCMCS[i] <- median(S_samples[i,])
}

# save 95% CIs
qs_2.5Ic5000N5c1.3e0.25fA0.5wd <- qs_2.5I
qs_97.5Ic5000N5c1.3e0.25fA0.5wd <- qs_97.5I
medMCMCIc5000N5c1.3e0.25fA0.5wd <- medMCMCI
qs_2.5Sc5000N5c1.3e0.25fA0.5wd <- qs_2.5S
qs_97.5Sc5000N5c1.3e0.25fA0.5wd <- qs_97.5S
medMCMCSc5000N5c1.3e0.25fA0.5wd <- medMCMCS

