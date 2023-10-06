# Analyze results of MCMC for estimating params of discrete case (p_A, p_B, f_A)
# Run SIR dynamics with estimated params and calculate 95% CIs

library(coda)
library(mcmcplots)

# read in MCMC results
fileName <- "2gpsMCContdata_k0.591715976331361_theta0.626097060979711_F5000_N5_N25_len3e+07_sims100seed.csv"
fileName <- "2gpsMC2gpsData_pA0.4999258614314_pB7.41385686000573e-05_fA0.5_F5000_N5_N25_len3e+07_sims100t2.csv"
fileName <- "2gpsMC2gpsData_pA0.74809046024617_pB0.125477384938458_fA0.2_F5000_N5_N25_len3e+07_sims100_N100subseed.csv"
fileName <- "2gpsMC2gpsDataHiT_pA0.74809046024617_pB0.125477384938458_fA0.2_F5000_N5_N25_len3e+07_sims100_HiTm0.5phi2seed.csv"
fileName <- "2gpsMCFNeg12adj_pA0.74809046024617_pB0.125477384938458_fA0.2_F1000_N5_N25_fneg0.2_sims1000_ninf825_finf116_seed.csv"
results <- read.csv(fileName, header=T, row.names=1)
results <- read.csv(fileName, header=T)
results <- results[,-1]
results <- results[-seq(16,21328,16),]
colnames(results) <- c("p_a", "p_b", "f_a", "L")
results$p_a <- as.numeric(results$p_a)
results$p_b <- as.numeric(results$p_b)
results$f_a <- as.numeric(results$f_a)
results$L <- as.numeric(results$L)

# find burn-in period
rmeanplot(results)
BurnIn <- 10000

# prepare data
results <- results[(BurnIn+1):nrow(results),]
resultsMCMC <- mcmc(results)

# plot and look at results to check for convergence
plot(resultsMCMC)

# save results data frame to more easily access later
results1000N5c1.3e0.25fA0.2fN0.2adjrd <- results

############################ SIR Dynamics ###############################
# run SIR model with different levels of heterogeneity in susceptibility
# function to run SIR model and return data frame of number of individuals in each class over time
SIR <- function(N,fA,I_init_a,I_init_b,beta_a,beta_b,gamma,mu,EndTime,deltat)
{
  Sa <- numeric(EndTime/deltat +1)
  Sb <- numeric(EndTime/deltat +1)
  Ia <- numeric(EndTime/deltat +1)
  Ib <- numeric(EndTime/deltat +1)
  Ra <- numeric(EndTime/deltat +1)
  Rb <- numeric(EndTime/deltat +1)
  
  print(N)
  
  Sa[1] <- round((N - I_init_a - I_init_b)*fA)
  Sb[1] <- (N - Sa[1] - I_init_a - I_init_b)
  
  Ia[1] <- I_init_a
  Ib[1] <- I_init_b
  
  Ra[1] <- 0
  Rb[1] <- 0
  
  
  t <- seq(0, EndTime, deltat)
  
  for(i in 1:(EndTime/deltat))
  {
    Sa[i+1] <- Sa[i] + (-beta_a*Sa[i]*(Ia[i] + Ib[i]))*deltat
    Sb[i+1] <- Sb[i] + (-beta_b*Sb[i]*(Ib[i] + Ia[i]))*deltat
    
    Ia[i+1] <- Ia[i] + (beta_a*Sa[i]*(Ia[i] + Ib[i]) - gamma*Ia[i])*deltat
    Ib[i+1] <- Ib[i] + (beta_b*Sb[i]*(Ib[i] + Ia[i]) - gamma*Ib[i])*deltat
    
    Ra[i+1] <- Ra[i] + (gamma*Ia[i])*deltat
    Rb[i+1] <- Rb[i] + (gamma*Ib[i])*deltat
  }
  
  return(data.frame(Sa,Sb,Ia,Ib,Ra,Rb))
}

# define parameters for SIR dynamics
N <- 20010  # initial population size (S+I)
I_init <- 10  # initial number infected
R0 <- 3.0
gam <- 0.1  # recovery rate
mu <- 0  # death rate
EndTime <- 150  # total time
deltat <- 0.001  # time step size

# true parameters dictating HiS distribution
pA <- 0.4999258614314
pB <- 7.41385686000573e-05
fA <- 0.5

pA <- 0.74809046024617
pB <- 0.125477384938458
fA <- 0.2

# average probability of infection
p <- pA*fA + pB*(1 - fA)

# contact rate
c <- (gam * R0) / (p * (N-I_init))

# run SIR dynamics for true parameters with HiS
I_init_a <- round(I_init * pA/(pA + pB))  # initial number type A infected
I_init_b <- I_init - I_init_a  # initial number type B infected
actual_results <- SIR(N,fA,I_init_a,I_init_b,beta_a=c*pA,beta_b=c*pB,gam,mu,EndTime,deltat)


# run SIR dynamics for true parameters with no HiS (homogeneity)
pA <- pA*fA + pB*(1-fA)
pB <- pA
I_init_a <- round(I_init * pA/(pA + pB))
I_init_b <- I_init - I_init_a
nohet_results <- SIR(N,fA,I_init_a,I_init_b,beta_a=c*pA,beta_b=c*pB,gam,mu,EndTime,deltat)


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
  pA <- results$p_a[param]
  pB <- results$p_b[param]
  fA <- results$f_a[param]
  p <- pA*fA + pB*(1 - fA)	# avg probability of infection
  c <- (gam * R0) / (p * (N-I_init))	# contact rate
  I_init_a <- round(I_init * pA/(pA + pB))
  I_init_b <- I_init - I_init_a

  # save results for S and I separately
  res <- SIR(N,fA,I_init_a,I_init_b,beta_a=c*pA,beta_b=c*pB,gam,mu,EndTime,deltat)
  S_samples[,n] <- (res$Sa + res$Sb) / (res$Sa[1] + res$Sb[1])
  I_samples[,n] <- (res$Ia + res$Ib) / N
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
medMCMCI1000N5c1.3e0.25fA0.2fN0.2adjrd <- medMCMCI
qs_2.5I1000N5c1.3e0.25fA0.2fN0.2adjrd <- qs_2.5I
qs_97.5I1000N5c1.3e0.25fA0.2fN0.2adjrd <- qs_97.5I
medMCMCS1000N5c1.3e0.25fA0.2fN0.2adjrd <- medMCMCS
qs_2.5S1000N5c1.3e0.25fA0.2fN0.2adjrd <- qs_2.5S
qs_97.5S1000N5c1.3e0.25fA0.2fN0.2adjrd <- qs_97.5S
