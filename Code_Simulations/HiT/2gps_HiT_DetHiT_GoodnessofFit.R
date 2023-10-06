# Calculate power to detect heterogeneity in transmission given contact (HiT) in the discrete case (with 2 types) across different levels of HiT in the presence of heterogeneity in susceptibility using simulated contact tracing data
# Test if number of naive individuals infected in each contact network follows a binomial distribution
# If so, there is not evidence for HiT.
# If not, there is evidence for HiT.


arg <- commandArgs(trailingOnly=T)

# number of focal individuals
Focal <- as.numeric(arg[1])

# set cv, exp_frac_inf, f_a for HiS
cv <- 1.3
f_inf <- 0.25
f_a_old <- 0.2

# values for shape of HiT gamma distribution to check
ms <- seq(0.1, 6, 0.1)

# total initial population size
N <- 5
N2 <- 5

# number of simulations to run
num_sims <- 10000

# initialize table to record results
binom_dist_table <- matrix(0,length(ms),5)

# function to find r_b with root finding
find_r_b <- function(r_b, f_a, cv, f_inf)
{
  abs(((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))) - r_b)*sqrt(f_a*(1 - f_a))) / ((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))))*f_a + r_b*(1 - f_a)) - cv)
}

# function to find avg force of inf Lamb with root finding
find_Lamb_disc <- function(Lamb, m, phi, f_A, r_A, r_B, f_inf)
{
  (1 - (1 + r_A*Lamb*phi)^(-m))*f_A + (1 - (1 + r_B*Lamb*phi)^(-m))*(1 - f_A) - f_inf
}


set.seed(3)

# function to find power to detect heterogeneity in transmission
find_power <- function(ra, rb, f_inf, m)
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
  

  # find avg force of infection
  Lamb <- uniroot(find_Lamb_disc, lower=0, upper=2000, m=m, phi=1/m, f_A=f_a_old, r_A=ra, r_B=rb, f_inf=f_inf)$root
  
  # collect focal indivs - simulate contact networks until get enough focal indivs (not inf)
  for(sim in 1:num_sims)
  {
    while(f[sim] < Focal)
    {
      ## first exposure events
      # draw transmission/force of infection parameter
      lamb <- rgamma(1, shape=m, scale=1/m)
      # calculate transmission * risk
      tr_a <- lamb * ra * Lamb
      tr_b <- lamb * rb * Lamb
      # calculate probs of being infected
      p_a <- 1 - exp(-tr_a)
      p_b <- 1 - exp(-tr_b)
      
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
  
  
  for(f in (1:Focal))
  {
    # second exposure events
    # draw transmission/force of infection parameter
    lamb2 <- rgamma(1, shape=m, scale=1/m)
    # calculate transmission * risk
    tr_a2 <- lamb2 * ra * Lamb
    tr_b2 <- lamb2 * rb * Lamb
    # calculate probs of being infected
    p_a2 <- 1 - exp(-tr_a2)
    p_b2 <- 1 - exp(-tr_b2)
    
    
    # number naive individuals of each type, add up to N2-1
    ### use negative binomial dist for number of contacts b/c most indivs have few contacts but some have many
    ### number of contacts typically modeled with neg binom!
    ### set mean (mu)=4, size=shape param for gamma dist so low values=lots of variability in number contacts
    ### add 1 to avoid N=0 contacts
    #N2 <- rnbinom(1, size=1, mu=4) + 2
    Na_new[f] <- rbinom(1, N2-1, f_a_old)
    Nb_new[f] <- N2 - 1 - Na_new[f]
    
    # infection of naive individuals
    Ia_new[f,] <- rbinom(n=num_sims, size=Na_new[f], prob=p_a2)
    Ib_new[f,] <- rbinom(n=num_sims, size=Nb_new[f], prob=p_b2)
    
    
    # infection of prev exposed indiv
    for(j in (1:num_sims))
    {
      if(!is.nan(f_a_new[j]))
      {
        # use unif RV to decide if choose type a or b, check for infection, record type
        if(runif(1,0,1) < f_a_new[j])
        {
          I[f,j] <- rbinom(1,1,p_a2)
        }
        else
        {
          I[f,j] <- rbinom(1,1,p_b2)
        }
      }
    }
    
  }
  
  
  # use naive indivs to decide if there is het in trans
  # simulate log-likelihoods assuming no HiT and check if data likelihood is within potential likelihoods
  
  ## calculate prob of being infected for naive indivs
  # total number of naive indivs, integer
  naive_pop <- sum(Na_new + Nb_new)
  # number of naive infected individuals in each focal group, matrix Focal x num_sims
  naive_inf <- Ia_new + Ib_new
  # prob of inf for each focal group, matrix Focal x num_sims
  p_naive <- apply(naive_inf, 2, sum) / naive_pop # vector, length=num_sims
  
  # Na_new + Nb_new = number naive indivs in each group, same for each sim, vector length=Focal
  N_new <- Na_new + Nb_new
  
  
  ### method to detect HiT
  not_binom <- rep(TRUE, num_sims) # initalize vector detailing whether data shows HiT (T) or not (F) and therefore is binom dist
  for(sim in 1:num_sims)
  {
    potential_Ls <- rep(0,1000) # initialize vector of likelihoods calculated from resimulating data with p_naive
    
    for(lhood in 1:1000)
    {
      # simulate one time number of naive indivs infected in each of Focal networks with prob of inf given by p_naive
      resim_data <- rbinom(Focal, N_new, p_naive[sim]) # vector of # inf in each network, length=Focal
      
      # calculate log-likelihood of re-simulated data, float
      potential_Ls[lhood] <- sum(dbinom(resim_data, N_new, p_naive[sim], log=T))
    }
    
    # likelihood of actual data seen with p_naive
    data_L <- sum(dbinom(naive_inf[,sim], N_new, p_naive[sim], log=T))
    
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
  
  
  return(c(m, f_inf, cv, not_binom, mean(p_naive)))
}



i <- 1

for(m in ms)
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
    
    # record NAs if error in finding ra, rb
    if(p_a < 0 || p_a > 1 || p_b < 0 || p_b > 1)
    {
      binom_dist_table[i,] <- c(m, f_inf, cv, NA, NA)
      i <- i + 1
    }
    else
    {
      results <- find_power(ra, rb, f_inf, m)
      
      # table of probability combos, power, and avg cv (across simulations) to plot
      binom_dist_table[i,] <- results[1:5]
      i <- i + 1
    }
}
#########################################################################
binom_dist_table <- as.data.frame(binom_dist_table)
colnames(binom_dist_table) <- c("m", "f_inf", "cv", "not_binom_dist", "p_naive")

# record in file
filename <- paste("2gpPowersDetHiT_fA",f_a_old,"_F",Focal,"_N",N,"_N2",N2,"_sims",num_sims,"c1.3e0.25m6seed.csv",sep="")
write.csv(binom_dist_table, file=filename)
