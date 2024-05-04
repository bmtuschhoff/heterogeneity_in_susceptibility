# look at ability to detect HiS with "naive" individuals that have been previously exposed so not actually naive
# generate number of exposures for naive indivs and determine their risk from gamma dist, for focal indivs add 1 exposure to each of generated exposures
# number of previous exposures is Poisson i.e., exposed at a constant rate

install.packages("mvtnorm", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
library("mvtnorm", lib="../R_libs")


num_focal <- 1000
num_naive <- num_focal*4

# set level of HiS
cv <- 1.3
f_inf <- 0.25

# determine k and theta to use for gamma dist, theta0=0 previous exposures
k <- 1 / (cv)^2
theta0 <- ((1 - f_inf)^(-1/k) - 1)

# avg number of previous exposures for "naive" indivs
avg_num_exps <- 1


set.seed(3)

# naive individuals
num_exps_naive <- rpois(num_naive, lambda=avg_num_exps)  # number of previous exposures
risks_naive <- rgamma(num_naive, shape=k, scale=theta0/(1 + num_exps_naive*theta0))  # risk of infection for that indiv
  
# simulate if naive individuals infected, vector len=num_naive
inf_naive <- rbinom(num_naive, 1, 1-exp(-risks_naive))


# focal individuals
num_exps_focal <- rpois(num_focal, lambda=avg_num_exps) + 1  # number of previous exposures
risks_focal <- rgamma(num_focal, shape=k, scale=theta0/(1 + num_exps_focal*theta0))  # risk of infection
  
# simulate if focal individuals infected, vector len=num_focal
inf_focal <- rbinom(num_focal, 1, 1-exp(-risks_focal))


###############################################
# write numbers infected to files along with their risk of infection and number of previous exposures
naiveFile <- paste("infMultiExps_naive_k",k,"_theta",theta0,"_F",num_focal,"_Naive",num_naive,"seedt2.csv", sep="")
prevFile <- paste("infMultiExps_prev_k",k,"_theta",theta0,"_F",num_focal,"_Naive",num_naive,"seedt2.csv", sep="")

naive_info <- data.frame(num_exps=num_exps_naive, risks=risks_naive, inf=inf_naive)
focal_info <- data.frame(num_exps=num_exps_focal, risks=risks_focal, inf=inf_focal)

write.table(naive_info, file=naiveFile, sep=",")
write.table(focal_info, file=prevFile, sep=",")


