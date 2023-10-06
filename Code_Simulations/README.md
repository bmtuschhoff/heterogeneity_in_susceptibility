Below are brief descriptions of the code in each directory.

All code file names that start with "2gps" are used to simulate the discrete case and all file names that start with "Cont" are used to simulate the continuous case.
All code file names used for MCMC that contain the phrase "_cont" use previously generated contact network data to continue running MCMC.

Detection runs simulations for the power to detect heterogeneity in susceptibility using the likelihood equations 1 and 2 in the main text.

Estimation runs simulations to estimate the parameters for the underlying susceptibility distributions (p_A, p_B, f_A, k, theta) using Metropolis-Hastings MCMC (and ABC in the discrete case).

ChangeN runs simulations to estimate the parameters for the underlying susceptibility distributions (p_A, p_B, f_A, k, theta) using Metropolis-Hastings MCMC (and ABC in the discrete case) to check the effect of changing N. These code input contact networks with N=100 and subset them to N=5 to minimize the stochasticity from simulating data while determining the effect of N.
	Contact network data input into this code should be simulated in 2gps_GenData.R in the Estimation directory.

DiscretePriors runs simulations to estimate the parameters for the underlying susceptibility distribution (p_A, p_B, f_A) in the discrete case using Metropolis-Hastings MCMC and ABC with informative priors for p_A, p_B, or f_A.
	Contact network data input into this code should be simulated in 2gps_GenData.R in the Estimation directory.

FalseNegatives runs simulations for both detection and estimation of heterogeneity in susceptibility in the presence of false negatives.
	Code file names with the phrase "adj" use the method adjusted to account for false negatives as described in the supplementary information. Files without this phrase use the unmodified method.
	Code file names with the phrase "Det" run simulations for the power to detect heterogeneity in susceptibility.
	Code file names with the phrase "MCMC" run simulations to estimate the parameters for the underlying susceptibility distributions (p_A, p_B, f_A, k, theta) using Metropolis-Hastings MCMC (and ABC in the discrete case).
	The ".cpp" files are called to calculate likelihoods for the adjusted method within other ".R" files.

HiT runs simulations for both detection and estimation of heterogeneity in susceptibility in the presence of heterogeneity in transmission.
	Code file names with the phrase "DetHiS" run simulations for the power to detect heterogeneity in susceptibility.
	Code file names with the phrase "MCMC" run simulations to estimate the parameters for the underlying susceptibility distributions (p_A, p_B, f_A, k, theta) using Metropolis-Hastings MCMC (and ABC in the discrete case). Note that "2gps_HiT_GenData.R" generates contact networks with HiT, but does not run MCMC. The output data from this code should be input into 2gps_MCMC_cont.R in the Estimation directory to run MCMC.
	Code file names with the phrase "DetHiT" run simulations for the power to detect heterogeneity in transmission using the simulation-based goodness of fit test described in the supplementary information.
