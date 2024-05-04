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
	Code file names with the phrase "GenData" generate contact networks with HiT. The output data from this code should be input into 2gps_MCMC_cont.R or Cont_MCMC_cont.R in the Estimation directory or Cont_HiT_MCMC_AccountForHiT.R in this directory to run MCMC.
	Code file names with the phrase "DetHiT" run simulations for the power to detect heterogeneity in transmission using the simulation-based goodness of fit test described in the supplementary information.
	The file Cont_HiT_MCMC_AccountForHiT.R runs MCMC to estimate the parameters for the underlying susceptibility distribution (k, theta) and forces of infection distribution (m).

DynamicEpidemic runs simulations for both detection and estimation of heterogeneity in susceptibility with the contact tracing data (i.e., contact networks) generated from a dynamic epidemic.
	The ".cpp" files simulate a stochastic, individual-based SIR model via the Gillespie algorithm. The output of these files is a file containing the contact network for each infected individual throughout the simulation, which is written as the number of naive and focal individuals exposed and infected by this individual. For the discrete case, the 2gps_FindrB.R file must be first used to calculate the risk of infection for type B individuals which can then be input into the IBM.
	Code file names with the phrase "MCMC" estimate the parameters for the underlying susceptibility distributions (p_A, p_B, f_A, k, theta) using Metropolis-Hastings MCMC (and ABC in the discrete case). The input for these files come from CleanData.R.
	The file DetHet_DynamEpi.R calculates the power to detect heterogeneity in susceptibility using the likelihood equations 1 and 2 in the main text. The input for this file is simulations from the ".cpp" files. We used 100 simulations for each parameter combination tested.
	CleanData.R separates the contact tracing data from each simulated epidemic into smaller datasets containing 50, 200, 1000, or 5000 focal individuals exposed that can then be input into the "MCMC" files to estimate parameters.

MissingContacts runs simulations for both detection and estimation of heterogeneity in susceptibility when there are missing contacts (individuals that have been previously exposed more times than observed).
	Cont_MissExps_GenData.R generates two files with a table of F focal individuals in the first file and F(N-1) naive individuals in the second containing their number of previous exposures, risk of infection, and whether they were infected (1) or not (0). The output data from this code should be input into Cont_MCMC_MissExps_cont.R to run MCMC.
	Cont_MCMC_MissExps_cont.R estimates the parameters for the underlying susceptibility distribution (k, theta) using Metropolis-Hastings MCMC. The input for this file comes from Cont_MissExps_GenData.R.
	Cont_DetHiS_MissExps.R runs simulations for the power to detect heterogeneity in susceptibility when there are missing contacts.
