This directory (Estimation) contains data used for parameter estimation and MCMC.

File names specify relevant parameters, such as the number of focal individuals (F), number of naive individuals (N), parameters dictating the underlying susceptibility distribution (p_A, p_B, f_A, k, theta), etc.

Contact_Networks contains data for the simulated contact networks. Each file consists of a vector of the number of individuals infected in each network. These were used to run all MCMC in this directory.
These networks were also used to run MCMC in the discrete case with informative priors for p_A, p_B, or f_A (Discrete Priors directory).

MCMC contains the parameter estimates from MCMC without any modifications. These data were used to generate Fig 7 in the main text.

MCMC_ChangingErrorTolerance contains the parameter estimates from MCMC for the discrete case where the error tolerance for ABC was changed. These data were used to generate Fig S7.
In the file names, "tol0" means 0% error tolerance, "tol0.1" means 10% error tolerance, and no specification means 1% error tolerance.

MCMC_WrongModel contains the parameter estimates from MCMC where the wrong underlying model was assumed. These data were used to generate Fig S8.
In the file names, "ContMC2gpsData" means the discrete case was the correct model with the continuous case incorrectly assumed and vice versa.