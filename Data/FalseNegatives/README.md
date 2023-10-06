This directory (FalseNegatives) contains data used for detecting heterogeneity in susceptibility, parameter estimation, and MCMC in the presence of false negatives with or without the adjusted method.

File names specify relevant parameters, such as the number of focal individuals (F), number of naive individuals (N), parameters dictating the underlying susceptibility distribution (p_A, p_B, f_A, k, theta), etc.
File names that contain the phrase "adj" mean that data was simulated with the method adjusted to account for false negatives.

Contact_Networks contains data for the simulated contact networks. Each file consists of a vector of the number of individuals infected in each network. These were used to run all MCMC in this directory.

Detection contains simulations for the power to detect heterogeneity in susceptibility in the presence of false negatives with or without the adjusted method. These data were used to generate Figs S11 and S12.

MCMC contains the parameter estimates from MCMC with or without the adjusted method for false negatives. These data were used to generate Fig S13.