This directory (HiT) contains data used for detecting heterogeneity in susceptibility, parameter estimation, and MCMC in the presence of heterogeneity in transmission. There is also data for the power to detect heterogeneity in transmission using the simulation-based goodness of fit method.

File names specify relevant parameters, such as the number of focal individuals (F), number of naive individuals (N), parameters dictating the underlying susceptibility distribution (p_A, p_B, f_A, k, theta), etc.

Contact_Networks contains data for the simulated contact networks used to check for the effect of HiT without accounting for it. Each file consists of a vector of the number of individuals infected in each network. These were used to run MCMC in this directory.

Contact_Networks_AccountForHiT contains data for the simulated contact networks used to run MCMC when accounting for HiT. Each file consists of a vector of the number of individuals infected in each network. "inf1" is data for individuals from the first exposure event, "inf2_naive" is for naive individuals from the second exposure event, and "inf2_prev" is for focal individuals from the second exposure event.

Detection_HiS contains simulations for the power to detect heterogeneity in susceptibility (HiS) in the presence of heterogeneity in transmission. These data were used to generate Figs S14 and S15.

Detection_HiT contains simulations for the power to detect heterogeneity in transmission (HiT) using the simulation-based goodness of fit test. These data were used to generate Fig S18.

MCMC contains the parameter estimates from MCMC for HiS (p_A, p_B, f_A, k, theta). These data were used to generate Fig S16.

