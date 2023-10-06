This directory (ChangeN) contains data used for parameter estimation and MCMC.

File names specify relevant parameters, such as the number of focal individuals (F), number of naive individuals (N), parameters dictating the underlying susceptibility distribution (p_A, p_B, f_A, k, theta), etc.

Contact_Networks contains data for the simulated contact networks. Each file consists of a vector of the number of individuals infected in each network. These were used to run all MCMC in this directory.

MCMC_ChangeN contains the parameter estimates from MCMC where N was changed between N=100 and N=5. These data were used to generate Fig S5.
In the file names, "N100sub" means that data came from contact networks with N=100 but subset to N=5 in order to better compare the effect of changing N while minimizing stochasticity from simulating the data.