Below are brief descriptions of the code. These scripts are used to analyze the simulated data and generate figures for the discrete (2gps) and continuous (Cont) cases.

2gps_AnalyzeMCMC_SIRdynams.R and Cont_AnalyzeMCMC_SIRdynams.R are used to analyze the MCMC data.
	MCMC output files are read in and checked for convergence.
	Then, 1,000 parameter combinations are sampled and used to run SIR dynamics and find the 95% CIs for the SIR dynamics.
	The saved 95% CIs from these scripts are then used to generate all figures in MakePlots_ParamEsts_SIRdynams.R.

2gps_DetHiSFigs.R and Cont_DetHiSFigs.R are used to generate all detection plot figures throughout the main text and supplementary information.

MakePlots_ParamEsts_SIRdynams.R is used to generate all figures related to parameter estimates from MCMC and SIR dynamics.
	The data used here are initially read in and analyzed in 2gps_AnalyzeMCMC_SIRdynams.R and Cont_AnalyzeMCMC_SIRdynams.R.

MakePlots_DetHit_GoodnessofFit.R is used to generate Figs S17 and S18 for the power to detect heterogeneity in transmission using the simulation-based goodness of fit method.

FinalEpidemicSize.R is used to simulate and generate Fig S9 determining the effect of CV of risk of infection and f_A on the final epidemic size in the discrete case.

LHoodeq_detHiS_pNvspF_FigS6.R is used to investigate the likelihood equations used for detection (Eqs 1, 2) and how they affect the power to detect heterogeneity in susceptibility.
	This script is used to generate Fig S6.