// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>

#include <iostream>
#include <cmath>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;


void log_like(double Lhom[], double Lhet[], int prev_inf[], int naive_inf[], int Focal, int N, int sims, double false_neg);


// [[Rcpp::export]]
double find_power(double k, double theta, double f_inf, double false_neg)
{
	int num_sims = 1000, Focal = 200, N = 5, N2 = 5;	// number of simulations, focal indivs, indivs per contact group
	int naive_inf[num_sims] = {0}, prev_inf[num_sims] = {0};	// number of infected individuals in 2nd round
	int sim, f, n, j;	// variables for iterating over
	double risk[N], riskp[Focal], prob_I[N];	// risks/probs of inf for temp storage
	int I, Iprev, check_fn, num_naiveI;	// temp variables for infected indiv and checking false negs
	double prev_indiv[1], risk_naive, prob_focal_inf, prob_naive_inf;	// risk/prob of inf for focal/naive indiv
	double Lhom[num_sims], Lhet[num_sims];	// log-likelihoods of data
	double ratio_test, power = 0.0;	// likelihood ratio test and power to detect het in susc
	int foc[num_sims] = {0};	// vector to track number of focal indivs

	// set up random number generator
	// GSL's Taus generator -- is this the best one? does it matter?
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	// initialize generator with seed 3
	gsl_rng_set(rng, 3);

	for(sim=0; sim<num_sims; sim++)
	{
		j = 0;		// initialize counter var to save risks
		while(foc[sim] < Focal)
		{
			// first exposure events
			// set up risk/prob of infection distributions
			for(n=0; n<N; n++)
			{
				if(k == -1.0)	// cv=0, no HiS=all indivs should have prob f_inf of inf
					risk[n] = -log(1.0 - f_inf);
				else
					risk[n] = gsl_ran_gamma(rng, k, theta);
				prob_I[n] = 1 - exp(-risk[n]);

				// infect individuals
				// use bernoulli for each indiv to check if infected, store risk for those not infected
				// check for infection
				I = gsl_ran_binomial(rng, prob_I[n], 1);

				// if indiv infected, check if false neg
				if(I == 1)
				{
					// determine if false neg
					check_fn = gsl_ran_binomial(rng, false_neg, 1);
					// if indiv is false neg, save as focal indiv with risk 0
					if(check_fn == 1)
					{
						if(foc[sim] < Focal)
						{
							riskp[j] = 0.0;
							j++;
							foc[sim]++;
						}
					}
				}
				if(I == 0)
				{
					if(foc[sim] < Focal)
					{
						riskp[j] = risk[n];
						j++;
						foc[sim]++;
					}
				}
			}
		}



		for(f=0; f<Focal; f++)
		{
 			// second exposure events
			// infection of focal indiv, riskp[f] is their risk
			prob_focal_inf = 1 - exp(-riskp[f]);
			Iprev = gsl_ran_binomial(rng, prob_focal_inf, 1);

			// determine risks/probs of inf for naive indivs and simulate infection
			num_naiveI = 0;
			for(n=0; n<N-1; n++)
			{
				if(k == -1.0)	// cv=0, no HiS=all indivs should have prob f_inf of inf
					risk_naive = -log(1.0 - f_inf);
				else
					risk_naive = gsl_ran_gamma(rng, k, theta);
				prob_naive_inf = 1 - exp(-risk_naive);
				I = gsl_ran_binomial(rng, prob_naive_inf, 1);

				// if naive indiv inf, add to number of naive indivs inf this round
				if(I == 1)
					num_naiveI++;
			}

			// store number of indivs infected in 2nd round, remove some indivs b/c false negs
			naive_inf[sim] += gsl_ran_binomial(rng, 1-false_neg, num_naiveI);
			prev_inf[sim] += gsl_ran_binomial(rng, 1-false_neg, Iprev);
		}
	}


	// calculate log-likelihoods of the data
	log_like(Lhom, Lhet, prev_inf, naive_inf, Focal, N, num_sims, false_neg);

	// calculate likelihood ratio tests and determine percentage (power) greater than critical value (3.84)
	for(sim=0; sim<num_sims; sim++)
	{
		ratio_test = -2.0 * (Lhom[sim] - Lhet[sim]);
		if(ratio_test > 3.84)
			power += 100.0 / num_sims;
	}

	// output back to R function
	return(power);
}








void log_like(double Lhom[], double Lhet[], int prev_inf[], int naive_inf[], int Focal, int N, int sims, double false_neg)
{
	int prev_inf_sim, naive_inf_sim;	// number of individuals inf observed with false negs in each simulation
	double p_avg_obs_sim, p_prev_obs_sim, p_naive_obs_sim;	// prob of inf that was observed each sim with false negs
	double p_avg_true, p_prev_true, p_naive_true;		// prob of inf that should have been observed without false negs
	double p_xI_hom, p_xI_het;		// prob of a focal indiv being prev inf (xI)
	int sim, xI, xF, xN;		// variables to track through for loops
	double LhomF, LhomN, LhetF, LhetN;	// variables for calculating likelihoods
	double temp;	// temporary storage for gsl calculation used for multiple likelihoods

	// calculate log-likelihoods correcting for false negs
	for(sim=0; sim<sims; sim++)
	{
		// initialize Lhom and Lhet to 0
		LhomF = 0.0;
		LhomN = 0.0;
		LhetF = 0.0;
		LhetN = 0.0;

		// set number of focal/naive indivs inf each simulation to reduce memory access in the for loops
		prev_inf_sim = prev_inf[sim];
		naive_inf_sim = naive_inf[sim];

		// calculate observed prob of being inf in each simulation, multiply by 1.0 to have type double instead of int
		p_prev_obs_sim = 1.0 * prev_inf_sim / Focal;
		p_naive_obs_sim = 1.0 * naive_inf_sim / (Focal * (N - 1));
		p_avg_obs_sim = 1.0 * (prev_inf_sim + naive_inf_sim) / (Focal * N);

		// calculate "true" prob of inf for naive indivs that should have been obs without false negs
		// cap probs at 1, b/c not accounting for randomness in obs values
		p_naive_true = std::min(p_naive_obs_sim / (1 - false_neg), 1.0);

		// calculate prob of a focal indiv being prev inf (xI) in heterogeneous case
		p_xI_het = (p_naive_true*false_neg) / (1 - p_naive_true + p_naive_true*false_neg);

		// calculate "true" probs of inf for avg indiv that should have been obs without false negs
		// cap probs at 1, b/c not accounting for randomness in obs values
		// approx number of focal indivs that should be removed b/c prev inf as an expected value
		p_avg_true = std::min((Focal*N*p_avg_obs_sim) / ((Focal*N - Focal*p_xI_het) * (1 - false_neg)), 1.0);

		// calculate prob of a focal indiv being prev inf (xI) in homogeneous case
		p_xI_hom = (p_avg_true*false_neg) / (1 - p_avg_true + p_avg_true*false_neg);


		// xI=number of prev inf indivs thought to be focal, between 0 and F but L=0 at F
		for(xI=0; xI<=Focal-1; xI++)
		{
			// calculate "true" probs of inf for focal indivs that should have been obs without false negs
			// cap probs at 1, b/c not accounting for randomness in obs values
			p_prev_true = std::min((Focal*p_prev_obs_sim) / ((Focal - xI) * (1 - false_neg)), 1.0);

			// xF=number of focal indivs that were inf but false neg
			for(xF=0; xF<=Focal-xI-prev_inf_sim; xF++)
			{
				temp = gsl_ran_binomial_pdf(xF, false_neg, prev_inf_sim+xF);	// calculate this value just once since used in both likelihoods

				LhomF += gsl_ran_binomial_pdf(prev_inf_sim+xF, p_avg_true, Focal-xI)
					* gsl_ran_binomial_pdf(xI, p_xI_hom, Focal)
					* temp;

				LhetF += gsl_ran_binomial_pdf(prev_inf_sim+xF, p_prev_true, Focal-xI)
					* gsl_ran_binomial_pdf(xI, p_xI_het, Focal)
					* temp;
			}
		}


		// xN=number of naive indivs that were inf but false neg
		for(xN=0; xN<=Focal*(N-1)-naive_inf_sim; xN++)
		{
			temp = gsl_ran_binomial_pdf(xN, false_neg, naive_inf_sim+xN);		// calculate this value just once since used in both likelihoods

			LhomN += gsl_ran_binomial_pdf(naive_inf_sim+xN, p_avg_true, Focal*(N-1))
				* temp;

			LhetN += gsl_ran_binomial_pdf(naive_inf_sim+xN, p_naive_true, Focal*(N-1))
				* temp;
		} 

		// sum log-likelihoods from focal and naive individuals for total likelihoods
		Lhom[sim] = log(LhomF) + log(LhomN);
		Lhet[sim] = log(LhetF) + log(LhetN);
	}
}



