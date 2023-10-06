// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>

#include <iostream>
#include <cmath>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;


void log_like(double Lhom[], double Lhet[], int prev_inf[], int naive_inf[], int Focal, int N, int sims, double false_neg);


// [[Rcpp::export]]
double find_power(double p_a, double p_b, double f_a, double false_neg)
{
	int Focal = 200, N = 5, N2 = 5;		// number of focal indivs and numbers of indivs per contact group in 1st and 2nd rounds
	int num_sims = 1000;	// number of simulations to run
	// initialize vectors and matrices to hold data
	int Na, Nb, Ia, Ib;	// number of type A/B indivs in first contact groups and number infected 
	int Na_new[Focal], Nb_new[Focal];	// number of type A/B indivs in second contact groups
	int Ia_new[Focal][num_sims], Ib_new[Focal][num_sims], I[Focal][num_sims];	// matrices for infected indivs in each group
	int f, sim;	// variable for iterations
	int fneg_sim;	// number of false negs in 1st round in each sim
	double f_a_new[num_sims], f_b_new[num_sims];	// fraction type A/B after 1st round of infection
	double rand_unif;	// uniform RV to decide which type indiv chosen as focal indiv
	int naive_inf[num_sims] = {0}, prev_inf[num_sims] = {0};	// number of inf indivs in 2nd round
	double Lhom[num_sims], Lhet[num_sims];	// log-likelihoods of data
	double ratio_test;	// likelihood ratio test
	double power = 0.0;	// power to detect het in susc
	int foc[num_sims] = {0};	// vector to track number of focal indivs
	int a_focal[num_sims] = {0}, b_focal[num_sims] = {0}, num_fneg[num_sims] = {0};	// vectors to track number of focal indivs that are A/B/FN

	// set up random number generator
	// GSL's Taus generator -- is this the best one? does it matter?
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	// initialize generator with seed 3
	gsl_rng_set(rng, 3);

	for(sim=0; sim<num_sims; sim++)
	{
		while(foc[sim] < Focal)
		{
			// first exposure events
			// number of indivs of each type, add up to total pop size N, based on fraction type A
			// gsl_ran_binomial(generator, prob, num trials)
			Na = gsl_ran_binomial(rng, f_a, N);
			Nb = N - Na;

			// infect individuals
			Ia = gsl_ran_binomial(rng, p_a, Na);
			Ib = gsl_ran_binomial(rng, p_b, Nb);

			// check if any indivs not inf -- if not, save as focal indivs
			if(Na - Ia > 0)
			{
				a_focal[sim] += Na - Ia;
				foc[sim] += Na - Ia;
			}
			if(Nb - Ib > 0)
			{
				b_focal[sim] += Nb - Ib;
				foc[sim] += Nb - Ib;
			}

			// number of indivs that appear as false negs, assume 100% immunity
			fneg_sim = gsl_ran_binomial(rng, false_neg, Ia + Ib);

			num_fneg[sim] += fneg_sim;
			foc[sim] += fneg_sim;
		}

		// calculate new fractions of type A and B after 1st infection round
		f_a_new[sim] = 1.0 * a_focal[sim] / (a_focal[sim] + b_focal[sim] + num_fneg[sim]);
		f_b_new[sim] = 1.0 * b_focal[sim] / (a_focal[sim] + b_focal[sim] + num_fneg[sim]);
	}


	// second exposure events
	for(f=0; f<Focal; f++)
	{
		// number naive indivs of each type, add up to N2-1
		Na_new[f] = gsl_ran_binomial(rng, f_a, N2-1);
		Nb_new[f] = N2 - 1 - Na_new[f];

		// determine number of infected individuals
		for(sim=0; sim<num_sims; sim++)
		{
			// infection of naive indivs
			Ia_new[f][sim] = gsl_ran_binomial(rng, p_a, Na_new[f]);
			Ib_new[f][sim] = gsl_ran_binomial(rng, p_b, Nb_new[f]);

			// infection of prev exposed individual
			// use uniform RV to decide if choose type A, B, or C
			rand_unif = gsl_ran_flat(rng, 0, 1);

			// check type and check for infection
			if(rand_unif < f_a_new[sim])	// type A
				I[f][sim] = gsl_ran_binomial(rng, p_a, 1);
			else if(rand_unif < f_a_new[sim] + f_b_new[sim])	// type B
				I[f][sim] = gsl_ran_binomial(rng, p_b, 1);
			else	// type C, prob of inf=0
				I[f][sim] = 0;


			// sum across focal indiv groups for number infected indivs in 2nd round total
			// remove some inf indivs b/c false neg tests
			naive_inf[sim] += gsl_ran_binomial(rng, 1-false_neg, Ia_new[f][sim] + Ib_new[f][sim]);
			prev_inf[sim] += gsl_ran_binomial(rng, 1-false_neg, I[f][sim]);
		}
	}

	// calculate log-likelihoods of data
	log_like(Lhom, Lhet, prev_inf, naive_inf, Focal, N, num_sims, false_neg);

	// calculate likelihood ratio tests and determine number greater than critical value (3.84)
	for(sim=0; sim<num_sims; sim++)
	{
		ratio_test = -2.0 * (Lhom[sim] - Lhet[sim]);
		if(ratio_test > 3.84)
			power++;
	}

	// calculate power as percentage
	power /= num_sims;
	power *= 100.0;

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



