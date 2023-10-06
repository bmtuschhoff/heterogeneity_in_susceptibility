// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>

#include <iostream>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;


// [[Rcpp::export]]
double likelihood(double kval, double thetaval, int Focal, int N, int prev_inf, int naive_inf, double false_neg)
{
	int num_sims = 100;	// number of simulations to reduce randomness in rbinom
	int FocalInf[num_sims], NonFocalInf[num_sims];	// number of infected indivs
	int sim, i;	// variables for iterating over
	NumericVector risk_focal_inf(Focal);	// risks of being infected for focal indivs
	double prob_focal_inf, prob_naive_inf;	// probabilities of being infected for each focal and naive indiv
	int num_focal;	// track number of focal indivs saved for each sim
	double prob_inf, risk;	// temp prob/risk of inf for an indiv
	int I, check_fn;	// temp variables for checking for inf in focal indivs
	double tol = 0.01;	// error tolerance allowed for ABC likelihood
	double LFocal = 0.0, LNonFocal = 0.0, L;	// likelihoods


	// set up random number generator
	// GSL's Taus generator -- is this the best one? does it matter?
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	// initialize generator with seed 3
	gsl_rng_set(rng, 3);


	for(sim=0; sim<num_sims; sim++)
	{
		num_focal = 0;	// initialize num focal indivs

		// determine probs of inf of focal indivs
		// draw original prob of inf, check if inf in 1st round and if not infected or if false neg save prob as focal indiv
		while(num_focal < Focal)
		{
			// draw probability of infection
			// gsl_ran_gamma(generator, k, theta)
			risk = gsl_ran_gamma(rng, kval, thetaval);
			prob_inf = 1 - exp(-risk);

			// simulate if individual is infected
			// gsl_ran_binomial(generator, prob, num trials)
			I = gsl_ran_binomial(rng, prob_inf, 1);

			// if individual is infected, check if false negative
			if(I == 1)
			{
				check_fn = gsl_ran_binomial(rng, false_neg, 1);

				// if individual is false neg, record as focal indiv with risk 0
				// if not false neg, do nothing
				if(check_fn == 1)
				{
					if(num_focal < Focal)
					{
						risk_focal_inf[num_focal] = 0.0;
						num_focal++;
					}	
				}
			}
			// if individual is not infected, record as focal indiv with that risk of inf
			if(I == 0)
			{
				if(num_focal < Focal)
				{
					risk_focal_inf[num_focal] = risk;
					num_focal++;
				}
			}
		}



		// second round of infections
		// get probabilities of being infected for naive indivs and check if infected
		for(i=0; i<Focal*(N-1); i++)
		{
			prob_naive_inf = 1 - exp(-gsl_ran_gamma(rng, kval, thetaval));
			NonFocalInf[sim] += gsl_ran_binomial(rng, prob_naive_inf, 1);
		}
		// calculate probs of infection for focal indivs from risks
		for(i=0; i<Focal; i++)
		{
			prob_focal_inf = 1 - exp(-risk_focal_inf[i]);
			FocalInf[sim] += gsl_ran_binomial(rng, prob_focal_inf, 1);
		}


		// remove false negs to get total num indivs infected in second round
		FocalInf[sim] = gsl_ran_binomial(rng, 1-false_neg, FocalInf[sim]);
		NonFocalInf[sim] = gsl_ran_binomial(rng, 1-false_neg, NonFocalInf[sim]);

		// check if FocalInf=prev_inf and NonFocalInf=naive_inf for this sim within some error tolerance
		// if so, add to likelihood calculation as a success
		if(1.0*abs(FocalInf[sim] - prev_inf) / prev_inf <= tol)
			LFocal += 1.0 / num_sims;
		if(1.0*abs(NonFocalInf[sim] - naive_inf) / naive_inf <= tol)
			LNonFocal += 1.0 / num_sims;
	}

	// calculate overall log-likelihood of the data
	L = log(LFocal) + log(LNonFocal);

	return(L);
}






















