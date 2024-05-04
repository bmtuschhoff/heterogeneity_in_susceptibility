// individual based SIR Gillespie model including different forms of heterogeneity, continuous HiS
// heterogeneity can be in contact rate or given contact and it can be HiT or HiS
// also set up different correlations between HiT and HiS (positive, negative, uncorrelated)
// output is contact networks containing number of naive and focal individuals exposed and infected by each infected individual

#include <iostream>
#include <cmath>		// sqrt and log
#include <math.h>
#include <gsl/gsl_rng.h>	// to set seeds for random dists
#include <gsl/gsl_randist.h>	// random dists
#include <fstream>		// file input/output
#include <random>		// set up distributions of rates to choose indiv infected/recovered
#include <cstring>		// to compare strings for inputting parameters from command line
#include <algorithm> 		// contains sort() function

using namespace std;


bool parseOptions(int argc, char *argv[]);
double find_tot_rate(double *expose_rate, double *recover_rate, char current_pop_status[], vector<double> c_rate_Snaive, vector<double> c_rate_Sprev, vector<double> c_rate_I, double beta_avg, int numS_naive, int numS_prev, int numI);
void expose(vector<double> &risk_Snaive, vector<double> &lamb_Snaive, vector<double> &c_rate_Snaive, vector<int> &index_Snaive, vector<double> &risk_Sprev, vector<double> &lamb_Sprev, vector<double> &c_rate_Sprev, vector<int> &index_Sprev, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<int> &expose_naive, vector<int> &inf_naive, vector<int> &expose_prev, vector<int> &inf_prev, double beta_avg, int *numS_naive, int *numS_prev, int *numI, gsl_rng *rng, mt19937 &rng2);
void infect(char current_pop_status[], vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, double beta_avg, int *numS, int *numI, mt19937 &rng2);
void recover(char current_pop_status[], vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<int> &expose_naive, vector<int> &inf_naive, vector<int> &expose_prev, vector<int> &inf_prev, int *numI, double t, mt19937 &rng2, ofstream &ofs);
double c_rate(int S_index, int I_index);


// default parameter settings, can change with flags
string output = "results.txt";	// output file name
char his = 'F';			// set whether to include HiS given contact as false
char hit = 'F';			// set whether to include HiT given contact as false
char hic = 'F';			// set whether to include het in contact as false
string cor = "none";		// set correlation between HiS and HiT -- could be pos, neg, or none
int pop_size = 4;		// population size
int init_inf = 20;		// initial number of infected indivs
double R0 = 3;			// R0
double gam = 0.1;		// recovery rate
double cv_s = 1.3;		// coefficient of variation of risk (HiS)
double cv_t = 1.0 / sqrt(0.5);		// coefficient of variation (HiT)
double cv_c = 1.0;		// coefficient of variation for contact rate
double f_inf = 0.25;		// expected fraction of naive individuals infected
int seed = 3;			// seed for random number generators


int main(int argc, char *argv[])
{
	parseOptions(argc, argv);	// set up parameters from command line input and flags

	double k, m, psi;	// shape params for HiS, HiT, HiC distributions
	double theta;		// scale param for HiS distribution
	double beta_avg, prob_inf_avg;	// avg transmission rate and probability of infection
	int time = 1000;	// total length of time to run the simulation
	int numS_naive, numS_prev, numI;		// number of S (naive and previously exposed) and I indivs in the population
	//int init_inf = 20;	// initial number of infected indivs
	vector<double> risk_Snaive(pop_size), lamb_Snaive(pop_size), c_rate_Snaive(pop_size);	// vectors of risks, forces of inf, and contact rates for each naive S indiv
	vector<double> risk_Sprev(pop_size), lamb_Sprev(pop_size), c_rate_Sprev(pop_size);	// vectors of risks, forces of inf, and contact rates for each previously exposed S indiv
	vector<double> risk_I(pop_size), lamb_I(pop_size), c_rate_I(pop_size);	// vectors of risks, forces of inf, and contact rates for each I indiv
	vector<int> index_Snaive(pop_size), index_Sprev(pop_size), index_I(pop_size);	// vectors of indices for each S, I indiv
	vector<int> expose_naive(pop_size), inf_naive(pop_size), expose_prev(pop_size), inf_prev(pop_size);	// vectors of num naive, focal indivs exposed and infected by each I indiv
	int indiv, i, j;	// variables for iteration
	int index;		// index of node in population
	char current_pop_status[pop_size];	// vector holding current status of each indiv (S,I,R)
	int init_inf_indivs[init_inf];		// vector holding the indices of the indivs initially infected
	int pop_indices[pop_size];		// vector of indices of all indivs in the pop
	double t = 0.0;		// track time throughout epidemic
	double expose_rate, recover_rate, tot_rate;	// rates of infection, recovery, and both across all indivs
	double event;	// uniform RV to determine which event occurs next

	// calculate other params from those read-in
	k = 1.0 / pow(cv_s, 2);
	theta = pow(1.0 - f_inf, -1.0/k) - 1.0;
	m = 1.0 / pow(cv_t, 2);
	psi = 1.0 / pow(cv_c, 2);


	// calculate avg transmission rate from R0, gam, pop size
	// if there is no heterogeneity in susceptibility (CV=0), avg probability is f_inf
	if(cv_s == 0.0)
	{
		his = 'F';
		prob_inf_avg = f_inf;
	}
	else
		prob_inf_avg = 1.0 - pow(1 + theta, -k);

	beta_avg = (gam * R0) / (prob_inf_avg * (pop_size - init_inf));


	// set up output file stream
	// check if file already exists
	bool file_exists = false;
	if(FILE *file = fopen(output.c_str(), "r"))
	{
		fclose(file);
		file_exists = true;
	}

	// declare file output stream
	ofstream ofs;
	ofs.open(output.c_str(), ofstream::app);	// set up to append new data to end of file

	// write header to file if it doesn't already exist
	if(file_exists == false)
	{
		ofs << "Time" << "\tID" << "\tNaive exposed" << "\tNaive infected" << "\tFocal exposed" << "\tFocal infected" << endl;
	}
	


	// set up random number generator for using with GSL functions - GSL's Taus generator
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	// initialize generator with seed
	gsl_rng_set(rng, seed);
	// set up second random number generator for using with discrete_distribution
	mt19937 rng2(seed);


	// draw risk value, force of infection, and contact rate for each indiv from gamma dist depending on whether including that het
	for(indiv = 0; indiv < pop_size; indiv++)
	{
		if(his == 'T')	// there is HiS
			risk_Snaive[indiv] = gsl_ran_gamma(rng, k, theta);
		else	// no HiS
			risk_Snaive[indiv] = -log(1.0 - f_inf);

		if(hit == 'T')	// there is HiT
			lamb_Snaive[indiv] = gsl_ran_gamma(rng, m, 1.0/m);
		else	// no HiT
			lamb_Snaive[indiv] = 1.0;

		if(hic == 'T')	// there is HiC
			c_rate_Snaive[indiv] = gsl_ran_gamma(rng, psi, 1.0/psi);
		else	// no HiC
			c_rate_Snaive[indiv] = 1.0;
	}



	// if command line param for correlation set to positive cor between HiS and HiT
	// order risks and lambs so that both ordered descending order
	if(strcmp(cor.c_str(), "pos") == 0)
	{
		sort(risk_Snaive.begin(), risk_Snaive.end(), greater<double>());
		sort(lamb_Snaive.begin(), lamb_Snaive.end(), greater<double>());
	}

	// if command line param for correlation set to negative cor between HiS and HiT
	// order risks and lambs so that ordered opposite each other - descending risks and ascending lambs
	if(strcmp(cor.c_str(), "neg") == 0)
	{
		sort(risk_Snaive.begin(), risk_Snaive.end(), greater<double>());
		sort(lamb_Snaive.begin(), lamb_Snaive.end());
	}




	// initialize vector with current status of each indiv
	// randomly select init_inf indiv(s) to be infected at the start
	for(i = 0; i < pop_size; i++)	// initialize vector of indices for each indiv in pop, for randomly choosing indivs
	{				// also initialize all indivs in pop as S
		pop_indices[i] = i;
		index_Snaive[i] = i;
		current_pop_status[i] = 'S';
	}	

	gsl_ran_choose(rng, init_inf_indivs, init_inf, pop_indices, pop_size, sizeof(int));	// randomly select the indices of init_inf indivs out of pop_size to start off as inf

	for(i = 0; i < init_inf; i++)	// set up indivs to be infected at start
		current_pop_status[init_inf_indivs[i]] = 'I';

	// move each initially infected indiv from S to I vectors
	for(i = 0; i < init_inf; i++)
	{
		// move each param from S to first available row in I (all rows available currently)
		risk_I[i] = risk_Snaive[init_inf_indivs[i]];
		lamb_I[i] = lamb_Snaive[init_inf_indivs[i]];
		c_rate_I[i] = c_rate_Snaive[init_inf_indivs[i]];
		index_I[i] = index_Snaive[init_inf_indivs[i]];

		// replace data for each param in S with data from last row of S
		risk_Snaive[init_inf_indivs[i]] = risk_Snaive[pop_size-1-i];
		lamb_Snaive[init_inf_indivs[i]] = lamb_Snaive[pop_size-1-i];
		c_rate_Snaive[init_inf_indivs[i]] = c_rate_Snaive[pop_size-1-i];
		index_Snaive[init_inf_indivs[i]] = index_Snaive[pop_size-1-i];
	}



	// set num S and num I indivs at start of epidemic
	numS_naive = pop_size - init_inf;
	numS_prev = 0;
	numI = init_inf;
	
	

	// start Gillespie algorithm - draw time, determine which event happened, carry out event, repeat
	while(t <= time)
	{

		// determine rate of each event and combined total rate -- events are infection or recovery
		tot_rate = find_tot_rate(&expose_rate, &recover_rate, current_pop_status, c_rate_Snaive, c_rate_Sprev, c_rate_I, beta_avg, numS_naive, numS_prev, numI);

		if(tot_rate == -1)	// no inf indivs, so epidemic is over - close output file and end simulation
		{
			ofs.close();
			cout << seed << "\t" << t << "\t" << pop_size - (numS_naive + numS_prev) << endl;	// write out seed, time, and FES to screen
			return(0);
		}

		// determine time of next event
		t += gsl_ran_exponential(rng, 1 / tot_rate);	// rate lambda for exp dist is tot_rate, this fntn takes the mean 1/lambda

		// determine next event and carry out
		// draw uniform RV and compare to probabilities of events
		event = gsl_ran_flat(rng, 0, 1);
		
		
		if(event <= expose_rate / tot_rate)	// event is exposure
		{
			expose(risk_Snaive, lamb_Snaive, c_rate_Snaive, index_Snaive, risk_Sprev, lamb_Sprev, c_rate_Sprev, index_Sprev, risk_I, lamb_I, c_rate_I, index_I, expose_naive, inf_naive, expose_prev, inf_prev, beta_avg, &numS_naive, &numS_prev, &numI, rng, rng2);
			//infect(current_pop_status, risk_S, lamb_S, c_rate_S, index_S, risk_I, lamb_I, c_rate_I, index_I, beta_avg, &numS, &numI, rng2);
		}
		else	// event is recovery
		{
			recover(current_pop_status, risk_I, lamb_I, c_rate_I, index_I, expose_naive, inf_naive, expose_prev, inf_prev, &numI, t, rng2, ofs);
		}

		// write out time and current compartment status of each indiv in population -- write out when an I indiv recovers, print FES at end???
		//ofs << t << "\t" << numS << "\t" << numI << "\t" << pop_size-numS-numI << endl;
	}

	ofs.close();	// close ouput file
	cout << seed << "\t" << t << "\t" << pop_size - (numS_naive + numS_prev) << endl;	// write out seed, time, and FES to screen
	return(0);
}




// find total rates of next event
double find_tot_rate(double *expose_rate, double *recover_rate, char current_pop_status[], vector<double> c_rate_Snaive, vector<double> c_rate_Sprev, vector<double> c_rate_I, double beta_avg, int numS_naive, int numS_prev, int numI)
{
	int i, j;	// variables for iteration
	
	*expose_rate = 0.0;
	*recover_rate = 0.0;
	for(i = 0; i < numI; i++)	// calculate infection and recovery rates across all infected indivs
	{

		*recover_rate += gam;	// each inf indiv could recover with rate gamma

		for(j = 0; j < numS_naive + numS_prev; j++)	// determine infection rate across all S indivs that could have contacted this indiv
		{
			*expose_rate += beta_avg * c_rate(j, i);
		}
	}

	if(numI == 0)	// no one is infected, epidemic is over, flag total rate
	{
		return(-1);
	}

	return(*expose_rate + *recover_rate);	// total rate of next event
}




// know event is that a S indiv is exposed, need to decide which one and whether inf
void expose(vector<double> &risk_Snaive, vector<double> &lamb_Snaive, vector<double> &c_rate_Snaive, vector<int> &index_Snaive, vector<double> &risk_Sprev, vector<double> &lamb_Sprev, vector<double> &c_rate_Sprev, vector<int> &index_Sprev, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<int> &expose_naive, vector<int> &inf_naive, vector<int> &expose_prev, vector<int> &inf_prev, double beta_avg, int *numS_naive, int *numS_prev, int *numI, gsl_rng *rng, mt19937 &rng2)
{
	int i, j;	// variables for iterations
	int numS = *numS_naive + *numS_prev;	// total number of S individuals
	int expose_index, S_expose, I_expose;	// indices denoting which S indiv is exposed and I indiv is the exposer
	vector<double> expose_rate(*numI * numS, 0);	// vector of exposure rate for each combination of I and S indiv
	double infect_rate;	// inf rate for S indiv exposed
	int is_inf;		// 1 if inf, 0 if not

	// determine exposure rates across all I and S indivs
	for(i = 0; i < *numI; i++)
	{
		for(j = 0; j < numS; j++)
		{
			expose_rate[i*numS + j] = beta_avg * c_rate(j, i);
		}
	}

	// have vector that contains exposure rate for each pair of I and S indivs
	// need to pick which indiv is exposed and by whom

	// choose index from exposure rates
	discrete_distribution<int> pick_exp(expose_rate.begin(), expose_rate.end());	// set up dist of exposure rates
	expose_index = pick_exp(rng2);		// pick index of I and S for who exposed

	// from expose_index, figure out which I indiv is the exposer and which S indiv is exposed
	I_expose = expose_index / numS;
	S_expose = expose_index % numS;

	// check if S indiv is naive or prev exposed and add to vector for that I indiv as exposed
	// also check if S indiv is infected -- if not move to Sprev, if so move to I and note inf by that I indiv
	if(S_expose < *numS_naive)	// S indiv is naive
	{
		expose_naive[index_I[I_expose]]++;	// add 1 exposed Snaive indiv for this I

		// check if S indiv is infected
		infect_rate = risk_Snaive[S_expose] * lamb_I[I_expose];
		is_inf = gsl_ran_binomial(rng, 1-exp(-infect_rate), 1);		// use bernoulli to check if indiv infected
		if(is_inf == 1)		// indiv is inf, move to I and note for that I indiv
		{
			inf_naive[index_I[I_expose]]++;		// note that I indiv inf a Snaive indiv

			// move each param from Snaive to first available row in I - since there are numI I indivs, the first row would be at index numI
			risk_I[*numI] = risk_Snaive[S_expose];
			lamb_I[*numI] = lamb_Snaive[S_expose];
			c_rate_I[*numI] = c_rate_Snaive[S_expose];
			index_I[*numI] = index_Snaive[S_expose];

			// replace data for each param in Snaive with data from last row of numS_naive Snaive indivs
			risk_Snaive[S_expose] = risk_Snaive[*numS_naive-1];
			lamb_Snaive[S_expose] = lamb_Snaive[*numS_naive-1];
			c_rate_Snaive[S_expose] = c_rate_Snaive[*numS_naive-1];
			index_Snaive[S_expose] = index_Snaive[*numS_naive-1];

			// keep track of number of Snaive and I indivs
			--*numS_naive;
			++*numI;
		}
		else		// indiv not inf, move to Sprev
		{
			// move each param from Snaive to first available row in Sprev - since there are numS_prev Sprev indivs, the first row would be at index numS_prev
			risk_Sprev[*numS_prev] = risk_Snaive[S_expose];
			lamb_Sprev[*numS_prev] = lamb_Snaive[S_expose];
			c_rate_Sprev[*numS_prev] = c_rate_Snaive[S_expose];
			index_Sprev[*numS_prev] = index_Snaive[S_expose];

			// replace data for each param in Snaive with data from last row of numS_naive Snaive indivs
			risk_Snaive[S_expose] = risk_Snaive[*numS_naive-1];
			lamb_Snaive[S_expose] = lamb_Snaive[*numS_naive-1];
			c_rate_Snaive[S_expose] = c_rate_Snaive[*numS_naive-1];
			index_Snaive[S_expose] = index_Snaive[*numS_naive-1];

			// keep track of number of Snaive and Sprev indivs
			--*numS_naive;
			++*numS_prev;
		}
	}
	else	// S indiv is previously exposed
	{
		expose_prev[index_I[I_expose]]++;	// add 1 exposed Sprev indiv for this I

		// check if S indiv is infected
		infect_rate = risk_Sprev[S_expose - *numS_naive] * lamb_I[I_expose];
		is_inf = gsl_ran_binomial(rng, 1-exp(-infect_rate), 1);		// use bernoulli to check if indiv infected
		if(is_inf == 1)		// indiv is inf, move to I and note for that I indiv
		{
			inf_prev[index_I[I_expose]]++;		// note that I indiv inf a Sprev indiv

			// move each param from Sprev to first available row in I - since there are numI I indivs, the first row would be at index numI
			risk_I[*numI] = risk_Sprev[S_expose - *numS_naive];
			lamb_I[*numI] = lamb_Sprev[S_expose - *numS_naive];
			c_rate_I[*numI] = c_rate_Sprev[S_expose - *numS_naive];
			index_I[*numI] = index_Sprev[S_expose - *numS_naive];

			// replace data for each param in Sprev with data from last row of numS_prev Sprev indivs
			risk_Sprev[S_expose - *numS_naive] = risk_Sprev[*numS_prev-1];
			lamb_Sprev[S_expose - *numS_naive] = lamb_Sprev[*numS_prev-1];
			c_rate_Sprev[S_expose - *numS_naive] = c_rate_Sprev[*numS_prev-1];
			index_Sprev[S_expose - *numS_naive] = index_Sprev[*numS_prev-1];

			// keep track of number of Sprev and I indivs
			--*numS_prev;
			++*numI;
		}
	}
}






// know event is that an indiv is inf, need to decide which one and implement
void infect(char current_pop_status[], vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, double beta_avg, int *numS, int *numI, mt19937 &rng2)
{
	int i, j;	// variables for iterations
	double prob_inf;	// probability of infection given contact
	int I_indiv;	// index of indiv that will be infected, index in terms of indivs currently S/S vector
	vector<double> infect_rate(*numS, 0);	// vector of infection rate for each contact
	
	for(i = 0; i < *numI; i++)
	{
		for(j = 0; j < *numS; j++)	// determine infection rate across all S indivs that could have contacted this indiv
		{
			infect_rate[j] += beta_avg * c_rate(j, i) * risk_S[j] * lamb_I[i];
		}
	}


	// have vector that contains infection rate for each contact
	// need to pick which indiv is inf
	
	// choose which contact is infected
	discrete_distribution<int> pick_inf(infect_rate.begin(), infect_rate.end());	// set up dist of infection rates to choose who is inf
	I_indiv = pick_inf(rng2);	// pick index of indiv who is infected

	current_pop_status[index_S[I_indiv]] = 'I';	// change status of chosen indiv to be inf

	// move each param from S to first available row in I - since there are numI I indivs, the first row would be at index numI
	risk_I[*numI] = risk_S[I_indiv];
	lamb_I[*numI] = lamb_S[I_indiv];
	c_rate_I[*numI] = c_rate_S[I_indiv];
	index_I[*numI] = index_S[I_indiv];

	// replace data for each param in S with data from last row of numS S indivs
	risk_S[I_indiv] = risk_S[*numS-1];
	lamb_S[I_indiv] = lamb_S[*numS-1];
	c_rate_S[I_indiv] = c_rate_S[*numS-1];
	index_S[I_indiv] = index_S[*numS-1];

	// keep track of number of S and I indivs
	--*numS;
	++*numI;
}




// know event is that an indiv recovers, need to decide which one and implement
void recover(char current_pop_status[], vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<int> &expose_naive, vector<int> &inf_naive, vector<int> &expose_prev, vector<int> &inf_prev, int *numI, double t, mt19937 &rng2, ofstream &ofs)
{
	int R_indiv;	// index of indiv that will recover, index in terms of indivs currently I/I vector

	// all indivs have same recovery rate gam so just have to pick number from 0 to numI-1
	uniform_int_distribution<int> pick_recover(0, *numI-1);
	R_indiv = pick_recover(rng2);	// pick index of indiv who recovers

	current_pop_status[index_I[R_indiv]] = 'R';	// change status of chosen indiv to be recovered

	// write out for this I indiv their ID number, number naive and prev exposed indivs exposed and infected, and time recover
	ofs << t << "\t" << index_I[R_indiv] << "\t" << expose_naive[index_I[R_indiv]] << "\t" << inf_naive[index_I[R_indiv]] << "\t" << expose_prev[index_I[R_indiv]] << "\t" << inf_prev[index_I[R_indiv]] << endl;

	// replace data for each param in I with data from last row of numI I indivs since indiv is no longer I
	risk_I[R_indiv] = risk_I[*numI-1];
	lamb_I[R_indiv] = lamb_I[*numI-1];
	c_rate_I[R_indiv] = c_rate_I[*numI-1];
	index_I[R_indiv] = index_I[*numI-1];

	// keep track of number of I indivs
	--*numI;
}



double c_rate(int S_index, int I_index)
{
	return(1.0);
}



// read in parameters and parse
bool parseOptions(int argc, char *argv[])
{
	for(int a = 1; a < argc; a++)
	{
		if(strcmp(argv[a], "-o") == 0 && a + 1 < argc)	// set output file name
			output = argv[++a];
		if(strcmp(argv[a], "-his") == 0)	// include heterogeneity in susceptibility given contact
			his = 'T';
		if(strcmp(argv[a], "-hit") == 0)	// include heterogeneity in transmission given contact
			hit = 'T';
		if(strcmp(argv[a], "-hic") == 0)	// include heterogeneity in contact rate
			hic = 'T';
		if(strcmp(argv[a], "-cor") == 0 && a + 1 < argc)	// set correlation between HiS and HiT
			cor = argv[++a];
		if(strcmp(argv[a], "-p") == 0 && a + 1 < argc)	// set population size
			pop_size = atoi(argv[++a]);
		if(strcmp(argv[a], "-i") == 0 && a + 1 < argc)	// set initial number of infected individuals
			init_inf = atoi(argv[++a]);
		if(strcmp(argv[a], "-r") == 0 && a + 1 < argc)	// set R0
			R0 = atof(argv[++a]);
		if(strcmp(argv[a], "-g") == 0 && a + 1 < argc)	// set recovery rate gamma
			gam = atof(argv[++a]);
		if(strcmp(argv[a], "-cs") == 0 && a + 1 < argc)	// set cv for HiS
			cv_s = atof(argv[++a]);
		if(strcmp(argv[a], "-ct") == 0 && a + 1 < argc)	// set cv for HiT
			cv_t = atof(argv[++a]);
		if(strcmp(argv[a], "-f") == 0 && a + 1 < argc)	// set scale theta for HiS
			f_inf = atof(argv[++a]);
		if(strcmp(argv[a], "-s") == 0 && a + 1 < argc)	// set seed for simulation
			seed = atoi(argv[++a]);
	}

	return(false);
}

