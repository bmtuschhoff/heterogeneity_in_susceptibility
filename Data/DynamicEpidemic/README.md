This directory (DynamicEpidemic) contains the contact tracing data output from the stochastic, individual-based SIR models and used for parameter estimation with MCMC.

File names specify relevant parameters, such as parameters dictating the underlying susceptibility distribution (his=C, finf=E, fA) and parameters important for running the SIR dynamics (g=gamma, p=pop size, R0, initI=I_0, s=seed).
The file with f_A specified is data from the discrete case, and the other file is data from the continuous case.

The contact networks from these files were used to run MCMC, but not for detection. To generate data for detection, we used a script such as below.

# compile C++ code
g++ -o ibmHiSTCTD -std=c++11 IBM_epi_model_HiST_contactTracing_discrete.cpp -lgsl -lgslcblas 

# run code to simulate epidemic
for j in $(seq 0.0 0.1 3);
do
	for k in $(seq 0.02 0.04 0.98);
	do
		rb=$(Rscript 2gps_FindrB.R $j $k $'0.5')
		for i in {1..100}
		do
			time ./ibmHiSTCTD -his -cs $j -f $k -rb $rb -fa 0.5 -r 3 -g 0.1 -cor none -p 500 -i 1 -s $i -o results_his${j}finf${k}fA0.5g0.1p500R03initI1s${i}.txt
		done
	done
done



