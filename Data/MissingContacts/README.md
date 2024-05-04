This directory (MissingContacts) contains data used for detecting heterogeneity in susceptibility, parameter estimation, and MCMC when there are missing contacts.

File names specify relevant parameters, such as the number of focal individuals (F), number of naive individuals (N), parameters dictating the underlying susceptibility distribution (k, theta), etc.

Contact_Networks contains data for the simulated contact networks. Each file consists of a table specifying each individual's number of previous exposures (num_exps), risk of infection (risks), and whether they were infected or not (inf). These were used to run MCMC.
"naive" is for individuals observed to be naive and "prev" is for individuals observed to be focal.

Detection contains simulations for the power to detect heterogeneity in susceptibility when there are missing contacts. These data were used to generate Figs S25 and S26.
