## Clean data from dynamic epidemic and output datasets to be used for parameter estimation with MCMC
## Remove contact networks that don't have both naive and focal individuals exposed
## Save datasets with different sample sizes (different numbers of focal individuals F exposed)

# read in contact tracing data from dynamic epidemic
fileName <- "results_his1.3finf0.25g0.1p10000R03initI10s3time.txt"
contact_data <- read.table(fileName, header=T, sep="\t")

# remove all contact networks where no one is exposed and there is not at least 1 naive and 1 focal individual
contact_data <- contact_data[-which(contact_data$Naive.exposed == 0 | contact_data$Focal.exposed == 0),]

# cumulative sum of focal individuals exposed - to find first point where there are 50, 200, 1000, 5000 focal indivs
Fexp_cumsum <- cumsum(contact_data$Focal.exposed)

# index of first instance where at least 50,200,... focal indivs exp (F=50,200,...)
which(Fexp_cumsum >= 50)[1]
Fexp_cumsum[which(Fexp_cumsum >= 50)[1]]
which(Fexp_cumsum >= 200)[1]
Fexp_cumsum[which(Fexp_cumsum >= 200)[1]]
which(Fexp_cumsum >= 1000)[1]
Fexp_cumsum[which(Fexp_cumsum >= 1000)[1]]
which(Fexp_cumsum >= 5000)[1]
Fexp_cumsum[which(Fexp_cumsum >= 5000)[1]]

# save datasets that contain contact networks up to first instance where at least 50, 200, 1000, or 5000 focal individuals are exposed
contact_data_F50 <- contact_data[1:which(Fexp_cumsum >= 50)[1],]
contact_data_F200 <- contact_data[1:which(Fexp_cumsum >= 200)[1],]
contact_data_F1000 <- contact_data[1:which(Fexp_cumsum >= 1000)[1],]
contact_data_F5000 <- contact_data[1:which(Fexp_cumsum >= 5000)[1],]

# write new datasets out to a file in order to run MCMC with them
write.csv(contact_data_F5000, "contact_data_dynamepi_pA0.74809046024617_pB0.125477384938458_c1.3e0.25fA0.2_F5000.csv")
