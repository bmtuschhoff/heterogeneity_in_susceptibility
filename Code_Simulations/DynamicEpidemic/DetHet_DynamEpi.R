# read in all files for dynamic epidemic to create power plot
# need to read in each file (100 sims per param combo), get number of naive and focal indivs exposed and infected at F=50
# save that info so 100 for that param combo and get power of detection for that param combo
# then save power for that param combo and then test for next param combo
# want to end up with table of cv, f_inf, power

# number of simulations using to calculate power
sims <- 100

# set working directory to get data from
setwd("C:/Users/bmt5507/Documents/Power_Cont_DynamEpi")

# function to read in datasets that are tab delimited
read_sav <- function(fileName) {read.table(fileName, header=T, sep="\t")}

# function to get numbers of indivs exposed and infected at F=50 for each data file
num_expinf <- function(contact_data)
{
  nets_notfull <- which(contact_data$Naive.exposed == 0 | contact_data$Focal.exposed == 0) # networks that don't have naive and focal indiv exposed

  if(!identical(nets_notfull, integer(0)))  # if there are networks that don't have naive and focal indivs, remove them
    contact_data <- contact_data[-nets_notfull,] # remove networks that don't have naive and focal indivs
  contact_data <- cumsum(contact_data[,-(1:2)])  # save the cumulative sum of all indivs exposed and infected
  
  Fexp_50 <- which(contact_data$Focal.exposed >= 50)  # cumulative sum of indivs up to F=50
  
  if(identical(Fexp_50, integer(0)))  # if no epidemic/data available, save as NAs so easier to deal with later
    return(rep(NA, 4))
  else
    return(contact_data[Fexp_50[1],])  # info wanted of number naive and focal indivs exp and inf with largest F, F<=50
}


cvs <- seq(0.0,3.0,0.1)
f_infs <- seq(0.02,0.98,0.04)

powers_table <- matrix(NA, nrow=length(cvs)*length(f_infs), ncol=5)
combo <- 1

for(cv in cvs)
{
  for(f_inf in f_infs)
  {
    # list of file names
    filePattern <- paste("results_his",format(cv,nsmall=1),"finf",format(f_inf,nsmall=2),"g0.1p500R03initI1s*", sep="")
    files <- list.files(pattern=filePattern)
    
    # make a list that contains all files
    read.all <- lapply(files, read_sav)
    
    # make data frame containing data from all 100 sims of this param combo for number indivs exposed and infected
    expinf_data <- data.frame(t(as.matrix(sapply(read.all, num_expinf))))
    
    # make data frame all numeric and make sure column names are correct
    expinf_data <- data.frame(sapply(expinf_data, as.numeric))
    colnames(expinf_data) <- c("Naive.exposed", "Naive.infected", "Focal.exposed", "Focal.infected")
    
    # remove all sims with NA values
    expinf_data <- na.omit(expinf_data)
    
    ##### find power for this param combo
    # MLE probabilities of infection
    pnaive <- expinf_data$Naive.infected / expinf_data$Naive.exposed
    pfocal <- expinf_data$Focal.infected / expinf_data$Focal.exposed
    pavg <- (expinf_data$Naive.infected + expinf_data$Focal.infected) / (expinf_data$Naive.exposed + expinf_data$Focal.exposed)
    
    # log-likelihoods for homogeneity vs heterogeneity
    Lhom <- dbinom(expinf_data$Naive.infected, expinf_data$Naive.exposed, pavg, log=T) + dbinom(expinf_data$Focal.infected, expinf_data$Focal.exposed, pavg, log=T)
    Lhet <- dbinom(expinf_data$Naive.infected, expinf_data$Naive.exposed, pnaive, log=T) + dbinom(expinf_data$Focal.infected, expinf_data$Focal.exposed, pfocal, log=T)
    
    # calculate likelihood ratio test statistics
    ratios <- -2*(Lhom - Lhet)
    
    # find power = percent of times the ratio test is greater than critical value 3.84 (df=1)
    power <- 100 * sum(ratios > 3.84) / nrow(expinf_data)
    
    # save cv, f_inf, power for this param combo
    powers_table[combo,] <- c(cv, f_inf, power, mean(expinf_data$Focal.exposed), nrow(expinf_data))
    combo <- combo + 1
  }
}



### make power plot and save powers table, paper figures actually made in 2gps_DetHiSFigs.R and Cont_DetHiSFigs.R
library(ggplot2)
library(ggpubr)
library(margins)

mybreaks <- seq(0,100,10)
mycolors <- function(x)
{
  colors <- colorRampPalette(c("darkblue","yellow"))( 10 )
  colors[1:x]
}

powers_table <- data.frame(powers_table)
colnames(powers_table) <- c("cv", "f_inf", "power", "avg_F", "num_sims")
powers_table$power[which(powers_table$power == 100)] <- 99.999

powers_table_cont_dynamepi <- powers_table

powers_table_cont_dynamepi[which(powers_table_cont_dynamepi$f_inf == 0.1)+1,]$f_inf <- rep(0.14, 31)

par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

ggplot(powers_table_cont_dynamepi, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.14,0.86), breaks=c(0.14,0.25,0.5,0.75,0.86)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(c)],")"))) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-10, b=8, l=0))


