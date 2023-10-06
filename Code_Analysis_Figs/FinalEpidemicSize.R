### Check final epidemic size for specific expected fraction infected across different CVs of risk and fAs
### Test does bigger CV (more het in susc) mean smaller epidemic?
# Fig S9

library(ggplot2)
library(dplyr)

# parameter grid to test
cvs <- seq(0.0,3,0.05)
fAs <- seq(0.05,0.95,0.05)
f_inf <- 0.25

# set the seed
set.seed(3)

# function to find r_b with root finding
find_r_b <- function(r_b, f_a, cv, f_inf)
{
  abs(((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))) - r_b)*sqrt(f_a*(1 - f_a))) / ((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))))*f_a + r_b*(1 - f_a)) - cv)
}


# function to run SIR model and return data frame of number of individuals in each class over time
SIR <- function(N,fA,I_init_a,I_init_b,beta_a,beta_b,gamma,mu,EndTime,deltat)
{
  Sa <- numeric(EndTime/deltat +1)
  Sb <- numeric(EndTime/deltat +1)
  Ia <- numeric(EndTime/deltat +1)
  Ib <- numeric(EndTime/deltat +1)
  Ra <- numeric(EndTime/deltat +1)
  Rb <- numeric(EndTime/deltat +1)
  
  print(N)
  
  Sa[1] <- round((N - I_init_a - I_init_b)*fA)
  Sb[1] <- (N - Sa[1] - I_init_a - I_init_b)
  
  Ia[1] <- I_init_a
  Ib[1] <- I_init_b
  
  Ra[1] <- 0
  Rb[1] <- 0
  
  
  t <- seq(0, EndTime, deltat)
  
  for(i in 1:(EndTime/deltat))
  {
    Sa[i+1] <- Sa[i] + (-beta_a*Sa[i]*(Ia[i] + Ib[i]))*deltat
    Sb[i+1] <- Sb[i] + (-beta_b*Sb[i]*(Ib[i] + Ia[i]))*deltat
    
    Ia[i+1] <- Ia[i] + (beta_a*Sa[i]*(Ia[i] + Ib[i]) - gamma*Ia[i])*deltat
    Ib[i+1] <- Ib[i] + (beta_b*Sb[i]*(Ib[i] + Ia[i]) - gamma*Ib[i])*deltat
    
    Ra[i+1] <- Ra[i] + (gamma*Ia[i])*deltat
    Rb[i+1] <- Rb[i] + (gamma*Ib[i])*deltat
  }
  
  return(data.frame(Sa,Sb,Ia,Ib,Ra,Rb))
}

# define parameters
N <- 20010  # pop size
I_init <- 10  # initial number infected
#c <- 0.00005  # contact rate
R0 <- 3.0
gam <- 0.1  # recovery rate
mu <- 0  # death rate
EndTime <- 150  # total time
deltat <- 0.001  # time step size


# find final epidemic size across various param combos
epi_table <- matrix(nrow=length(cvs)*length(fAs), ncol=5)
n <- 1

for(cv in cvs)
{
  for(f_a_old in fAs)
  {
    if(cv == 0)
    {
      pA <- f_inf
      pB <- f_inf
    } else
    {
      # find risk values numerically solving for root
      if(f_inf < f_a_old)
      {
        if(f_inf < 1-f_a_old)
        {
          rb <- optimize(find_r_b, interval=c(0,-log(1 - (f_inf / (1 - f_a_old)))), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        } else
        {
          rb <- optimize(find_r_b, interval=c(0,50), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        }
      } else
      {
        if(f_inf < 1-f_a_old)
        {
          rb <- optimize(find_r_b, interval=c(-log(1 - ((f_inf - f_a_old) / (1 - f_a_old))),-log(1 - (f_inf / (1 - f_a_old)))), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        } else
        {
          rb <- optimize(find_r_b, interval=c(-log(1 - ((f_inf - f_a_old) / (1 - f_a_old))),50), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
        }
      }
      ra <- -log(1 - (1/f_a_old)*(f_inf - (1 - exp(-rb))*(1 - f_a_old)))
      # back calculate for p_a and p_b
      pA <- 1 - exp(-ra)
      pB <- 1 - exp(-rb)
    }

    # run SIR dynamics to get final epidemic size
    # parameters for SIR model
    p <- pA*f_a_old + pB*(1 - f_a_old)
    c <- (gam * R0) / (p * (N-I_init))

    # run SIR model and calculate final epidemic size
    res <- SIR(N,f_a_old,I_init_a,I_init_b,beta_a=c*pA,beta_b=c*pB,gam,mu,EndTime,deltat)
    final_size <- 1 - (res$Sa[EndTime/deltat+1] + res$Sb[EndTime/deltat+1]) / (res$Sa[1] + res$Sb[1])
    
    # save results
    epi_table[n,] <- c(pA, pB, f_a_old, cv, final_size)
    n <- n + 1
  }
}

##############################################################################
# plot results
epi_table <- data.frame(epi_table)
colnames(epi_table) <- c("pA", "pB", "fA", "CV", "final_size")

mybreaks <- seq(0.2,0.95,length.out=11)
mycolors <- function(x)
{
  colors <- colorRampPalette(c("darkblue","red"))( 10 )
  colors[1:x]
}


#### remove all duplicated values for p_a and p_b (some vals of efi and cv give same probs)
epi_table_uniq <- epi_table %>%
  mutate(Index=row_number()) %>%
  filter(!duplicated(pA))

## add back cv=0 to epi_table_uniq
epi_table_uniq <- rbind(epi_table[(1:19),], epi_table_uniq[-1,-6])

### find maximum cv for each f_a of unique values in order to plot line on plots
epi_maxC <- epi_table_uniq %>%
  group_by(fA) %>%
  summarise_at(vars(CV),
               list(max = max))

par(mar=c(5.1,5.1,4.1,2.1), pty='s') # change margin sizes

# create plot
ggplot(epi_table, aes(x=fA, y=CV, z=final_size)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.05,0.95), breaks=c(0.05,0.25,0.5,0.75,0.95)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F) +
  theme_minimal() +
  labs(fill="Final Epidemic Size") +
  xlab(expression(paste("Fraction of the population that is more susceptible (",italic(f)[italic(A)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_text(), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=0)) +
  geom_line(data=epi_maxC, aes(x=fA, y=max), col="gray", lty="dashed")
  #geom_hline(yintercept=1.3, color="black", cex=1.5) +
  #geom_point(data=results, aes(x=f_a, y=cv, color=power), size=2) + ## next 3 lines=add real data pts to check with generic binom
  #scale_color_gradientn(breaks=mybreaks, colors=mycolors(10)) +
  #geom_point(data=results[sample((1:nrow(results)), 500, replace=F),], aes(x=f_a, y=cv), shape=1, size=2, color="gray")


# 10.53 x 7.73 pdf
