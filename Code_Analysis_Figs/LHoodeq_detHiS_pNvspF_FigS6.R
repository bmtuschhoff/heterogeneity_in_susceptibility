## use generic binomial to investigate the likelihood equations used for detection (Eqs 1, 2 in main text) and how they affect the power to detect heterogeneity in susceptibility
# Fig S6

library(ggplot2)
library(dplyr)

# pop size
n <- 50
# number of simulations
sims <- 10000
# seq of 1st probs
p1s <- seq(0,1,0.01)

prob_L_power_table <- matrix(NA, nrow=length(p1s)^2, ncol=13)

set.seed(3)
i <- 1
# calc log likelihoods for generic binom with dif probs p1, p2
for(p1 in p1s)
{
  # set 2nd prob p2 to always be <= to 1st p1
  p2s <- seq(0,p1,0.01)
  
  for(p2 in p2s)
  {
    # simulate rand nums for p1, p2 with rbinom
    rand_nums1 <- rbinom(sims, n*4, p1)
    rand_nums2 <- rbinom(sims, n, p2)
    
    # log likelihoods for dif probs p1, p2  
    L1 <- dbinom(rand_nums1, n*4, (rand_nums1+rand_nums2)/(n*5), log=T) + dbinom(rand_nums2, n, (rand_nums1+rand_nums2)/(n*5), log=T)
    L2 <- dbinom(rand_nums1, n*4, rand_nums1/(n*4), log=T) + dbinom(rand_nums2, n, rand_nums2/n, log=T)
      
    # log likelihood ratio test statistics, vector len=sims
    ratios <- -2 * (L1 - L2)
    
    # find power = percent of times the ratio test is greater than crit val 3.84
    power <- sum(ratios > 3.84, na.rm=T)/length(ratios) * 100
    
    prob_L_power_table[i,] <- c(p1, p2, sd(L1), mean(L1), sd(L2), mean(L2), sd(ratios), mean(ratios), sd(rand_nums1), mean(rand_nums1), sd(rand_nums2), mean(rand_nums2), power)
    i <- i + 1
  }
}

prob_L_power_table <- data.frame(prob_L_power_table)
colnames(prob_L_power_table) <- c("p1", "p2", "sd_L1", "m_L1", "sd_L2", "m_L2", "sd_ratio", "m_ratio", "sd_nums1", "m_nums1", "sd_nums2", "m_nums2", "power")
prob_L_power_table <- prob_L_power_table[1:5151,]


# plot heat map
library(ggplot2)
library(dplyr)

mybreaks <- seq(0,100,10)
mycolors <- function(x)
{
  colors <- colorRampPalette(c("darkblue","yellow"))( 10 )
  colors[1:x]
}


prob_L_power_table$power[which(prob_L_power_table$power == 100)] <- 99.9

labels <- data.frame(p1=seq(0.5,1,0.05), p2=seq(0.5,0,-0.05))

# read in data to set boundaries for what params are possible in discrete case
fileName <- "2gpPowers_fA0.5_F50_N5_N25_sims1000sdseed.csv"
powers_table <- read.csv(fileName, header=T)
powers_table <- powers_table[-seq(50,7550,50),]
powers_table <- powers_table[,-1]
colnames(powers_table) <- c("f_inf", "cv", "power", "cv_ratio", "p_avg", "sd_p_avg", "p_naive", "sd_p_naive", "p_prev", "sd_p_prev")
powers_table$power[which(powers_table$power == 100)] <- 99.999
powers_table$cv <- as.numeric(powers_table$cv)
powers_table$f_inf <- as.numeric(powers_table$f_inf)
powers_table$power <- as.numeric(powers_table$power)
powers_table$p_naive <- as.numeric(powers_table$p_naive)
powers_table$p_prev <- as.numeric(powers_table$p_prev)
powers_table$p_avg <- as.numeric(powers_table$p_avg)

# find minimum vals for real data to set boundary for possible params
powers_min <- powers_table %>%
  mutate(across("p_naive", round, 2)) %>%
  group_by(p_naive) %>%
  summarise_at(vars(p_prev),
               list(min = min))

# plot
par(mar=c(5.1,5.1,4.1,2.1), pty="s")
ggplot(prob_L_power_table, aes(x=p1, y=p2, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  scale_fill_manual(values=mycolors(10), drop=F) +
  theme_minimal() +
  labs(fill="Power") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        text=element_text(size=20), axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_text(hjust=0.5), panel.border = element_rect(fill = NA, size = 1),
        axis.ticks=element_line()) +
  xlab(expression(paste("Probability of infection for naive individuals (",italic(p)[italic(n)],")"))) +
  ylab(expression(paste("Probability of infection for focal individuals (",italic(p)[italic(f)],")"))) +
  stat_smooth(data=powers_min, aes(x=p_naive, y=min), color="gray", method="lm", formula=y~poly(x,9), se=F) +
  geom_abline(intercept=-seq(0,1,0.1), slope=1, linetype="dashed", color="red") + ## lines for deltap=p1-p2=0,0.1,0.2,...
  geom_text(data=labels, aes(x=p1, y=p2, label=paste("\u0394 p=", abs(p1-p2))), angle=45, col="red", size=5) +
  geom_segment(aes(x = 0.4, y = 0.45, xend = 0.1, yend = 0.15), arrow = arrow(length = unit(0.5, "cm"))) + # add arrows to plot to describe how power increases
  geom_text(aes(x=0.25, y=0.33, label="Power increases"), angle=45, size=5, check_overlap=T) +
  geom_segment(aes(x = 0.6, y = 0.65, xend = 0.9, yend = 0.95), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_text(aes(x=0.75, y=0.83, label="Power increases"), angle=45, size=5, check_overlap=T) +
  geom_segment(aes(x = 0.1, y = 0.9, xend = 0.4, yend = 0.6), arrow = arrow(length = unit(0.5, "cm")), cex=3) +
  geom_text(aes(x=0.25, y=0.78, label="Power increases"), angle=-45, size=5, check_overlap=T)
  
# 9.6 x 7.73 pdf
