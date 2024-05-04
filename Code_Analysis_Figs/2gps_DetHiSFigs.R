# Plot power to detect het in susc for discrete case (2 types), CV vs expected fraction inf
# Figs 3, S1, S2, S11, S14, S21

library(ggplot2)
library(ggpubr)
library(margins)
library(dplyr)

# set up color palette
mybreaks <- seq(0,100,10)
mycolors <- function(x)
{
  colors <- colorRampPalette(c("darkblue","yellow"))( 10 )
  colors[1:x]
}


# read in power table
fileName01 <- "2gpPowers_fA0.5_F50_N5_N25_sims1000sdseed.csv"
fileName01 <- "2gpPowersFalseNeg12adj_fA0.5_F200_N5_N25_fneg0.1_sims1000seed.csv"
fileName01 <- "2gpPowersFalseNeg12_fA0.5_F200_N5_N25_fneg0.1_sims1000seed.csv"
fileName01 <- "2gpPowersHiT_fA0.5_F200_N5_N25_sims1000HiTm0.5phi2cv03e0.98Lambseed.csv"
powers_table <- read.csv(fileName01, header=T)
powers_table <- powers_table[-seq(50,7550,50),]
powers_table <- powers_table[,-1]
colnames(powers_table) <- c("cv", "f_inf", "power")
colnames(powers_table) <- c("f_inf", "cv", "power")
colnames(powers_table) <- c("f_inf", "cv", "power", "p_avg", "p_naive", "p_prev")
colnames(powers_table) <- c("f_inf", "cv", "power", "cv_ratio", "p_avg", "sd_p_avg", "p_naive", "sd_p_naive", "p_prev", "sd_p_prev")
powers_table$power[which(powers_table$power == 100)] <- 99.999
powers_table$cv <- as.numeric(powers_table$cv)
powers_table$f_inf <- as.numeric(powers_table$f_inf)
powers_table$power <- as.numeric(powers_table$power)
powers_table$p_naive <- as.numeric(powers_table$p_naive)
powers_table$p_prev <- as.numeric(powers_table$p_prev)
powers_table$p_avg <- as.numeric(powers_table$p_avg)
powers_table$cv_ratio <- as.numeric(powers_table$cv_ratio)


###################### find area where get computationally indistinguishable p_A and p_B
######## calculate p_A and p_B from table
# function to find r_b with root finding
find_r_b <- function(r_b, f_a, cv, f_inf)
{
  abs(((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))) - r_b)*sqrt(f_a*(1 - f_a))) / ((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))))*f_a + r_b*(1 - f_a)) - cv)
}

# function to find p_A, p_B
find_p_a_p_b <- function(cv, f_inf, f_a)
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
  p_a <- 1 - exp(-ra)
  p_b <- 1 - exp(-rb)
  
  return(c(p_a, p_b))
}

# set f_A
f_a_old <- 0.5

# calculate p_A and p_B
for(r in 1:(nrow(powers_table)))
{
  pApB_res <- find_p_a_p_b(powers_table$cv[r], powers_table$f_inf[r], f_a_old)
  powers_table$p_a[r] <- pApB_res[1]
  powers_table$p_b[r] <- pApB_res[2]
}

##### create table with unique values by removing repeated p_A's
powers_table_uniq <- powers_table %>%
  mutate(Index=row_number()) %>%
  filter(!duplicated(p_a))

### find maximum cv for each f_inf of unique values in order to plot line on plots for where computationally indistinguishable
powers_maxC <- powers_table_uniq %>%
  group_by(f_inf) %>%
  summarise_at(vars(cv),
               list(max = max))

#####################################################################
## make plot for paper
## main text, Fig 3, changing F and f_A

par(mar=c(5.1,5.1,4.1,2.1)) # change margin sizes

fA0.1_F50 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98), position="top") +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab("") +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_blank(), axis.ticks.length.x=unit(1.3,"cm"), axis.ticks.y=element_line(),
        plot.margin=margin(t=8, r=-16, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")
  #annotate("text", x=c(0.5,-0.5), y=c(3.5, 1.5), label=c("0.1","50"), size=8) +
  #coord_cartesian(xlim=c(0,1), ylim=c(0,3), clip="off")

fA0.5_F50 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), position="top", limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab("") +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_blank(), axis.ticks.length.x=unit(1.3,"cm"),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"),
        plot.margin=margin(t=8, r=-8, b=8, l=-8)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.9_F50 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), position="top", limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab("") +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_blank(), axis.ticks.length.x=unit(1.3,"cm"),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"),
        plot.margin=margin(t=8, r=0, b=8, l=-16)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.1_F200 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-16, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.5_F200 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=-8, b=8, l=-8)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.9_F200 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=-16)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")


annF50 <- ggplot() +
  geom_text(aes(x=0, y=0, label="50"),
            parse=TRUE, size=8, vjust=0.5, hjust=1.65, angle=90) +
  theme_void()

annF200 <- ggplot() +
  geom_text(aes(x=0, y=0, label="200"),
            parse=TRUE, size=8, vjust=0.5, hjust=-0.4, angle=90) +
  theme_void()

annfA0.1 <- ggplot() +
  geom_text(aes(x=0, y=0, label="0.1"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=90, r=-16, b=-90, l=16)) #hjust=-0.6

annfA0.5 <- ggplot() +
  geom_text(aes(x=0, y=0, label="0.5"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=90, r=-8, b=-90, l=8)) #hjust=-0.4

annfA0.9 <- ggplot() +
  geom_text(aes(x=0, y=0, label="0.9"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=90, r=0, b=-90, l=0)) #hjust=-0.3

blank <- ggplot(hjust=-0.4) + theme_void()

ggarrange(blank, annfA0.1, annfA0.5, annfA0.9,
          annF50, fA0.1_F50, fA0.5_F50, fA0.9_F50, 
          annF200, fA0.1_F200, fA0.5_F200, fA0.9_F200, 
          nrow=3, ncol=4, common.legend=T, legend="right",
          widths=c(0.1,0.3,0.3,0.3),
          heights=c(0.1,0.45,0.45))

plt <- ggarrange(blank, annfA0.1, annfA0.5, annfA0.9,
                 annF50, fA0.1_F50, fA0.5_F50, fA0.9_F50, 
                 annF200, fA0.1_F200, fA0.5_F200, fA0.9_F200, 
                 nrow=3, ncol=4, common.legend=T, legend="right",
                 widths=c(0.1,0.3,0.3,0.3),
                 heights=c(0.1,0.45,0.45))

annotate_figure(plt, top=text_grob(expression(paste("Fraction of population that is more susceptible (",italic(f)[italic(A)],")")), size=25, hjust=0.47, vjust=3.4), 
                left=text_grob(expression(paste("Number of focal individuals (",italic(F),")")), size=25, vjust=1, hjust=0.65, rot=90))
# 16 x 10.5 pdf


#################################################################################
### make plot for paper, effect of changing N, Figs S1 and S2
# N=5 vs N=100
# F=200, f_A=0.5

par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

fA0.5_F200N5 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-10, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.5_F200N100 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=-10)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

annN5 <- ggplot() +
  geom_text(aes(x=0, y=0, label="5"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=10, r=-24, b=-45, l=8)) #hjust=-0.6

annN100 <- ggplot() +
  geom_text(aes(x=0, y=0, label="100"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=10, r=0, b=-45, l=8)) #hjust=-0.4

ggarrange(annN5, annN100,
          fA0.5_F200N5, fA0.5_F200N100,
          nrow=2, ncol=2, common.legend=T, legend="right",
          heights=c(0.1,0.9))
  
plt <- ggarrange(annN5, annN100,
                 fA0.5_F200N5, fA0.5_F200N100,
                 nrow=2, ncol=2, common.legend=T, legend="right",
                 heights=c(0.1,0.9))

annotate_figure(plt, top=text_grob(expression(paste("Number of contacts per network (",italic(N),")")), size=25, hjust=0.55, vjust=1.6))

# 11 x 7 pdf



## find dif in power w/ N=5, N=100
powers_tableN5 <- powers_table
#powers_tableN100 <- powers_table
powers_table_dif <- powers_tableN5[,(1:3)]
powers_table_dif$power <- powers_tableN100$power - powers_tableN5$power

# change color scale, red=negative=less power w/ N=100, blue=positive=more power w/ N=100
mybreaks <- c(seq(-15,-5,5),seq(5,30,5))
mycolors_neg <- function(x)
{
  colors <- colorRampPalette(c("red","white"))( 3 )
  colors[1:x]
}
mycolors_pos <- function(x)
{
  colors <- colorRampPalette(c("white","blue"))( 6 )
  colors[2:x]
}

par(mar=c(5.1,5.1,4.1,4.1))

ggplot(powers_table_dif, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=c(mycolors_neg(3), mycolors_pos(6)), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=8, b=8, l=8), panel.border=element_rect(fill = NA, size = 1), 
        legend.key=element_rect()) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

# 8.39 x 7 pdf


#################################################################################
# make plot for paper
# effect of false negatives, Fig S11
# 1 x 3 with normal, false negs, false negs adjusted method
# F=200, N=5, fA=0.5, beta=0.1
par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

fA0.5_Norm <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=30), legend.text=element_text(size=20), legend.title = element_text(size=30),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-8, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.5_FN <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=30), legend.text=element_text(size=20), legend.title = element_text(size=30),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks.x=element_line(),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"),
        plot.margin=margin(t=8, r=0, b=8, l=2)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

fA0.5_FNadj <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=30), legend.text=element_text(size=20), legend.title = element_text(size=30),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks.x=element_line(),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"),
        plot.margin=margin(t=8, r=0, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

annNorm <- ggplot() +
  geom_text(aes(x=0, y=0, label="'No false negatives'"),
            parse=TRUE, size=10) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=30), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=0, r=-16, b=-60, l=16)) #hjust=-0.6

annFN <- ggplot() +
  geom_text(aes(x=0, y=0, label="'False negatives'"),
            parse=TRUE, size=10) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=30), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=0, r=-8, b=-60, l=8)) #hjust=-0.4

annFNadj <- ggplot() +
  geom_text(aes(x=0, y=0, label="'False negatives adjusted'"),
            parse=TRUE, size=10) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=30), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=0, r=-8, b=-60, l=8)) #hjust=-0.4


ggarrange(annNorm, annFN, annFNadj,
          fA0.5_Norm, fA0.5_FN, fA0.5_FNadj, 
          nrow=2, ncol=3, common.legend=T, legend="right",
          widths=c(1/3,1/3,1/3),
          heights=c(0.1,0.9))

#25in x 10in pdf


#################################################################################
## make plot for paper
## effect of HiT, Fig S14
## F=200, N=5, fA=0.5, m=0.5, phi=2
par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

fA0.5_F200NHiT <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-16, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")


fA0.5_F200HiT <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=-16)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")

annNoHiT <- ggplot() +
  geom_text(aes(x=0, y=0, label="'Homogeneity in transmission'"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=10, r=-16, b=-45, l=10)) #hjust=-0.6

annHiT <- ggplot() +
  geom_text(aes(x=0, y=0, label="'Heterogeneity in transmission'"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=10, r=0, b=-45, l=0)) #hjust=-0.4


ggarrange(annNoHiT, annHiT,
          fA0.5_F200NHiT, fA0.5_F200HiT,
          nrow=2, ncol=2, common.legend=T, legend="right",
          heights=c(0.1,0.9))

# 12 x 7, pdf


##################################################################################
## make plot for paper
## power with contact networks from dynamic epidemic vs static, Fig S21
## F=50, N=5, fA=0.5
par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

fA0.5_F50_stat <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.14,0.86), breaks=c(0.14,0.3,0.5,0.7,0.86)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(d)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-10, b=8, l=0)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")


fA0.5_F50_dynam <- ggplot(powers_table_disc_dynamepi, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.14,0.86), breaks=c(0.14,0.3,0.5,0.7,0.86)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(d)],")"))) +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=-10)) +
  geom_line(data=powers_maxC, aes(x=f_inf, y=max), col="gray", lty="dashed")



annStat <- ggplot() +
  geom_text(aes(x=0, y=0, label="'Static'"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=10, r=-16, b=-45, l=10)) #hjust=-0.6

annDynam <- ggplot() +
  geom_text(aes(x=0, y=0, label="'Dynamic'"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=20), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=10, r=0, b=-45, l=0)) #hjust=-0.4


ggarrange(annStat, annDynam,
          fA0.5_F50_stat, fA0.5_F50_dynam,
          nrow=2, ncol=2, common.legend=T, legend="right",
          heights=c(0.1,0.9))

# 12 x 7, pdf