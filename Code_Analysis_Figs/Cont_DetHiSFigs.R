# Plot power to detect het in susc for continuous case, CV vs expected fraction inf
# Figs 3, S3, S4, S12, S15

library(ggplot2)
library(ggpubr)
library(margins)

# set up color palette
mybreaks <- seq(0,100,10)
mycolors <- function(x)
{
  colors <- colorRampPalette(c("darkblue","yellow"))( 10 )
  colors[1:x]
}

# read in power table
fileName01 <- "ContPowers_F50_N5_N25_sims1000_c03e0.98v5.csv"
fileName01 <- "ContPowers_FI1_FI21_F50_N100_N25_sims100_cv1.483_efi0.020.9v5.csv"
fileName01 <- "ContPowersFNeg12adj_F200_N5_N25_sims1000_fneg0.1c3e0.98seed.csv"
fileName01 <- "ContPowersFNeg12_F200_N5_N25_sims1000_fneg0.1_c3e0.98seed.csv"
fileName01 <- "ContPowers_FI1_FI21_F1000_N100_N25_sims100_cv03e0.9v5.csv"
fileName01 <- "ContPowersHiT_F200_N5_N25_sims1000_c03_e0.98_HiTm0.5phi2Lambseed.csv"
powers_table <- read.csv(fileName01, header=T)
powers_table <- powers_table[-seq(50,7550,50),]
powers_table <- powers_table[,-1]
colnames(powers_table) <- c("cv", "f_inf", "power")
colnames(powers_table) <- c("cv", "f_inf", "power", "p_naive", "p_prev")
colnames(powers_table) <- c("cv", "f_inf", "power", "p_naive", "p_prev", "mean_risk_pre", "mean_risk_post")
powers_table$power[which(powers_table$power == 100)] <- 99.999
powers_table$cv <- as.numeric(powers_table$cv)
powers_table$f_inf <- as.numeric(powers_table$f_inf)
powers_table$power <- as.numeric(powers_table$power)
powers_table$p_naive <- as.numeric(powers_table$p_naive)
powers_table$p_prev <- as.numeric(powers_table$p_prev)


powers_tableN100 <- powers_table

###############################################################################
## make plot for paper
## main text, Fig 3, changing F

par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

cont_F50 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(c)],")"))) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=20), legend.title = element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-8, b=8, l=0))

cont_F200 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=20), legend.title = element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=2))

cont_F1000 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=20), legend.title = element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=0))

annF50 <- ggplot() +
  geom_text(aes(x=0, y=0, label="50"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=30), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=0, r=-8, b=-50, l=16)) #hjust=-0.6

annF200 <- ggplot() +
  geom_text(aes(x=0, y=0, label="200"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=30), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=0, r=-8, b=-50, l=8)) #hjust=-0.4

annF1000 <- ggplot() +
  geom_text(aes(x=0, y=0, label="1000"),
            parse=TRUE, size=8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        text=element_text(size=30), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        plot.margin=margin(t=0, r=-6, b=-50, l=8)) #hjust=-0.4

ggarrange(annF50, annF200, annF1000,
          cont_F50, cont_F200, cont_F1000, 
          nrow=2, ncol=3, common.legend=T, legend="right",
          widths=c(1/3,1/3,1/3),
          heights=c(0.1,0.9))

plt <- ggarrange(annF50, annF200, annF1000,
                 cont_F50, cont_F200, cont_F1000, 
                 nrow=2, ncol=3, common.legend=T, legend="right",
                 widths=c(1/3,1/3,1/3),
                 heights=c(0.1,0.9))


annotate_figure(plt, top=text_grob(expression(paste("Number of focal individuals (",italic(F),")")), size=25, hjust=0.6, vjust=0.8))

# 15 x 7 pdf

###############################################################################
## make plot for paper, changing N, Figs S3 and S4
## F=200, N=5,100
par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

cont_N5 <- ggplot(powers_tableN5, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
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

cont_N100 <- ggplot(powers_tableN100, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=-10))


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
          cont_N5, cont_N100, 
          nrow=2, ncol=2, common.legend=T, legend="right",
          heights=c(0.1,0.9))

plt <- ggarrange(annN5, annN100,
                 cont_N5, cont_N100, 
                 nrow=2, ncol=2, common.legend=T, legend="right",
                 heights=c(0.1,0.9))

annotate_figure(plt, top=text_grob(expression(paste("Number of contacts per network (",italic(N),")")), size=25, hjust=0.55, vjust=1.6))

# 11 x 7 pdf



## find dif in power w/ N=5, N=100
powers_table_dif <- powers_tableN5[,(1:3)]
powers_table_dif$power <- powers_tableN100$power - powers_tableN5$power
# 1.24 for N1=10

# change color scale, red=negative=less power w/ N=100, blue=positive=more power w/ N=100
mybreaks <- c(-5,seq(5,20,5))
mycolors_neg <- function(x)
{
  colors <- colorRampPalette(c("white"))( 1 )
  colors[1:x]
}
mycolors_pos <- function(x)
{
  colors <- colorRampPalette(c("white","blue"))( 4 )
  colors[2:x]
}

ggplot(powers_table_dif, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=c(mycolors_neg(1), mycolors_pos(4)), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(c)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=0), legend.key=element_rect())

# 8.39 x 7 pdf

#################################################################################
# make plot for paper
# effect of false negatives, Fig S13
# 1 x 3 with normal, false negs, false negs adj
# F=200, N=5
par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

cont_F200 <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(c)],")"))) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=30), legend.text=element_text(size=20), legend.title = element_text(size=30),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-8, b=8, l=0))

cont_F200_FN <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=30), legend.text=element_text(size=20), legend.title = element_text(size=30),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks.x=element_line(),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"),
        plot.margin=margin(t=8, r=0, b=8, l=2))

cont_F200_FNadj <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(" ") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=30), legend.text=element_text(size=20), legend.title = element_text(size=30),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks.x=element_line(),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"),
        plot.margin=margin(t=8, r=0, b=8, l=0))

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
          cont_F200, cont_F200_FN, cont_F200_FNadj, 
          nrow=2, ncol=3, common.legend=T, legend="right",
          widths=c(1/3,1/3,1/3),
          heights=c(0.1,0.9))

#25in x 10in pdf

#################################################################################
# make plot for paper
# effect of HiT, Fig S15
# F=200, N=5
par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

cont_F200NHiT <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab(expression(paste("Coefficient of variation (",italic(C)[italic(c)],")"))) +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust=0.5), axis.ticks=element_line(),
        plot.margin=margin(t=8, r=-16, b=8, l=0))


cont_F200HiT <- ggplot(powers_table, aes(x=f_inf, y=cv, z=power)) +
  geom_contour_filled(breaks=mybreaks, show.legend=T) +
  scale_x_continuous(expand=c(0,0), limits=c(0.02,0.98), breaks=c(0.02,0.25,0.5,0.75,0.98)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3)) +
  scale_fill_manual(values=mycolors(10), drop=F, na.value="white") +
  theme_minimal() +
  labs(fill="Power") +
  xlab(expression(paste("Expected fraction infected (",italic(E)[italic(c)],")"))) +
  ylab("") +
  #xlim(0,1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), text=element_text(size=20),
        axis.text.x=element_text(angle=90),
        axis.text.y=element_blank(), axis.ticks.length.y=unit(0.71, "cm"), axis.ticks.x=element_line(),
        plot.margin=margin(t=8, r=0, b=8, l=-16))

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
          cont_F200NHiT, cont_F200HiT,
          nrow=2, ncol=2, common.legend=T, legend="right",
          heights=c(0.1,0.9))



# 12 x 7, pdf


##################################################################################



