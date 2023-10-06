## Make plots for the supplementary information
## Goodness of fit method to detect heterogeneity in transmission given contact in the presence of HiS
## Figs S17 and S18


#################################################################
## Example of goodness of fit method to detect heterogeneity in transmission given contact in the presence of HiS in the discrete case
## Fig S17, plotted while running simulation
## Plotting distribution of simulated likelihoods, 95% CIs, and likelihood of real data
## C_d=1.3, E_d=0.25, f_A=0.2, F=200, N=5, m=0.5, phi=2

# potential_Ls=simulated likelihoods, quantLs=95% CI, data_L=real likelihood
plot(density(potential_Ls), yaxs="i", main="", xlab="Log-likelihood", zero.line=F, 
     xlim=c(-320,-200), ylim=c(0,0.055), cex.lab=1.5, cex.axis=1.5)
abline(v=quantLs, lty="dashed", col="gray", ylim=c(0,0.05))
abline(v=data_L, col="blue")
legend("topright", bty="n", legend=c("Simulated","95% CI","Real"),
       col=c("black","gray","blue"),
       lty=c("solid","dashed","solid"), lwd=2, cex=1.2)

# 9.6 x 7.73 pdf


#####################################################################
## Power to detect heterogeneity in transmission given contact in the presence of HiS across varying levels of HiT - goodness of fit method
## discrete and continuous, Fig S18
## C=1.3, E=0.25, fA=0.2, m in [0.1,6], F=200 or 1000, N=5 or 10

# read in data
fileName <- "2gpPowersDetHiT_fA0.2_F1000_N5_N25_sims10000c1.3e0.25m6seed.csv"
fileName <- "ContPowersDetHiT_F1000_N5_N25_sims1000_cv1.3_e0.25m6seed.csv"
binom_dist_table_F1000c <- read.csv(fileName, header=T, row.names=1)

par(mar=c(6.1,6.1,5.1,2.1)) # change margin sizes

### discrete
HiTpow2 <- ggplot() +
  geom_point(data=binom_dist_table_F200, aes(x=m, y=not_binom_dist, col="F=200"), shape=17, cex=3.5) +
  geom_point(data=binom_dist_table_F200N10, aes(x=m, y=not_binom_dist, col="F=200, N=10"), shape=6, cex=3) +
  geom_point(data=binom_dist_table_F1000, aes(x=m, y=not_binom_dist, col="F=1000"), shape=19, cex=3) +
  xlab(expression(paste("Dispersion parameter (",italic(m),")"))) +
  ylab("Power to detect heterogeneity in transmission (%)") +
  theme_minimal() +
  theme(text=element_text(size=20), legend.text=element_text(size=15), legend.position=c(0.87,0.96),
        legend.title.align=0.5, panel.border = element_rect(fill = NA, size = 1), axis.ticks=element_line()) +
  scale_x_continuous(limits=c(0,6), breaks=seq(0,6,1)) +
  scale_y_continuous(limits=c(0,100)) +
  scale_color_manual(name="", limits=c("F=200","F=200, N=10","F=1000"), values = c("gray85","gray50","black")) +
  guides(color=guide_legend(override.aes=list(shape=c(17,6,19)))) +
  geom_text(aes(x=0, y=100, label="a"), size=7, vjust=-0.3, hjust=1.9)

### continuous
HiTpowc <- ggplot() +
  geom_point(data=binom_dist_table_F200c, aes(x=m, y=not_binom_dist, col="F=200"), shape=17, cex=3.5) +
  geom_point(data=binom_dist_table_F200cN10, aes(x=m, y=not_binom_dist, col="F=200, N=10"), shape=6, cex=3) +
  geom_point(data=binom_dist_table_F1000c, aes(x=m, y=not_binom_dist, col="F=1000"), shape=19, cex=3) +
  xlab(expression(paste("Dispersion parameter (",italic(m),")"))) +
  ylab("Power to detect heterogeneity in transmission (%)") +
  theme_minimal() +
  theme(text=element_text(size=20), legend.text=element_text(size=15), legend.position=c(0.87,0.96),
        legend.title.align=0.5, panel.border = element_rect(fill = NA, size = 1), axis.ticks=element_line()) +
  scale_x_continuous(limits=c(0,6), breaks=seq(0,6,1)) +
  scale_y_continuous(limits=c(0,100)) +
  scale_color_manual(name="", limits=c("F=200","F=200, N=10","F=1000"), values = c("gray85","gray50","black")) +
  guides(color=guide_legend(override.aes=list(shape=c(17,6,19)))) +
  geom_text(aes(x=0, y=100, label="b"), size=7, vjust=-0.3, hjust=1.9)

ggarrange(HiTpow2, HiTpowc,
          nrow=1, ncol=2, common.legend=F,
          widths=c(1/2,1/2))

# 17 x 8 pdf
