## Make plots for the paper and supplementary information
## Parameter estimates and SIR dynamics
## Figs 5, 6, 7, 8, S5, S7, S8, S10, S13, S16, S19, S20, S23, S24, S27, S28


#################################################################################
## Parameter estimates plot for the discrete case (p_A vs p_B vs f_A), main text Fig 5
## C_d=1.3, E_d=0.25, f_A=0.2, N=5, p_A = 0.74809046024617, p_B = 0.125477384938458

par(mfrow=c(1,3))

## only plot 1000 points for each param so easier to see
set.seed(5)
sub_results50N5c1.3e0.25fA0.2rd <- results50N5c1.3e0.25fA0.2rd[sample((1:nrow(results50N5c1.3e0.25fA0.2rd)), 1000, replace=F),]
sub_results200N5c1.3e0.25fA0.2rd <- results200N5c1.3e0.25fA0.2rd[sample((1:nrow(results200N5c1.3e0.25fA0.2rd)), 1000, replace=F),]
sub_results1000N5c1.3e0.25fA0.2rd <- results1000N5c1.3e0.25fA0.2rd[sample((1:nrow(results1000N5c1.3e0.25fA0.2rd)), 1000, replace=F),]
sub_results5000N5c1.3e0.25fA0.2rd <- results5000N5c1.3e0.25fA0.2rd[sample((1:nrow(results5000N5c1.3e0.25fA0.2rd)), 1000, replace=F),]

## p_B vs p_A
plot(sub_results50N5c1.3e0.25fA0.2rd$p_b, sub_results50N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="gray85", xlim=c(0,0.35), ylim=c(0.1,1.04), 
     xlab=expression(paste("Probability of infection, less susceptible (",italic(p)[italic(B)],")")), 
     ylab=expression(paste("Probability of infection, more susceptible (",italic(p)[italic(A)],")")), 
     cex.axis=2, cex.lab=2)
points(sub_results200N5c1.3e0.25fA0.2rd$p_b, sub_results200N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="gray70")
points(sub_results1000N5c1.3e0.25fA0.2rd$p_b, sub_results1000N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="gray50")
points(sub_results5000N5c1.3e0.25fA0.2rd$p_b, sub_results5000N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="black")
points(0.125477384938458, 0.74809046024617, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=50","F=200","F=1000","F=5000"),
       col=c("red","gray85","gray70","gray50","black"), pch=19, cex=1.5)
text(0, 1.04, labels="a", cex=2)


## f_A vs p_A
plot(sub_results50N5c1.3e0.25fA0.2rd$f_a, sub_results50N5c1.3e0.25fA0.2rd$p_a,
     pch=19, col="gray85", xlim=c(0,1), ylim=c(0.1,1.04), 
     xlab=expression(paste("Fraction of population more susceptible (",italic(f)[italic(A)],")")), 
     ylab=expression(paste("Probability of infection, more susceptible (",italic(p)[italic(A)],")")), 
     cex.axis=2, cex.lab=2)
points(sub_results200N5c1.3e0.25fA0.2rd$f_a, sub_results200N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="gray70")
points(sub_results1000N5c1.3e0.25fA0.2rd$f_a, sub_results1000N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="gray50")
points(sub_results5000N5c1.3e0.25fA0.2rd$f_a, sub_results5000N5c1.3e0.25fA0.2rd$p_a,
       pch=19, col="black")
points(0.2, 0.74809046024617, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=50","F=200","F=1000","F=5000"),
       col=c("red","gray85","gray70","gray50","black"), pch=19, cex=1.5)
text(0, 1.04, labels="b", cex=2)


## f_A vs p_B
plot(sub_results50N5c1.3e0.25fA0.2rd$f_a, sub_results50N5c1.3e0.25fA0.2rd$p_b,
     pch=19, col="gray85", xlim=c(0,1), ylim=c(0,0.35), 
     xlab=expression(paste("Fraction of population more susceptible (",italic(f)[italic(A)],")")), 
     ylab=expression(paste("Probability of infection, less susceptible (",italic(p)[italic(B)],")")), 
     cex.axis=2, cex.lab=2)
points(sub_results200N5c1.3e0.25fA0.2rd$f_a, sub_results200N5c1.3e0.25fA0.2rd$p_b,
       pch=19, col="gray70")
points(sub_results1000N5c1.3e0.25fA0.2rd$f_a, sub_results1000N5c1.3e0.25fA0.2rd$p_b,
       pch=19, col="gray50")
points(sub_results5000N5c1.3e0.25fA0.2rd$f_a, sub_results5000N5c1.3e0.25fA0.2rd$p_b,
       pch=19, col="black")
points(0.2, 0.125477384938458, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=50","F=200","F=1000","F=5000"),
       col=c("red","gray85","gray70","gray50","black"), pch=19, cex=1.5)
text(0, 0.35, labels="c", cex=2)

par(mfrow=c(1,1))

# 14 x 7 pdf


####################################################################
## Parameter estimates plot for the continuous case (k vs theta), main text Fig 6
## C_c=1.3, E_c=0.25, N=5, k = 0.591715976331361, theta = 0.626097060979711

## only plot 1000 points for each param so easier to see
set.seed(5)
sub_resultsc50N5c1.3e0.25rd <- resultsc50N5c1.3e0.25rd[sample((1:nrow(resultsc50N5c1.3e0.25rd)), 1000, replace=F),]
sub_resultsc200N5c1.3e0.25rd <- resultsc200N5c1.3e0.25rd[sample((1:nrow(resultsc200N5c1.3e0.25rd)), 1000, replace=F),]
sub_resultsc1000N5c1.3e0.25rd <- resultsc1000N5c1.3e0.25rd[sample((1:nrow(resultsc1000N5c1.3e0.25rd)), 1000, replace=F),]
sub_resultsc5000N5c1.3e0.25rd <- resultsc5000N5c1.3e0.25rd[sample((1:nrow(resultsc5000N5c1.3e0.25rd)), 1000, replace=F),]

## k vs theta
plot(sub_resultsc50N5c1.3e0.25rd$k, sub_resultsc50N5c1.3e0.25rd$theta, 
     pch=19, col="gray85", xlim=c(0,2), ylim=c(0,10), xlab=expression(paste("Shape parameter (",italic(k),")")), 
     ylab=expression(paste("Scale parameter (",italic(theta),")")), 
     cex.axis=2, cex.lab=2)
points(sub_resultsc200N5c1.3e0.25rd$k, sub_resultsc200N5c1.3e0.25rd$theta,
       pch=19, col="gray70")
points(sub_resultsc1000N5c1.3e0.25rd$k, sub_resultsc1000N5c1.3e0.25rd$theta,
       pch=19, col="gray50")
points(sub_resultsc5000N5c1.3e0.25rd$k, sub_resultsc5000N5c1.3e0.25rd$theta,
       pch=19, col="black")
points(0.591715976331361, 0.626097060979711, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=50","F=200","F=1000","F=5000"),
       col=c("red","gray85","gray70","gray50","black"), pch=19, cex=1.5)

# 9.5 x 7 pdf


#########################################################
## Parameter estimates plot for the CV of risk of infection and expected fraction of naive individuals infected
## discrete and continuous, main text Fig 7
## C=1.3, E=0.25, f_A=0.2, N=5

par(mar=c(5.1,5.1,4.1,2.1), mfrow=c(1,2)) # change margin sizes

### discrete
# find CV of risk from p_A, p_B, f_A
cv_risk <- function(pA,pB,fA)
{
  rA <- -log(1 - pA)
  rB <- -log(1 - pB)
  
  return(((rA - rB)*sqrt(fA*(1 - fA))) / (rA*fA + rB*(1 - fA)))
}

# find exp frac inf from p_A, p_B, f_A
exp_frac_inf <- function(pA,pB,fA)
{
  pA*fA + pB*(1-fA)
}

## only plot 1000 points for each param so easier to see
set.seed(5)
sub_results50N5c1.3e0.25fA0.2rd <- results50N5c1.3e0.25fA0.2rd[sample((1:nrow(results50N5c1.3e0.25fA0.2rd)), 1000, replace=F),]
sub_results200N5c1.3e0.25fA0.2rd <- results200N5c1.3e0.25fA0.2rd[sample((1:nrow(results200N5c1.3e0.25fA0.2rd)), 1000, replace=F),]
sub_results1000N5c1.3e0.25fA0.2rd <- results1000N5c1.3e0.25fA0.2rd[sample((1:nrow(results1000N5c1.3e0.25fA0.2rd)), 1000, replace=F),]
sub_results5000N5c1.3e0.25fA0.2rd <- results5000N5c1.3e0.25fA0.2rd[sample((1:nrow(results5000N5c1.3e0.25fA0.2rd)), 1000, replace=F),]

## C_d vs E_d
plot(data.frame(efi=exp_frac_inf(sub_results50N5c1.3e0.25fA0.2rd$p_a, sub_results50N5c1.3e0.25fA0.2rd$p_b, sub_results50N5c1.3e0.25fA0.2rd$f_a), 
                cv=cv_risk(sub_results50N5c1.3e0.25fA0.2rd$p_a, sub_results50N5c1.3e0.25fA0.2rd$p_b, sub_results50N5c1.3e0.25fA0.2rd$f_a)), 
     pch=19, col="gray85", xlim=c(0.1,0.5), ylim=c(0,3), xlab=expression(paste("Expected fraction infected (",italic(E)[italic(d)],")")), 
     ylab=expression(paste("Coefficient of variation (",italic(C)[italic(d)],")")), cex.axis=2, cex.lab=2)
points(data.frame(efi=exp_frac_inf(sub_results200N5c1.3e0.25fA0.2rd$p_a, sub_results200N5c1.3e0.25fA0.2rd$p_b, sub_results200N5c1.3e0.25fA0.2rd$f_a), 
                  cv=cv_risk(sub_results200N5c1.3e0.25fA0.2rd$p_a, sub_results200N5c1.3e0.25fA0.2rd$p_b, sub_results200N5c1.3e0.25fA0.2rd$f_a)), 
       pch=19, col="gray70")
points(data.frame(efi=exp_frac_inf(sub_results1000N5c1.3e0.25fA0.2rd$p_a, sub_results1000N5c1.3e0.25fA0.2rd$p_b, sub_results1000N5c1.3e0.25fA0.2rd$f_a), 
                  cv=cv_risk(sub_results1000N5c1.3e0.25fA0.2rd$p_a, sub_results1000N5c1.3e0.25fA0.2rd$p_b, sub_results1000N5c1.3e0.25fA0.2rd$f_a)), 
       pch=19, col="gray50")
points(data.frame(efi=exp_frac_inf(sub_results5000N5c1.3e0.25fA0.2rd$p_a, sub_results5000N5c1.3e0.25fA0.2rd$p_b, sub_results5000N5c1.3e0.25fA0.2rd$f_a), 
                  cv=cv_risk(sub_results5000N5c1.3e0.25fA0.2rd$p_a, sub_results5000N5c1.3e0.25fA0.2rd$p_b, sub_results5000N5c1.3e0.25fA0.2rd$f_a)), 
       pch=19, col="black")
points(0.25,1.3,pch=19,col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=50","F=200","F=1000","F=5000"),
       col=c("red","gray85","gray70","gray50","black"), pch=19, cex=1.3)
text(0.1, 3, labels="a", cex=2)



### continuous
# find CV of risk from k
cv_risk <- function(k)
{
  1/sqrt(k)
}

# find exp frac inf from k and theta
exp_frac_inf <- function(k, theta)
{
  1 - (1 + theta)^(-k)
}

## only plot 1000 points for each param so easier to see
set.seed(5)
sub_resultsc50N5c1.3e0.25rd <- resultsc50N5c1.3e0.25rd[sample((1:nrow(resultsc50N5c1.3e0.25rd)), 1000, replace=F),]
sub_resultsc200N5c1.3e0.25rd <- resultsc200N5c1.3e0.25rd[sample((1:nrow(resultsc200N5c1.3e0.25rd)), 1000, replace=F),]
sub_resultsc1000N5c1.3e0.25rd <- resultsc1000N5c1.3e0.25rd[sample((1:nrow(resultsc1000N5c1.3e0.25rd)), 1000, replace=F),]
sub_resultsc5000N5c1.3e0.25rd <- resultsc5000N5c1.3e0.25rd[sample((1:nrow(resultsc5000N5c1.3e0.25rd)), 1000, replace=F),]


## C_c vs E_c
plot(data.frame(efi=exp_frac_inf(sub_resultsc50N5c1.3e0.25rd$k, sub_resultsc50N5c1.3e0.25rd$theta), 
                cv=cv_risk(sub_resultsc50N5c1.3e0.25rd$k)), 
     pch=19, col="gray85", xlim=c(0.1,0.5), ylim=c(0,3), xlab=expression(paste("Expected fraction infected (",italic(E)[italic(c)],")")), 
     ylab=expression(paste("Coefficient of variation (",italic(C)[italic(c)],")")), cex.axis=2, cex.lab=2)
points(data.frame(efi=exp_frac_inf(sub_resultsc200N5c1.3e0.25rd$k, sub_resultsc200N5c1.3e0.25rd$theta), 
                  cv=cv_risk(sub_resultsc200N5c1.3e0.25rd$k)), 
       pch=19, col="gray70")
points(data.frame(efi=exp_frac_inf(sub_resultsc1000N5c1.3e0.25rd$k, sub_resultsc1000N5c1.3e0.25rd$theta), 
                  cv=cv_risk(sub_resultsc1000N5c1.3e0.25rd$k)), 
       pch=19, col="gray50")
points(data.frame(efi=exp_frac_inf(sub_resultsc5000N5c1.3e0.25rd$k, sub_resultsc5000N5c1.3e0.25rd$theta), 
                  cv=cv_risk(sub_resultsc5000N5c1.3e0.25rd$k)), 
       pch=19, col="black")
points(0.25,1.3,pch=19,col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=50","F=200","F=1000","F=5000"),
       col=c("red","gray85","gray70","gray50","black"), pch=19, cex=1.3)
text(0.1, 3, labels="b", cex=2)

# 10 x 7 pdf

par(mfrow=c(1,1))


##################################################################################
## SIR dynamics, main text, Fig 8
## discrete and continuous
## C=1.3, E=0.25, f_A=0.2, N=5
## make sure true (actual) and no HiS (nohet) dynamics are right for each case

par(mfrow=c(1,2))

### discrete
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)


lines(t, qs_2.5S50N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50N5c1.3e0.25fA0.2rd, rev(qs_97.5S50N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S200N5c1.3e0.25fA0.2rd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5S200N5c1.3e0.25fA0.2rd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S200N5c1.3e0.25fA0.2rd, rev(qs_97.5S200N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="twodash")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="twodash")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2rd, rev(qs_97.5S1000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2rd, rev(qs_97.5S5000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=200","F=1000","F=5000"),
       col=c("maroon","blue","gray85","gray70","gray50","black"),
       lty=c("solid","solid","dotted","dashed","twodash","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="a", cex=1.5)


### continuous
par(mar=c(5.1,5.1,4.1,2.1)) # change margin sizes

t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc50N5c1.3e0.25rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5Sc50N5c1.3e0.25rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5Sc50N5c1.3e0.25rd, rev(qs_97.5Sc50N5c1.3e0.25rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc200N5c1.3e0.25rd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc200N5c1.3e0.25rd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc200N5c1.3e0.25rd, rev(qs_97.5Sc200N5c1.3e0.25rd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc1000N5c1.3e0.25rd, col="gray50", lwd=1.5, lty="twodash")
lines(t, qs_97.5Sc1000N5c1.3e0.25rd, col="gray50", lwd=1.5, lty="twodash")
polygon(c(t, rev(t)), c(qs_2.5Sc1000N5c1.3e0.25rd, rev(qs_97.5Sc1000N5c1.3e0.25rd)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc5000N5c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000N5c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25rd, rev(qs_97.5Sc5000N5c1.3e0.25rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=200","F=1000","F=5000"),
       col=c("maroon","blue","gray85","gray70","gray50","black"),
       lty=c("solid","solid","dotted","dashed","twodash","solid"), lwd=2, cex=1.2)

par(mfrow=c(1,1))


##################################################################################
## SIR dynamics showing effect of changing N
## discrete and continuous, Fig S5
## C=1.3, E=0.25, f_A=0.2, F=200 or 5000, N=5 or 100
## make sure true (actual) and no HiS (nohet) dynamics are right for each case

layout(matrix(c(1,1,1,
                2,3,4,
                5,6,7,
                8,9,10), ncol=3, byrow=T),
       heights=c(0.05,0.05,0.45,0.45), widths=c(0.1,0.45,0.45))

# top labels for F
par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")
text(1.065, 1, labels=expression(paste("Number of focal individuals (",italic(F),")")), cex=2.5)

plot(1, type="n", axes=F, xlab="", ylab="")

plot(1, type="n", axes=F, xlab="", ylab="")
text(1.06, 1, labels="200", cex=2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1.06, 1, labels="5000", cex=2)


# Discrete label
par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")
text(1.3, 1.06, labels="Discrete", srt=90, cex=2)

# change plot margins
par(mar=c(4.0,5.1,1.0,1.0))


### discrete
# F=200
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S200N5c1.3e0.25fA0.2N100subrd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5S200N5c1.3e0.25fA0.2N100subrd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S200N5c1.3e0.25fA0.2N100subrd, rev(qs_97.5S200N5c1.3e0.25fA0.2N100subrd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S200N100c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S200N100c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S200N100c1.3e0.25fA0.2rd, rev(qs_97.5S200N100c1.3e0.25fA0.2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","N=5","N=100"),
       col=c("maroon","blue","gray70","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)


# F=5000
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="", ylab="", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2N100subrd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2N100subrd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2N100subrd, rev(qs_97.5S5000N5c1.3e0.25fA0.2N100subrd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N100c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N100c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N100c1.3e0.25fA0.2rd, rev(qs_97.5S5000N100c1.3e0.25fA0.2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","N=5","N=100"),
       col=c("maroon","blue","gray70","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)


# Continuous label
par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")
text(1.3, 1.06, labels="Continuous", srt=90, cex=2)

### continuous
# F=200
t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc200N5c1.3e0.25N100subrd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc200N5c1.3e0.25N100subrd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc200N5c1.3e0.25N100subrd, rev(qs_97.5Sc200N5c1.3e0.25N100subrd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc200N100c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc200N100c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc200N100c1.3e0.25rd, rev(qs_97.5Sc200N100c1.3e0.25rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","N=5","N=100"),
       col=c("maroon","blue","gray70","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)


# F=5000
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc5000N5c1.3e0.25N100subrd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc5000N5c1.3e0.25N100subrd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25N100subrd, rev(qs_97.5Sc5000N5c1.3e0.25N100subrd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc5000N100c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000N100c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N100c1.3e0.25rd, rev(qs_97.5Sc5000N100c1.3e0.25rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","N=5","N=100"),
       col=c("maroon","blue","gray70","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)

# 10 x 8 pdf

par(mfrow=c(1,1))


######################################################################################
## SIR dynamics showing effect of changing error tolerance for ABC in the discrete case
## Fig S7
## C_d=1.3, E_d=0.25, f_A=0.2, N=5, F=200 or 1000, error tolerance = 10%, 1%, or 0%

par(mfrow=c(1,2))


# F=200
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S200N5c1.3e0.25fA0.2tol0.1rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S200N5c1.3e0.25fA0.2tol0.1rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S200N5c1.3e0.25fA0.2tol0.1rd, rev(qs_97.5S200N5c1.3e0.25fA0.2tol0.1rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S200N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="dashed")
lines(t, qs_97.5S200N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S200N5c1.3e0.25fA0.2rd, rev(qs_97.5S200N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S200N5c1.3e0.25fA0.2tol0.0rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S200N5c1.3e0.25fA0.2tol0.0rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S200N5c1.3e0.25fA0.2tol0.0rd, rev(qs_97.5S200N5c1.3e0.25fA0.2tol0.0rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","10%","1%","0%"),
       col=c("maroon","blue","gray85","gray50","black"),
       lty=c("solid","solid","dotted","dashed","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="a", cex=1.5)


#F=1000
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2tol0.1rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2tol0.1rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2tol0.1rd, rev(qs_97.5S1000N5c1.3e0.25fA0.2tol0.1rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="dashed")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2rd, rev(qs_97.5S1000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2tol0.0rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2tol0.0rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2tol0.0rd, rev(qs_97.5S1000N5c1.3e0.25fA0.2tol0.0rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","10%","1%","0%"),
       col=c("maroon","blue","gray85","gray50","black"), 
       lty=c("solid","solid","dotted","dashed","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="b", cex=1.5)

par(mfrow=c(1,1))

# 10 x 7 pdf


##################################################################################
## SIR dynamics showing effect of assuming wrong underlying model for estimation
## discrete and continuous, Fig S8
## C=1.3, E=0.25, f_A=0.2 or 0.5, F=5000, N=5
## make sure true (actual) and no HiS (nohet) dynamics are right for each case

par(mfrow=c(1,3))

### discrete is correct model and continuous is incorrectly assumed, f_A=0.2
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5Sc5000N5c1.3e0.25wd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000N5c1.3e0.25wd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25wd, rev(qs_97.5Sc5000N5c1.3e0.25wd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=5000"),
       col=c("maroon","blue","black"),
       lty=c("solid","solid","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="a", cex=1.5)


### discrete is correct model and continuous is incorrectly assumed, f_A=0.5
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5Sc5000N5c1.3e0.25fA0.5wd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000N5c1.3e0.25fA0.5wd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25fA0.5wd, rev(qs_97.5Sc5000N5c1.3e0.25fA0.5wd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=5000"),
       col=c("maroon","blue","black"),
       lty=c("solid","solid","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="b", cex=1.5)


### continuous is correct model and discrete is incorrectly assumed
t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5S5000N5c1.3e0.25wd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25wd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25wd, rev(qs_97.5S5000N5c1.3e0.25wd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=5000"),
       col=c("maroon","blue","black"),
       lty=c("solid","solid","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="c", cex=1.5)

par(mfrow=c(1,1))

# 14 x 7 for pdf


#################################################################################################
## SIR dynamics showing effect of having informative priors for p_A, p_B, or f_A in the discrete case, Fig S10
## C_d=1.3, E_d=0.25, f_A=0.2, N=5
## first plot shows dynamics with uninformative priors

par(mar=c(4.0,5.1,1.0,1.0), mfrow=c(2,2))

### uninformative priors
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S50N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50N5c1.3e0.25fA0.2rd, rev(qs_97.5S50N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2rd, rev(qs_97.5S5000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=5000"),
      col=c("maroon","blue","gray85","black"),
      lty=c("solid","solid","dotted","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="a", cex=1.5)


### p_A informative prior
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="", ylab="", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S50N5c1.3e0.25fA0.2pApriorrd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50N5c1.3e0.25fA0.2pApriorrd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50N5c1.3e0.25fA0.2pApriorrd, rev(qs_97.5S50N5c1.3e0.25fA0.2pApriorrd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2pApriorrd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2pApriorrd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2pApriorrd, rev(qs_97.5S5000N5c1.3e0.25fA0.2pApriorrd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=5000"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dotted","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="b", cex=1.5)


### p_B informative prior
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S50N5c1.3e0.25fA0.2pBpriorrd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50N5c1.3e0.25fA0.2pBpriorrd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50N5c1.3e0.25fA0.2pBpriorrd, rev(qs_97.5S50N5c1.3e0.25fA0.2pBpriorrd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2pBpriorrd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2pBpriorrd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2pBpriorrd, rev(qs_97.5S5000N5c1.3e0.25fA0.2pBpriorrd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=5000"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dotted","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="c", cex=1.5)


### f_A informative prior
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S50N5c1.3e0.25fA0.2fApriorrd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50N5c1.3e0.25fA0.2fApriorrd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50N5c1.3e0.25fA0.2fApriorrd, rev(qs_97.5S50N5c1.3e0.25fA0.2fApriorrd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2fApriorrd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2fApriorrd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2fApriorrd, rev(qs_97.5S5000N5c1.3e0.25fA0.2fApriorrd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=5000"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dotted","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="d", cex=1.5)

# 10.5 x 8 pdf


###################################################################################
## SIR dynamics showing effect of false negatives and method adjusted for false negatives
## discrete and continuous, Fig S13
## F=1000, N=5
## C_d=1.3, E_d=0.25, f_A=0.2, beta_d=0.2, C_c=1, E_c=0.75, beta_c=0.1
## make sure true (actual) and no HiS (nohet) dynamics are right for each case

par(mar=c(5.1,5.1,4.1,2.1), mfrow=c(1,2)) # change margin sizes

### discrete
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2fN0.2rd, col="gray85", lwd=1.5, lty="dashed")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2fN0.2rd, col="gray85", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2fN0.2rd, rev(qs_97.5S1000N5c1.3e0.25fA0.2fN0.2rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2fN0.2adjrd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2fN0.2adjrd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2fN0.2adjrd, rev(qs_97.5S1000N5c1.3e0.25fA0.2fN0.2adjrd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","False negs","False negs adj"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="a", cex=1.5)


### continuous
t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc1000N5c1e0.75fN0.1rd, col="gray85", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc1000N5c1e0.75fN0.1rd, col="gray85", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc1000N5c1e0.75fN0.1rd, rev(qs_97.5Sc1000N5c1e0.75fN0.1rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc1000N5c1e0.75fN0.1adjrd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc1000N5c1e0.75fN0.1adjrd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc1000N5c1e0.75fN0.1adjrd, rev(qs_97.5Sc1000N5c1e0.75fN0.1adjrd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","False negs","False negs adj"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="b", cex=1.5)

par(mfrow=c(1,1))

# 10 x 7 pdf


############################################################################################
## SIR dynamics showing effect of heterogeneity in transmission given contact - plotting dynamics without and with HiT
## discrete and continuous, Fig S16
## F=5000, N=5
## C=1.3, E=0.25, fA=0.2, m=0.5, phi=2
## make sure true (actual) and no HiS (nohet) dynamics are right for each case

par(mar=c(5.1,5.1,4.1,2.1)) # change margin sizes
par(mfrow=c(1,2))

### discrete
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dashed")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2rd, rev(qs_97.5S5000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2m0.5phi2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2m0.5phi2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2m0.5phi2rd, rev(qs_97.5S5000N5c1.3e0.25fA0.2m0.5phi2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","No HiT","HiT"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="a", cex=1.5)



### continuous
t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc5000N5c1.3e0.25rd, col="gray85", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc5000N5c1.3e0.25rd, col="gray85", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25rd, rev(qs_97.5Sc5000N5c1.3e0.25rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc5000N5c1.3e0.25m0.5phi2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000N5c1.3e0.25m0.5phi2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25m0.5phi2rd, rev(qs_97.5Sc5000N5c1.3e0.25m0.5phi2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","No HiT","HiT"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)

text(0, 1.04, labels="b", cex=1.5)

par(mfrow=c(1,1))

# 10 x 7 pdf


#################################################################################
## Parameter estimates plot for accounting for HiT (C_c vs E_c vs m), Fig S19
## C_c=0.8, E_c=0.75, N=5, k=1.5625, theta=1.42838976879009, m=1
## note Lamb is included in theta, Lamb=2.26209772417053 so theta=3.231157

## find Lamb
# function to numerically integrate to find Lamb
integ_over_lamb <- function(lamb, m, phi, k, theta, Lamb)
{
  (1 + lamb*Lamb*theta)^(-k) * (1 / (gamma(m) * phi^m)) * lamb^(m-1) * exp(-lamb/phi)
}

# function to find avg force of inf Lamb with root finding
find_Lamb_cont <- function(Lamb, m, phi, k, theta, f_inf)
{
  1 - f_inf - integrate(integ_over_lamb, lower=0, upper=Inf, m=m, phi=phi, k=k, theta=theta, Lamb=Lamb)$value
}

## Lamb
uniroot(find_Lamb_cont, lower=0, upper=300, m=1, phi=1, k=1.5625, theta=1.42838976879009, f_inf=0.75)$root

# find CV of risk from k
cv_risk <- function(k)
{
  1/sqrt(k)
}

# find exp frac inf from k and theta
# function to numerically integrate to get exp_frac_inf w/ HiT
integ_over_lamb <- function(lamb, m, phi, k, theta)
{
  (1 + lamb*theta)^(-k) * (1 / (gamma(m) * phi^m)) * lamb^(m-1) * exp(-lamb/phi)
}

exp_frac_inf <- function(k, theta, m)
{
  1 - integrate(integ_over_lamb, lower=0, upper=Inf, m=m, phi=1/m, k=k, theta=theta)$value
}

set.seed(5)
sub_resultsc1000N5c0.8e0.75m1HiT <- resultsc1000N5c0.8e0.75m1HiT[sample((1:nrow(resultsc1000N5c0.8e0.75m1HiT)), 1000, replace=F),]

### calculate all f_infs with HiT
f_infs_mcmc <- rep(NA, nrow(sub_resultsc1000N5c0.8e0.75m1HiT))
for(row in 1:nrow(sub_resultsc1000N5c0.8e0.75m1HiT))
{
  f_infs_mcmc[row] <- exp_frac_inf(sub_resultsc1000N5c0.8e0.75m1HiT$k[row], sub_resultsc1000N5c0.8e0.75m1HiT$theta[row], sub_resultsc1000N5c0.8e0.75m1HiT$m[row])
}

par(mar=c(5.1,5.1,4.1,2.1), mfrow=c(1,3))

## E_c vs C_c
plot(f_infs_mcmc, 
     cv_risk(sub_resultsc1000N5c0.8e0.75m1HiT$k),
     pch=19, col="gray85", xlim=c(0.65,0.85), ylim=c(0.4,1.2), 
     xlab=expression(paste("Expected fraction infected (",italic(E)[italic(c)],")")), 
     ylab=expression(paste("Coefficient of variation, HiS (",italic(C)[italic(c)],")")), 
     cex.axis=2, cex.lab=2)
points(0.75, 0.8, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.5)
text(0.65, 1.2, labels="a", cex=2)


## m vs C_c
plot(1/sqrt(sub_resultsc1000N5c0.8e0.75m1HiT$m),
     cv_risk(sub_resultsc1000N5c0.8e0.75m1HiT$k),
     pch=19, col="gray85", xlim=c(0.6,1.4), ylim=c(0.4,1.2), 
     xlab=(expression(paste("Coefficient of variation, HiT (",1/sqrt(italic(m)),")"))), 
     ylab=expression(paste("Coefficient of variation, HiS (",italic(C)[italic(c)],")")), 
     cex.axis=2, cex.lab=2)
points(1, 0.8, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.5)
text(0.6, 1.2, labels="b", cex=2)


## E_c vs m
plot(f_infs_mcmc, 
     1/sqrt(sub_resultsc1000N5c0.8e0.75m1HiT$m),
     pch=19, col="gray85", xlim=c(0.65,0.85), ylim=c(0.6,1.4), 
     xlab=expression(paste("Expected fraction infected (",italic(E)[italic(c)],")")), 
     ylab=(expression(paste("Coefficient of variation, HiT (",1/sqrt(italic(m)),")"))), 
     cex.axis=2, cex.lab=2)
points(0.75, 1, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.5)
text(0.65, 1.4, labels="c", cex=2)

par(mfrow=c(1,1))

# 12 x 5.5 pdf


#################################################################################
## Parameter estimates plot for accounting for HiT (C_c vs E_c vs m), Fig S20
## C_c=1.3, E_c=0.25, N=5, k=0.591715976331361, theta=0.626097060979711, m=0.5
## note Lamb is included in theta, Lamb=1.739117 so theta=1.088856

## find Lamb
# function to numerically integrate to find Lamb
integ_over_lamb <- function(lamb, m, phi, k, theta, Lamb)
{
  (1 + lamb*Lamb*theta)^(-k) * (1 / (gamma(m) * phi^m)) * lamb^(m-1) * exp(-lamb/phi)
}

# function to find avg force of inf Lamb with root finding
find_Lamb_cont <- function(Lamb, m, phi, k, theta, f_inf)
{
  1 - f_inf - integrate(integ_over_lamb, lower=0, upper=Inf, m=m, phi=phi, k=k, theta=theta, Lamb=Lamb)$value
}

## Lamb
uniroot(find_Lamb_cont, lower=0, upper=300, m=0.5, phi=2, k=0.591715976331361, theta=0.626097060979711, f_inf=0.25)$root

# find CV of risk from k
cv_risk <- function(k)
{
  1/sqrt(k)
}

# find exp frac inf from k and theta
# function to numerically integrate to get exp_frac_inf w/ HiT
integ_over_lamb <- function(lamb, m, phi, k, theta)
{
  (1 + lamb*theta)^(-k) * (1 / (gamma(m) * phi^m)) * lamb^(m-1) * exp(-lamb/phi)
}

exp_frac_inf <- function(k, theta, m)
{
  1 - integrate(integ_over_lamb, lower=0, upper=Inf, m=m, phi=1/m, k=k, theta=theta)$value
}

set.seed(5)
sub_resultsc1000N5c1.3e0.25m0.5HiT <- resultsc1000N5c1.3e0.25m0.5HiT[sample((1:nrow(resultsc1000N5c1.3e0.25m0.5HiT)), 1000, replace=F),]

### calculate all f_infs with HiT
f_infs_mcmc <- rep(NA, nrow(sub_resultsc1000N5c1.3e0.25m0.5HiT))
for(row in 1:nrow(sub_resultsc1000N5c1.3e0.25m0.5HiT))
{
  f_infs_mcmc[row] <- exp_frac_inf(sub_resultsc1000N5c1.3e0.25m0.5HiT$k[row], sub_resultsc1000N5c1.3e0.25m0.5HiT$theta[row], sub_resultsc1000N5c1.3e0.25m0.5HiT$m[row])
}

par(mar=c(5.1,5.1,4.1,2.1), mfrow=c(1,3))

## E_c vs C_c
plot(f_infs_mcmc, 
     cv_risk(sub_resultsc1000N5c1.3e0.25m0.5HiT$k),
     pch=19, col="gray85", xlim=c(0.15,0.35), ylim=c(0.8,2.4), yaxt="n", 
     xlab=expression(paste("Expected fraction infected (",italic(E)[italic(c)],")")), 
     ylab=expression(paste("Coefficient of variation, HiS (",italic(C)[italic(c)],")")), 
     cex.axis=2, cex.lab=2)
axis(side=2, at=seq(0.8,2.4,0.2), cex.axis=2)
points(0.25, 1.3, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.5)
text(0.15, 2.4, labels="a", cex=2)


## m vs C_c
plot(1/sqrt(sub_resultsc1000N5c1.3e0.25m0.5HiT$m),
     cv_risk(sub_resultsc1000N5c1.3e0.25m0.5HiT$k),
     pch=19, col="gray85", xlim=c(1,2), ylim=c(0.8,2.4), yaxt="n", 
     xlab=(expression(paste("Coefficient of variation, HiT (",1/sqrt(italic(m)),")"))), 
     ylab=expression(paste("Coefficient of variation, HiS (",italic(C)[italic(c)],")")), 
     cex.axis=2, cex.lab=2)
axis(side=2, at=seq(0.8,2.4,0.2), cex.axis=2)
points(1/sqrt(0.5), 1.3, pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.5)
text(1, 2.4, labels="b", cex=2)


## E_c vs m
plot(f_infs_mcmc, 
     1/sqrt(sub_resultsc1000N5c1.3e0.25m0.5HiT$m),
     pch=19, col="gray85", xlim=c(0.15,0.35), ylim=c(1,2), 
     xlab=expression(paste("Expected fraction infected (",italic(E)[italic(c)],")")), 
     ylab=(expression(paste("Coefficient of variation, HiT (",1/sqrt(italic(m)),")"))), 
     cex.axis=2, cex.lab=2)
points(0.25, 1/sqrt(0.5), pch=19, col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.5)
text(0.15, 2, labels="c", cex=2)

par(mfrow=c(1,1))

# 12 x 5.5 pdf


##################################################################################
### make plot for paper/proposal of SIR dynamics
### 2 type static vs dynamic contact networks, Fig S23
### N=5 for static
### make sure actual/nohet results are right for each type of SIR model

par(mfrow=c(1,2))

### static
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S50N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50N5c1.3e0.25fA0.2rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50N5c1.3e0.25fA0.2rd, rev(qs_97.5S50N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S200N5c1.3e0.25fA0.2rd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5S200N5c1.3e0.25fA0.2rd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S200N5c1.3e0.25fA0.2rd, rev(qs_97.5S200N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S1000N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="twodash")
lines(t, qs_97.5S1000N5c1.3e0.25fA0.2rd, col="gray50", lwd=1.5, lty="twodash")
polygon(c(t, rev(t)), c(qs_2.5S1000N5c1.3e0.25fA0.2rd, rev(qs_97.5S1000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000N5c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000N5c1.3e0.25fA0.2rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000N5c1.3e0.25fA0.2rd, rev(qs_97.5S5000N5c1.3e0.25fA0.2rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=200","F=1000","F=5000"),
       col=c("maroon","blue","gray85","gray70","gray50","black"),
       lty=c("solid","solid","dotted","dashed","twodash","solid"), lwd=2, cex=1.2)
text(0, 1.04, labels="a", cex=1.5)


### dynamic
t <- seq(0, EndTime, deltat)
plot(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.axis=1.5, cex.lab=1.5)
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

lines(t, qs_2.5S50c1.3e0.25fA0.2depi, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5S50c1.3e0.25fA0.2depi, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5S50c1.3e0.25fA0.2depi, rev(qs_97.5S50c1.3e0.25fA0.2depi)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S201c1.3e0.25fA0.2depi, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5S201c1.3e0.25fA0.2depi, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5S201c1.3e0.25fA0.2depi, rev(qs_97.5S201c1.3e0.25fA0.2depi)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S1002c1.3e0.25fA0.2depi, col="gray50", lwd=1.5, lty="twodash")
lines(t, qs_97.5S1002c1.3e0.25fA0.2depi, col="gray50", lwd=1.5, lty="twodash")
polygon(c(t, rev(t)), c(qs_2.5S1002c1.3e0.25fA0.2depi, rev(qs_97.5S1002c1.3e0.25fA0.2depi)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5S5000c1.3e0.25fA0.2depi, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5S5000c1.3e0.25fA0.2depi, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5S5000c1.3e0.25fA0.2depi, rev(qs_97.5S5000c1.3e0.25fA0.2depi)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, (actual_results$Sa + actual_results$Sb)/(actual_results$Sa[1] + actual_results$Sb[1]), col="blue", lwd=1, lty="solid")
lines(t, (nohet_results$Sa + nohet_results$Sb) / (nohet_results$Sa[1] + nohet_results$Sb[1]), col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=200","F=1000","F=5000"),
       col=c("maroon","blue","gray85","gray70","gray50","black"),
       lty=c("solid","solid","dotted","dashed","twodash","solid"), lwd=2, cex=1.2)
text(0, 1.04, labels="b", cex=1.5)

# 10 x 7 pdf


##################################################################################
### make plot for paper/proposal of SIR dynamics
### continuous static vs dynamic contact networks, Fig S24
### N=5 for static
### make sure actual/nohet results are right for each type of SIR model

par(mar=c(5.1,5.1,4.1,2.1)) # change margin sizes
par(mfrow=c(1,2))

### static
t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc50N5c1.3e0.25rd, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5Sc50N5c1.3e0.25rd, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5Sc50N5c1.3e0.25rd, rev(qs_97.5Sc50N5c1.3e0.25rd)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc200N5c1.3e0.25rd, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc200N5c1.3e0.25rd, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc200N5c1.3e0.25rd, rev(qs_97.5Sc200N5c1.3e0.25rd)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc1000N5c1.3e0.25rd, col="gray50", lwd=1.5, lty="twodash")
lines(t, qs_97.5Sc1000N5c1.3e0.25rd, col="gray50", lwd=1.5, lty="twodash")
polygon(c(t, rev(t)), c(qs_2.5Sc1000N5c1.3e0.25rd, rev(qs_97.5Sc1000N5c1.3e0.25rd)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc5000N5c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000N5c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000N5c1.3e0.25rd, rev(qs_97.5Sc5000N5c1.3e0.25rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=200","F=1000","F=5000"),
       col=c("maroon","blue","gray85","gray70","gray50","black"),
       lty=c("solid","solid","dotted","dashed","twodash","solid"), lwd=2, cex=1.2)
text(0, 1.04, labels="a", cex=1.5)


### dynamic
t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1.04),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

lines(t, qs_2.5Sc50c1.3e0.25dynamepi_tprior, col="gray85", lwd=1.5, lty="dotted")
lines(t, qs_97.5Sc50c1.3e0.25dynamepi_tprior, col="gray85", lwd=1.5, lty="dotted")
polygon(c(t, rev(t)), c(qs_2.5Sc50c1.3e0.25dynamepi_tprior, rev(qs_97.5Sc50c1.3e0.25dynamepi_tprior)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc200c1.3e0.25dynamepi, col="gray70", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc200c1.3e0.25dynamepi, col="gray70", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc200c1.3e0.25dynamepi, rev(qs_97.5Sc200c1.3e0.25dynamepi)),
        col=adjustcolor("gray70", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc1000c1.3e0.25dynamepi, col="gray50", lwd=1.5, lty="twodash")
lines(t, qs_97.5Sc1000c1.3e0.25dynamepi, col="gray50", lwd=1.5, lty="twodash")
polygon(c(t, rev(t)), c(qs_2.5Sc1000c1.3e0.25dynamepi, rev(qs_97.5Sc1000c1.3e0.25dynamepi)),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc5000c1.3e0.25dynamepi, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc5000c1.3e0.25dynamepi, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc5000c1.3e0.25dynamepi, rev(qs_97.5Sc5000c1.3e0.25dynamepi)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

lines(t, actual_results$
S/actual_results$S[1], col="blue", lwd=1, lty="solid")
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)

legend("topright", bty="n", legend=c("No het","True","F=50","F=200","F=1000","F=5000"),
       col=c("maroon","blue","gray85","gray70","gray50","black"),
       lty=c("solid","solid","dotted","dashed","twodash","solid"), lwd=2, cex=1.2)
text(0, 1.04, labels="b", cex=1.5)
# 10 x 7 pdf


#########################################################
## Parameter estimates plot for the CV of risk of infection and expected fraction of naive individuals infected
## Effect of missing contacts
## Continuous, Fig S27
## C=1.3, E=0.25, F=1000, N=5

# find CV of risk from k
cv_risk <- function(k)
{
  1/sqrt(k)
}

# find exp frac inf from k and theta
exp_frac_inf <- function(k, theta)
{
  1 - (1 + theta)^(-k)
}

## only plot 1000 points for each param so easier to see
set.seed(5)
sub_resultsc1000N5c1.3e0.25missexps <- resultsc1000N5c1.3e0.25missexps[sample((1:nrow(resultsc1000N5c1.3e0.25missexps)), 1000, replace=F),]


## C_c vs E_c
plot(data.frame(efi=exp_frac_inf(sub_resultsc1000N5c1.3e0.25missexps$k, sub_resultsc1000N5c1.3e0.25missexps$theta), 
                cv=cv_risk(sub_resultsc1000N5c1.3e0.25missexps$k)), 
     pch=19, col="gray85", xlim=c(0.15,0.3), ylim=c(0.5,2.5), xlab=expression(paste("Expected fraction infected (",italic(E)[italic(c)],")")), 
     ylab=expression(paste("Coefficient of variation (",italic(C)[italic(c)],")")), cex.axis=2, cex.lab=2)

points(0.25,1.3,pch=19,col="red", cex=1.2)
legend("topright", bty="n", legend=c("True", "F=1000"),
       col=c("red","gray85"), pch=19, cex=1.3)

# 9 x 7 pdf


############################################################################################
## SIR dynamics showing effect of missing contacts
## continuous, Fig S28
## F=1000, N=5
## C=1.3, E=0.25

par(mar=c(5.1,5.1,4.1,2.1)) # change margin sizes

t <- seq(0, EndTime, deltat)
plot(t, actual_results$S/actual_results$S[1], type="l", col="blue", ylim=c(0,1),
     xlab="Time", ylab="Fraction of susceptible individuals", lwd=1, cex.lab=1.5, cex.axis=1.5)
lines(t, nohet_results$S/nohet_results$S[1], col="maroon", lwd=1)
#90,70,50,20

lines(t, qs_2.5Sc1000c1.3e0.25multiexp, col="gray85", lwd=1.5, lty="dashed")
lines(t, qs_97.5Sc1000c1.3e0.25multiexp, col="gray85", lwd=1.5, lty="dashed")
polygon(c(t, rev(t)), c(qs_2.5Sc1000c1.3e0.25multiexp, rev(qs_97.5Sc1000c1.3e0.25multiexp)),
        col=adjustcolor("gray85", alpha.f=0.2) , lty = 0)

lines(t, qs_2.5Sc1000N5c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
lines(t, qs_97.5Sc1000N5c1.3e0.25rd, col="black", lwd=1.5, lty="solid")
polygon(c(t, rev(t)), c(qs_2.5Sc1000N5c1.3e0.25rd, rev(qs_97.5Sc1000N5c1.3e0.25rd)),
        col=adjustcolor("black", alpha.f=0.2) , lty = 0)

legend("topright", bty="n", legend=c("No het","True","Missing Contacts","No Missing Contacts"),
       col=c("maroon","blue","gray85","black"),
       lty=c("solid","solid","dashed","solid"), lwd=2, cex=1.2)

# 9.7 x 7.6 pdf

