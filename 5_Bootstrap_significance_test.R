
################################################################################
## R-Script - 5_Bootstrap_significance_test.R                                 ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Mapping plant species in mixed grassland communities using     ##
##             close range imaging spectroscopy                               ##
##                                                                            ##
## description: One-sided bootstrap significance test                         ## 
##                                                                            ##
################################################################################


dat_potVal <- read.table("D:/Sp_Images/BootsClass_out/Fits_all_pot.txt", header = T)
dat_rfVal  <- read.table("D:/Sp_Images/BootsClass_out/Fits_all_rf.txt", header = T)

sig_canopy <- significanceTest_CanopyLevel(dat_rfVal, dat_potVal)

# Hist do not have to overlap with 0, otherwise is not significant
# The black "zero-line" needs to be left of the blue "alpha-line". The green line is just the upper quantile.
par(mfrow=c(3,2), mar=c(2,3,3,1))
main <- names(sig_canopy)
for(i in 1:length(sig_canopy)){
  hist(unlist(sig_canopy[i]), main=main[i], col="grey", border="white", xlab="", ylab="")
  abline(v=quantile(unlist(sig_canopy[i]), probs=c(0.05, 0.95)), col=c("blue", "green"))
  abline(v=0, col=c("black"))
  box()
}
