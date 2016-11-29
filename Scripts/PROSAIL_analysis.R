
################################################################################
## R-Script - PROSAIL_analysis.R                                              ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Hyperspectral classification of grassland species: towards an  ##
## application for semi-automatic field surveys                               ##
##                                                                            ##
## description: PROSAIL sensitivity analysis to help in the variable          ## 
##              importance discussion. Results presented in Supplementary data## 
##                                                                            ##
################################################################################

library(hsdar)

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
setwd(home)

# lai
parameter <- data.frame(N = rep(1.7, 100),
                        LAI = seq(from = 0.1, to = 8, length.out = 100),
                        Cab = rep(35, 100),
                        Car = rep(12, 100),
                        Cm = rep(0.007, 100),
                        Cw = rep(0.025, 100))

spec_lai = spectra(PROSAIL(parameterList = parameter))      
sens_lai = apply(as.matrix(spec_lai), 2, FUN = mad)

plot(400:2500, sens_lai, ylim=c(0,1), type="l")                               
for(i in 1:100){
  lines(400:2500,spec_lai[i,])
}
lines(400:2500, sens_lai, col="red")

# Cab
parameter <- data.frame(N = rep(1.7, 100),
                        LAI = rep(4, 100),
                        Cab = seq(from = 15, to = 45, length.out = 100),
                        Car = rep(12, 100),
                        Cm = rep(0.007, 100),
                        Cw = rep(0.025, 100))

spec_cab = spectra(PROSAIL(parameterList = parameter))      
sens_cab = apply(as.matrix(spec_cab), 2, FUN = mad)

plot(400:2500, sens_cab, ylim=c(0,1), type="l")                               
for(i in 1:100){
  lines(400:2500,spec_cab[i,])
}
lines(400:2500, sens_cab, col="red") 

# Car
parameter <- data.frame(N = rep(1.7, 100),
                        LAI = rep(4, 100),
                        Cab = rep(35, 100),
                        Car = seq(from = 5, to = 20, length.out = 100),
                        Cm = rep(0.007, 100),
                        Cw = rep(0.025, 100))

spec_car = spectra(PROSAIL(parameterList = parameter))      
sens_car = apply(as.matrix(spec_car), 2, FUN = mad)

plot(400:2500, sens_car, ylim=c(0,1), type="l")                               
for(i in 1:100){
  lines(400:2500,spec_car[i,])
}
lines(400:2500, sens_car, col="red")

# Cm
parameter <- data.frame(N = rep(1.7, 100),
                        LAI = rep(4, 100),
                        Cab = rep(35, 100),
                        Car = rep(12, 100),
                        Cm = seq(from = 0.0025, to = 0.01, length.out = 100),
                        Cw = rep(0.025, 100))

spec_cm = spectra(PROSAIL(parameterList = parameter))      
sens_cm = apply(as.matrix(spec_cm), 2, FUN = mad)

plot(400:2500, sens_cm, ylim=c(0,1), type="l")                               
for(i in 1:no_sim){
  lines(400:2500,spec_cm[i,])
}
lines(400:2500, sens_cm, col="red") 

# Cw
parameter <- data.frame(N = rep(1.7, 100),
                        LAI = rep(4, 100),
                        Cab = rep(35, 100),
                        Car = rep(12, 100),
                        Cm = rep(0.007, 100),
                        Cw = seq(from = 0.01, to = 0.04, length.out = 100))

spec_cw = spectra(PROSAIL(parameterList = parameter))      
sens_cw = apply(as.matrix(spec_cw), 2, FUN = mad)

plot(400:2500, sens_cw, ylim=c(0,1), type="l")                               
for(i in 1:100){
  lines(400:2500,spec_cw[i,])
}
lines(400:2500, sens_cw, col="red")


# Leaf angles
parameter <- data.frame(N = rep(1.7, 100),
                        LAI = rep(4, 100),
                        Cab = rep(35, 100),
                        Car = rep(12, 100),
                        Cm = rep(0.007, 100),
                        Cw = rep(0.025, 100),
                        lidfa = rnorm(100, mean = 0.5, 0.7),
                        lidfb = rnorm(100, mean = -0.048, 0.07))

spec_ang = spectra(PROSAIL(parameterList = parameter))      
sens_ang = apply(as.matrix(spec_ang), 2, FUN = mad)

plot(400:2500, sens_ang, ylim=c(0,1), type="l")                               
for(i in 1:100){
  lines(400:2500,spec_ang[i,])
}
lines(400:2500, sens_ang, col="red")s

### Plot median absolute deviation results
pdf(file = "Figures/PROSAIL.pdf", width=7, height=4)
par(mai=c(1, 1, 0.3, 0.3))
plot(400:2500, ylim=c(0,0.1), axes = FALSE, ylab="Median absolute deviation",
     xlab=expression(lambda(mu*m)))
axis(1, pos=0, at = seq(0, length(sens_lai), 200), labels = seq(0.4, 2.5, 0.2))
axis(2, pos=0, las=1)
#lines(sens_N, col = "yellow")
lines(sens_lai, col="black", lwd=1.5)
lines(sens_cab, col="orange", lwd=1.5)
lines(sens_car, col="red", lwd=1.5)
lines(sens_cm, col="green", lwd=1.5)
lines(sens_cw, col="blue", lwd=1.5)
lines(sens_ang, col="purple",lwd=1.5)
legend("topright", c("LAI", expression(C[ab]), expression(C[ar]), expression(C[m]), expression(C[w]), "ALA"), 
       col=c("black", "orange", "red", "green", "blue", "purple"), lwd=1,
       bty = "n", text.font = 3)
dev.off()
