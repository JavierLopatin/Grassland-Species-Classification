## R-Script - Analysis
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
DDir = "D:/Sp_Images"

setwd(home)


## load the data
fit_potVal <-  read.table("Data/Fits_potVal.csv", sep = ",", header = T)
fit_rf     <-  read.table("Data/Fits_rf.csv", sep = ",", header = T)

## Obtain cover fits
setwd(DDir)

potVal_cover <- coverSummary("potVal")
rf_cover     <- coverSummary("rf")

potVal_cover$Species <- factor(potVal_cover$Species)
setwd(home)

# obtain R2, RMSE and bias per specie 
spList = factor( unique(potVal_cover$Species) )

outputGOF <- data.frame(Species=character(), r2=double(), RMSE=double(),
                        bias=double(), Models=character(), Normalization=character(),
                        Validation=character())

for (i in 1:length(spList)){
  x = grep( spList[i], potVal_cover$Species)
  pv = potVal_cover[x, ]
  
  x = grep(spList[i], rf_cover$Species)
  rf = rf_cover[x, ]
  
  #### Pot validation method ####
  ## Spectral values
  x = grep("spect", pv$Normalization)
  pv_spect_all <- pv[x, ]
  # take care of BN datas
  x = grep("BN", pv_spect_all$Normalization)
  pv_spect_BN <- pv_spect_all[x, ]
  pv_spect <- pv_spect_all[-x, ]
  
  ## MNF values
  x = grep("MNF", pv$Normalization)
  pv_MNF_all <- pv[x, ]
  # take care of BN datas
  x = grep("BN", pv_MNF_all$Normalization)
  pv_MNF_BN <- pv_MNF_all[x, ]
  pv_MNF <- pv_MNF_all[-x, ]
  
  # PLS
  x = grep("PLS", pv_spect$Model)
  pv_spect_PLS <- pv_spect[x, ]
  
  x = grep("PLS", pv_spect_BN$Model)
  pv_spect_BN_PLS <- pv_spect_BN[x, ]
  
  x = grep("PLS", pv_MNF$Model)
  pv_MNF_PLS <- pv_MNF[x, ]
  
  x = grep("PLS", pv_MNF_BN$Model)
  pv_MNF_BN_PLS <- pv_MNF_BN[x, ]
  
  # RF
  x = grep("RF", pv_spect$Model)
  pv_spect_RF <- pv_spect[x, ]
  
  x = grep("RF", pv_spect_BN$Model)
  pv_spect_BN_RF <- pv_spect_BN[x, ]
  
  x = grep("RF", pv_MNF$Model)
  pv_MNF_RF <- pv_MNF[x, ]
  
  x = grep("RF", pv_MNF_BN$Model)
  pv_MNF_BN_RF <- pv_MNF_BN[x, ]
  
  # SVM
  x = grep("SVM", pv_spect$Model)
  pv_spect_SVM <- pv_spect[x, ]
  
  x = grep("SVM", pv_spect_BN$Model)
  pv_spect_BN_SVM <- pv_spect_BN[x, ]
  
  x = grep("SVM", pv_MNF$Model)
  pv_MNF_SVM <- pv_MNF[x, ]
  
  x = grep("SVM", pv_MNF_BN$Model)
  pv_MNF_BN_SVM <- pv_MNF_BN[x, ]
  
  #### rip-it-off validation method ####
  ## Spectral values
  x = grep("spect", rf$Normalization)
  rf_spect_all <- rf[x, ]
  # take care of BN datas
  x = grep("BN", rf_spect_all$Normalization)
  rf_spect_BN <- rf_spect_all[x, ]
  rf_spect <- rf_spect_all[-x, ]
  
  ## MNF values
  x = grep("MNF", rf$Normalization)
  rf_MNF_all <- rf[x, ]
  # take care of BN datas
  x = grep("BN", rf_MNF_all$Normalization)
  rf_MNF_BN <- rf_MNF_all[x, ]
  rf_MNF <- rf_MNF_all[-x, ]
  
  # PLS
  x = grep("PLS", rf_spect$Model)
  rf_spect_PLS <- rf_spect[x, ]
  
  x = grep("PLS", rf_spect_BN$Model)
  rf_spect_BN_PLS <- rf_spect_BN[x, ]
  
  x = grep("PLS", rf_MNF$Model)
  rf_MNF_PLS <- rf_MNF[x, ]
  
  x = grep("PLS", rf_MNF_BN$Model)
  rf_MNF_BN_PLS <- rf_MNF_BN[x, ]
  
  # RF
  x = grep("RF", rf_spect$Model)
  rf_spect_RF <- rf_spect[x, ]
  
  x = grep("RF", rf_spect_BN$Model)
  rf_spect_BN_RF <- rf_spect_BN[x, ]
  
  x = grep("RF", rf_MNF$Model)
  rf_MNF_RF <- rf_MNF[x, ]
  
  x = grep("RF", rf_MNF_BN$Model)
  rf_MNF_BN_RF <- rf_MNF_BN[x, ]
  
  # SVM
  x = grep("SVM", rf_spect$Model)
  rf_spect_SVM <- rf_spect[x, ]
  
  x = grep("SVM", rf_spect_BN$Model)
  rf_spect_BN_SVM <- rf_spect_BN[x, ]
  
  x = grep("SVM", rf_MNF$Model)
  rf_MNF_SVM <- rf_MNF[x, ]
  
  x = grep("SVM", rf_MNF_BN$Model)
  rf_MNF_BN_SVM <- rf_MNF_BN[x, ]
  
  
  ## Goodness-of-fits
  
  out <- matrix(nrow = 24, ncol = 7)
  colnames(out) <- c("Species", "r2", "RMSE", "bias", "Models", "Normalization", "Validation")
  
  out[,1] <- as.character(spList[i])
  out[,5] <- rep(c("PLS", "RF", "SVM"), 8)
  out[,6] <- rep( c(rep("spectra", 3), rep("spectra_BN",3), rep("MNF", 3), rep("MNF_BN", 3 ) ), 2 )
  out[,7] <- c( rep("potVal", 12), rep("rf", 12) )
  
  # r²
  # pot Val
  out[1,2] <- (cor(pv_spect_PLS$Predicted, pv_spect_PLS$Observed, method="pearson"))^2
  out[2,2]  <- (cor(pv_spect_RF$Predicted,  pv_spect_RF$Observed, method="pearson"))^2
  out[3,2] <- (cor(pv_spect_SVM$Predicted, pv_spect_SVM$Observed, method="pearson"))^2
  
  out[4,2] <- (cor(pv_spect_BN_PLS$Predicted, pv_spect_BN_PLS$Observed, method="pearson"))^2
  out[5,2] <- (cor(pv_spect_BN_RF$Predicted, pv_spect_BN_RF$Observed, method="pearson"))^2
  out[6,2] <- (cor(pv_spect_BN_SVM$Predicted, pv_spect_BN_SVM$Observed, method="pearson"))^2
 
  out[7,2] <- (cor(pv_MNF_PLS$Predicted, pv_MNF_PLS$Observed, method="pearson"))^2
  out[8,2] <- (cor(pv_MNF_RF$Predicted,  pv_MNF_RF$Observed, method="pearson"))^2
  out[9,2] <- (cor(pv_MNF_SVM$Predicted, pv_MNF_SVM$Observed, method="pearson"))^2
  
  out[10,2] <- (cor(pv_MNF_BN_PLS$Predicted, pv_MNF_BN_PLS$Observed, method="pearson"))^2
  out[11,2] <- (cor(pv_MNF_BN_RF$Predicted,  pv_MNF_BN_RF$Observed, method="pearson"))^2
  out[12,2] <- (cor(pv_MNF_BN_SVM$Predicted, pv_MNF_BN_SVM$Observed, method="pearson"))^2
  
  # rf
  out[13,2] <- (cor(rf_spect_PLS$Predicted, rf_spect_PLS$Observed, method="pearson"))^2
  out[14,2]  <- (cor(rf_spect_RF$Predicted, rf_spect_RF$Observed, method="pearson"))^2
  out[15,2] <- (cor(rf_spect_SVM$Predicted, rf_spect_SVM$Observed, method="pearson"))^2
  
  out[16,2] <- (cor(rf_spect_BN_PLS$Predicted, rf_spect_BN_PLS$Observed, method="pearson"))^2
  out[17,2]  <- (cor(rf_spect_BN_RF$Predicted, rf_spect_BN_RF$Observed, method="pearson"))^2
  out[18,2] <- (cor(rf_spect_BN_SVM$Predicted, rf_spect_BN_SVM$Observed, method="pearson"))^2
  
  out[19,2] <- (cor(rf_MNF_PLS$Predicted, rf_MNF_PLS$Observed, method="pearson"))^2
  out[20,2]  <- (cor(rf_MNF_RF$Predicted, rf_MNF_RF$Observed, method="pearson"))^2
  out[21,2] <- (cor(rf_MNF_SVM$Predicted, rf_MNF_SVM$Observed, method="pearson"))^2
  
  out[22,2] <- (cor(rf_spect_PLS$Predicted, rf_spect_PLS$Observed, method="pearson"))^2
  out[23,2]  <- (cor(rf_spect_RF$Predicted, rf_spect_RF$Observed, method="pearson"))^2
  out[24,2] <- (cor(rf_spect_SVM$Predicted, rf_spect_SVM$Observed, method="pearson"))^2
  
  # RMSE
  # pot Val
  out[1,3] <- sqrt(mean((pv_spect_PLS$Observed - pv_spect_PLS$Predicted)^2))
  out[2,3]  <- sqrt(mean((pv_spect_RF$Observed- pv_spect_RF$Predicted)^2))
  out[3,3] <- sqrt(mean((pv_spect_SVM$Observed- pv_spect_SVM$Predicted)^2))
  
  out[4,3] <- sqrt(mean((pv_spect_BN_PLS$Observed- pv_spect_BN_PLS$Predicted)^2))
  out[5,3]  <- sqrt(mean((pv_spect_BN_RF$Observed- pv_spect_BN_RF$Predicted)^2))
  out[6,3] <- sqrt(mean((pv_spect_BN_SVM$Observed- pv_spect_BN_SVM$Predicted)^2))
  
  out[7,3] <- sqrt(mean((pv_MNF_PLS$Observed- pv_MNF_PLS$Predicted)^2))
  out[8,3]  <- sqrt(mean((pv_MNF_RF$Observed- pv_MNF_RF$Predicted)^2))
  out[9,3] <- sqrt(mean((pv_MNF_SVM$Observed- pv_MNF_SVM$Predicted)^2))
  
  out[10,3] <- sqrt(mean((pv_MNF_BN_PLS$Observed- pv_MNF_BN_PLS$Predicted)^2))
  out[11,3]  <- sqrt(mean((pv_MNF_BN_RF$Observed- pv_MNF_BN_RF$Predicted)^2))
  out[12,3] <- sqrt(mean((pv_MNF_BN_SVM$Observed- pv_MNF_BN_SVM$Predicted)^2))
  
  # rf
  out[13,3] <- sqrt(mean((rf_spect_PLS$Observed- rf_spect_PLS$Predicted)^2))
  out[14,3]  <- sqrt(mean((rf_spect_RF$Observed- rf_spect_RF$Predicted)^2))
  out[15,3] <- sqrt(mean((rf_spect_SVM$Observed- rf_spect_SVM$Predicted)^2))

  out[16,3] <- sqrt(mean((rf_spect_BN_PLS$Observed- rf_spect_BN_PLS$Predicted)^2))
  out[17,3]  <- sqrt(mean((rf_spect_BN_RF$Observed- rf_spect_BN_RF$Predicted)^2))
  out[18,3] <- sqrt(mean((rf_spect_BN_SVM$Observed- rf_spect_BN_SVM$Predicted)^2))
  
  out[19,3] <- sqrt(mean((rf_MNF_PLS$Observed- rf_MNF_PLS$Predicted)^2))
  out[20,3]  <- sqrt(mean((rf_MNF_RF$Observed- rf_MNF_RF$Predicted)^2))
  out[21,3] <- sqrt(mean((rf_MNF_SVM$Observed- rf_MNF_SVM$Predicted)^2))
  
  out[22,3] <- sqrt(mean((rf_MNF_BN_PLS$Observed- rf_MNF_BN_PLS$Predicted)^2))
  out[23,3]  <- sqrt(mean((rf_MNF_BN_RF$Observed- rf_MNF_BN_RF$Predicted)^2))
  out[24,3] <- sqrt(mean((rf_MNF_BN_SVM$Observed- rf_MNF_BN_SVM$Predicted)^2))
  
  # bias
  # pot Val
  out[1,4] <- 1 - coef( lm(pv_spect_PLS$Predicted~ pv_spect_PLS$Observed - 1) )
  out[2,4]  <- 1 - coef( lm(pv_spect_RF$Predicted~ pv_spect_RF$Observed - 1) )
  out[3,4] <- 1 - coef( lm(pv_spect_SVM$Predicted~ pv_spect_SVM$Observed - 1) )
  
  out[4,4] <- 1 - coef( lm(pv_spect_BN_PLS$Predicted~ pv_spect_BN_PLS$Observed - 1) )
  out[5,4]  <- 1 - coef( lm(pv_spect_BN_RF$Predicted~ pv_spect_BN_RF$Observed - 1) )
  out[6,4] <- 1 - coef( lm(pv_spect_BN_SVM$Predicted~ pv_spect_BN_SVM$Observed - 1) )
  
  out[7,4] <- 1 - coef( lm(pv_MNF_PLS$Predicted~ pv_MNF_PLS$Observed - 1) )
  out[8,4]  <- 1 - coef( lm(pv_MNF_RF$Predicted~ pv_MNF_RF$Observed - 1) )
  out[9,4] <- 1 - coef( lm(pv_MNF_SVM$Predicted~ pv_MNF_SVM$Observed - 1) )
  
  out[10,4] <- 1 - coef( lm(pv_MNF_BN_PLS$Predicted~ pv_MNF_BN_PLS$Observed - 1) )
  out[11,4] <- 1 - coef( lm(pv_MNF_BN_RF$Predicted~ pv_MNF_BN_RF$Observed - 1) )
  out[12,4] <- 1 - coef( lm(pv_MNF_BN_SVM$Predicted~ pv_MNF_BN_SVM$Observed - 1) )
  
  # rf
  out[13,4] <- 1 - coef( lm(rf_spect_PLS$Predicted~ rf_spect_PLS$Observed - 1) )
  out[14,4]  <- 1 - coef( lm(rf_spect_RF$Predicted~ rf_spect_RF$Observed - 1) )
  out[15,4] <- 1 - coef( lm(rf_spect_SVM$Predicted~ rf_spect_SVM$Observed - 1) )

  out[16,4] <- 1 - coef( lm(rf_spect_BN_PLS$Predicted~ rf_spect_BN_PLS$Observed - 1) )
  out[17,4]  <- 1 - coef( lm(rf_spect_BN_RF$Predicted~ rf_spect_BN_RF$Observed - 1) )
  out[18,4] <- 1 - coef( lm(rf_spect_BN_SVM$Predicted~ rf_spect_BN_SVM$Observed - 1) )
  
  out[19,4] <- 1 - coef( lm(rf_MNF_PLS$Predicted~ rf_MNF_PLS$Observed - 1) )
  out[20,4]  <- 1 - coef( lm(rf_MNF_RF$Predicted~ rf_MNF_RF$Observed - 1) )
  out[21,4] <- 1 - coef( lm(rf_MNF_SVM$Predicted~ rf_MNF_SVM$Observed - 1) )
  
  out[22,4] <- 1 - coef( lm(rf_MNF_BN_PLS$Predicted~ rf_MNF_BN_PLS$Observed - 1) )
  out[23,4]  <- 1 - coef( lm(rf_MNF_BN_RF$Predicted~ rf_MNF_BN_RF$Observed - 1) )
  out[24,4] <- 1 - coef( lm(rf_MNF_BN_SVM$Predicted~ rf_MNF_BN_SVM$Observed - 1) )
  
  if ( length(outputGOF[,1])==0 ){
    outputGOF <- as.data.frame(out)
  }
  if ( length(outputGOF[,1])!=0 ){
    outputGOF <- merge(outputGOF, as.data.frame(out), by = intersect(colnames(outputGOF), colnames(out)), all = TRUE)
  }
}

# erase empty "Species"
outputGOF <- outputGOF[- seq(1,24,1),]

save(outputGOF, file="outputGOF.RData")

# get GOF per model/tag
#######
x = grep("potVal", outputGOF$Validation)
gof_pv <- outputGOF[x, ]

x = grep("spect", gof_pv$Normalization)
gof_pv_spect <- gof_pv[x, ]
x = grep("BN", gof_pv_spect$Normalization)
gof_pv_spect_BN <- gof_pv_spect[x, ]
gof_pv_spect <- gof_pv_spect[-x, ]

x = grep("MNF", gof_pv$Normalization)
gof_pv_MNF <- gof_pv[x, ]
x = grep("BN", gof_pv_MNF$Normalization)
gof_pv_MNF_BN <- gof_pv_MNF[x, ]
gof_pv_MNF <- gof_pv_MNF[-x, ]

######
x = grep("rf", outputGOF$Validation)
gof_rf <- outputGOF[x, ]

x = grep("spect", gof_rf$Normalization)
gof_rf_spect <- gof_rf[x, ]
x = grep("BN", gof_rf_spect$Normalization)
gof_rf_spect_BN <- gof_rf_spect[x, ]
gof_rf_spect <- gof_rf_spect[-x, ]

x = grep("MNF", gof_rf$Normalization)
gof_rf_MNF <- gof_rf[x, ]
x = grep("BN", gof_rf_MNF$Normalization)
gof_rf_MNF_BN <- gof_rf_MNF[x, ]
gof_rf_MNF <- gof_rf_MNF[-x, ]




## plot
library(ggplot2)
library(gridExtra)

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

## Overall accuracy
plot_potVal <- ggplot(data = fit_potVal, mapping = aes(x = Models, y = OA, fill = factor(Normalization))) +
   geom_violin(mapping=aes(ymin = OA, ymax = OA), position = position_dodge(width = 0.6), size = 0.5) +
   stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + # fixt ylim
   ylab("Accuracy [0-1]") + xlab("Models") + guides(fill=FALSE) + # legend off
   theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(), 
                                  axis.title.x = element_blank()) +
   geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
   ggtitle("Pot validation method")

plot_rf <- ggplot(data = fit_rf, mapping = aes(x = Models, y = OA, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = OA, ymax = OA), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
  ylab("") + xlab("Models") + labs(fill="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
  ggtitle("Rip-it-of validation method")

### R²
ggmodelR2 <- na.omit( data.frame(r2=as.numeric(as.character(gof_pv$r2)), Models=gof_pv$Models, Normalization=gof_pv$Normalization) )

r2_potVal <- ggplot(data = ggmodelR2, mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + # fixt ylim
  ylab(expression(r^2)) + xlab("Models") + guides(fill=FALSE) + # legend off
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray")

ggmodelR2 <- na.omit( data.frame(r2=as.numeric(as.character(gof_rf$r2)), Models=gof_rf$Models, Normalization=gof_rf$Normalization) )

r2_rf <- ggplot(data = ggmodelR2, mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
  ylab("") + xlab("Models") + labs(fill="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray")

### RMSE
ggmodelRMSE <- na.omit( data.frame(RMSE=as.numeric(as.character(gof_pv$RMSE)), Models=gof_pv$Models, Normalization=gof_pv$Normalization) )

RMSE_potVal <- ggplot(data = ggmodelRMSE, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) + # fixt ylim
  ylab("RMSE [%]") + xlab("Models") + guides(fill=FALSE) + # legend off
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") 

ggmodelRMSE <- na.omit( data.frame(RMSE=as.numeric(as.character(gof_rf$RMSE)), Models=gof_rf$Models, Normalization=gof_rf$Normalization) )

RMSE_rf <- ggplot(data = ggmodelRMSE, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) + 
  ylab("RMSE [%]") + xlab("Models") + labs(fill="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray")

### bias
ggmodelbias <- na.omit( data.frame(bias=as.numeric(as.character(gof_pv$bias)), Models=gof_pv$Models, Normalization=gof_pv$Normalization) )

bias_potVal <- ggplot(data = ggmodelbias, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.4)) + # fixt ylim
  ylab("Bias") + xlab("Models") + guides(fill=FALSE) + # legend off
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
  geom_hline(yintercept = 0, lty = 2)

ggmodelbias <- na.omit( data.frame(bias=as.numeric(as.character(gof_rf$bias)), Models=gof_rf$Models, Normalization=gof_rf$Normalization) )

bias_rf <- ggplot(data = ggmodelbias, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.4)) + 
  ylab("") + xlab("Models") + labs(fill="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
  geom_hline(yintercept = 0, lty = 2)


plot_OA <- grid.arrange(plot_potVal, plot_rf, r2_potVal, r2_rf, 
                        RMSE_potVal, RMSE_rf, ncol=2, bias_potVal, bias_rf, 
                        nrow=4#, 
                        #widths = c(5, 6), heights = c(5,6) 
                        )

 