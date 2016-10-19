## R-Script - Analysis
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
DDir = "D:/Sp_Images"

setwd(home)

load("outputGOF.RData")

#### Source Functions from GitHub
source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
} 
source_github("https://raw.githubusercontent.com/JavierLopatin/Herbaceous-Species-Classification/master/Scripts/0_Functions.R")

## load the data
fit_potVal <-  read.table("Data/Fits_potVal.csv", sep = ",", header = T)
fit_rf     <-  read.table("Data/Fits_rf.csv", sep = ",", header = T)

## Obtain cover fits
setwd(DDir)

potVal_cover <- coverSummary("potVal")
rf_cover     <- coverSummary("rf")

potVal_cover$Species <- factor(potVal_cover$Species)
rf_cover$Species     <- factor(rf_cover$Species)

setwd(home)

save(potVal_cover, file="potVal_cover.RData")
save(rf_cover,     file="rf_cover.RData")

# Obtain R2, RMSE and Bias per species
outputGOF <- GOF(potVal_cover, rf_cover)
# erase empty "Species"
# <- outputGOF[- seq(1,24,1),]
outputGOF$r2 <- as.numeric(as.character(outputGOF$r2))
outputGOF$RMSE <- as.numeric(as.character(outputGOF$RMSE))
outputGOF$bias <- as.numeric(as.character(outputGOF$bias))

save(outputGOF, file="outputGOF.RData")


# get GOF per model/tag
x = grep("potVal", outputGOF$Validation)
gof_pv <- outputGOF[x, ]

x = grep("rf", outputGOF$Validation)
gof_rf <- outputGOF[x, ]

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
# Add a common legend for multiple ggplot2 graphs
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## Overall accuracy
plot_potVal <- ggplot(data = fit_potVal, mapping = aes(x = Models, y = Kappa, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = Kappa, ymax = Kappa), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +  ylab("Kappa [0-1]") + 
  ggtitle("Pot validation method")

plot_potVal <- plot_potVal+ theme(legend.position="none")

plot_rf <-ggplot(data = fit_rf, mapping = aes(x = Models, y = Kappa, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = Kappa, ymax = Kappa), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +   ylab("") + 
  ggtitle("Rip-it-of validation method")

plot_rf <- plot_rf + theme(legend.position="none")

### R2
ggmodelR2 <- na.omit( gof_pv )

r2_potVal <- ggplot(data = gof_pv, mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-0.03,1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") + ylab(expression(r^2)) 

r2_potVal <- r2_potVal+ theme(legend.position="none")

###
ggmodelR2 <- na.omit( data.frame(r2=as.numeric(as.character(gof_rf$r2)), Models=gof_rf$Models, Normalization=gof_rf$Normalization) )

r2_rf <- ggplot(data = ggmodelR2, mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-0.03,1), breaks = seq(0,1,0.2)) + 
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") + ylab("")

r2_rf <- r2_rf+ theme(legend.position="none")

### RMSE
ggmodelRMSE <- na.omit( data.frame(RMSE=as.numeric(as.character(gof_pv$RMSE)), Models=gof_pv$Models, Normalization=gof_pv$Normalization) )

RMSE_potVal <- ggplot(data = ggmodelRMSE, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-3.5,60), breaks = seq(0,60,20)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") + ylab("RMSE [%]")

RMSE_potVal <- RMSE_potVal + theme(legend.position="none")

ggmodelRMSE <- na.omit( data.frame(RMSE=as.numeric(as.character(gof_rf$RMSE)), Models=gof_rf$Models, Normalization=gof_rf$Normalization) )

RMSE_rf <- ggplot(data = ggmodelRMSE, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-3,60), breaks = seq(0,60,20)) + 
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray")

RMSE_rf <- RMSE_rf + theme(legend.position="none")

### bias
ggmodelbias <- na.omit( data.frame(bias=as.numeric(as.character(gof_pv$bias)), Models=gof_pv$Models, Normalization=gof_pv$Normalization) )

bias_potVal <- ggplot(data = ggmodelbias, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-0.5,1.05), breaks = seq(-0.4,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) + ylab("Bias") + xlab("Models") +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
  geom_hline(yintercept = 0, lty = 2) 

bias_potVal <- bias_potVal + theme(legend.position="none")

ggmodelbias <- na.omit( data.frame(bias=as.numeric(as.character(gof_rf$bias)), Models=gof_rf$Models, Normalization=gof_rf$Normalization) )

bias_rf <- ggplot(data = ggmodelbias, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) + aes(shape = factor(Normalization)) + 
  scale_y_continuous(limits = c(-0.5,1.05), breaks = seq(-0.4,1,0.2)) +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +  
  geom_hline(yintercept = 0, lty = 2) 

Legend <- get_legend(bias_rf)

# 3. Remove the legend from the box plot
bias_rf <- bias_rf + theme(legend.position="none") 

plot_fits <- grid.arrange(plot_potVal, plot_rf, r2_potVal, r2_rf, 
                        RMSE_potVal, RMSE_rf, bias_potVal, bias_rf, Legend,
                        nrow=4, ncol=3, layout_matrix = rbind(c(1,2,9), 
                                                              c(3,4,9),
                                                              c(5,6,9),
                                                              c(7,8,9)), 
                        heights = c(5.8,5,5,6), widths = c(5,5,1.5))
# Save
ggsave("Figures/Fits_all.pdf", plot_fits, width = 10, height = 8)
