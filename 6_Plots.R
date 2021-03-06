
################################################################################
## R-Script - 5_Plots.R                                                       ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Mapping plant species in mixed grassland communities using     ##
##             close range imaging spectroscopy                               ##
##                                                                            ##
## description: Plots presented in the paper                                  ## 
##                                                                            ##
################################################################################

################
### Figure 3 ###
################

dat_rfVal  <- read.table("D:/Sp_Images/BootsClass_out/Fits_all_rf.txt", header = T)

library(beanplot)

pdf(file = "Figures/Figure3.pdf", width=7, height=6)
mat <- layout(rbind(c(1,1,1),c(2,3,4)), heights=c(1,1), TRUE) 
par(mai=c(0.6,0.7,0.3,0.3))
# Classification
beanplot(dat_rfVal$OA_SVM, dat_rfVal$Kappa_SVM, col = "gray", 
         border = NA, innerboerder=NA, beanlines="median", ll = 0, log="",  main = "", 
         names=c("Overall Accuracy", "Kappa"), ylab = "Accuracy [0-1]",
         ylim = c(0.4,1), yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1) 
text(x = c(0.95,0.95), labels = "*", cex = 3)
mtext("A", side=3, line=0.5, adj=0, cex=1.3)

# Cover prediction
# r2
par(mai=c(0.6,0.7,0,0.3))
beanplot( subset(gof_rf, Models=="SVM" & Normalization=="spectra")$r2,
          col = "gray", border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, log="",  main = "", names="", 
          ylab = expression(r^2), yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
mtext("B", side=3, line=0.5, adj=0, cex=1.3)
# RMSE
par(mai=c(0.6,0.7,0,0.3))
beanplot( subset(gof_rf, Models=="SVM" & Normalization=="spectra")$RMSE,
          col = "gray", border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, log="",  main = "", names="", 
          ylab = "RMSE [%]", yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
mtext("C", side=3, line=0.5, adj=0, cex=1.3)
# bias
par(mai=c(0.6,0.7,0,0.3))
beanplot( subset(gof_rf, Models=="SVM" & Normalization=="spectra")$bias,
          col = "gray", border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, log="",  main = "", names="", 
          ylab = "Bias", yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
mtext("D", side=3, line=0.5, adj=0, cex=1.3)

dev.off()

################
### Figure 4 ###
################

library(ggplot2)
library(gridExtra)

well <- subset(complex_all, ClassPresence == "Well")
well <- well[order(well$complex, decreasing = T), ]
not_well <- subset(complex_all, ClassPresence == "Over" | ClassPresence == "Miss")
not_well <- not_well[order(not_well$complex, decreasing = T), ]

# pallete
mycolors <- (c("blue", "green", "orange", "red"))


p1 <- ggplot(data = well, aes(x = Predicted, y = Observed, color = factor(complex))) + 
  geom_point(alpha=0.5, aes(size = factor(complex))) + 
  scale_colour_manual(values=mycolors) +
  scale_y_continuous(limits = c(0,120), breaks = seq(0, 100, 020)) +
  scale_x_continuous(limits = c(0,120), breaks = seq(0, 100, 020)) +
  theme_bw(base_size=10) + xlab("Predicted [%]") + ylab("Observed [%]") +
  geom_abline(intercept = 0, lty = 2, lwd = 1.2) +
  geom_smooth(method = "lm", se=F, fullrange=T) +
  # add false positives and negatives
  geom_point(data = not_well2, alpha=0.5, aes(size = factor(complex))) +
  # add text
  geom_segment(aes(x=3, y=0, xend=3, yend=95), color = "black", lwd = 1.3) +
  geom_segment(aes(x=3, y=95, xend=0, yend=95), color = "black", lwd = 1.3) +
  geom_segment(aes(x=3, y=0, xend=0, yend=0), color = "black", lwd = 1.3) +
  geom_segment(aes(x=0, y=3, xend=115, yend=3), color = "black", lwd = 1.3) + 
  geom_segment(aes(x=0, y=3, xend=0, yend=0), color = "black", lwd = 1.3) +
  geom_segment(aes(x=115, y=3, xend=115, yend=0), color = "black", lwd = 1.3) +
  geom_segment(aes(x=100, y=20, xend=80, yend=0), color = "black", lwd = 0.8, 
               arrow = arrow()) +
  geom_segment(aes(x=20, y=100, xend=0, yend=80), color = "black", lwd = 0.8, 
               arrow = arrow()) +
  geom_text(label="False negatives", x=30, y=108, color="black", check_overlap = TRUE) +
  geom_text(label="False positives", x=108, y= 25, color="black", check_overlap = TRUE) +
  theme(panel.grid.major.x = element_blank(), legend.key = element_blank())+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=14)) +
  theme(panel.border = element_blank(), axis.line.x = element_line(),axis.line.y = element_line()) + 
  labs(colour="Complexity\ngradient", size="Complexity\ngradient") +
  theme(legend.text=element_text(size=10)) +
  ggtitle("A") + theme(plot.title = element_text(size=18, hjust = 0))

p2 <- ggplot(data = complex_all2, aes(x = complex, fill = factor(ClassPresence))) +
  geom_bar() + 
  theme(legend.title=element_blank(), panel.grid.major.x = element_blank(), legend.key = element_blank()) +
  theme_bw(base_size=10) + xlab("Complexity gradient") + ylab("Frequency") + labs(fill="") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=14)) +
  theme(panel.border = element_blank(), axis.line.x = element_line(),axis.line.y = element_line()) +  
  scale_fill_grey(start = 0, end = .9, labels=c("False negatives", "False positives", "Correctly classified")) +
  theme(legend.text=element_text(size=9)) +
  ggtitle("B") + theme(plot.title = element_text(size=18, hjust = 0))

p3  <- grid.arrange(p1, p2, nrow=1, ncol=2) 
                                                                
ggsave("Figures/Figure4.pdf", p3, width = 12, height = 4.5)

 

################
### Figure 7 ###
################
  
### Resolution 
## r2
p1 <- ggplot(data = na.omit(gof_resolution), mapping = aes(x = factor(Resolution), y = r2)) +
   geom_violin(mapping=aes(ymin = r2, ymax = r2), scale = "width", fill="gray90", trim=T) + 
   stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
   stat_summary(fun.data=data_summary) + labs(y=expression(r^2), x="Resolution [cm]", title="Square Pearson´s \n correlation coefficient \n") +
   scale_y_continuous(limits = c(-0.05,1), breaks = seq(0,1,0.2)) + # fixt ylim
   theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                  axis.text.y=element_text(size=14), text = element_text(size=15), 
                                  plot.title = element_text(hjust = 0.5))
## RMSE
p2 <- ggplot(data = na.omit(gof_resolution), mapping = aes(x = factor(Resolution), y = RMSE)) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=T) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="RMSE [%]",  x="Resolution [cm]", title="Root mean square error \n \n A: Spatial analysis") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## Bias
p3 <- ggplot(data = na.omit(gof_resolution), mapping = aes(x = factor(Resolution), y = bias)) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=T) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="Bias", x="Resolution [cm]", title="Bias \n \n ") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2) 

### Analysis of architectural complexities
##r2
p4 <- ggplot(data = na.omit(gof_complex), mapping = aes(x = factor(complex), y = r2)) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2),  fill="gray90", trim=T) + 
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) +  labs(y=expression(r^2), x="Complexity gradient", title="") + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## RMSE
p5 <- ggplot(data = na.omit(gof_complex), mapping = aes(x = factor(complex), y = RMSE)) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=T) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="RMSE [%]", x="Complexity gradient", title="B: Complexity gradient") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## Bias
p6 <- ggplot(data = na.omit(gof_complex), mapping = aes(x = factor(complex), y = bias)) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=T) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="Bias", x="Complexity gradient", title="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))+
geom_hline(yintercept = 0, lty = 2) 

### analysis per cover percentage
##r2
p7 <- ggplot(data = na.omit(gof_cover), mapping = aes(x = factor(CoverRange), y = r2)) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2),  fill="gray90", trim=T) + 
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y=expression(r^2), x="Cover [%]", title="") +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## RMSE
p8 <- ggplot(data = na.omit(gof_cover), mapping = aes(x = factor(CoverRange), y = RMSE)) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=T) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="RMSE [%]", x="Cover [%]", title="C: Species cover") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## Bias
p9 <- ggplot(data = na.omit(gof_cover), mapping = aes(x = factor(CoverRange), y = bias)) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=T) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) +  labs(y="Bias", x="Cover [%]", title="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2) 

### Analysis per Growth forms 
##r2
p10 <- ggplot(data = na.omit(PFT), mapping = aes(x = factor(PFT), y = r2)) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2),  fill="gray90", trim=F) + 
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y=expression(r^2), x="Functional types", title="") +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## RMSE
p11 <- ggplot(data = na.omit(PFT), mapping = aes(x = factor(PFT), y = RMSE)) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=F) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="RMSE [%]", x="Functional types", title="D: Functional types") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15))
## Bias
p12 <- ggplot(data = na.omit(PFT), mapping = aes(x = factor(PFT), y = bias)) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=F) +
  stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
  stat_summary(fun.data=data_summary) + labs(y="Bias", x="Functional types", title="") +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14), text = element_text(size=15)) +
  geom_hline(yintercept = 0, lty = 2) 

plot_fits <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                          nrow=4, ncol=3 ,heights = c(5,4.5,5,5))

# Save
ggsave("Figures/Figure7.pdf", plot_fits, width = 13, height = 10)


#################
### Figure S1 ###
#################

dat_potVal <- read.table("D:/Sp_Images/BootsClass_out/Fits_all_pot.txt", header = T)
dat_rfVal  <- read.table("D:/Sp_Images/BootsClass_out/Fits_all_rf.txt", header = T)

library(beanplot)

pdf(file = "Figures/Figure2.pdf", width=7, height=6)
mat <- layout(rbind(c(1,1,1),c(2,3,4)), heights=c(1,1), TRUE) 
par(mai=c(0.6,0.7,0.3,0.3))
# Classification
beanplot(dat_potVal$OA_SVM, dat_rfVal$OA_SVM, dat_potVal$Kappa_SVM, dat_rfVal$Kappa_SVM, col = list("black", "gray"), 
         border = NA, innerboerder=NA, beanlines="median", ll = 0, side = "b", log="",  main = "", 
         names=c("Overall Accuracy", "Kappa"), ylab = "Accuracy [0-1]",
         ylim = c(0.4,1), yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1) 
legend("bottomright", legend=c("Pot method", "Isolation method"), fill=c("black", "gray"), bty="n", cex=1.3)
text(x = c(0.95,0.95), labels = "*", cex = 3)
mtext("A", side=3, line=0.5, adj=0, cex=1.3)

# Cover prediction
# r2
par(mai=c(0.6,0.7,0,0.3))
beanplot( subset(gof_pv, Models=="SVM" & Normalization=="spectra")$r2 , 
          subset(gof_rf, Models=="SVM" & Normalization=="spectra")$r2,
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, side = "b", log="",  main = "", names="", 
          ylab = expression(r^2), yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
mtext("B", side=3, line=0.5, adj=0, cex=1.3)
# RMSE
par(mai=c(0.6,0.7,0,0.3))
beanplot( subset(gof_pv, Models=="SVM" & Normalization=="spectra")$RMSE , 
          subset(gof_rf, Models=="SVM" & Normalization=="spectra")$RMSE,
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, side = "b", log="",  main = "", names="", 
          ylab = "RMSE [%]", yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
mtext("C", side=3, line=0.5, adj=0, cex=1.3)
# bias
par(mai=c(0.6,0.7,0,0.3))
beanplot( subset(gof_pv, Models=="SVM" & Normalization=="spectra")$bias, 
          subset(gof_rf, Models=="SVM" & Normalization=="spectra")$bias,
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, side = "b", log="",  main = "", names="", 
          ylab = "Bias", yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
mtext("D", side=3, line=0.5, adj=0, cex=1.3)

dev.off()

#################
### Figure S2 ###
#################

load("mat_mrpp.RData")
load("Gramm_Imp.RData")
load("Fobs_Imp.RData")
#load("results_leaf/ClassLeafLevel.Rdata")

pdf(file = "Figures/VarImport.pdf", width=8, height=6)
mat <- layout(matrix(1:4,ncol=1), heights=c(1,0.3,0.3,0.5), TRUE) 
#layout.show(mat)
#
par(mai=c(0,1,0.1,0.1))
plot.spectra(spectra = rf_spec[, 3:length(rf_spec)]/10000,
             wl = wl,
             xaxis=F, ylab = T, ymax=0.6)

par(mai=c(0.1,1,0,0.1))
plot.importance(gram_mrpp, wl, FALSE)

par(mai=c(0.1,1,0,0.1))
plot.importance(frobs_mrpp, wl, FALSE)

par(mai=c(0.1,1,0,0.1))
plot.importance(mat_mrpp, wl, TRUE)

dev.off()

# legend 
colfunc <- c( colorRampPalette(brewer.pal(9,"Blues"))(100) )
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

library(fields)

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
colorbar.plot( 1, 0, 1:255, col = colfunc, strip.length = 2, strip.width = 0.3)


#################
### Figure S3 ###
#################

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

## Overall accuracy all models (Supplementary data)
p1 <- ggplot(data = fit_potVal, mapping = aes(x = Models, y = Kappa, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = Kappa, ymax = Kappa), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +  ylab("Kappa [0-1]") + 
  ggtitle("Pot validation method")

p1 <- p1+ theme(legend.position="none")

p2 <- ggplot(data = fit_rf, mapping = aes(x = Models, y = Kappa, fill = factor(Normalization))) +
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

p2 <- p2 + theme(legend.position="none")

### R2
p3 <- ggplot(data = na.omit(gof_pv), mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = r2, ymax = r2), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-0.03,1), breaks = seq(0,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") + ylab(expression(r^2)) 

p3 <- p3+ theme(legend.position="none")

###
p4 <- ggplot(data = na.omit(gof_rf), mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
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

p4 <- p4+ theme(legend.position="none")

### RMSE
p5 <- ggplot(data = gof_pv, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-3.5,60), breaks = seq(0,60,20)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") + ylab("RMSE [%]")

p5 <- p5 + theme(legend.position="none")

p6 <- ggplot(data = gof_rf, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
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

p6 <- p6 + theme(legend.position="none")

### bias
p7 <- ggplot(data = gof_pv, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-0.5,1.05), breaks = seq(-0.4,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) + ylab("Bias") + xlab("Models") +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
  geom_hline(yintercept = 0, lty = 2) 

p7 <- p7 + theme(legend.position="none")

p8 <- ggplot(data = gof_rf, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6)) + aes(shape = factor(Normalization)) + 
  scale_y_continuous(limits = c(-0.5,1.05), breaks = seq(-0.4,1,0.2)) +
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), legend.key = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +  
  geom_hline(yintercept = 0, lty = 2) 

Legend <- get_legend(p8)

# 3. Remove the legend from the box plot
p8 <- p8 + theme(legend.position="none") 

plot_fits <- grid.arrange(p1, p2, p3, p4, 
                          p5, p6, p7, p8, Legend,
                          nrow=4, ncol=3, layout_matrix = rbind(c(1,2,9), 
                                                                c(3,4,9),
                                                                c(5,6,9),
                                                                c(7,8,9)), 
                          heights = c(5.8,5,5,6), widths = c(5,5,1.5))
# Save
ggsave("Figures/Fits_all.pdf", plot_fits, width = 10, height = 8)

