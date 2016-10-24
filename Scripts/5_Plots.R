
##############################################
### plot variable importance at leaf level ###
##############################################

pdf(file = "Figures/VarImplot_LeafLevel.pdf", width=8, height=5)
plot.classificationEnsemble(hyperASD$spc, fitASD)
dev.off()

############################
### plot general results ###
############################

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

#############################
### scaterplot best model ###
#############################

well <- subset(complex_all, ClassPresence == "Well")
not_well <- subset(complex_all, ClassPresence == "Over" | ClassPresence == "Miss")

lm <- coef(lm(Observed ~ Predicted, data = well))

p1 <- ggplot(data = well, aes(x = Predicted, y = Observed, shape = factor(complex)
                                        , color = factor(complex))) +
  geom_point(alpha = 0.5, size = 2) + scale_colour_hue(l=50) +
  scale_y_continuous(limits = c(0,120), breaks = seq(0, 100, 020)) +
  scale_x_continuous(limits = c(0,120), breaks = seq(0, 100, 020)) +
  theme_bw(base_size=10) + xlab("Predicted [%]") + ylab("Observed [%]") +
  geom_abline(intercept = 0, lty = 2, lwd = 1.2) +
  geom_abline(intercept = 0, slope = lm[2], lwd=1.2) +
  geom_point(data = not_well, alpha = 0.5, size = 2) +
  theme(legend.title=element_blank()) + 
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
  geom_text(label="Misclassified", x=30, y=108, color="black", check_overlap = TRUE) +
  geom_text(label="Omitted", x=108, y= 30, color="black", check_overlap = TRUE)

scale_shape_discrete(breaks=c(1, 2, 3, 4),
                     labels=c("Complexity 1", "Complexity 2", "Complexity 3", "Complexity 4")) +
  
  
###########################
### Best model analysis ###
###########################

### Resolution 
## r2
p1 <- ggplot(data = na.omit(gof_resolution), mapping = aes(x = factor(Resolution), y = r2)) +
      geom_violin(mapping=aes(ymin = r2, ymax = r2), scale = "width", fill="gray90", trim=F) + 
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab(expression(r^2)) + xlab("Resolution [cm]") +
      scale_y_continuous(limits = c(-0.05,1), breaks = seq(0,1,0.2)) + # fixt ylim
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank())
## RMSE
p2 <- ggplot(data = na.omit(gof_resolution), mapping = aes(x = factor(Resolution), y = RMSE)) +
      geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("RMSE [%]") + xlab("Resolution [cm]") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank())
## Bias
p3 <- ggplot(data = na.omit(gof_resolution), mapping = aes(x = factor(Resolution), y = bias)) +
      geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("Bias") + xlab("Resolution [cm]") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) +
      geom_hline(yintercept = 0, lty = 2) 

### Analysis of architectural complexities
##r2
p4 <- ggplot(data = na.omit(gof_complex), mapping = aes(x = factor(complex), y = r2)) +
      geom_violin(mapping=aes(ymin = r2, ymax = r2),  fill="gray90", trim=F) + 
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab(expression(r^2)) + xlab("complex [cm]") +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + # fixt ylim
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(), 
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(), 
                                     axis.title.x = element_blank())
## RMSE
p5 <- ggplot(data = na.omit(gof_complex), mapping = aes(x = factor(complex), y = RMSE)) +
      geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("RMSE [%]") + xlab("complex [cm]") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(), 
                                     axis.title.x = element_blank())
## Bias
p6 <- ggplot(data = na.omit(gof_complex), mapping = aes(x = factor(complex), y = bias)) +
      geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("Bias") + xlab("complex [cm]") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(), 
                                     axis.title.x = element_blank()) +
      geom_hline(yintercept = 0, lty = 2) 

### analysis per cover percentage
##r2
p7 <- ggplot(data = na.omit(gof_cover), mapping = aes(x = factor(CoverRange), y = r2)) +
      geom_violin(mapping=aes(ymin = r2, ymax = r2),  fill="gray90", trim=F) + 
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab(expression(r^2)) + xlab("Cover [%]") +
      scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.2)) + # fixt ylim
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank())
## RMSE
p8 <- ggplot(data = na.omit(gof_cover), mapping = aes(x = factor(CoverRange), y = RMSE)) +
      geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("RMSE [%]") + xlab("Cover [%]") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank())
## Bias
p9 <- ggplot(data = na.omit(gof_cover), mapping = aes(x = factor(CoverRange), y = bias)) +
      geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("Bias") + xlab("Cover [%]") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) +
      geom_hline(yintercept = 0, lty = 2) 

### Analysis per Growth forms 
##r2
p10 <- ggplot(data = na.omit(PFT), mapping = aes(x = factor(PFT), y = r2)) +
      geom_violin(mapping=aes(ymin = r2, ymax = r2),  fill="gray90", trim=F) + 
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab(expression(r^2)) + xlab("Growth forms") +
      scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.2)) + # fixt ylim
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank())
## RMSE
p11 <- ggplot(data = na.omit(PFT), mapping = aes(x = factor(PFT), y = RMSE)) +
      geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("RMSE [%]") + xlab("Growth forms") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank())
## Bias
p12 <- ggplot(data = na.omit(PFT), mapping = aes(x = factor(PFT), y = bias)) +
      geom_violin(mapping=aes(ymin = bias, ymax = bias), scale = "width", fill="gray90", trim=F) +
      stat_summary(fun.y=mean, geom="line", aes(group=1), lty=2, col = "gray50", lwd = 1.5) +
      stat_summary(fun.data=data_summary) + ylab("Bias") + xlab("Growth forms") +
      theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) +
      geom_hline(yintercept = 0, lty = 2) 

plot_fits <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                          nrow=4, ncol=3 ,heights = c(5,4.5,5,5))
# Save
ggsave("Figures/Best_analysis.pdf", plot_fits, width = 10, height = 8)



###########################
### Per species results ###
###########################



