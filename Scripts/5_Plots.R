


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
r2_potVal <- ggplot(data = na.omit(gof_pv), mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
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
r2_rf <- ggplot(data = na.omit(gof_rf), mapping = aes(x = Models, y = r2, fill = factor(Normalization))) +
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
RMSE_potVal <- ggplot(data = gof_pv, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = RMSE, ymax = RMSE), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-3.5,60), breaks = seq(0,60,20)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(), 
                                 axis.title.x = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") + ylab("RMSE [%]")

RMSE_potVal <- RMSE_potVal + theme(legend.position="none")

RMSE_rf <- ggplot(data = gof_rf, mapping = aes(x = Models, y = RMSE, fill = factor(Normalization))) +
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
bias_potVal <- ggplot(data = gof_pv, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
  geom_violin(mapping=aes(ymin = bias, ymax = bias), position = position_dodge(width = 0.6), size = 0.5) +
  stat_summary(fun.data=data_summary, position = position_dodge(width = 0.6), aes(shape = factor(Normalization))) +
  scale_y_continuous(limits = c(-0.5,1.05), breaks = seq(-0.4,1,0.2)) + # fixt ylim
  theme_bw(base_size=10) + theme(panel.grid.major.x = element_blank()) + ylab("Bias") + xlab("Models") +
  geom_vline(xintercept = c(1.5, 2.5), colour = "gray") +
  geom_hline(yintercept = 0, lty = 2) 

bias_potVal <- bias_potVal + theme(legend.position="none")

bias_rf <- ggplot(data = gof_rf, mapping = aes(x = Models, y = bias, fill = factor(Normalization))) +
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

