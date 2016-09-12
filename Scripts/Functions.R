
## Functions used in the paper


################################################################################
##                                                                            ##
## Multi-method ensemble selection of spectral bands                          ##
##                                                                            ##
## This function performs a band selection based on a multi-method ensemble   ##
## assessment of the variable importance and classification coefficients of   ##
## three different model types: Partial Least Squares Discriminant Analysis,  ##
## Random Forest and Support Vector Machine classifications                   ## 
##                                                                            ## 
## Arguments:                                                                 ##
## x        Numeric matrix containing the spectra (samples as rows)           ##
## y        Numeric vector containing the response variable                   ##
## wl       Numeric vector containing the wavelength information of the bands ##
##                                                                            ##
## function based on the paper:                                               ##
## Feilhauer, H., Asner, G. P., & Martin, R. E. (2015). Multi-method ensemble ##
## selection of spectral bands related to leaf biochemistry. Remote Sensing of## 
## Environment, 164(November), 57-65. http://doi.org/10.1016/j.rse.2015.03.033##                                           ##
##                                                                            ##
################################################################################

classificationEnsemble <- function(classes, spec, wl=NA){
  
  ## load required libraries
  library (caret)
  library(e1071)
  library(doParallel)
  
  # set data
  data <- data.frame(classes, spec)
  data <- na.omit(data)
  
  # Set the random number seed so we can reproduce the results
  set.seed(123)
  # Split data in training and test
  forTraining <- createDataPartition(data$classes, p = 0.6, list=F)
  train <- data [ forTraining,]
  test<- data [-forTraining,]
  
  # Each model used 5 repeated 10-fold cross-validation. Use AUC to pick the best model
  controlObject <- trainControl(method = "cv", number = 10,  repeats=10, classProbs=TRUE, allowParallel = TRUE)
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  #############################
  ### PLS-DA classification ###
  #############################
  
  # apply classification
  set.seed(123)
  plsClas <- train(x=train[, 2:length(train)], y=train$classes, method = "pls", tuneLength=20, 
                   preProc = c("center", "scale"), trControl = controlObject)
  
  # predict
  pls.pred <- predict(plsClas, test[, 2:length(train)])
  
  # confusion matix
  preMatrix <- function(pred, test){ # functionn to prevent caret error for different length
    u = union(pred, test)
    t = table(factor(pred, u), factor(test, u))
    return(t)
  }  
  conf.pls <- confusionMatrix(preMatrix(pls.pred, test$classes))
  
  # get accuracies
  PA.pls    <- conf.pls$byClass[,3] 
  UA.pls    <- conf.pls$byClass[,4]
  OA.pls    <- conf.pls$overall["Accuracy"]
  kappa.pls <- conf.pls$overall["Kappa"]
  
  ### variable importance
  plscf <- as.vector(rowMeans(plsClas$finalModel$coefficients)) ## extract coeff.
  plscf <- plscf / sd (plscf) ## scale regression coefficients
  
  #########################
  ### RF classification ###
  #########################
  
  set.seed(123)
  rfClas <- train(x=train[, 2:length(train)], y=train$classes, method = "rf", tuneLength=15, trControl = controlObject)
  
  # predict
  rf.pred <- predict(rfClas, test[,2:length(train)])
  
  # confusion matix
  conf.rf <- confusionMatrix(preMatrix(rf.pred, test$classes))
  
  # get accuracies
  PA.rf    <- conf.rf$byClass[,3] 
  UA.rf    <- conf.rf$byClass[,4]
  OA.rf    <- conf.rf$overall["Accuracy"]
  kappa.rf <- conf.rf$overall["Kappa"]
  
  ### variable importance
  rfcf <- varImp(rfClas$finalModel)[[1]]
  rfcf <- as.vector (rfcf / sd(rfcf))
  
  ##########################
  ### SVM classification ###
  ##########################
  
  set.seed(123)
  svmClas <- train(x=train[, 2:length(train)], y=train$classes, method = "svmLinear2", tuneLength=15, 
                   preProc = c("center", "scale"), trControl = controlObject)
  
  # predict
  svm.pred <- predict(svmClas, test[,2:length(train)])
  
  # confusion matix
  conf.svm <- confusionMatrix(preMatrix(svm.pred, test$classes))
  
  # get accuracies
  PA.svm    <- conf.svm$byClass[,3] 
  UA.svm    <- conf.svm$byClass[,4]
  OA.svm    <- conf.svm$overall["Accuracy"]
  kappa.svm <- conf.svm$overall["Kappa"]
  
  ### variable importance
  vr.alpha <- t(svmClas$finalModel$coefs) ## extract alpha vector
  svr.index <-  svmClas$finalModel$index ## extract alpha index
  ## calculate pseudo-regression coefficients from the alpha vector
  svrcf <- numeric (ncol (spec))
  for(i in 1:ncol(spec)){
    svrcf[i] <- svr.alpha %*% spec[svr.index, i]
  } 
  svrcf <- svrcf / sd (svrcf) ## scale pseudo-coefficients
  
  #####################################################################    
  ### get ensemble from all models and identify important variables ###
  #####################################################################
  
  ## get ensemble from all models and identify important variables
  ensemblecf <- abs(plscf) * OA.pls + abs(rfcf) * OA.rf + abs(svrcf) * OA.svm 
  th <- mean(ensemblecf) + sd(ensemblecf) ## calculate threshold
  r <- ensemblecf > th ## apply threshold
  
  # stop parallel process
  stopCluster(cl)
  
  ######################
  ### prepare output ###
  ######################
  
  cf <- rbind (wl, plscf, rfcf, svrcf, ensemblecf, selbands)
  colnames(cf) <- colnames(spec)
  
  fit <- c (OA.pls, OA.rf, OA.svm)
  names (fit) <- c ("PLS-DA OA", "RF OA", "SVR OA")
  output <- list (cf, fit, th, plsClas, rfClas, svmClas)
  names (output) <- c ("selection", "fits", "threshold", "PLS", "RF", "SVM")
  class (output) <- "ensemble"
  output
  
}

################################################################################
#################### END function classificationEnsemble #######################
################################################################################


################################################################################
##                                                                            ##
## plot.ensemble: visualization of ensemble objects                           ##
##                                                                            ##
## Arguments:                                                                 ##
## spec     spectral information. Used to create the quantiles of spectra     ##
## en       classificationEnsemble object                                     ##
##                                                                            ##
################################################################################

plot.classificationEnsemble <- function (spec, en, label=TRUE, ...) {
  # extract the data from the classification Ensamble function
  wl <- en[[1]][1,]
  cf <- en[[1]][2:4,]
  fit <- en[[2]]
  cf <- cf * fit
  fit <- round (fit, 2)
  encf <- en[[1]][5,]
  z1 <- matrix (rep (encf, 100), ncol=100)
  sel <- as.logical (en[[1]][6,])
  z2 <- matrix (rep (sel, 100), ncol=100)
  z2[,11:100] <- NA
  z2[z2==0] <- NA
  # obtain the quantiles of the spectras
  quant <- apply(spec, 2, quantile, probs =c(0.05, 0.25, 0.5, 0.75, 0.95))
  # Plot the spectra
  par (mar=c (5, 4, 4, 7) + 0.1, xpd=NA)
  plot(wl, quant[1,], type="l", ylim = c(0,0.8), axes=F, ylab = NA, xlab=NA, las=1)
  lines(wl, quant[2,], type="l")
  lines(wl, quant[3,], type="l")
  lines(wl, quant[4,], type="l")
  lines(wl, quant[5,], type="l")
  polygon(c(wl, rev(wl)), c(quant[2,], rev(quant[1,])), col = "grey70")
  polygon(c(wl, rev(wl)), c(quant[3,], rev(quant[2,])), col = "grey50")
  polygon(c(wl, rev(wl)), c(quant[4,], rev(quant[3,])), col = "grey50")
  polygon(c(wl, rev(wl)), c(quant[5,], rev(quant[4,])), col = "grey70")
  axis(side = 4, las=1)
  mtext(side = 4, line = 3, 'Reflectance')
  # add coefficients
  par(new = T)
  plot(wl,  cf[1,], type = "l", col=2, xlab=expression(lambda(nm)), ylab="Weighted coefficients", las=1,
       ylim=c(min(cf), max(cf)))
  lines(wl,  cf[2,], type = "l", col=3)
  lines(wl,  cf[3,], type = "l", col=4) 
  # add spectral bands selected by the ensemble method
  image (wl, seq(min (cf) * 1.1, max (cf) * 1.1, length.out=100), z2, col=7, 
         xlab="", ylab="", add=T)
  # Labels
  if(label==TRUE){
    labels <- c (paste (c ("PLS ", "RF ", "SVM "), "OA", "=", fit, sep=""), NA, 
                 "Ensemble selection")
    legend ("topright", bty="n", col=c (2, 3, 4, NA, NA, rep (1, 3)), 
            pt.bg=c(rep (NA, 4), 7), lwd=c(rep (2, 3), rep (NA, 2)),
            pch=c (rep (NA, 4), 22), cex=0.7, pt.cex=1, legend=labels)
   }
 }

################################################################################
################### END function plot.classificationEnsemble ###################
################################################################################


###

rasterListNames <- function(fileExtantion, folder){
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # delete the ".dat" from the name
  rast_list = gsub('.{4}$', '', rast_list)
  return(rast_list)
}

rasterList <- function(fileExtantion, folder, rasterNames=NULL){
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # import rasters
  setwd(file.path(home, folder))
  rasterlist <- list()
  for(i in 1:length(rast_list)){
    rast <- stack(rast_list[i])
    if (is.null(rasterNames)){
      names(rast) <- rasterNames
      }
    rasterlist[[i]] <- rast
  }
  setwd(home)
  return(rasterlist)
}


GLCM <- function(img){
  # Estimate the GLCM values of one band
  texture <- glcm(img, n_grey = 32, window = c(3, 3), shift = c(1, 1), 
                  statistics =c("homogeneity", "contrast", "dissimilarity", "entropy","second_moment", "correlation"),
                  min_x=NULL, max_x=NULL, na_opt="any",na_val=NA, scale_factor=1, asinteger=FALSE)
  stats <- c("homogeneity", "contrast", "dissimilarity", "entropy","second_moment", "correlation")
  Names <- paste(names(img), "_", stats, sep="")
  names(texture) <- Names
  return (texture)
}


