
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
  ensemblecf <- abs (plscf) * OA.pls + abs (svrcf) * OA.svm + abs (rfcf) * OA.rf
  th <- mean (ensemblecf) + sd (ensemblecf) ## calculate threshold
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
############################# END function ensemble ############################
################################################################################


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


