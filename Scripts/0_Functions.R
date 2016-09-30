################################################################################
## R-Script - 0_Functions.R                                                   ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## description: 
##
## Manuscript: 
##
## last changes: 
##                                                                            ##
################################################################################

##----------------------------------------------------------------------------##
##                                                                            ##
## Multi-method ensemble selection of spectral bands                          ##
##                                                                            ##
## This function performs a band selection based on a multi-method ensemble   ##
## assessment of the variable importance and classification coefficients of   ##
## three different model types: Partial Least Squares Discriminant Analysis,  ##
## Random Forest and Support Vector Machine classifications                   ## 
##                                                                            ## 
## Arguments:                                                                 ##
## - x      Numeric matrix containing the spectra (samples as rows)           ##
## - y      Numeric vector containing the response variable                   ##
## - wl     Numeric vector containing the wavelength information of the bands ##
##                                                                            ##
## function based on the paper:                                               ##
## Feilhauer, H., Asner, G.P., & Martin, R.E. (2015). Multi-method ensemble   ##
## selection of spectral bands related to leaf biochemistry. Remote Sensing   ## 
## of Environment, 164, 57-65. http://doi.org/10.1016/j.rse.2015.03.033       ##
##                                                                            ##
##----------------------------------------------------------------------------##

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
  svr.alpha <- t(svmClas$finalModel$coefs) ## extract alpha vector
  svr.alpha <- colMeans(svr.alpha)
  svr.index <-  svmClas$finalModel$index ## extract alpha index
  ## calculate pseudo-regression coefficients from the alpha vector
  svrcf <- numeric (ncol (spec))
  for(i in 1:ncol(spec))
    svrcf[i] <- svr.alpha %*% spec[svr.index, i]
  svrcf <- svrcf / sd (svrcf) ## scale pseudo-coefficients
  
  #####################################################################    
  ### get ensemble from all models and identify important variables ###
  #####################################################################
  
  ## get ensemble from all models and identify important variables
  ensemblecf <- abs(plscf) * OA.pls + abs(rfcf) * OA.rf + abs(svrcf) * OA.svm 
  th <- mean(ensemblecf) + sd(ensemblecf) ## calculate threshold
  selbands <- ensemblecf > th ## apply threshold
  
  # stop parallel process
  stopCluster(cl)
  
  ######################
  ### prepare output ###
  ######################
  
  cf <- rbind (wl, plscf, rfcf, svrcf, ensemblecf, selbands)
  colnames(cf) <- colnames(spec)
  
  fit <- c (OA.pls, OA.rf, OA.svm)
  names (fit) <- c ("PLS-DA OA", "RF OA", "SVR OA")
  output <- list (cf, fit, th, plsClas, rfClas, svmClas, conf.pls, conf.rf, conf.svm)
  names (output) <- c ("selection", "fits", "threshold", "PLS", "RF", "SVM", "confusionMatrix.PLS", "confusionMatrix.RF", "confusionMatrix.SVM")
  class (output) <- "ensemble"
  output
  
}


##----------------------------------------------------------------------------##
##                                                                            ##
## ApplyBestClassification:                                                   ##
##                                                                            ##
## Apply the best model from classificationEnsemble to the plots using        ##
## a bootstrapping procidure. Exported results are in shapefile format        ## 
##                                                                            ##
## Arguments:                                                                 ##
## - spec     spectral information. Used to create the quantiles of spectra   ##
## - en       classificationEnsemble object                                   ##
##                                                                            ##
##----------------------------------------------------------------------------##

ApplyBootsClassification <- function(classes, spec, en, rasterPlots, boots=100, outDir, modelTag){  
  
  library(raster)
  library(rgdal)
  library(caret)
  library(e1071)
  library(randomForest)
  library(pls)
  library(doParallel)
  
  # extract the data from the classification Ensamble function
  wl <- en[[1]][1,]
  
  ncomp = en$PLS$finalModel$ncomp
  probMethod = en$PLS$finalModel$probMethod
  
  bestNtree = en$RF$finalModel$ntree
  bestMtry = en$RF$finalModel$mtry
  
  bestCost <- en$SVM$finalModel$cost
  bestGamma <- en$SVM$finalModel$gamma
  
  data <- data.frame(classes, spec)
  
  ## apply funtion to predict species cover with SVM
  
  # list of accuracies
  PA.PLS    <- list()
  UA.PLS    <- list()
  OA.PLS    <- list()
  kappa.PLS <- list()
  model.PLS <- list()
  predict.PLS <- list()
  rasterList.PLS <- list()
  shpList.PLS <- list()
  
  PA.RF    <- list()
  UA.RF    <- list()
  OA.RF    <- list()
  kappa.RF <- list()
  model.RF <- list()
  predict.RF <- list()
  rasterList.RF <- list()
  shpList.RF <- list()
  
  PA.SVM    <- list()
  UA.SVM    <- list()
  OA.SVM    <- list()
  kappa.SVM <- list()
  model.SVM <- list()
  predict.SVM <- list()
  rasterList.SVM <- list()
  shpList.SVM <- list()
  
  OBS <- list()
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  for (i in 1:boots){
    
    N = length(data[,1])
    
    # create random numbers with replacement to select samples from each group
    idx = sample(1:N, N, replace=TRUE)
    
    train <- data[idx,]
    val <- data[-idx,]
    
    # store and select the observations
    obs <- val$classes
    OBS[[i]]<-obs
    
    #################
    ### Apply PLS ###
    #################
    
    PLS  <- plsda(x =  train[,2:length(train)], y = as.factor( train$classes ), ncomp = ncomp, 
                  probMethod = probMethod)
    # use all data to apply to image
    PLS2 <- plsda(x =  data[,2:length(data)], y = as.factor( data$classes ), ncomp = ncomp, 
                  probMethod = probMethod)
    
    # sotore tunning model
    model.PLS[[i]] <- PLS
    
    # predict
    pred    <- predict(PLS, val[,2:length(val)])
    predict.PLS[[i]] <- pred
    
    # confusion matix
    conf   <- confusionMatrix(pred, val$classes)
    
    # get accuracies
    PA.PLS[[i]]       <- conf$byClass[,3] 
    UA.PLS[[i]]       <- conf$byClass[,4]
    OA.PLS[[i]]       <- conf$overall["Accuracy"]
    kappa.PLS[[i]]    <- conf$overall["Kappa"]
    
    #################
    ### Apply SVM ###
    #################
    
    RF  <- randomForest( y = as.factor( train$classes ), x = train[,2:length(train)],
                         ntree= bestNtree, mtry = bestMtry)
    # use all data to apply to image
    RF2 <- randomForest( y = as.factor( data$classes ), x = data[,2:length(data)],
                         ntree= bestNtree, mtry = bestMtry)
    
    # sotore tunning model
    model.RF[[i]] <- RF
    
    # predict
    pred    <- predict(RF, val[,2:length(val)])
    predict.RF[[i]] <- pred
    
    # confusion matix
    conf   <- confusionMatrix(pred, val$classes)
    
    # get accuracies
    PA.RF[[i]]       <- conf$byClass[,3] 
    UA.RF[[i]]       <- conf$byClass[,4]
    OA.RF[[i]]       <- conf$overall["Accuracy"]
    kappa.RF[[i]]    <- conf$overall["Kappa"]

    #################
    ### Apply SVM ###
    #################
    
    SVM  <- svm(train[,2:length(train)], train$classes, kernel = "linear",
                gamma = bestGamma, cost = bestCost, probability = TRUE)
    # use all data to apply to image
    SVM2 <- svm(data[,2:length(data)], data$classes, kernel = "linear",
                gamma = bestGamma, cost = bestCost, probability = TRUE)
    
    # sotore tunning model
    model.SVM[[i]] <- SVM
    
    # predict
    pred    <- predict(SVM, val[,2:length(val)])
    predict.SVM[[i]] <- pred
    
    # confusion matix
    conf   <- confusionMatrix(pred, val$classes)
    
    # get accuracies
    PA.SVM[[i]]       <- conf$byClass[,3] 
    UA.SVM[[i]]       <- conf$byClass[,4]
    OA.SVM[[i]]       <- conf$overall["Accuracy"]
    kappa.SVM[[i]]    <- conf$overall["Kappa"]
    
    
    ##############################
    ### Apply models to raster ### 
    ##############################
    
    dummyList.PLS <- list()
    dummyList.RF <- list()
    dummyList.SVM <- list()
    
    for (j in 1:length(rasterPlots)){
      raster = rasterPlots[[j]]
      names(raster) <- paste( rep("B", 61), seq(1,61,1),  sep="" )
      # mask out raster zones with NDVI below 0.3
      red <- raster[[31]]
      Ired <- raster[[43]]
      NDVI <- (Ired-red)/(Ired+red)
      NDVI[NDVI<0.3]<- NA
      # apply mask
      raster <- mask(raster, NDVI)
      
      #### Predict PLS DA
      r_PLS  <- predict(raster, PLS2, type="class")
      dummyList.PLS[[j]] <- r_PLS
      #### Predict RF
      r_RF <- predict(raster, RF2, type="class")
      dummyList.RF[[j]] <- r_RF
      #### Predict SVM
      r_SVM  <- predict(raster, SVM2, type="class")
      dummyList.SVM[[j]] <- r_SVM
      
      ## poliginize to keep the species name label
      poly_PLS <- rasterToPolygons(r_PLS, dissolve = T)
      poly_RF  <- rasterToPolygons(r_RF,  dissolve = T)
      poly_SVM <- rasterToPolygons(r_SVM, dissolve = T)
      
      ## add spp manes to the shapefiles
      # create a function to finde classes that were not predicted
      # lexicographic vector
      not.pred <- function(pred, poly){
        a = c(1, seq(10, length(levels(pred)[[1]][,1]), 1), seq(2,9,1))
        b = setdiff(levels(pred[[i]])[[1]][,1], poly@data$layer)  
        c = which(a == b)
        return(c)
      }
      
      ### PLS-DA
      if ( length(levels(r_PLS)[[1]][,1]) < 10 ){
        a = c(1,2,3,4,5,6,7,8,9)
      } else if  (length(levels(r_PLS)[[1]][,1]) >= 10 & length(levels(r_PLS)[[1]][,1]) < 20){
        a = c(1, seq(10, length(levels(r_PLS)[[1]][,1]), 1), seq(2,9,1)) 
      } else {
        a = c(1, seq(10, 19, 1), 2, seq(20,length(levels(r_PLS)[[1]][,1]),1), seq(3,9,1)) 
      }
      
      if ( !is.null( setdiff( levels(r_PLS)[[1]][,1], poly@data$layer ) ) ){
        poly$sp <- levels(r_PLS)[[1]][,2][poly@data$layer]
      } else { 
        b = as.vector_PLS(levels(r_PLS)[[1]][,1]) 
        b[a]
        b[not.pred(r_PLS, poly)] <- NA
        poly$sp <- na.omit(b)
      }
      
      ### RF
      if ( length(levels(r_RF)[[1]][,1]) < 10 ){
        a = c(1,2,3,4,5,6,7,8,9)
      } else if  (length(levels(r_RF)[[1]][,1]) >= 10 & length(levels(r_RF)[[1]][,1]) < 20){
        a = c(1, seq(10, length(levels(r_RF)[[1]][,1]), 1), seq(2,9,1)) 
      } else {
        a = c(1, seq(10, 19, 1), 2, seq(20,length(levels(r_RF)[[1]][,1]),1), seq(3,9,1)) 
      }
      
      if ( !is.null( setdiff( levels(r_RF)[[1]][,1], poly@data$layer ) ) ){
        poly$sp <- levels(r_RF)[[1]][,2][poly@data$layer]
      } else { 
        b = as.vector_RF(levels(r_RF)[[1]][,1]) 
        b[a]
        b[not.pred(r_RF, poly)] <- NA
        poly$sp <- na.omit(b)
      }
      
      ### SVM-DA
      if ( length(levels(r_SVM)[[1]][,1]) < 10 ){
        a = c(1,2,3,4,5,6,7,8,9)
      } else if  (length(levels(r_SVM)[[1]][,1]) >= 10 & length(levels(r_SVM)[[1]][,1]) < 20){
        a = c(1, seq(10, length(levels(r_SVM)[[1]][,1]), 1), seq(2,9,1)) 
      } else {
        a = c(1, seq(10, 19, 1), 2, seq(20,length(levels(r_SVM)[[1]][,1]),1), seq(3,9,1)) 
      }
      
      if ( !is.null( setdiff( levels(r_SVM)[[1]][,1], poly@data$layer ) ) ){
        poly$sp <- levels(r_SVM)[[1]][,2][poly@data$layer]
      } else { 
        b = as.vector_SVM(levels(r_SVM)[[1]][,1]) 
        b[a]
        b[not.pred(r_SVM, poly)] <- NA
        poly$sp <- na.omit(b)
      }
      
      
      ### export shapefiles
      # create a folder per plot to store results
      plotName = paste( names(rasterPlots)[j], "_PLS_", modelTag, sep=""  )
      dir.create(file.path(outDir, plotName), showWarnings = FALSE)
      outPolydir_PLS = file.path(outDir, plotName)
      
      plotName = paste( names(rasterPlots)[j], "_RF_", modelTag, sep=""  )
      dir.create(file.path(outDir, plotName), showWarnings = FALSE)
      outPolydir_RF = file.path(outDir, plotName)
      
      plotName = paste( names(rasterPlots)[j], "_SVM_", modelTag, sep=""  )
      dir.create(file.path(outDir, plotName), showWarnings = FALSE)
      outPolydir_SVM = file.path(outDir, plotName)
      
      outPoly_PLS = paste( names(rasterPlots)[j], "_PLS_", i, sep="" )
      outPoly_RF  = paste( names(rasterPlots)[j], "_RF_", i, sep="" )
      outPoly_SVM = paste( names(rasterPlots)[j], "_SVM_", i, sep="" )
      
      shpList.PLS[[i]] <- poly_PLS
      shpList.RF[[i]]  <- poly_RF
      shpList.SVM[[i]] <- poly_SVM
      
      writeOGR(poly_PLS, outPolydir, outPoly_PLS, driver="ESRI Shapefile", overwrite_layer = T)
      writeOGR(poly_RF,  outPolydir, outPoly_RF,  driver="ESRI Shapefile", overwrite_layer = T)
      writeOGR(poly_SVM, outPolydir, outPoly_SVM, driver="ESRI Shapefile", overwrite_layer = T)
      
    }
    
    rasterList.SVM[[i]] <- dummyList.SVM
  }
  
  # stop parallel process
  stopCluster(cl) 
  
  ######################
  ### prepare output ###
  ######################
  
  fit <- c (PA.PLS, UA.PLS, OA.PLS, kappa.PLS, PA.RF, UA.RF, OA.RF, kappa.RF, PA.SVM, UA.SVM, OA.SVM, kappa.SVM)
  names (fit) <- c ("PLS PA", "PLS UA", "PLS OA", "PLS Kappa", 
                    "RF PA", "RF UA", "RF OA", "RF Kappa", 
                    "SVM PA", "SVM UA", "SVM OA", "SVM Kappa")
  output <- list (fit, model.PLS, predict.PLS, rasterList.PLS, shpList.PLS,
                  model.RF, predict.RF, rasterList.RF, shpList.RF,
                  model.SVM, predict.SVM, rasterList.SVM, shpList.SVM)
  names (output) <- c ("fits", "PLS Model", "PLS Predictions", "PLS rasterPredictions", "PLS shpPredictions",
                       "RF Model", "RF Predictions", "RF rasterPredictions", "RF shpPredictions",
                       "SVM Model", "SVM Predictions", "SVM rasterPredictions", "SVM shpPredictions")
  class (output) <- "BootsClassification"
  output
  
}


##----------------------------------------------------------------------------##
##                                                                            ##
## plot.ensemble: visualization of ensemble objects                           ##
##                                                                            ##
## Arguments:                                                                 ##
## - spec     spectral information. Used to create the quantiles of spectra   ##
## - en       classificationEnsemble object                                   ##
##                                                                            ##
##----------------------------------------------------------------------------##

plot.classificationEnsemble <- function (spec, en, xlab_tag,label=TRUE, ...) {
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
  plot(wl, quant[1,], type="l", ylim = c(0,max(spec)), axes=F, ylab = NA, xlab=NA, las=1)
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
  plot(wl,  cf[1,], type = "l", col=2, ylab="Weighted coefficients", las=1, xlab=xlab_tag,
       ylim=c(min(cf), max(cf)), lty=1, lwd=2)
  lines(wl,  cf[2,], type = "l", col=3, lty=2, lwd=2)
  lines(wl,  cf[3,], type = "l", col=4, lty=3, lwd=2) 
  par(xpd = F) 
  abline(h=0, lty=2)
  # add spectral bands selected by the ensemble method
  image (wl, seq(min (cf) * 1.1, max (cf) * 1.1, length.out=100), z2, col=7, 
         xlab="", ylab="", add=T)
  # Labels
  if(label==TRUE){
    labels <- c (paste (c ("PLS ", "RF ", "SVM "), "OA", "=", fit, sep=""), NA, 
                 "Ensemble selection")
    legend ("topright", bty="n", col=c (2, 3, 4, NA, NA, rep (1, 3)), 
            pt.bg=c(rep (NA, 4), 7), lwd=c(rep (2, 3), rep (NA, 2)),
            pch=c (rep (NA, 4), 22), cex=0.7, pt.cex=1, legend=labels, lty=c(1,2,3,rep (NA, 2)))
   }
 }


##----------------------------------------------------------------------------##
##                                                                            ##
## significanceTest_LeafLevel: Apply bootstrap significance test to leaf data ##
##                                                                            ##
## Arguments:                                                                 ##
## -  data:    input dataset                                                  ##
## - fitASD:   classificationEnsemble ASD object                              ##
## - fitAISA:   classificationEnsemble AISA object                            ##
##                                                                            ##
## Function based on:                                                         ##
## Lopatin, J., Dolos, K., Hernández, H. J., Galleguillos, M., & Fassnacht,   ##
## F. E. (2016). Comparing Generalized Linear Models and random forest to     ##
## model vascular plant species richness using LiDAR data in a natural forest ##
## in central Chile. Remote Sensing of Environment, 173, 200-210.             ##
## http://doi.org/10.1016/j.rse.2015.11.029                                   ##
##                                                                            ##
##----------------------------------------------------------------------------##

significanceTest_LeafLevel <- function(data, fitASD, fitAISA, B=500){
  
  library(hyperSpec)
  library(boot)
  library(e1071)
  library(doParallel)
  library(caret)
      
  ####################
  ### prepare data ###
  ####################
  classes <- data$Species 
  # spectral bands
  spectra <- data[, 2:(length(data)-1)]

  new("hyperSpec")
  hyperASD <- new("hyperSpec", spc=spectra, data=data, wavelength = seq(350, 2500, 1),
                 label=list(spc="Reflectance", .wavelength =  expression(lambda(nm))))
  
  hyperASD <- hyperASD[,, c(390~max)]
  hyperASD <- spc.bin (hyperASD, 10)
  
  hyperAISA <- hyperASD[,, c(390~990)]
  
  #########################
  ### prepare bootstrap ###
  #########################
  set.seed(123)
  # set the bootstrap parameters
  N = length(data[,1]) # N° of observations
  
  # list to store
  diff_OA <- list()
  diff_kappa <- list()
  
  OA.ASD <- list()
  kappa.ASD <- list()
  PA.ASD <- list()
  UA.ASD <- list()
  
  OA.AISA <- list()
  kappa.AISA <- list()
  PA.AISA <- list()
  UA.AISA <- list()
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  # start loop
  for(i in 1:B){
   # create random numbers with replacement to select samples from each group
   idx = sample(1:N, N, replace=TRUE)
    
   ###########
   ### ASD ###
   ###########
   m1 <-  svm(hyperASD$spc[idx,], hyperASD@data$Species[idx], gamma = fitASD$SVM$finalModel$gamma, cost = fitASD$SVM$finalModel$cost, probability = TRUE)
   
   m1.pred <- predict(m1, hyperASD$spc[-idx,])
   
   # confusion matix
   preMatrix <- function(pred, test){ # functionn to prevent caret error for different length
     u = union(pred, test)
     t = table(factor(pred, u), factor(test, u))
     return(t)
   }  
   m1.conf <- confusionMatrix(preMatrix(m1.pred, hyperASD@data$Species[-idx]))
   
   # get accuracies
   m1.OA    <- m1.conf$overall["Accuracy"]
   m1.kappa <- m1.conf$overall["Kappa"]
   m1.PA    <- m1.conf$byClass[,3] 
   m1.UA    <- m1.conf$byClass[,4]
   
   OA.ASD[[i]] <- m1.OA
   kappa.ASD[[i]] <- m1.kappa
   PA.ASD[[i]] <- m1.PA
   UA.ASD[[i]] <- m1.UA
   
   ############
   ### AISA ###
   ############
   m2 <-  svm(hyperAISA$spc[idx,], hyperAISA@data$Species[idx], gamma = fitAISA$SVM$finalModel$gamma, cost = fitAISA$SVM$finalModel$cost, probability = TRUE)
   
   m2.pred <- predict(m2, hyperAISA$spc[-idx,])
   
   # confusion matix
   m2.conf <- confusionMatrix(preMatrix(m2.pred, hyperASD@data$Species[-idx]))
   # get accuracies
   m2.OA    <- m2.conf$overall["Accuracy"]
   m2.kappa <- m2.conf$overall["Kappa"]
   m2.PA    <- m2.conf$byClass[,3] 
   m2.UA    <- m2.conf$byClass[,4]
   
   OA.AISA[[i]] <- m2.OA
   kappa.AISA[[i]] <- m2.kappa
   PA.AISA[[i]] <- m2.PA
   UA.AISA[[i]] <- m2.UA
   
   ### compute the differences between ASD and AISA band settings
   ### OA of ASD should be larger. So, if OA(m1) - OA(m2) is positive, ASD is significantly better.
   diff_OA[[i]] <- m1.OA - m2.OA
   
   ### Kappa of ASD should be larger. So, if kappa(m1) - kappa(m2) is positive, ASD is significantly better.
   diff_kappa[[i]] <- m1.kappa - m2.kappa
    
  }
  
  # stop parallel process
  stopCluster(cl)
  
  # prepare output
  fit <- list(OA.ASD, OA.AISA, kappa.ASD, kappa.AISA, PA.ASD, PA.AISA, UA.ASD, UA.AISA)
  names(fit) <- c("OA.ASD", "OA.AISA", "kappa.ASD", "kappa.AISA", "PA.ASD", "PA.AISA", "UA.ASD", "UA.AISA")
  boot_test <- list(diff_OA, diff_kappa)
  names(boot_test) <- c("OA", "Kappa")
  output <- list(boot_test, fit)
  names(output) <- c("boot_test", "fit")
  class(output) <- "boot_test"
  output
 
  }


##----------------------------------------------------------------------------##
##                                                                            ##
## GLCM: Apply Gray Level Covariance Matrix to an image                       ##
##                                                                            ##
## Arguments:                                                                 ##
## -  img       input image                                                   ##
##                                                                            ##
##----------------------------------------------------------------------------##

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

##----------------------------------------------------------------------------##
##                                                                            ##
## Small functions to help in repetitive tasks                                ##
##                                                                            ##
##                                                                            ##
##----------------------------------------------------------------------------##

## List the names of the rasters in a folder
rasterListNames <- function(fileExtantion, folder){
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # delete the ".dat" from the name
  rast_list = gsub('.{4}$', '', rast_list)
  return(rast_list)
}

## List and load the rasters contained in a folder
rasterList <- function(fileExtantion, folder, dir=NULL, select=NULL){
  # if dir = NULL, set it to "home" by default
  if (is.null(dir)){
    dir = home
  }
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # select only rasters with a especific pattern
  if (!is.null(select)){
    rast_list <- rast_list[ grep(select, rast_list) ]
  }
  # raster names
  rasterNames = gsub('.{4}$', '', rast_list)
  # import rasters
  setwd(file.path(dir, folder))
  rasterlist <- list()
  for(i in 1:length(rast_list)){
    rast <- stack(rast_list[i])
    rasterlist[[i]] <- rast
  }
  names(rasterlist) <- rasterNames
  setwd(dir)
  return(rasterlist)
}
