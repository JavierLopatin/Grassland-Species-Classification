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
## tunningModels and  Multi-method ensemble selection of spectral bands       ##
##                                                                            ##
## This function performs a tunning procidure on the models and               ##
## a band selection based on a multi-method ensemble                          ##
## assessment of the variable importance and classification coefficients of   ##
## three different model types: Partial Least Squares Discriminant Analysis,  ##
## Random Forest and Support Vector Machine classifications                   ## 
##                                                                            ## 
## Arguments:                                                                 ##
## - x      Numeric matrix containing the spectra (samples as rows)           ##
## - y      Numeric vector containing the response variable                   ##
## - wl     Numeric vector containing the wavelength information of the bands ##
##                                                                            ##
## variable importance based on the paper:                                    ##
## Feilhauer, H., Asner, G.P., & Martin, R.E. (2015). Multi-method ensemble   ##
## selection of spectral bands related to leaf biochemistry. Remote Sensing   ## 
## of Environment, 164, 57-65. http://doi.org/10.1016/j.rse.2015.03.033       ##
##                                                                            ##
##----------------------------------------------------------------------------##

tunningModels <- function(data, Site, wl=NA){
  
  ## load required libraries
  library (caret)
  library(e1071)
  library(doParallel)
  
  # set data
  x = grep(Site, data$Site) 
  x1 <- data[x, ]
  x1$Species <- factor(x1$Species)
  
  spec <- x1[, 3:length(x1)]
  data2 <- data.frame(classes = x1$Species, x1[, 3:length(x1)])
  data2 <- na.omit(data2)
  
  # Set the random number seed so we can reproduce the results
  set.seed(123)
  # Split data in training and test
  forTraining <- createDataPartition(data2$classes, p = 0.6, list=F)
  train <- data2 [ forTraining,]
  test<- data2 [-forTraining,]
  
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

ApplyBootsClassification <- function(data, en, Site, rasterPlots, boots=100, outDir, modelTag){  
  
  library(raster)
  library(rgdal)
  library(caret)
  library(e1071)
  library(randomForest)
  library(pls)
  library(doParallel)
  
  # extract the data from the classification Ensamble function
  x = grep(Site, data$Site) 
  x1 <- data[x, ]
  x1$Species <- factor(x1$Species)
  
  data2 <- data.frame(classes = x1$Species, x1[, 3:length(x1)])
  
  wl <- en[[1]][1,]
  
  ncomp = en$PLS$finalModel$ncomp
  probMethod = en$PLS$finalModel$probMethod
  
  bestNtree = en$RF$finalModel$ntree
  bestMtry = en$RF$finalModel$mtry
  
  bestCost <- en$SVM$finalModel$cost
  bestGamma <- en$SVM$finalModel$gamma
  

  ## apply funtion to predict species cover with SVM
  
  # list of accuracies
  PA.PLS    <- list()
  UA.PLS    <- list()
  OA.PLS    <- list()
  kappa.PLS <- list()
  predict.PLS <- list()
 
  PA.RF    <- list()
  UA.RF    <- list()
  OA.RF    <- list()
  kappa.RF <- list()
  predict.RF <- list()

  PA.SVM    <- list()
  UA.SVM    <- list()
  OA.SVM    <- list()
  kappa.SVM <- list()
  predict.SVM <- list()

  OBS <- list()
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  for (i in 1:boots){
    
    N = length(data2[,1])
    
    # create random numbers with replacement to select samples from each group
    idx = sample(1:N, N, replace=TRUE)
    
    train <- data2[idx,]
    val <- data2[-idx,]
    
    # store and select the observations
    obs <- val$classes
    OBS[[i]]<-obs
    
    #################
    ### Apply PLS ###
    #################
    
    PLS  <- plsda(x =  train[,2:length(train)], y = as.factor( train$classes ), ncomp = ncomp, 
                  probMethod = probMethod)
    
    # predict
    pred_pls    <- predict(PLS, val[,2:length(val)])
    predict.PLS[[i]] <- pred_pls
    
    # confusion matix
    conf   <- confusionMatrix(pred_pls, val$classes)
    
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

    # predict
    pred_rf    <- predict(RF, val[,2:length(val)])
    predict.RF[[i]] <- pred_rf
    
    # confusion matix
    conf   <- confusionMatrix(pred_rf, val$classes)
    
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

    # predict
    pred_svm    <- predict(SVM, val[,2:length(val)])
    predict.SVM[[i]] <- pred_svm
    
    # confusion matix
    conf   <- confusionMatrix(pred_svm, val$classes)
    
    # get accuracies
    PA.SVM[[i]]       <- conf$byClass[,3] 
    UA.SVM[[i]]       <- conf$byClass[,4]
    OA.SVM[[i]]       <- conf$overall["Accuracy"]
    kappa.SVM[[i]]    <- conf$overall["Kappa"]
    
    
    ##############################
    ### Apply models to raster ### 
    ##############################
    
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
      
      ### Predict PLS DA
      r_PLS  <- predict(raster, PLS, type="class")
      ### Predict RF
      r_RF <- predict(raster, RF, type="class")
      #### Predict SVM
      r_SVM  <- predict(raster, SVM, type="class")
      
      ### export rasters
      # create a folder per plot to store results
      plotName = paste( names(rasterPlots)[j], "_PLS_", modelTag, sep=""  )
      dir.create(file.path(outDir, plotName), showWarnings = FALSE)
      outdir_PLS = file.path(outDir, plotName)
      
      plotName = paste( names(rasterPlots)[j], "_RF_", modelTag, sep=""  )
      dir.create(file.path(outDir, plotName), showWarnings = FALSE)
      outdir_RF = file.path(outDir, plotName)
      
      plotName = paste( names(rasterPlots)[j], "_SVM_", modelTag, sep=""  )
      dir.create(file.path(outDir, plotName), showWarnings = FALSE)
      outdir_SVM = file.path(outDir, plotName)
      
      out_PLS = paste( names(rasterPlots)[j], "_PLS_", i, ".tif", sep="" )
      out_RF  = paste( names(rasterPlots)[j], "_RF_", i, ".tif", sep="" )
      out_SVM = paste( names(rasterPlots)[j], "_SVM_", i, ".tif", sep="" )
      
      out_PLS = file.path(outdir_PLS, out_PLS)
      out_RF  = file.path(outdir_RF, out_RF)
      out_SVM = file.path(outdir_PLS, out_SVM)
      
      # Export rasters
      writeRaster(r_PLS, filename=out_PLS, format="GTiff", overwrite = T)
      writeRaster(r_RF,  filename=out_RF,  format="GTiff", overwrite = T)
      writeRaster(r_SVM, filename=out_SVM, format="GTiff", overwrite = T)
      
    }
  }
  
  # stop parallel process
  stopCluster(cl) 
  
  ######################
  ### prepare output ###
  ######################
  

  fits <- data.frame( unlist(OA.PLS), unlist(OA.RF), unlist(OA.SVM),
                     unlist(kappa.PLS), unlist(kappa.RF), unlist(kappa.SVM))
  names(fits) <- c("OA_PLS", "OA_RF", "OA_SVM", "Kappa_PLS", "Kappa_RF", "Kapa_SVM")
  
  fits_2 <- data.frame( PLS$obsLevels , unlist(PA.PLS), unlist(PA.RF), unlist(PA.SVM), 
                        unlist(OA.PLS), unlist(OA.RF), unlist(OA.SVM))
  names(fits_2) <- c("Species", "PA_PLS", "PA_RF", "PA_SVM", "OA_PLS", "OA_RF", "OA_SVM")
  
  predict_all <- data.frame( unlist(OBS), unlist(predict.PLS), unlist(predict.RF), unlist(predict.SVM) )
  names(predict_all) <- c("Observed", "PLS", "RF", "SVM")
  
  ### Export 
  
  write.table(fits, file = file.path(outDir, paste("Fits_Site_", Site, "_", modelTag, ".txt", sep="")), row.names = F, col.names = T)
  write.table(fits_2, file = file.path(outDir, paste("FitsPA_OA_Site_", Site, "_", modelTag, ".txt", sep="")), row.names = F, col.names = T)
  write.table(predict_all, file = file.path(outDir, paste("Predicts_Site_", Site, "_", modelTag, ".txt", sep="")), row.names = F, col.names = T)
  
}


##----------------------------------------------------------------------------##
##                                                                            ##
## obstainCovers: Obtain the covers prediction valiuos per plot               ##
##                                                                            ##
## Arguments:                                                                 ##
## - rasterDir:   folder of the predicted plots of ApplyBootsClassification   ##
## - shpDir:      folder of the shapefiles used to mask the image             ##
## - shpMaskName: Name of the plot to use as a mask                           ## 
##                                                                            ##
##----------------------------------------------------------------------------##

obstainCovers <- function(rasterDir, shpDir, shpMaskName){ 
  
  library(raster)
  library(rgdal)
  
  
  #######################################
  ### obtain observed covers per plot ###
  #######################################
  
  x = grep(18, species$Plot) 
  y <- species[x,]
  y$Species = factor(y$Species) 
  
  store_plot <-  matrix(nrow = 1, ncol = length( levels(y$Species) ))
  colnames(store_plot) <- levels(y$Species)
  
  for (i in 1:length( levels(y$Species))){
    sp = levels(y$Species)[i]
    z =  grep(sp, y$Species)
    z <- y[z,]
    sumCov = sum(z$Cover)
    cov = (sumCov*100)/1600
    store_plot[[1,i]] <- cov
  }
  
  #########################################
  ### estimate predicted cover per plot ###
  #########################################
  setwd(rasterDir)
  
  rstLisr <- rasterList(fileExtantion = ".tif", folder = ".", dir = rasterDir, select=NULL)
  shp <- readOGR(dsn = shpDir, layer = shpMaskName)
  
  levs = levels(rstLisr[[1]])
  levels = levs[[1]][,2]
  levelsNumber = length(levels)
  store_areas <- matrix(nrow = 100, ncol = levelsNumber)
  colnames(store_areas) <- levels
  for (i in 1:100){
    clip <- crop(rstLisr[[i]], shp)
    
    # count number of pixels
    area = ncell(clip)
    
    # cover per class
    for (j in 1:levelsNumber){ 
      count_class = freq(clip, value=j)
      percent = (count_class*100)/area
      # add value to store_areas matrix
      store_areas[[i,j]] <- percent
    }
  }
  
  ################################
  ### Export results to tables ###
  ################################
  
  setwd(home)
  
  dir.create("Covers_results", showWarnings = FALSE)
  ObsCov = paste( "ObsCov_", shpMaskName, ".txt", sep=""  )
  ObsCov = file.path( home, "Covers_results", ObsCov )
  PredCov = paste( "PredCov_", shpMaskName, ".txt", sep=""  )
  PredCov = file.path( home, "Covers_results", PredCov )
  
  write.table(store_plot, file = ObsCov, sep = " ", row.names = F, col.names = T )
  write.table(store_areas, file = PredCov, sep = " ", row.names = F, col.names = T )
  
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
## Lopatin, J., Dolos, K., Hernández, H. J., Galleguillos, M., & Fassnacht,  ##
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
  x = grep(".tif.aux.xml", rast_list)
  if ( length(x) > 0 ){ rast_list <- rast_list[-x] }
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
  x = grep(".tif.aux.xml", rast_list) 
  if ( length(x) > 0 ){ rast_list <- rast_list[-x] }
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
