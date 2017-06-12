################################################################################
## R-Script - 0_Functions.R                                                   ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Mapping plant species in mixed grassland communities using     ##
##             close range imaging spectroscopy                               ##
##                                                                            ##
## description: This R-code provide most of the functions applied in the paper## 
##                                                                            ##
################################################################################

##----------------------------------------------------------------------------##
##                                                                            ##
## Apply the functions tuningModels, BootsClassification and obstainCovers    ##
## in a row to obtain the classification results per dataset                  ## 
##                                                                            ##
## Arguments:                                                                 ##
## - valData     classes to classify. I this case the species list per subplot##
## - potVa-      reflectance obtained with the pot validation method          ##
## - rf          reflectance obtained with the rip-it-off validation method   ##
## - raster_List list of the plots rasters                                    ##
## - wl          wavelength vector                                            ##
## - modelTag    tag to assign to the exported results                        ##
## - boots       number of bootstraps to perform                              ##
##                                                                            ##
##----------------------------------------------------------------------------##

ApplyModels <- function(valData, potVal, rf, raster_List, wl, modelTag, boots){
  
  for (i in 1:length(raster_List)){ 
    
    # obtain the validation data per plot
    raster = raster_List[[i]]
    names(raster) <- paste0( rep("B", nlayers(raster)), seq(1, nlayers(raster), 1) )
    plot = unique(na.omit(as.numeric(unlist(strsplit( names( raster_List )[[i]], "[^0-9]+")))))
    plot_name = paste0("plot_", plot)  
    
    x = grep( plot, valData$Plot )  
    data = valData[x, ] 
    data$Species <- factor( data$Species ) # reset species Levels
    # get species to classify in the plot
    classes = unique( data$Species )
    
    # obtain subset of data of potVal and rf that include these species
    x = grep( paste(classes, collapse = "|") , potVal$Species )
    data_potVal = potVal[x, ]
    data_potVal$Species <- factor(data_potVal$Species) 
    
    x = grep( paste(classes, collapse = "|") , rf$Species )
    data_rf = rf[x, ]
    data_rf$Species <- factor(data_rf$Species) 
    
    ####################################
    ### Apply tuningModels function ###
    ####################################
    
    print(paste("### Tuning", plot_name, modelTag, "###"))
    print("")
    
    print("Tuning pot Validation model")
    print("")
    
    fit_potVal <- tuningModels(classes = data_potVal$Species, 
                                spectra = data_potVal[, 3:length( data_potVal )],
                                wl = wl)
    
    print("Tuning rip-it-off model")
    print("")
    
    fit_rf     <- tuningModels(classes = data_rf$Species, 
                                spectra = data_rf[, 3:length( data_rf )],
                                wl = wl)
    # save tuning models
    dir.create(file.path(home, "tuningOutputs"), showWarnings = FALSE)
    
    save(fit_potVal, file=paste0(home,  "/tuningOutputs/", "potVal_", 
                                 plot_name, "_", modelTag, ".RData"))
    save(fit_rf,     file=paste0(home,  "/tuningOutputs/", "rf_",
                                 plot_name,  "_", modelTag, ".RData"))
    
    print("Done!")
    print("")
    
    #################################
    ### Apply BootsClassification ###
    #################################
    
    dir.create(file.path(home, "BootsClass_out"), showWarnings = FALSE)
    
    print("### Starting bootstrap predictions ###")
    print("")
    
    BootsClassification(classes = data_potVal$Species, 
                        spectra = data_potVal[, 3:length( data_potVal )],
                        en = fit_potVal, 
                        raster = raster, 
                        boots = boots, 
                        outDir = file.path(home, "BootsClass_out"), 
                        modelTag = paste0("potVal_", modelTag),
                        plotName = plot_name)
    
    BootsClassification(classes = data_rf$Species, 
                        spectra = data_rf[, 3:length( data_potVal )],
                        en = fit_rf, 
                        raster = raster, 
                        boots = boots, 
                        outDir = file.path(home, "BootsClass_out"), 
                        modelTag = paste0("rf_", modelTag),
                        plotName = plot_name)
    
    print("Done!")
    print("")
    
    #####################################
    ### Apply obstainCovers function ####
    #####################################
    
    print("### Starting Cover estimation ###")
    print("")
    
    # for PLS-DA potVal
    obstainCovers(ObservedSpecies = valData, 
                  rasterDir =  paste0( home, "/BootsClass_out/", plot_name, "_PLS_",  
                                       paste0("potVal_", modelTag) ), 
                  subplotDir = subplotDir, 
                  shpMaskName = plot_name, 
                  plotNumber = plot, 
                  Iter = boots,
                  algorithm = "PLS_potVal")
    
    # for PLS-DA rf
    obstainCovers(ObservedSpecies = valData, 
                  rasterDir =  paste0( home, "/BootsClass_out/", plot_name, "_PLS_", 
                                       paste0("rf_", modelTag) ), 
                  subplotDir = subplotDir, 
                  shpMaskName = plot_name, 
                  plotNumber = plot, 
                  Iter = boots,
                  algorithm = "PLS_rf")
    
    print("PLS-DA done!")
    
    # for RF potVal
    obstainCovers(ObservedSpecies = valData, 
                  rasterDir =  paste0( home, "/BootsClass_out/", plot_name, "_RF_",  
                                       paste0("potVal_", modelTag) ), 
                  subplotDir = subplotDir, 
                  shpMaskName = plot_name, 
                  plotNumber = plot, 
                  Iter = boots,
                  algorithm = "RF_potVal")
    
    # for RF rf
    obstainCovers(ObservedSpecies = valData, 
                  rasterDir =  paste0( home, "/BootsClass_out/", plot_name, "_RF_",  
                                       paste0("rf_", modelTag) ), 
                  subplotDir = subplotDir, 
                  shpMaskName = plot_name, 
                  plotNumber = plot, 
                  Iter = boots,
                  algorithm = "RF_rf")
    
    print("RF done!")
    
    # for SVM potVal
    obstainCovers(ObservedSpecies = valData, 
                  rasterDir =  paste0( home, "/BootsClass_out/", plot_name, "_SVM_",  
                                       paste0("potVal_", modelTag) ), 
                  subplotDir = subplotDir, 
                  shpMaskName = plot_name, 
                  plotNumber = plot, 
                  Iter = boots,
                  algorithm = "SVM_potVal")
    
    # for SVM rf potVal
    obstainCovers(ObservedSpecies = valData, 
                  rasterDir =  paste0( home, "/BootsClass_out/", plot_name, "_SVM_",  
                                       paste0("rf_", modelTag) ), 
                  subplotDir = subplotDir, 
                  shpMaskName = plot_name, 
                  plotNumber = plot, 
                  Iter = boots,
                  algorithm = "SVM_rf")
    
    print("SVM done!")
    print("")
    
  }
}


##----------------------------------------------------------------------------##
##                                                                            ##
## tuningModels function                                                      ##
##                                                                            ##
## This function performs a tuning procidure on the models and                ##
## a band selection based on a multi-method ensemble                          ##
## assessment of the variable importance and classification coefficients of   ##
## three different model types: Partial Least Squares Discriminant Analysis,  ##
## Random Forest and Support Vector Machine classifications (not presented in ##
## the paper)                                                                 ## 
##                                                                            ## 
## Arguments:                                                                 ##
## - x      Numeric matrix containing the spectra (samples as rows)           ##
## - y      Numeric vector containing the response variable                   ##
## - wl     Numeric vector containing the wavelength information of the bands ##
##                                                                            ##
##----------------------------------------------------------------------------##

tuningModels <- function(classes, spectra, wl=NA){
  
  ## load required libraries
  library(caret)
  library(e1071)
  library(doParallel)
  
  # set data
  data2 <- data.frame(classes = classes, spectra)
  data2 <- na.omit(data2)
  
  # Set the random number seed so we can reproduce the results
  set.seed(123)
  # Split data in training and test
  forTraining <- createDataPartition(data2$classes, p = 0.6, list=F)
  train <- data2 [ forTraining,]
  test<- data2 [-forTraining,]
  
  # Each model used 5 repeated 10-fold cross-validation. Use AUC to pick the best model
  controlObject <- trainControl(method = "cv", number = 5, classProbs=TRUE, 
                                allowParallel = TRUE, seeds = set.seed(123))
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  #############################
  ### PLS-DA classification ###
  #############################
  
  print("Tuning PLS-DA Model...")

  # apply classification
  set.seed(123)
  plsClas <- train(x=train[, 2:length(train)], y=make.names( train$classes ), method = "pls", 
                   tuneLength=20, preProc = c("center", "scale"), trControl = controlObject)
  
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
  
  print("Done!")
  
  #########################
  ### RF classification ###
  #########################
  
  print("Tuning RF model...")

  set.seed(123)
  rfClas <- train(x=train[, 2:length(train)], y=make.names( train$classes ), method = "rf", 
                  tuneLength=15, trControl = controlObject)
  
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
  
  print("Done!")
  
  ##########################
  ### SVM classification ###
  ##########################
  
  print("Tuning SVM Model...")

  set.seed(123)
  svmClas <- train(x=train[, 2:length(train)], y=make.names( train$classes ), method = "svmLinear2", 
                   tuneLength=10, preProc = c("center", "scale"), trControl = controlObject)
  
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
  svrcf <- numeric (ncol (spectra))
  for(i in 1:ncol(spectra))
    svrcf[i] <- svr.alpha %*% spectra[svr.index, i]
  svrcf <- svrcf / sd (svrcf) ## scale pseudo-coefficients
  
  print("Done!")
  
  # stop parallel process
  stopCluster(cl) 
  
  #####################################################################    
  ### get ensemble from all models and identify important variables ###
  #####################################################################
  
  ## get ensemble from all models and identify important variables
  ensemblecf <- abs(plscf) * OA.pls + abs(rfcf) * OA.rf + abs(svrcf) * OA.svm 
  th <- mean(ensemblecf) + sd(ensemblecf) ## calculate threshold
  selbands <- ensemblecf > th ## apply threshold

  ######################
  ### prepare output ###
  ######################
  
  cf <- rbind (wl, plscf, rfcf, svrcf, ensemblecf, selbands)
  colnames(cf) <- colnames(spectra)
  
  fit <- c (OA.pls, OA.rf, OA.svm)
  names (fit) <- c ("PLS-DA OA", "RF OA", "SVR OA")
  output <- list (cf, fit, th, plsClas, rfClas, svmClas, conf.pls, conf.rf, conf.svm)
  names (output) <- c ("selection", "fits", "threshold", "PLS", "RF", "SVM", "confusionMatrix.PLS", "confusionMatrix.RF", "confusionMatrix.SVM")
  class (output) <- "ensemble"
  output
  
}


##----------------------------------------------------------------------------##
##                                                                            ##
## Apply the best model from classificationEnsemble to the plots using        ##
## a bootstrapping procidure.                                                 ## 
##                                                                            ##
## Arguments:                                                                 ##
## - spectra  spectral information. Used to create the quantiles of spectra   ##
## - en       classificationEnsemble object                                   ##
##                                                                            ##
##----------------------------------------------------------------------------##

BootsClassification <- function(classes, spectra, en, raster, boots, 
                                outDir, modelTag, plotName){  
  
  library(raster)
  library(rgdal)
  library(caret)
  library(e1071)
  library(randomForest)
  library(pls)
  library(doParallel)
  
  # extract the data from the classification Ensamble function
  data2 <- data.frame(classes = classes, spectra)
  data2 <- na.omit(data2)
  data2$classes <- factor(data2$classes)
  
  ncomp      = en$PLS$finalModel$ncomp
  probMethod = en$PLS$finalModel$probMethod
  
  bestNtree = en$RF$finalModel$ntree
  bestMtry  = en$RF$finalModel$mtry
  
  bestCost  = en$SVM$finalModel$cost
  bestGamma = en$SVM$finalModel$gamma
  
  
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
 
  # progress bar
  print(paste0(plotName, " ", modelTag))
  pb <- txtProgressBar(min = 0, max = boots, style = 3)

  for (i in 1:boots){
    
    # progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    # stratify samplig. All species get selected at least once
    samp <- stratifySampling(data2, classes)
    
    train <- samp$train
    val <- samp$validation
    
    # store and select the observations
    obs <- val$classes
    OBS[[i]]<-obs
    
    #################
    ### Apply PLS ###
    #################
    
    PLS  <- plsda(x =  train[, 2:length(train)], y = train$classes, ncomp = ncomp, 
                  probMethod = probMethod)
    
    # predict
    pred_pls    <- predict(PLS, val[, 2:length(val)])
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
    
    RF  <- randomForest( y = train$classes, x = train[, 2:length(train)],
                         ntree= bestNtree, mtry = bestMtry)

    # predict
    pred_rf    <- predict(RF, val[, 2:length(val)])
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
    
    SVM  <- svm(train[, 2:length(train)], train$classes, kernel = "linear",
                gamma = bestGamma, cost = bestCost, probability = TRUE)

    # predict
    pred_svm    <- predict(SVM, val[, 2:length(val)])
    predict.SVM[[i]] <- pred_svm
    
    # confusion matix
    conf <- confusionMatrix(pred_svm, val$classes)
    
    # get accuracies
    PA.SVM[[i]]       <- conf$byClass[,3] 
    UA.SVM[[i]]       <- conf$byClass[,4]
    OA.SVM[[i]]       <- conf$overall["Accuracy"]
    kappa.SVM[[i]]    <- conf$overall["Kappa"]
    
    
    ##############################
    ### Apply models to raster ### 
    ##############################
    
    # NDVI mask to spectral images
    if ( nlayers(raster)==61 ){ 
      names(raster) <- paste( rep("B", 61), seq(1,61,1),  sep="" )
      # mask out raster zones with NDVI below 0.3
      red <- raster[[31]]
      Ired <- raster[[43]]
      NDVI <- (Ired-red)/(Ired+red)
      NDVI[NDVI < 0.3]<- NA
      # apply mask
      raster <- mask(raster, NDVI)
    } else { # for MNF components
      names(raster) <- paste( rep("B", 10), seq(1,10,1),  sep="" )
    }
    
    ### Predict PLS DA
    r_PLS  <- predict(raster, PLS, type="class")
    ### Predict RF
    r_RF <- predict(raster, RF, type="class")
    #### Predict SVM
    r_SVM  <- predict(raster, SVM, type="class")
      
    ### export rasters
    # create a folder per plot to store results
    plotName_PLS = paste( plotName, "_PLS_", modelTag, sep=""  )
    dir.create(file.path(outDir, plotName_PLS), showWarnings = FALSE)
    outdir_PLS = file.path(outDir, plotName_PLS)
      
    plotName_RF = paste(plotName, "_RF_", modelTag, sep=""  )
    dir.create(file.path(outDir, plotName_RF), showWarnings = FALSE)
    outdir_RF = file.path(outDir, plotName_RF)
      
    plotName_SVM = paste( plotName, "_SVM_", modelTag, sep=""  )
    dir.create(file.path(outDir, plotName_SVM), showWarnings = FALSE)
    outdir_SVM = file.path(outDir, plotName_SVM)
      
    out_PLS = paste( plotName, "_PLS_", i, ".tif", sep="" )
    out_RF  = paste( plotName, "_RF_", i, ".tif", sep="" )
    out_SVM = paste( plotName, "_SVM_", i, ".tif", sep="" )
    
    out_PLS = file.path(outdir_PLS, out_PLS)
    out_RF  = file.path(outdir_RF, out_RF)
    out_SVM = file.path(outdir_SVM, out_SVM)
      
    # Export rasters
    writeRaster(r_PLS, filename=out_PLS, format="GTiff", overwrite = T)
    writeRaster(r_RF,  filename=out_RF,  format="GTiff", overwrite = T)
    writeRaster(r_SVM, filename=out_SVM, format="GTiff", overwrite = T)
      
  }

  # close progress bar
  close(pb)
  
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
  
  write.table(fits, file = file.path(outDir, paste("Fits_", plotName, "_", modelTag, ".txt", sep="")),
              row.names = F, col.names = T)
  write.table(fits_2, file = file.path(outDir, paste("FitsPA_OA_", plotName, "_", modelTag, ".txt", sep="")), 
              row.names = F, col.names = T)
  write.table(predict_all, file = file.path(outDir, paste("Predicts_", plotName, "_", modelTag, ".txt", sep="")), 
              row.names = F, col.names = T)
  
}


##----------------------------------------------------------------------------##
##                                                                            ##
## obstainCovers: Obtain the covers prediction values per plot                ##
##                                                                            ##
## Arguments:                                                                 ##
## - rasterDir    folder of the predicted plots of ApplyBootsClassification   ##
## - shpDir       folder of the shapefiles used to mask the image             ##
## - maskName     Name of the plot to use as a mask                           ## 
##                                                                            ##
##----------------------------------------------------------------------------##

obstainCovers <- function(ObservedSpecies, rasterDir, subplotDir, maskName,
                          plotNumber, Iter, algorithm){ 
  
  library(raster)
  library(rgdal)
  
  #######################################
  ### obtain observed covers per plot ###
  #######################################
  
  x = grep(plotNumber, ObservedSpecies$Plot) 
  obs_plot <- ObservedSpecies[x,]
  obs_plot$Species = factor(obs_plot$Species) 
  
  store_plot <-  matrix(nrow = length( levels(obs_plot$Species)), ncol = 2)
  store_plot[,1] <- levels(obs_plot$Species)
  colnames(store_plot) <- c("Species", "Covers")
  
  for (i in 1:length( levels(obs_plot$Species))){
    sp = levels(obs_plot$Species)[i]
    z =  grep(sp, obs_plot$Species)
    z <- obs_plot[z,]
    sumCov = sum(z$Cover)
    cov = (sumCov*100)/1600 # 1600: ~ N° of pixels per subplot
    store_plot[[i,2]] <- cov
  }
  
  rownames(store_plot) <- store_plot[,1]
  
  ############################################
  ### estimate predicted cover per subplot ###
  ############################################
  
  # make a rasterlist from the rasterDir folder
  rstLisr <- rasterList(fileExtantion = ".tif", folder = ".", dir = rasterDir, select=NULL)
  levs <- levels(rstLisr[[1]])
  levelsNumber <- length( levs[[1]][,1] ) 
  levels_plot <- levs[[1]][2:length(levs[[1]][,1]), 2]
  
   
  # load subplots from to copy the extation
  subplots <- rasterList(fileExtantion = ".tif", folder = "subplots", dir = home)
  x = grep(plotNumber, names( subplots ))  
  subplots <- subplots[x]
  
  # subplot names
  subplotNames <-  c("A1","A2","A3","A4","B1","B2","B3","B4","C1", "C2","C3","C4","D1","D2","D3","D4")
  
  ### start loop to obtain covers ###
  areasList <- list()
  
  # loop through the prediction maps 
  for (i in 1:length( rstLisr )){
    raster = rstLisr[[i]]
    #store_areas <- as.data.frame(levs[[1]][2:length(levs[[1]][,1]),])
    store_areas <- as.data.frame(levs)
    colnames(store_areas) <- c("ID", "Species")
    rownames( store_areas ) <-  store_areas$category
    
    # loop through the subplots
    for (i2 in 1:16){
      # get subplot extent
      ext = extent( subplots[[i2]] )
      clip <- crop(raster, ext)
      
      # count number of pixels
      area = ncell(clip)
      
      # count number os of pixels per class
      count_class = freq(clip)
      # estimate area of per class
      percent = (count_class[[1]][,2]*100)/area
      percent = data.frame(ID = count_class[[1]][,1], Cover = percent)
      #percent = data.frame(ID = count_class[[1]][,1]-1, Cover = percent)
      ## correct the shift
      #if ( is.na( percent$ID[length(percent[,1])] ) ){
      #  percent[length(percent[,1]), 1] <- levelsNumber
      #}
      
      # add value to store_areas matrix
      store_areas <- merge(store_areas, percent, by.x = "ID", all = TRUE)
      
      ### Merge species with flowers
      store_areas <- mergeSpecies(store_areas)
      
      # aggregate merge
      store_areas <- aggregate( . ~ Species, data=store_areas, sum, na.action = na.pass)
      colnames( store_areas )[2+i2] <- subplotNames[i2] 
      
    }
    
    store_areas$Species <- factor(store_areas$Species) 
    areasList [[i]] <- store_areas
  }
  
  # unlist 
  areasPlots <- do.call( "rbind", areasList )
  
  # get the median value per specie
  dummy_matrix <- matrix( ncol = ncol(areasList[[1]]), nrow = nrow(areasList[[1]]) )
  colnames(dummy_matrix) <- colnames(areasList[[1]])
  dummy_matrix[,1] <- as.character( areasList[[1]][,1] )
  dummy_matrix[,2] <- areasList[[1]][,2]
  
  for (i in 1:length( unique(areasList[[1]]$Species) )){#(levelsNumber-1)
    x = grep(areasList[[1]]$Species[i], areasPlots$Species)
    sp = areasPlots[x, ]
    for (i2 in 1:15){
      y = sp[,i2+2]
      y = na.omit(y)
      med = median(y) 
      dummy_matrix[i,i2+2] <- as.numeric(med)
    }
  }
  
  MedianCover = as.data.frame(dummy_matrix)
  
  # observed subplot covers
  Species_plot <- as.data.frame(levs)
  colnames(Species_plot) <- c("ID", "Species")
  for (i in 1:16){
    x = grep( subplotNames[i], obs_plot$Subplot )
    obs_subplot <- obs_plot[x, ]
    obs_subplot <- obs_subplot[, 4:5]
    Species_plot <- merge(Species_plot, obs_subplot, by.x = "Species", all = TRUE)
    colnames( Species_plot )[2+i] <- subplotNames[i2] 
  }
  
  x =grep("flowers", Species_plot$Species)
  Species_plot <- Species_plot[-x, ]
  x = duplicated(Species_plot$Species)
  Species_plot <- Species_plot[!x,]
  
  ################################
  ### Export results to tables ###
  ################################
  
  # create folder to store results
  dir.create( file.path(home, "Covers_results"), showWarnings = FALSE)
  
  # create output names
  boot_covers = paste0(algorithm, "_Boot_Covers_", maskName, "_", modelTag, ".txt")
  boot_covers = file.path( home, "Covers_results", boot_covers )
  
  median_covers = paste0( algorithm, "_Median_Covers_", maskName, "_", modelTag, ".txt")
  median_covers = file.path( home, "Covers_results", median_covers )
  
  observed_covers = paste0( algorithm, "_Obs_covers", maskName, "_", modelTag, ".txt")
  observed_covers = file.path( home, "Covers_results", observed_covers )
  
  # write results
  write.table(areasPlots,   file = boot_covers, sep = " ", row.names = F, col.names = T )
  write.table(MedianCover,  file = median_covers, sep = " ", row.names = F, col.names = T )
  write.table(Species_plot, file = observed_covers, sep = " ", row.names = F, col.names = T )
  
}

##------------------------------------------------------------------------------##
##                                                                              ##
## plots: visualization of spectras with the 5, 25, 50, 75 and 95 percentiles   ##
##                                                                              ##
## Arguments:                                                                   ##
## - spectra    spectral information. Used to create the quantiles of spectra   ##
## - wl         classificationEnsemble object                                   ##
## - xaxis      if TRUE the X axis is ploted                                    ##
## - ylab       if TRUE the Y axis label is ploted                              ##
## - ylabside   if TRUE the Y axis label is ploted at the left side of the plot ##
## - ymax       maximum value of the Y axis                                     ##
##                                                                              ##
##------------------------------------------------------------------------------##

plot.spectra <- function(spectra, wl, xaxis=TRUE, ylab = TRUE, ymax, ...){
  
  # obtain the quantiles of the spectras
  quant <- apply(spectra, 2, quantile, probs =c(0.05, 0.25, 0.5, 0.75, 0.95))
  
  if (xaxis == TRUE){
    
    plot(wl, quant[1,], type="l", ylim = c(0,ymax), xlim = c(min(wl)-10, max(wl)+10), 
         xaxs = "i", ylab= NA, las=1, xlab=expression(lambda(nm)), axes=F, cex.lab = 1.3)
    lines(wl, quant[2,], type="l")
    lines(wl, quant[3,], type="l")
    lines(wl, quant[4,], type="l")
    lines(wl, quant[5,], type="l")
    polygon(c(wl, rev(wl)), c(quant[2,], rev(quant[1,])), col = "grey70")
    polygon(c(wl, rev(wl)), c(quant[3,], rev(quant[2,])), col = "grey50")
    polygon(c(wl, rev(wl)), c(quant[4,], rev(quant[3,])), col = "grey50")
    polygon(c(wl, rev(wl)), c(quant[5,], rev(quant[4,])), col = "grey70")
    axis(side = 1, pos = 0, las=1, cex.axis = 1.3)
    if (ylab == TRUE){ 
      axis(side = 2, las=1, pos=min(wl), cex.axis = 1.3)
      mtext(side = 2, line = 3, 'Reflectance', cex=1.3)
    }
  } else {
    
    plot(wl, quant[1,], type="l", ylim = c(0,ymax), xlim = c(min(wl)-10, max(wl)+10), 
         xaxs = "i", axes=F, ylab = NA, xlab=NA, las=1, cex.lab = 1.3)
    lines(wl, quant[2,], type="l")
    lines(wl, quant[3,], type="l")
    lines(wl, quant[4,], type="l")
    lines(wl, quant[5,], type="l")
    polygon(c(wl, rev(wl)), c(quant[2,], rev(quant[1,])), col = "grey70")
    polygon(c(wl, rev(wl)), c(quant[3,], rev(quant[2,])), col = "grey50")
    polygon(c(wl, rev(wl)), c(quant[4,], rev(quant[3,])), col = "grey50")
    polygon(c(wl, rev(wl)), c(quant[5,], rev(quant[4,])), col = "grey70")
    axis(side = 2, pos=min(wl), labels=F, lwd.ticks=0)
    if (ylab == TRUE){ 
      axis(side = 2, las=1, pos=min(wl), cex.axis = 1.3)
      mtext(side = 2, line = 3, 'Reflectance', cex=1.3)
    }
  }
 }

##------------------------------------------------------------------------------##
##                                                                              ##
## plots: the MRPP variable importance and select the values over 0.3           ##
##                                                                              ##
## Arguments:                                                                   ##
## - spectra    spectral information. Used to create the quantiles of spectra   ##
## - wl         classificationEnsemble object                                   ##
##                                                                              ##
##------------------------------------------------------------------------------##
plot.importance <- function(varImport, wl, xaxis=TRUE, ...){
  
  library(grDevices)
  library(RColorBrewer)
  
  # variables from ensemble
  A <- varImport[1, ]
  imp = A
  imp[imp > 0.3] = 1 ; imp[imp < 0.3] = 0 
  imp <- as.logical(imp)

  # matices of varImport
  z1 <- matrix (rep (A, 100), ncol=100)
  z1[,0:50] <- NA
  # selection
  z2 <- matrix (rep (imp, 100), ncol=100)
  z2[z2==0] <- NA
  z2[,50:100] <- NA
  
  # add coefficients
  blueish <- colorRampPalette(brewer.pal(9,"Blues"))(100) 
  
  # MRFF coefficients
  if (xaxis == TRUE){
    image(wl, seq(0, 100, 1), z1, xlim = c(min(wl)-10, max(wl)+10), xlab=expression(lambda(nm)), col=blueish, ylab="", axes=F, cex.lab = 1.3)
  } else {
    image(wl, seq(0, 100, 1), z1, xlim = c(min(wl)-10, max(wl)+10), xlab="", col=blueish, ylab="", axes=F, cex.lab = 1.3)
  }
  # selection
  image(wl, seq(0, 100, 1), z2, xlim = c(min(wl)-10, max(wl)+10), col=1, ylab="", axes=F, add=T, cex.lab = 1.3)
  if (xaxis == TRUE){
    axis(side = 1, las=1, cex.axis = 1.3)
  }
  abline(h = 49.5)
  box()
  
}


##----------------------------------------------------------------------------##
##                                                                            ##
## significanceTest_CanopyLevel: Apply one-side bootstrap significance test to##
##                             canopy-level data                              ##
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
  
significanceTest_CanopyLevel <- function(model1, model2){
  
  ### PLS
  
  ### compute the differences between ASD and AISA band settings
  ### OA of rip-it-off should be larger. So, if OA(m1) - OA(m2) is positive, rip-it-off is significantly better.
  PLS_OA <- model1$OA_PLS - model2$OA_PL
  
  ### Kappa of rip-it-off should be larger. So, if OA(m1) - OA(m2) is positive, rip-it-off is significantly better.
  PLS_kappa <- model1$Kappa_PLS - model2$Kappa_PLS
  
  ### RF
  RF_OA <- model1$OA_RF - model2$OA_RF
  RF_kappa <- model1$Kappa_RF - model2$Kappa_RF
  
  ### SVM
  SVM_OA <- model1$OA_SVM - model2$OA_SVM
  SVM_kappa <- model1$Kappa_SVM - model2$Kappa_SVM
  
  # prepare output
  output <- list(PLS_OA, RF_OA, SVM_OA, PLS_kappa, RF_kappa, SVM_kappa)
  names(output) <- c("PLS_OA", "RF_OA", "SVM_OA", "PLS_kappa", "RF_kappa", "SVM_kappa")
  class(output) <- "boot_test"
  output
}

##----------------------------------------------------------------------------##
##                                                                            ##
## Small functions to help in repetitive tasks                                ##
##                                                                            ##
##                                                                            ##
##----------------------------------------------------------------------------##

###############################################
## stratified bootstrap selection of samples ##
###############################################

# So each species get selected for training and validation 
stratifySampling <- function(data, classes){
  TRAIN <- list()
  VAL <- list()
  
 for (i in 1:length( levels(data$classes ))){
    x = grep( levels(data$classes)[i], data$classes  )
    x <- data[x, ]
    # random sampling
    if ( length(x[, 1])==1 ){
      TRAIN[[i]] <- x
      VAL[[i]] <- x
    }
    else {
      idx = sample(1:length(x[,1]), length(x[,1]), replace=TRUE)
      TRAIN[[i]] <- x[idx, ]
      VAL[[i]] <- x[-idx, ]
    }
  }
  
  # unlist
  TRAIN2 <- do.call("rbind", TRAIN)
  VAL2   <- do.call("rbind", VAL)
  
  # prepare exit
  output <- list(TRAIN2, VAL2)
  names(output) <- c("train", "validation")
  output
  
}

###############################################
## List the names of the rasters in a folder ##
###############################################

rasterListNames <- function(fileExtantion, folder, dir=NULL){
  # if dir = NULL, set it to "home" by default
  if (is.null(dir)){
    dir = home
  }
  
  setwd(dir)
  
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # delete the ".dat" from the name
  x = grep(".tif.aux.xml", rast_list)
  if ( length(x) > 0 ){ rast_list <- rast_list[-x] }
  rast_list = gsub('.{4}$', '', rast_list)
  return(rast_list)
}

#####################################################
## List and load the rasters contained in a folder ##
#####################################################

rasterList <- function(fileExtantion, folder, dir=NULL, select=NULL){
  # if dir = NULL, set it to "home" by default
  if (is.null(dir)){
    dir = home
  }
  
  setwd(dir)
  
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

###################################################
## Merge species when flower classes are present ##
###################################################

mergeSpecies <- function(data){
  
  if ( length( grep("Achillea_millefolium_flowers", data$Species) ) > 0 ){
    x = grep( "Achillea_millefolium_flowers", data$Species )
    data$Species[x] <- "Achillea_millefolium"
  }
  
  if ( length( grep("Anagallis_arvensis_flowers", data$Species) ) > 0 ){
    x = grep( "Anagallis_arvensis_flowers", data$Species )
    data$Species[x] <- "Anagallis_arvensis"
  }
  
  if ( length( grep("Cichorium_intybus_flowers", data$Species) ) > 0 ){
    x = grep( "Cichorium_intybus_flowers", data$Species )
    data$Species[x] <- "Cichorium_intybus"
  }
  
  if ( length( grep("Daucum_carota_flowers", data$Species) ) > 0 ){
    x = grep( "Daucum_carota_flowers", data$Species )
    data$Species[x] <- "Daucum_carota"
  }
  
  if ( length( grep("Echium_vulgare_flowers", data$Species) ) > 0 ){
    x = grep( "Echium_vulgare_flowers", data$Species )
    data$Species[x] <- "Echium_vulgare"
  }
  
  if ( length( grep("Erigoron_annuus_flowers", data$Species) ) > 0 ){
    x = grep( "Erigoron_annuus_flowers", data$Species )
    data$Species[x] <- "Erigoron_annuus"
  }
  
  if ( length( grep("Medicago_lupulina_flowers", data$Species) ) > 0 ){
    x = grep( "Medicago_lupulina_flowers", data$Species )
    data$Species[x] <- "Medicago_lupulina"
  }
  
  if ( length( grep("Sp_13_flowers", data$Species) ) > 0 ){
    x = grep( "Sp_13_flowers", data$Species )
    data$Species[x] <- "Sp_13"
  }
  
  if ( length( grep("Sp_53_flowers", data$Species) ) > 0 ){
    x = grep( "Sp_53_flowers", data$Species )
    data$Species[x] <- "Sp_53"
  }
  
  if ( length( grep("Trifolium_pratense_flowers", data$Species) ) > 0 ){
    x = grep( "Trifolium_pratense_flowers", data$Species )
    data$Species[x] <- "Trifolium_pratense"
  } 
  
  if ( length( grep("Verbena_officinalis_flowers", data$Species) ) > 0 ){
    x = grep( "Verbena_officinalis_flowers", data$Species )
    data$Species[x] <- "Verbena_officinalis"
  } 
  
  if ( length( grep("Crepis_capillaris_flowers", data$Species) ) > 0 ){
    x = grep( "Crepis_capillaris_flowers", data$Species )
    data$Species[x] <- "Crepis_capillaris"
  } 
  
  data
}

###########################
## Delate flower classes ##
###########################

delFlowers <- function(data){
  
  if ( length( grep("Achillea_millefolium_flowers", data$Species) ) > 0 ){
    x = grep( "Achillea_millefolium_flowers", data$Species )
    data[-x, ] <- "Achillea_millefolium"
  }
  
  if ( length( grep("Anagallis_arvensis_flowers", data$Species) ) > 0 ){
    x = grep( "Anagallis_arvensis_flowers", data$Species )
    data[-x, ] <- "Anagallis_arvensis"
  }
  
  if ( length( grep("Cichorium_intybus_flowers", data$Species) ) > 0 ){
    x = grep( "Cichorium_intybus_flowers", data$Species )
    data[-x, ] <- "Cichorium_intybus"
  }
  
  if ( length( grep("Daucum_carota_flowers", data$Species) ) > 0 ){
    x = grep( "Daucum_carota_flowers", data$Species )
    data[-x, ] 
  }
  
  if ( length( grep("Echium_vulgare_flowers", data$Species) ) > 0 ){
    x = grep( "Echium_vulgare_flowers", data$Species )
    data[-x, ]
  }
  
  if ( length( grep("Erigoron_annuus_flowers", data$Species) ) > 0 ){
    x = grep( "Erigoron_annuus_flowers", data$Species )
    data[-x, ] 
  }
  
  if ( length( grep("Medicago_lupulina_flowers", data$Species) ) > 0 ){
    x = grep( "Medicago_lupulina_flowers", data$Species )
    data[-x, ] 
  }
  
  if ( length( grep("Sp_13_flowers", data$Species) ) > 0 ){
    x = grep( "Sp_13_flowers", data$Species )
    data[-x, ] 
  }
  
  if ( length( grep("Sp_53_flowers", data$Species) ) > 0 ){
    x = grep( "Sp_53_flowers", data$Species )
    data[-x, ] 
  }
  
  if ( length( grep("Trifolium_pratense_flowers", data$Species) ) > 0 ){
    x = grep( "Trifolium_pratense_flowers", data$Species )
    data[-x, ]
  } 
  
  if ( length( grep("Verbena_officinalis_flowers", data$Species) ) > 0 ){
    x = grep( "Verbena_officinalis_flowers", data$Species )
    data[-x, ] 
  } 
  
  if ( length( grep("Crepis_capillaris_flowers", data$Species) ) > 0 ){
    x = grep( "Crepis_capillaris_flowers", data$Species )
    data[-x, ] 
  } 
  
  data
}

###############################################
## Obtain cover summay per validation method ##
###############################################

coverSummary <- function(validation, na.replace = TRUE){ 
  
  modelTag = c("PLS", "RF", "SVM")
  Normalization = c("spect", "spect_BN", "MNF", "MNF_BN")
  
  output <- data.frame( Species=character(), Obserced=double(),Predicted=double(),	
                        Model=character(), Normalization=character(), plot=integer() )
  
  for (i in 1:4){ # Normalization 
    
    for (i2 in 1:3){ # modelTag
      
      for (i3 in 1:11){ # plots 
        
        obs <- read.table( paste0( "Covers_results/", modelTag[i2], "_", validation, "_Obs_coversplot_", 
                                  seq(9,19,1)[i3], "_", Normalization[i], ".txt"), header = T)
        pred <-  read.table( paste0( "Covers_results/", modelTag[i2], "_", validation, "_Median_Covers_plot_", 
                                    seq(9,19,1)[i3], "_", Normalization[i], ".txt"), header = T)
        
        for (i4 in 1:16){ #subplots
          
          sp <- obs$Species
          OBS <- obs[, i4+2]
          PRED <- pred[, i4+2]  
          
          if (nrow(pred) != nrow(obs)){
            minim = min(c(nrow(pred), nrow(obs)))
            df <- data.frame(Species=sp[1:minim], Observed=OBS[1:minim], Predicted=PRED[1:minim],
                             Model=modelTag[i2], Normalization=Normalization[i], Plot=seq(9,19,1)[i3])
          }
          if (nrow(pred) == nrow(obs)){
            df <- data.frame(Species=sp, Observed=OBS, Predicted=PRED,
                             Model=modelTag[i2], Normalization=Normalization[i], Plot=seq(9,19,1)[i3])
          }
          
          # eliminate rows when OBS= NA and PRED = NA 
          df <- df[rowSums(is.na(df)) != 2, ]
          # replace NA for zeros
          if (na.replace == TRUE){ 
            df[is.na(df)] <- 0
          }
          if (length(output[,1])==0){
            output = df
          } else { 
            output <-  merge(output, df, by = intersect(names(output), names(df)), all = TRUE)
          }
          
        }
        
      }
      
    }
    
  }
  
  output
  
}

#################################################
### Obtain R2, RMSE and Bias per specie       ###
### Inputs: is a results from coverSummary    ###
#################################################

GOF <- function(data){ 
  
  spList = factor( unique(data$Species) )
  
  outputGOF <- data.frame(Species=character(), r2=double(), RMSE=double(),
                          bias=double(), Models=character(), Normalization=character())
  
  for (i in 1:length(spList)){
    tryCatch({
  
      x = grep( spList[i], data$Species)
      pv = data[x, ]

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
      
      
      ## Goodness-of-fits
      out <- matrix(nrow = 24, ncol = 6)
      colnames(out) <- c("Species", "r2", "RMSE", "bias", "Models", "Normalization")
      
      out[,1] <- as.character(spList[i])
      out[,5] <- rep(c("PLS", "RF", "SVM"), 8)
      out[,6] <- rep( c(rep("spectra", 3), rep("spectra_BN",3), rep("MNF", 3), rep("MNF_BN", 3)), 2 )

      # r2
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
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
  }
  
  outputGOF$r2   <- as.numeric(as.character(outputGOF$r2))
  outputGOF$RMSE <- as.numeric(as.character(outputGOF$RMSE))
  outputGOF$bias <- as.numeric(as.character(outputGOF$bias))
  
  outputGOF
  
}

#################################################
### same as the GOF function, but only for    ###
### one model                                 ###
### Inputs: is a results from coverSummary    ###
###                                           ###
#################################################

GOFbest <- function(data, modelTag){ 
  
  spList = factor( unique(data$Species) )
  
  outputGOF <- data.frame(Species=character(), r2=double(), RMSE=double(),
                          bias=double())
  
  for (i in 1:length(spList)){
    tryCatch({
      
      x = grep( spList[i], data$Species)
      model = data[x, ]
      
      ## Goodness-of-fits
      out <- matrix(nrow = 1, ncol = 4)
      colnames(out) <- c("Species", "r2", "RMSE", "bias")
      
      out[,1] <- as.character(spList[i])
       
      # r2
      out[1,2] <- (cor(model$Predicted, model$Observed, method="pearson"))^2
      # RMSE
      out[1,3] <- sqrt(mean((model$Observed - model$Predicted)^2))
      # bias
      out[1,4] <- 1 - coef( lm(model$Predicted ~ model$Observed - 1) )
 
      if ( length(outputGOF[,1])==0 ){
        outputGOF <- as.data.frame(out)
      }
      if ( length(outputGOF[,1])!=0 ){
        outputGOF <- merge(outputGOF, as.data.frame(out), by = intersect(colnames(outputGOF), colnames(out)), all = TRUE)
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
  }
  
  outputGOF$r2   <- as.numeric(as.character(outputGOF$r2))
  outputGOF$RMSE <- as.numeric(as.character(outputGOF$RMSE))
  outputGOF$bias <- as.numeric(as.character(outputGOF$bias))
  
  outputGOF
  
}
  
#####################################################
### Count for well, miss and over classifications ###
#####################################################
  
ClassPresence <- function(data){
  data$ClassPresence <- NA
  
  for ( i in 1:nrow(data) ){
    if(data$Observed[i]==0 & data$Predicted[i]!=0){
      data$ClassPresence[i] <- "Over"
    }
    if(data$Observed[i]!=0 & data$Predicted[i]==0){
      data$ClassPresence[i] <- "Miss"
    }
    if(data$Observed[i]!=0 & data$Predicted[i]!=0){
      data$ClassPresence[i] <- "Well"
    }
  }
  data
}

###############################################
### function to calculate Camargo's eveness ###
###############################################
  
camargo <- function(n_spec, include_zeros = T)
{
  if (include_zeros) n <- n_spec else n <- n_spec[n_spec > 0]
  S <- length(n)
  camar <- 1
  for (i in 1:(S - 1))
  {
    for (j in (i + 1):S)
    {
      p_i <- n[i]/sum(n)
      p_j <- n[j]/sum(n)
      camar <- camar - abs(p_i - p_j)/S
    }
  }
  return(camar)
}
