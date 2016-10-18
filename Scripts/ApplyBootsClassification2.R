

ApplyBootsClassification2 <- function(data, en, Site, rasterPlots, boots=100, outDir, modelTag){  
  
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
  data2$classes <- factor(data2$classes)
  data2 <- na.omit(data2)
  
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
  
  # progress bar
  print(paste0(modelTag, "_", Site))
  pb <- txtProgressBar(min = 0, max = boots, style = 3)
  
  for (i in 1:boots){
    
    # progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    N = length(data2[,1])
    
    # create random numbers with replacement to select samples from each group
    idx = sample(1:N, N, replace=TRUE)
    
    train <- data2[idx,]
   # train$classes <- factor(train$classes)
    val <- data2[-idx,]
   # val$classes <- factor(val$classes)
    
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
    
    RF  <- randomForest( y = factor( train$classes ), x = train[,2:length(train)],
                         ntree= bestNtree, mtry = bestMtry )
    
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
  
  write.table(fits, file = file.path(outDir, paste("Fits_Site_", Site, "_", modelTag, ".txt", sep="")), row.names = F, col.names = T)
  write.table(fits_2, file = file.path(outDir, paste("FitsPA_OA_Site_", Site, "_", modelTag, ".txt", sep="")), row.names = F, col.names = T)
  write.table(predict_all, file = file.path(outDir, paste("Predicts_Site_", Site, "_", modelTag, ".txt", sep="")), row.names = F, col.names = T)
  
 }
