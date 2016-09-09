## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

#### run Clasification!

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
# home = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
polyOut = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp/shp_outPixel"
# polyOut = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp/shp_outPixel"

setwd(home)

pkgs<-c("caret", "raster", "rgdal", "e1071", "gtools", "doParallel", "autopls", "maptools")
lapply(pkgs, require, character.only=T)

load("ClassPixel.RData")

# load species table
species <- read.table("Species.csv", header = T, sep=",")

#### tune and perform SVM cassifier to all plots
# images of only plots 
imgIn <- rasterlist[9:length(rasterlist)-1]

# list of accuracies
PA.spec    <- list()
PA.MNF     <- list()
PA.GLCM    <- list()

UA.spec    <- list()
UA.MNF     <- list()
UA.GLCM    <- list()

OA.spec    <- list()
OA.MNF     <- list()
OA.GLCM    <- list()

kappa.spec <- list()
kappa.MNF  <- list()
kappa.GLCM <- list()

Pred_spc  <- list()
Pred_MNF  <- list()
Pred_GLCM <- list()
for (i in 2:length(imgIn)){
      # obtain training data for the plot
      x = grep(i, species$Plot) # which observations are in this plot?
      y <- species[x, ]
      z <- list()
      for (i2 in 1:length(unique(y$Species))){
        xx = grep(unique(y$Species)[i2], species.all$spp)
        yy <- species.all[xx, ]
        z[[i2]] <- yy
      }
      
      # data frame for the plot
      plot.all <- do.call("rbind", z)
      
      # Set the random number seed so we can reproduce the results
      set.seed(123)
      # Split data in training and test
      forTraining <- createDataPartition(plot.all$spp, p = 0.6, list=F)
      train <- plot.all [ forTraining,]
      test<- plot.all [-forTraining,]
    
      # Each model used 5 repeated 10-fold cross-validation. Use AUC to pick the best model
      controlObject <- trainControl(method = "repeatedcv", number = 10, repeats = 5,  classProbs=TRUE, allowParallel = TRUE)
      
      ## set tuning paramiters
      grid <- expand.grid(cost  = seq(.001, .1, by = .05), gamma  = seq(1,100, by = 10))
      
      # initialize parallel processing
      cl <- makeCluster(detectCores())
      registerDoParallel(cl)
      
      # tune
      set.seed(123)
      tune.spec <- train(x=train[,1:61], y=train$spp, method = "svmLinear2", tuneLength=10, preProc = c("center", "scale"), trControl = controlObject)
      #set.seed(123)
      #tune.MNF <- train(x=train[,62:71], y=train$spp, method = "svmLinear2", tuneLength=10, preProc = c("center", "scale"), trControl = controlObject)
      #set.seed(123)
      #tune.GLCM <- train(x=train[,72:131], y=train$spp, method = "svmLinear2", tuneLength=10, preProc = c("center", "scale"), trControl = controlObject)
      
      # predict
      pred.spec    <- predict(tune.spec, test[,1:61])
      #pred.MNF     <- predict(tune.MNF, test[,62:71])
      #pred.GLCM    <- predict(tune.GLCM, test[,72:131])
     
      # confusion matix
      preMatrix <- function(pred, test){ # functionn to prevent caret error for different length
        u = union(pred, test)
        t = table(factor(pred, u), factor(test, u))
        return(t)
      }  
      conf.spec    <- confusionMatrix(preMatrix(pred.spec, test$spp))
      #conf.MNF     <- confusionMatrix(preMatrix(pred.MNF, test$spp))
      #conf.GLCM    <- confusionMatrix(preMatrix(pred.GLCM, test$spp))
      
      # get accuracies
      PA.spec[[i]]    <- conf.spec$byClass[,3] 
      #PA.MNF[[i]]     <- conf.MNF$byClass[,3]
      #PA.GLCM[[i]]    <- conf.GLCM$byClass[,3]
      
      UA.spec[[i]]    <- conf.spec$byClass[,4]
      #UA.MNF[[i]]     <- conf.MNF$byClass[,4]
      #UA.GLCM[[i]]    <- conf.GLCM$byClass[,4]
      
      OA.spec[[i]]    <- conf.spec$overall["Accuracy"]
      #OA.MNF[[i]]     <- conf.MNF$overall["Accuracy"]
      #OA.GLCM[[i]]    <- conf.GLCM$overall["Accuracy"]
        
      kappa.spec[[i]]    <- conf.spec$overall["Kappa"]
      #kappa.MNF[[i]]     <- conf.MNF$overall["Kappa"]
      #kappa.GLCM[[i]]    <- conf.GLCM$overall["Kappa"]
       
      # Apply model using all data
      model.spec    <- e1071:::svm(plot.all[,1:61], as.factor(plot.all$spp), gamma = tune.spec$bestTune$gamma, cost = tune.spec$bestTune$cost, probability = TRUE, kernel="linear")
      #model.MNF     <- e1071:::svm(plot.all[,62:71], as.factor(plot.all$spp), gamma = tune.MNF$bestTune$gamma, cost = tune.MNF$bestTune$cost, probability = TRUE, kernel="linear")
      #model.GLCM    <- e1071:::svm(plot.all[,72:131], as.factor(plot.all$spp), gamma = tune.GLCM $bestTune$gamma, cost = tune.GLCM $bestTune$cost, probability = TRUE, kernel="linear")
      
      #apply model to image
      Pred_spc[[i]]  <- predict(imgIn[[i]][[1:61]], model.spec, type="class")
      #Pred_MNF[[i]]  <- predict(imgIn[[i]][[62:71]], model.MNF, type="class")
      #Pred_GLCM[[i]] <- predict(imgIn[[i]][[72:131]], model.GLCM, type="class")
     
      # stop parallel process
      stopCluster(cl)
    }

save.image("ClassPixel.RData")

## convert prediction to shepfiles with the classes
for (i in 3:length(Pred_spc)){ 
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  ## poliginize to keep the species name label
  poly_spec <- rasterToPolygons(Pred_spc[[i]], dissolve = T)
  #poly_MNF  <- rasterToPolygons(Pred_MNF[[i]], dissolve = T)
  #poly_GLCM <- rasterToPolygons(Pred_GLCM[[i]], dissolve = T)

  ## add spp manes to the shapefiles
  # create a function to finde classes that were not predicted
  # lexicographic vector
  not.pred <- function(pred, poly){
    a = c(1, seq(10, length(levels(pred)[[1]][,1]), 1), seq(2,9,1))
    b = setdiff(levels(pred[[i]])[[1]][,1], poly@data$layer)  
    c = which(a == b)
    return(c)
  }
  
  if (length(levels(Pred_spc[[i]])[[1]][,1]) < 10){
    a = c(1,2,3,4,5,6,7,8,9)
  } else if  (length(levels(Pred_spc[[i]])[[1]][,1]) >= 10 & length(levels(Pred_spc[[i]])[[1]][,1]) < 20){
    a = c(1, seq(10, length(levels(Pred_spc[[i]])[[1]][,1]), 1), seq(2,9,1)) 
  } else {
    a = c(1, seq(10, 19, 1), 2, seq(20,length(levels(Pred_spc[[i]])[[1]][,1]),1), seq(3,9,1)) 
  }
  
  if (!is.null(setdiff(levels(Pred_spc[[i]])[[1]][,1], poly_spec@data$layer))){
    poly_spec$spp <- levels(Pred_spc[[i]])[[1]][,2][poly_spec@data$layer]
  } else { 
    b = as.vector(levels(Pred_spc[[i]])[[1]][,1]) 
    b[a]
    b[not.pred(Pred_spc[[i]], poly_spec)] <- NA
    poly_spec$spp <- na.omit(b)
    }
  
  #if (!is.null(setdiff(levels(Pred_MNF[[i]])[[1]][,1], poly_MNF@data$layer))){
 #   poly_MNF$spp <- levels(Pred_MNF[[i]])[[1]][,2][poly_MNF@data$layer]
 # } else { 
  #  b = as.vector(levels(Pred_MNF[[i]])[[1]][,1]) 
  #  b[a]
 #   b[not.pred(Pred_MNF[[i]], poly_MNF)] <- NA
 #   poly_MNF$spp <- na.omit(b) 
 # }
  
 # if (!is.null(setdiff(levels(Pred_GLCM[[i]])[[1]][,1], poly_GLCM@data$layer))){
#    poly_GLCM$spp <- levels(Pred_GLCM[[i]])[[1]][,2][poly_GLCM@data$layer]
#  } else { 
#    b = as.vector(levels(Pred_GLCM[[i]])[[1]][,1]) 
#    b[a]
#    b[not.pred(Pred_GLCM[[i]], poly_GLCM)] <- NA
#    poly_GLCM$spp <- na.omit(b) 
#  }
  
  # export shapefiles
  writeOGR(poly_spec, polyOut, paste("plot", i, "_spec", sep=""), driver="ESRI Shapefile", overwrite_layer = T)
  #writeOGR(poly_MNF,  polyOut, paste("plot", i, "_MNF", sep=""),  driver="ESRI Shapefile", overwrite_layer = T)
  #writeOGR(poly_GLCM, polyOut, paste("plot", i, "_GLCM", sep=""), driver="ESRI Shapefile", overwrite_layer = T)
  
  
  # stop parallel process
  stopCluster(cl)
}

save.image("ClassPixel.RData")
