## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 


home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
# home = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

library(hyperSpec)

# load the data
data <- read.table("LeafHerbaceous.txt", sep = "", header = T)

# add Species
SpNames <- read.table("SpNamesLeafClip.csv", sep = "", header = T)
# add to data
data$Species <- SpNames$Species
# erase bad data
data <- subset(data, Species != "Erase")
data$Species <- factor(data$Species)
# spectral bands
spectra <- data[, 2:(length(data)-1)]

# create an hyperSpec object
new("hyperSpec")
hyperASD <- new("hyperSpec", spc=spectra, data=data, wavelength = seq(350, 2500, 1),
                label=list(spc="Reflectance", .wavelength =  expression(lambda(nm))))

plot(hyperASD)
plot(sample(hyperASD, 3))
plot(hyperASD[2,, 390~1000])
plot(hyperASD, "spcprctl5")

# eliminate the noisy 350-390 bands
hyperASD <- hyperASD[,, c(390~max)]
plot(sample(hyperASD, 10))
plot(hyperASD, "spcprctl5")

# set a data with the AISA+ spectral characteristics
hyperAISA <- hyperASD[,, c(390~990)]
hyperAISA <- spc.bin (hyperAISA, 10)
nwl(hyperAISA)
plot(hyperAISA, "spcprctl5")

par(mfrow=c(1,2),lend = 1, mai = c(1.2, 1.2, 0.5, 0.5))
plot(hyperASD, "spcprctl5")
plot(hyperAISA, "spcprctl5")
lines(svrcf)

#### Source Functions from GitHub
source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
} 
source_github("https://raw.githubusercontent.com/JavierLopatin/Herbaceous-Species-Classification/master/Scripts/Functions.R")

#########
subdata = data[1:100,]
subdata$Species <- factor(subdata$Species)
classes = subdata$Species
wl = hyperAISA@wavelength
spec = hyperAISA$spc[1:100,]

classificationLeafLevel <- function(classes, spec, wl=NA){
  
  if (is.na (wl[1]) | length (wl)!=ncol (classes))
    wl <- 1:ncol (classes)
  
  ## load required libraries
  #library (e1071)
  library (caret)
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
  controlObject <- trainControl(method = "repeatedcv", number = 10, repeats = 5,  classProbs=TRUE, allowParallel = TRUE)
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  ##########################
  ### PLS classification ###
  ##########################
  
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
  plscf <- as.vector (coef (plsClas$finalModel, ncomp=plsClas$finalModel$ncomp, intercept=F)) ## extract coeff.
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
  for(i in 1:ncol(spec)) 
    svrcf[i] <- svr.alpha %*% spec[svr.index, i]
  svrcf <- svrcf / sd (svrcf) ## scale pseudo-coefficients
  
  #####################################################################    
  ### get ensemble from all models and identify important variables ###
  #####################################################################
  
  ## get ensemble from all models and identify important variables
  ensemblecf <- abs (plscf) * plsrsq + abs (svrcf) * svrrsq + abs (rfcf) * rfrsq
  th <- mean (ensemblecf) + sd (ensemblecf) ## calculate threshold
  selbands <- ensemblecf > th ## apply threshold
  
  ######################
  ### prepare output ###
  ######################
  
  cf <- rbind (wl, plscf, rfcf, svrcf, ensemblecf, selbands)
  colnames (cf) <- colnames (x)
  
  fit <- c (plsrsq, rfrsq, svrrsq)
  names (fit) <- c ("PLS R2", "RF R2", "SVR R2")
  output <- list (cf, fit, th, pls, rf, svr)
  names (output) <- c ("selection", "fits", "threshold", "PLS", "RF", "SVM")
  class (output) <- "ensemble"
  output
  
  # stop parallel process
  stopCluster(cl)
 }

plot(output)

save.image("ClassLeafLevel.RData")
