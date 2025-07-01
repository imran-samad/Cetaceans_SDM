#This script takes in species occurrence and environmental data and obtains tuned hyperparameters for the six models. These hyperparameters are used in the subsequent model fitting script.

library(raster)
library(sdm)
library(ENMeval)
library(tidyverse)
library(sf)
library(ecospat)
library(caret)
library(earth)
library(doParallel)
library(dismo)
library(gbm)
library(rJava)

#Set seed for replicating results
set.seed(1243)

rasterlist =  list.files("./Data/Rasters/Processed/", pattern= "\\.tif$", full.names=T) 
rasterlist

#Keeping 8 uncorrelated predictors
rasterlist = rasterlist[c(1,4,9,13,22,25,26,27)]
rasterlist = rasterlist[c(1,2,7,3,4,5,6,8)]
rasterlist_all = rasterlist

#Function to unregister cores used in parallel processing
unregister_dopar = function() {
  env = foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#Loading presence data
sightings = read.csv("./Data/Sightings_processed.csv")

#create data frames to store tuning parameters
tune_glm_par = tune_gam_par = tune_maxent_par = tune_mars_par = tune_brt_par = tune_rf_par = tune_svm_par = tune_maxent_table = data.frame()

#for each species
for(spc in 1:length(unique(sightings$name_combined)))
{
  if(unique(sightings$name_combined)[spc] == "Humpback_dolphin" | unique(sightings$name_combined)[spc] == "Indo-Pacific_finless_porpoise" | unique(sightings$name_combined)[spc] == "Irrawaddy_dolphin"){
    spcs = 1 #Coastal
  }else{
    spcs = 2 #Oceanic
  }
  
  #Choose rasters based on species type
  if(spcs == 1)
  {
    rasterlist = rasterlist_all[-8] # because TempF don't act near coasts. Not removing this leads to a lot of NA values and drops data points
    pred_ras = stack(rasterlist)
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR")
  } else
  {
    pred_ras = stack(rasterlist_all)
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR", "TempF")
  }
  
  
  species = sightings |> filter(name_combined == unique(sightings$name_combined)[spc])
  
  #Thin species records to reduce bias
  species_thinned = spThin::thin(species, lat.col = "lat", long.col = "lon", spec.col = "name_combined", thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = F, write.log.file = FALSE)[[1]] #thin.par = 100 (for thinning at 100 km scale)
  
  #Extract covariate data for the thinned locations
  species_thinned = cbind(species_thinned, raster::extract(pred_ras, species_thinned))
  colnames(species_thinned)[1:2] = c("x","y")
  species_thinned$Depth = species_thinned$Depth * -1
  species_thinned$PA = 1

  #sampling bg from buffer
  spdf = st_as_sf(species_thinned, coords = c(1:2))
  spdf = st_buffer(spdf, ifelse(spcs == 1, 1, 2)) # buffer radius in degrees
  buff_ras = raster::mask(pred_ras[[1]], spdf)
  bg = data.frame(sampleRandom(buff_ras, 10000, xy = TRUE)) # number of bg points
  bg = bg[,-3]
  bg = cbind(bg, raster::extract(pred_ras, bg))
  colnames(bg)[1:2] = c("x","y")
  bg$Depth = bg$Depth * -1
  bg$PA = 0
  
  #model tuning gam
  # tuneGrid <- expand.grid(
  #   select = c(FALSE),
  #   method = c("REML", "GCV.Cp")
  # )
  # fitControl <- trainControl(method = "cv", number = 3, verboseIter = TRUE)
  # fit.crt = caret::train(PA~.,
  #                 data = training, method = "gam",
  #                 trControl = fitControl, tuneGrid = tuneGrid, family = binomial) #tune.length = 3
  # fit.crt$call
  # 
  # tune_gam_par = rbind(tune_gam_par, cbind("species" = unique(sightings$name_combined)[spc], fit.crt$bestTune))

  #Model tuning Maxent
  unregister_dopar() #Keeping cores registered for parallel processing causes issues with some functions, so it is ideal to unregister cores at the start of each model tuning
  tune_maxent = ENMevaluate(occs = species_thinned[,1:2], envs = pred_ras, bg = bg[,1:2], tune.args = list(fc = c("L","Q","LQ","LQT"), rm = c(0.5,1,2)), partitions = "randomkfold", algorithm = "maxent.jar", doClamp = T, parallel = T, numCores = parallel::detectCores())

  tune_maxent_table_temp = eval.results(tune_maxent)
  tune_maxent_table_temp = arrange(tune_maxent_table_temp, desc(auc.val.avg))
  tune_maxent_par = rbind(tune_maxent_par, cbind("species" = unique(sightings$name_combined)[spc], tune_maxent_table_temp))

  # tune_maxent_par = rbind(tune_maxent_par, data.frame("species" = unique(sightings$name_combined)[spc], "fc" =tune_maxent_table_temp$fc[1:3], "rm" = tune_maxent_table_temp$rm[1:3], "auc.val.avg" = tune_maxent_table_temp$auc.val.avg[1:3]))
  #eval.variable.importance(tune_maxent)

  #Model tuning MARS
  unregister_dopar()
  training = rbind(species_thinned, bg)
  training$PA = as.factor(training$PA)
  levels(training$PA) = c("c0","c1")
  training = training[,-c(1,2)]
  training = training |> dplyr::select("PA", everything()) |> na.omit()

  mytuneGrid = expand.grid(nprune = 2:20,
                            degree = 1) # no interaction
  mycontrol = trainControl(method = "cv",
                            number = 5, # 5-fold cross-validation
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary,
                            allowParallel = TRUE,
                            p = 0.7)

  tune_mars = caret::train(form = PA ~ .,
                    data = training,
                    method = "earth",
                    metric = "ROC",
                    trControl = mycontrol,
                    tuneGrid = mytuneGrid,
                    thresh = 0.00001)

  tune_mars_par = rbind(tune_mars_par, cbind("species" = unique(sightings$name_combined)[spc], tune_mars$bestTune))

  #SVM
  unregister_dopar()
  #levels(training$PA) = c("c0","c1")
  mycontrol = trainControl(method="cv",
                           number = 5,      # do 5 repetitions of cv
                      summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                      classProbs=TRUE,
                      allowParallel = TRUE,
                      p = 0.7)

  tune_svm = caret::train(PA~., data=training,
                   method = "svmRadial",   # Radial kernel
                   tuneLength = 5,                   # 5 values of the cost function
                   #preProc = c("center","scale"),  # Center and scale data
                   metric="ROC",
                   trControl=mycontrol)

  tune_svm_par = rbind(tune_svm_par, cbind("species" = unique(sightings$name_combined)[spc], tune_svm$bestTune))
  
  #model tuning BRT
  unregister_dopar()
  levels(training$PA) = c(0,1)
  prNum = as.numeric(table(training$PA)["1"]) # number of presences
  bgNum = as.numeric(table(training$PA)["0"]) # number of backgrounds
  wt = ifelse(training$PA == 1, 1, prNum / bgNum)
  training$PA = as.numeric(training$PA)
  if(max(as.numeric(training$PA))==2)
    training$PA = as.numeric(training$PA - 1)
  
  tune_settings = expand.grid(lrt = c(0.001, 0.005, 0.01), tc = 1:3, bag_rate = c(0.7))
  
  cl = makeCluster(detectCores())
  registerDoParallel(cl, cores = detectCores())
  clusterExport(cl = cl, varlist = 'training')
  registerDoSEQ()
  tune_brt = foreach(i = 1:nrow(tune_settings), .combine = rbind, .errorhandling = "remove") %dopar% {
    test = dismo::gbm.step(data = training,
                           gbm.x = 2:ncol(training),
                           gbm.y = 1, #response
                           family = "bernoulli",
                           tree.complexity = tune_settings$tc[i],
                           learning.rate = tune_settings$lrt[i],
                           bag.fraction = tune_settings$bag_rate[i],
                           n.folds = 10,
                           #n.trees = prNum,
                           site.weights = wt)
    
    print( c(i, test$n.trees, test$cv.statistics$deviance.mean) )
  }
  stopCluster(cl)
  tune_brt = data.frame(tune_brt)
  if(ncol(tune_brt)>1){
    colnames(tune_brt) = c("row","n.trees", "dev_mean")
    tune_settings[c(tune_brt$row),c('n_trees', 'dev_mean')] = tune_brt
  
  tune_brt_par = rbind(tune_brt_par, cbind("species" = unique(sightings$name_combined)[spc], tune_settings))
}else
  tune_brt_par = rbind(tune_brt_par, cbind("species" = unique(sightings$name_combined)[spc], "lrt" = NA, "tc" = NA, "bag_rate" = NA, "n_trees" = NA, "dev_mean" = NA))
  
  #model tuning RF
  unregister_dopar()
  training$PA = as.factor(training$PA)
  control = trainControl(method="cv", number=10, search="grid", allowParallel = TRUE, p=0.7)
  tunegrid = expand.grid(.mtry=c(1:8))

  tune_rf = caret::train(PA~., data=training, method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=control, weights = wt)

  tune_rf_par = rbind(tune_rf_par, cbind("species" = unique(sightings$name_combined)[spc], tune_rf$bestTune))
}

#Processing final tuning parameters
tune_brt_par_processed = tune_brt_par |> group_by(species) |> slice(which.min(dev_mean))
tune_maxent_par_processed = tune_maxent_par |> group_by(species) |> slice(which.max(auc.val.avg))

write.csv(tune_brt_par_processed, "./Results/Tuning/tune_brt_par.csv", row.names = F)
write.csv(tune_mars_par, "./Results/Tuning/tune_mars_par.csv", row.names = F)
write.csv(tune_rf_par, "./Results/Tuning/tune_rf_par.csv", row.names = F)
write.csv(tune_svm_par, "./Results/Tuning/tune_svm_par.csv", row.names = F)
write.csv(tune_maxent_par_processed, "./Results/Tuning/tune_maxent_par.csv", row.names = F)
