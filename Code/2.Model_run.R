#This code runs hyperparameter optimised models with spatial cross validation for each species. Predictions are made for all models (including ensembles for all and the top 3 models) for all species and saved as rasters. Model summaries and relevant parameters are also saved as csv files.

library(raster)
library(sdm)
library(blockCV)
library(tidyverse)
library(sf)

#Set seed for replicating results
set.seed(1243)

rasterlist =  list.files("./Data/Rasters/Processed/", pattern= "\\.tif$", full.names=T) 

#Keeping 8 uncorrelated predictors
rasterlist = rasterlist[c(1,4,9,13,22,25,26,27)]
rasterlist = rasterlist[c(1,2,7,3,4,5,6,8)]
rasterlist_all = rasterlist

#Processing presence data
sightings = read.csv("./Data/Sightings_processed.csv")

#Loading tuning parameters
metadata = data.frame()
tune_brt_par = read.csv("./Results/Tuning/tune_brt_par.csv")
tune_mars_par = read.csv("./Results/Tuning/tune_mars_par.csv")
tune_rf_par = read.csv("./Results/Tuning/tune_rf_par.csv")
tune_svm_par = read.csv("./Results/Tuning/tune_svm_par.csv")
tune_maxent_par = read.csv("./Results/Tuning/tune_maxent_par.csv")

#For each species
for(spc in 1:length(unique(sightings$name_combined)))
{
  if(unique(sightings$name_combined)[spc] == "Humpback_dolphin" | unique(sightings$name_combined)[spc] == "Indo-Pacific_finless_porpoise" | unique(sightings$name_combined)[spc] == "Irrawaddy_dolphin"){
    spcs = 1 #Coastal
  }else{
    spcs = 2 #Oceanic
  }
  
  if(spcs == 1)
  {
    rasterlist = rasterlist_all[-8] # because TempF don't act near coasts. Not removing this leads to a lot of NA values and drops data points (rows)
    pred_ras = stack(rasterlist)
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR")
  } else
  {
    pred_ras = stack(rasterlist_all)
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR", "TempF")
  }
  
  if(unique(sightings$name_combined)[spc] == "Irrawaddy_dolphin") # Irrawaddy dolphin only occurs in a small subset of the study area
  {
    rasterlist = rasterlist_all[-8] # because TempF don't act near coasts
    pred_ras = stack(rasterlist)
    pred_ras = stack(crop(pred_ras, extent(82,100,17,26)))
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR")
  }
  
  species = sightings |> filter(name_combined == unique(sightings$name_combined)[spc])
  species[,c("lat","lon")] = sapply(species[,c("lat","lon")], as.numeric)
  
  #Thinning species occurrence to reduce bias
  species_thinned = spThin::thin(species, lat.col = "lat", long.col = "lon", spec.col = "name_combined", thin.par = 100, reps = 10, locs.thinned.list.return = TRUE, write.files = F, write.log.file = FALSE)[[1]]
  
  species_thinned = cbind(species_thinned, raster::extract(pred_ras, species_thinned))
  colnames(species_thinned)[1:2] = c("x","y")
  species_thinned$Depth = species_thinned$Depth * -1
  species_thinned$PA = 1
  
  #sampling bg from buffer
  spdf = st_as_sf(species_thinned, coords = c(1:2))
  spdf = st_buffer(spdf, ifelse(spcs==1,1,2)) # buffer radius in degrees; 1 degree for coastal and 2 degrees for oceanic species
  buff_ras = raster::mask(pred_ras[[1]], spdf)
  bg = data.frame(sampleRandom(buff_ras, 10000, xy = TRUE)) # number of bg points
  bg = bg[,-3]
  bg = cbind(bg, raster::extract(pred_ras, bg))
  colnames(bg)[1:2] = c("x","y")
  bg$Depth = bg$Depth * -1
  bg$PA = 0
  pb = rbind(species_thinned,bg)
  
  #Estimating spatial autocorrelation range across rasters for blocking
  sac = spatialAutoRange(rasterLayer = pred_ras, sampleNumber = 2000, doParallel = F, showPlots = F)
  pb_data = st_as_sf(pb, coords = c("x", "y"), crs = crs(pred_ras))
  
  # spatial blocking by specified range with random assignment
  sb = spatialBlock(speciesData = pb_data,
                     species = "PA",
                     rasterLayer = pred_ras,
                     theRange = round(sac$range, digits = -(floor(log10(sac$range))-1)), # size of the blocks
                     k = ifelse(unique(sightings$name_combined)[spc] == "Irrawaddy_dolphin",3,
                         ifelse(unique(sightings$name_combined)[spc] == "Humpback_dolphin",4,5)),
                     selection = "random",
                     iteration = 100, # find evenly dispersed folds
                     biomod2Format = TRUE,
                     xOffset = 0, # shift the blocks horizontally
                     yOffset = 0)
  pb$folds = sb$foldID
  
  
  row_brt = which(tune_brt_par$species == unique(sightings$name_combined)[spc])
  row_mars = which(tune_mars_par$species == unique(sightings$name_combined)[spc])
  row_rf = which(tune_rf_par$species == unique(sightings$name_combined)[spc])
  row_svm = which(tune_svm_par$species == unique(sightings$name_combined)[spc])
  row_maxent = which(tune_maxent_par$species == unique(sightings$name_combined)[spc])
  
  zz=1 #Counter for the loop below
  #loop across each spatial fold for each species; needs to be changed for Irrawaddy and Humpback dolphins
  for(i in 1:5)
  {
    tryCatch({
    if(sum(pb[which(pb$folds==i),-c(1,2,ncol(pb))][complete.cases(pb[which(pb$folds==i),-c(1,2,ncol(pb))]),]$PA) == 0)
      next
    
    sdm_data = sdmData(train = pb[which(pb$PA==1 & pb$folds!=i),-c(1,2,ncol(pb))],
                       bg = pb[which(pb$PA==0 & pb$folds!=i),-c(1,2,ncol(pb))],
                       test = pb[which(pb$folds==i),-c(1,2,ncol(pb))])
    
    t1 = Sys.time()
    assign(paste0("m_",i),
           sdm(data = sdm_data,
               methods = c("maxent", "gam", "brt", "mars", "rf", "svm"), 
               #replicatin='cv',
               #test.percent=30,
               #n=1,
               parallelSettings = list(ncore = 12, method = 'parallel'),
               modelSettings=list(gam=list(method="REML"),
                                  brt=list(n.trees=tune_brt_par$n_trees[row_brt], bag.fraction = tune_brt_par$bag_rate[row_brt], interaction.depth = tune_brt_par$tc[row_brt], shrinkage = tune_brt_par$lrt[row_brt]),
                                  mars=list(nprune=tune_mars_par$nprune[row_mars]),
                                  rf=list(mtry=tune_rf_par$mtry[row_rf], sampsize=nrow(species_thinned)),
                                  svm=list(sigma=tune_svm_par$sigma[row_svm], C=tune_svm_par$C[row_svm]),
                                  maxent=list(fc=tune_maxent_par$fc[row_maxent], rm=tune_maxent_par$rm[row_maxent]))))

    # assign(paste0("p",i), predict(get(paste0("m",i)), newdata=pred_ras, nc = 12, method = "gam", mean = T, overwrite = T))
    
    metadata = rbind(metadata, cbind("Sl" = spc,
                                     "Species" = unique(sightings$name_combined)[spc],
                                     "rep" = i,
                                     "nPresence" = nrow(species),
                                     "nPresenceThinned" = nrow(species_thinned),
                                     "nbg" = length(which(complete.cases(bg)==T)),
                                     "repPres" = nrow(pb[which(pb$PA==1 & pb$folds!=i & complete.cases(pb) == T),]),
                                     "repBg" = nrow(pb[which(pb$PA==0 & pb$folds!=i & complete.cases(pb) == T),]),
                                     "time_elapsed" = round(as.numeric(difftime(Sys.time(), t1, units='mins')),2)))
    
    if(zz==1)
    {
      assign(paste0("m",spc), get(paste0("m_",i)))
      rm(list = paste0("m_",i))
      zz = zz+1
    } else
    {
      assign(paste0("m",spc), get(paste0("m",spc))+get(paste0("m_",i)))
      rm(list = paste0("m_",i))
    }
    }, error=function(e){})
  }
}


#summarising results
thresholds = data.frame()
stats = data.frame()
var_imp = data.frame()

#loop to extract thresholds, stats and variable importance for each species and model
for(spc in 1:length(unique(sightings$name_combined)))
{
  sp = unique(sightings$name_combined)[spc]
  tm = get(paste0("m",spc))
  for(i in 1:length(tm@models$PA))
  {
    model_names = names(tm@models$PA)
    for(j in 1:length(tm@models$PA[[i]]))
    {
      if(length(tm@models$PA[[i]][[j]]@varImportance$test.indep) != 0)
      {
        thresholds = rbind(thresholds, cbind("Species" = sp,
                                             "Model" = model_names[i],
                                             "Model_no" = i,
                                             "rep" = j,
                                             "Train" = tm@models$PA[[i]][[j]]@evaluation$training@threshold_based,
                                             "Test" = if(length(tm@models$PA[[i]][[j]]@evaluation$test.indep) != 0){
                                               tm@models$PA[[i]][[j]]@evaluation$test.indep@threshold_based
                                             }else{
                                               data.frame("criteria" = NA, "threshold" = NA, "sensitivity" = NA, "specificity" = NA, "TSS" = NA, "Kappa" = NA, "NMI" = NA, "phi" = NA, "ppv" = NA, "npv" = NA, "ccr" = NA, "prevalence" = NA)}))
        
        stats = rbind(stats, cbind("Species" = sp,
                                   "Model" = model_names[i],
                                   "Model_no" = i,
                                   "rep" = j,
                                   "Train.AUC" = tm@models$PA[[i]][[j]]@evaluation$training@statistics$AUC,
                                   "Test.AUC" = ifelse(length(tm@models$PA[[i]][[j]]@evaluation$test.indep) != 0, tm@models$PA[[i]][[j]]@evaluation$test.indep@statistics$AUC, NA),
                                   "Training.COR" = tm@models$PA[[i]][[j]]@evaluation$training@statistics$COR[1],
                                   "Test.COR" = ifelse(length(tm@models$PA[[i]][[j]]@evaluation$test.indep) != 0, tm@models$PA[[i]][[j]]@evaluation$test.indep@statistics$COR[1], NA),
                                   "Training.Prevalance" = tm@models$PA[[i]][[j]]@evaluation$training@statistics$Prevalence,
                                   "Test.Prevalance" = ifelse(length(tm@models$PA[[i]][[j]]@evaluation$test.indep) != 0, tm@models$PA[[i]][[j]]@evaluation$test.indep@statistics$Prevalence, NA),
                                   "Training.Deviance" = tm@models$PA[[i]][[j]]@evaluation$training@statistics$Deviance,
                                   "Test.Deviance" = ifelse(length(tm@models$PA[[i]][[j]]@evaluation$test.indep) != 0, tm@models$PA[[i]][[j]]@evaluation$test.indep@statistics$Deviance, NA)))
        
        var_imp = rbind(var_imp, cbind("Species" = sp,
                                       "Model" = model_names[i],
                                       "Model_no" = i,
                                       "rep" = j,
                                       "Train" = tm@models$PA[[i]][[j]]@varImportance$training@varImportance,
                                       "Test" = tm@models$PA[[i]][[j]]@varImportance$test.indep@varImportance
        ))
      }
    }
  }
}

stats[,-c(1,2)] = sapply(stats[,-c(1,2)], as.numeric)
stats$Diff.AUC = abs(stats$Train.AUC - stats$Test.AUC)

#Summarising and visualising stats
stats_mean = stats |> group_by(Species, Model) |> summarise_all(mean, na.rm = T)

thresholds = thresholds[-which(colnames(thresholds) == "Test.criteria"),]
thresholds_mean = thresholds |> group_by(Species, Model, Train.criteria) |> summarise_all(mean, na.rm = T)

var_imp_mean = var_imp |> group_by(Species, Model, Train.variables) |> summarise_all(mean, na.rm = T)

stats_mean = thresholds_mean[,c("Species", "Model", "Train.criteria", "Test.TSS")] |> filter(Train.criteria == "max(se+sp)") |> merge(stats_mean, by = c("Species", "Model"))
stats_mean = stats_mean[,-3]
ggplot(stats_mean, aes(y = Test.AUC, x = Test.TSS, col = Species)) + geom_point() + geom_smooth(method = "lm")
stats_mean |> ggplot(aes(x = Model, y = Test.TSS, col = "TSS")) + geom_point() + facet_wrap(~Species) + theme(axis.text.x = element_text(angle = 45)) + geom_point(aes(x = Model, y = Test.AUC, col = "AUC"))

write.csv(stats_mean, "./Results/ModelStats.csv", row.names = F)
write.csv(thresholds_mean, "./Results/ModelThresholds.csv", row.names = F)
write.csv(var_imp_mean, "./Results/ModelVarImpMean.csv", row.names = F)

#getting the top 3 models for creating an ensemble
best_models = stats_mean[, c("Species", "Model", "Test.TSS", "Test.AUC", "Diff.AUC")] |> arrange(desc(Test.AUC)) |> group_by(Species) |> slice(1:3)

#need to predict all models for CBI evaluation
all_models = stats_mean[, c("Species", "Model", "Test.TSS", "Test.AUC", "Diff.AUC")]#predicting all models

#loop across each species for generating predictions
for(i in 1:length(unique(sightings$name_combined)))
{
  all_models_sp = all_models$Model[which(all_models$Species == unique(sightings$name_combined)[i])]
  if(unique(sightings$name_combined)[i] == "Humpback_dolphin" | unique(sightings$name_combined)[i] == "Indo-Pacific_finless_porpoise")
  {
    rasterlist = rasterlist_all[-8] # because TempF don't act near coasts
    pred_ras = stack(rasterlist)
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR")
  } else
  {
    pred_ras = stack(rasterlist_all)
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR", "TempF")
  }
  if(unique(sightings$name_combined)[spc] == "Irrawaddy_dolphin")
  {
    rasterlist = rasterlist_all[-8] # because TempF don't act near coasts
    pred_ras = stack(rasterlist)
    pred_ras = stack(crop(pred_ras, extent(82,100,17,26)))
    names(pred_ras) = c("Dist_land", "Depth", "Slope", "Cur_vel", "PP", "TempM", "TempR")
  }
  
  for(j in 1:length(all_models_sp))
  {
    assign(paste0("p",i,j), predict(get(paste0("m",i)), newdata=pred_ras, method = all_models_sp[j], mean = T, overwrite = T, nc = 12))#
    assign(paste0("p",i,j), mean(stack(get(paste0("p",i,j)))))
    writeRaster(get(paste0("p",i,j)), paste0("./Results/Rasters/", unique(sightings$name_combined)[i],"_",all_models_sp[j],".tif"))
    rm(list = paste0("p",i,j))
  }
  
  #Ensemble:All models
  assign(paste0("e",i), ensemble(get(paste0("m",i)), newdata=pred_ras, setting=list(method='weighted',stat='AUC')))
  writeRaster(get(paste0("e",i)), paste0("./Results/Rasters/", unique(sightings$name_combined)[i],"_ensemble_all",".tif"))
  rm(list = paste0("e",i))
  
  #For top3 models
  best_models_sp = best_models$Model[which(best_models$Species == unique(sightings$name_combined)[i])]
  for(j in 1:length(best_models_sp))
  {
    mod_pos = c(which(names(get(paste0("m",i))@models$PA) %in% best_models_sp))
      if(!exists("ss")){
          ss = sdm::subset(get(paste0("m",i)),mod_pos[j],drop=F)
      }else{
          ss = ss + subset(get(paste0("m",i)),mod_pos[j],drop=F)
        }
  }
  assign(paste0("e",i), ensemble(ss, newdata=pred_ras, setting=list(method='weighted',stat='AUC')))
  writeRaster(get(paste0("e",i)), paste0("./Results/Rasters/", unique(sightings$name_combined)[i],"_ensemble_t3",".tif"))
  rm(list = paste0("e",i))
  rm(ss)
}