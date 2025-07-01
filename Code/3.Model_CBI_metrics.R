#This code generates evaluation metrics for ensemble models. CBI evaluation metrics and thresholds are generated for all models
library(raster)
library(ecospat)
library(dplyr)
library(sf)

#This functions generates key parameters like AUC, prevalance, etc. for ensemble models since they were not generated in the previous script
calc_stats = function(ras, spdf, species_thinned, sp, mod_name, thresholds, stats)
{
  buff_ras = raster::mask(ras, spdf)
  bg = data.frame(sampleRandom(buff_ras, 10000, xy = TRUE)) # number of bg points
  bg$PA = 0
  bg = bg[,-c(1,2)]
  colnames(bg) = c("HS","PA")
  pr = raster::extract(ras, species_thinned[,1:2])
  vals = rbind(cbind("HS"=pr, "PA"=rep(1,length(pr))), bg)
  ens_eval = sdm::evaluates(x = vals$PA, p = vals$HS)
  thresholds = rbind(thresholds, cbind("Species" = sp,
                                       "Model" = mod_name,
                                       "Test" = ens_eval@threshold_based))
  
  stats = rbind(stats, cbind("Species" = sp,
                             "Model" = mod_name,
                             "Test.AUC" = ens_eval@statistics$AUC,
                             "Test.COR" = ens_eval@statistics$COR[1],
                             "Test.Prevalance" = ens_eval@statistics$Prevalence,
                             "Test.Deviance" = ens_eval@statistics$Deviance))
  return(list(stats,thresholds))
}

#Evaluation using CBI
sightings = read.csv("./Data/Sightings_processed.csv")

thresholds = stats = cbi = cbi_thresholds = data.frame() # To store CBI metrics for each species and model
methods = c("gam","mars","rf","maxent","brt","svm")

for(spc in 1:length(unique(sightings$name_combined)))
{
  if(unique(sightings$name_combined)[spc] == "Humpback_dolphin" | unique(sightings$name_combined)[spc] == "Indo-Pacific_finless_porpoise" | unique(sightings$name_combined)[spc] == "Irrawaddy_dolphin"){
    spcs = 1
  }else{
    spcs = 2
  }
  set.seed(1243)
  species = sightings |> filter(name_combined == unique(sightings$name_combined)[spc])
  species_thinned = spThin::thin(species, lat.col = "lat", long.col = "lon", spec.col = "name_combined", thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = F, write.log.file = FALSE)[[1]]
  colnames(species_thinned)[1:2] = c("x","y")
  
  spdf = st_as_sf(species_thinned, coords = c(1:2))
  spdf = st_buffer(spdf, ifelse(spcs==1,1,2))
  
  for(j in methods)
  {
    if(!file.exists(paste0("./Results/Rasters/", unique(sightings$name_combined)[spc],"_",j,".tif"))) # incase some models failed to predict rasters
      next
    ras = raster(paste0("./Results/Rasters/", unique(sightings$name_combined)[spc],"_",j,".tif"))
    #ras = crop(ras, raster_extent)
    
    HS = raster::extract(ras, species_thinned) #Habitat suitability
    
    classification=data.frame(ecospat.boyce(na.omit(getValues(ras)), na.omit(HS), nclass=0, window.w="default", res=100, PEplot = F, rm.duplicate = TRUE, method = 'spearman'))
    cbi = rbind(cbi, cbind("Species" = unique(sightings$name_combined)[spc], "Model" = j, "cbi" = mean(classification$cor)))
    
    cbi_thresholds = rbind(cbi_thresholds, cbind("Species" = spc, "model" = method, "x1" = classification$HS[min(which(classification$F.ratio>=1))], "x2" = classification$HS[min(which(classification$F.ratio>=2))], "x3" = classification$HS[min(which(classification$F.ratio>=3))], "x5" = classification$HS[min(which(classification$F.ratio>=5))], "x10" = classification$HS[min(which(classification$F.ratio>=10))]))
    
    #Predicted to expected plot
    #ggplot(classification, aes(x=HS, y=F.ratio))+geom_line(cex=.6)+geom_hline(yintercept=20)+labs(x="Habitat Suitability", y="P/E ratio") +  theme(axis.text.x = element_text(size = 14,family="sans"), axis.title.x = element_text(size = 16,family="sans"), axis.text.y = element_text(size = 14,family="sans"), axis.title.y = element_text(size = 16,family="sans"))
  }
  
  if(!file.exists(paste0("./Results/Rasters/", unique(sightings$name_combined)[spc],"_maxent.tif")))
    next
  #Ensemble_all
  ras = raster(paste0("./Results/Rasters/", unique(sightings$name_combined)[spc],"_ensemble_all.tif"))
  HS = raster::extract(ras, species_thinned)
  classification=data.frame(ecospat.boyce(na.omit(getValues(ras)), na.omit(HS), nclass=0, window.w="default", res=100, PEplot = F, rm.duplicate = TRUE, method = 'spearman'))
  cbi = rbind(cbi, cbind("Species" = unique(sightings$name_combined)[spc], "Model" = "Ensemble_all", "cbi" = mean(classification$cor)))
  
  cbi_thresholds = rbind(cbi_thresholds, cbind("Species" = spc, "model" = "Ensemble_all", "x1" = classification$HS[min(which(classification$F.ratio>=1))], "x2" = classification$HS[min(which(classification$F.ratio>=2))], "x3" = classification$HS[min(which(classification$F.ratio>=3))], "x5" = classification$HS[min(which(classification$F.ratio>=5))], "x10" = classification$HS[min(which(classification$F.ratio>=10))]))
  
  temp_op = calc_stats(ras, spdf, species_thinned, sp = unique(sightings$name_combined)[spc], mod_name = "Ensemble_all", thresholds, stats)
  stats = temp_op[[1]]
  thresholds = temp_op[[2]]
  
  #Ensemble_t3
  ras = raster(paste0("./Results/Rasters/", unique(sightings$name_combined)[spc],"_ensemble_t3.tif"))
  HS = raster::extract(ras, species_thinned)
  classification=data.frame(ecospat.boyce(na.omit(getValues(ras)), na.omit(HS), nclass=0, window.w="default", res=100, PEplot = F, rm.duplicate = TRUE, method = 'spearman'))
  cbi = rbind(cbi, cbind("Species" = unique(sightings$name_combined)[spc], "Model" = "Ensemble_t3", "cbi" = mean(classification$cor)))
  
  cbi_thresholds = rbind(cbi_thresholds, cbind("Species" = spc, "model" = "Ensemble_t3", "x1" = classification$HS[min(which(classification$F.ratio>=1))], "x2" = classification$HS[min(which(classification$F.ratio>=2))], "x3" = classification$HS[min(which(classification$F.ratio>=3))], "x5" = classification$HS[min(which(classification$F.ratio>=5))], "x10" = classification$HS[min(which(classification$F.ratio>=10))]))
  
  temp_op = calc_stats(ras, spdf, species_thinned, sp = unique(sightings$name_combined)[spc], mod_name = "Ensemble_t3", thresholds, stats)
  stats = temp_op[[1]]
  thresholds = temp_op[[2]]
} 

write.csv(cbi, paste0("./Results/Evaluation_metrics_all.csv"), row.names = F)
write.csv(thresholds, paste0("./Results/ModelThresholds_Ensemble.csv"), row.names = F)
write.csv(cbi_thresholds, paste0("./Results/Evaluation_thresholds_CBI.csv"), row.names = F)
write.csv(stats, paste0("./Results/ModelStats_Ensemble.csv"), row.names = F)