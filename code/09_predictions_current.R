## predicting distribution

library(sp)
library(rgdal)
library(raster)
library(randomForest)
library(kernlab)
library(rgl)
library(ks)
library(sm)
library(sf)
library(tidyverse)
library(mapview)
library(caret)
library(rpart)


# upload data -------------------------------------------------------------

## gridded env df with comids
load(file = "ignore/00_RB9_grdded_data.RData") 
head(data_hyd_sf)

## upload nhd shape
nhd <- st_read("/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
## simplify
nhd <- nhd %>%
  st_as_sf %>%
  st_simplify(dTolerance = 0.5, preserveTopology = T)

crs(nhd) 

## define folders
gridFile <- paste0("ignore/ModelResults/Gridded/", models[m],"/")
COMIDFile <- paste0("ignore/ModelResults/COMID/", models[m],"/")

# Predict on all rb9 region-------------------------------------------------------

#Index refers to the right column of probabilities - in this model the second column, which is probs of "1"
all_data <- na.omit(data_hyd_sf)
# all_data <- all_data[,c("COMID", sel.vars)] 
head(all_data)
str(all_data)

## define models
models <- paste0("model",seq(1, 10,1))
m=2

## start loop

for(m in models) {
  ## define folder (model) to pull from
  gridFile <- paste0("ignore/ModelResults/Gridded/", models[m],"/")
  
  ## upload model
  load(file = paste0(gridFile,"rf_model.RData"))
  
  ##  upload presence absence
  load(file=paste0(gridFile, "all_presAbs_env_data.RData"))
  
  
  pred <- predict(rf.final, all_data, filename= paste0(gridFile, "SppProbs_no_clim_gridded.img"), type="prob",  index=2, 
                  na.rm=TRUE, overwrite=TRUE, progress="window")
  
  
  pred_df <- as.data.frame(predict(rf.final, all_data, filename= paste0(gridFile, "SppProbs_no_clim_gridded_df.img"), type="prob",  index=2, 
                                   na.rm=TRUE, overwrite=TRUE, progress="window"))
  ## add comids and cells
  pred_df$cells <- all_data$cells
  pred_df$COMID <- all_data$COMID
  

  
  ## get probability, known occs and env data all together - join by cells and comid
  pred_env <- pred_df %>%
    full_join(all_data, by =  c("cells", "COMID")) %>%
    rename(probOcc = 1) %>%
    full_join(obs, by = c("cells", "COMID")) %>%
    select(probOcc, cells, COMID, Latitude, Longitude, X,Y)

  ## change nas in obs to -999
  # pred_env$PresAbs[is.na(pred_env$PresAbs)] <- -999
  
  write.csv(pred_env, paste0(gridFile,"probOccs_gridded.csv"))
  
  head(pred_env)
}

# Combining probabilities --------------------------------------------------------

## define models
models <- paste0("model",seq(1, 2,1))
# m=2
## empty dataframe
probsx <- NULL

## combine all probabilities
for(m in models) {
  
  ## upload data and format, add model number
  gridFile <- paste0("ignore/ModelResults/Gridded/", m,"/")
  probs <- read.csv(paste0(gridFile,"probOccs_gridded.csv")) %>% dplyr::select(-X.1) %>%
    mutate(Model = m)
  
  ## combine
  probsx <- bind_rows(probsx, probs)

}

### summarise probability over models
probsx_mean <- probsx %>%
 group_by(cells, COMID, Latitude, Longitude, X, Y) %>%
  summarise(MeanProb = mean(probOcc))

head(probsx_mean)

## save out
write.csv(probsx_mean, "ignore/ModelResults/Gridded/09_Av_Probs_Current_RB9.csv")

# Plotting probability of occurrence --------------------------------------

head(nhd)
crs(nhd)

names(probsx_mean)

## get obs and format
obs <- NewDataObsSub %>% as.data.frame() %>%
  dplyr::select(PresAbs, cells, COMID)

obs

nhd_lines_rb9 <- nhd %>%
  dplyr::select(Shape_Leng, COMID) %>%
  filter(COMID %in% pred_df$COMID) %>%
  mutate(COMID = as.integer(COMID))

head(nhd_lines_prob)
dim(nhd_lines_rb9)
str(nhd_lines_rb9)

nhd_lines_prob <- full_join(nhd_lines_rb9, pred_df, by = "COMID")
nhd_lines_prob <- st_zm(nhd_lines_prob)

# plot(nhd_lines_prob[5])

### round prob values - change later
nhd_lines_prob <- nhd_lines_prob %>%
  rename(probOcc = 3) %>%
  mutate(probRound = round(probOcc, digits=1))

## join in obs

nhd_lines_obs <- right_join(nhd_lines_rb9, obs, by = "COMID")
nhd_lines_obs <- st_zm(nhd_lines_obs) %>% mutate(nhd_lines_obs, NewObs = as.factor(NewObs)) %>%
  drop_na(Shape_Leng) %>% st_centroid()

nhd_lines_obs
nhd_lines_prob

library(viridis)

my_breaks <- c(0, 0.25, 0.5, 0.75, 1)

map1 <- ggplot() +
  geom_sf(data = nhd_lines_prob, aes(color = probOcc)) + 
  scale_colour_gradientn(colours=inferno(6), name="Probability of Occurrence",
                         labels = my_breaks) +
  geom_sf(data = subset(nhd_lines_obs, NewObs == 1), shape = 1, size = 2) 


map1

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/01_full_model_prob_occs_map_no_clim.jpg"
ggsave(map1, filename=file.name1, dpi=300, height=5, width=8)

library(mapview)
library(sf)
library(RColorBrewer)
webshot::install_phantomjs()

## map
# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, vector.palette = colorRampPalette(c(  "red", "green")) , fgb = FALSE)


m1 <- mapview(nhd_lines_prob, zcol = "probOcc",  legend = TRUE, layer.name = "Probability of Occurrence") +
  mapview(subset(nhd_lines_obs, NewObs == 1), col.regions = "black",cex = 2, layer.name = "Observations")

m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

mapshot(m1, url = paste0(getwd(), "/ignore/map_no_clim_gridded.html"),
        file = paste0(getwd(), "/ignore/map_no_clim_gridded.png"))
getwd()


