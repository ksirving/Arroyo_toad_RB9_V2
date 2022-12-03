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
library(tidylog)


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
  summarise(MeanProb = mean(na.omit(probOcc)))

head(probsx_mean)

## save out
write.csv(probsx_mean, "ignore/ModelResults/Gridded/09_Av_Probs_Current_RB9.csv")

# Plotting probability of occurrence --------------------------------------

names(probsx_mean)
sum(is.na(probsx_mean$Latitude))
ind <- which(is.na(probsx_mean$Latitude))
probsx_mean[ind,] ### fix these nas, might be the same issue as above, remove for now
probsx_mean <- na.omit(probsx_mean)

## make probabilities spatial

probs_sf <- probsx_mean %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=crs(nhd), remove=F) 

## save out
st_write(probs_sf, "ignore/ModelResults/Gridded/09_Av_Probs_Current_RB9.shp")

## make raster of probs

## make spatial and transformCRS
coordinates(probsx_mean) <- ~X+Y
projection(probsx_mean) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"

## create template raster

## dims etc from CurrentGridFeb14.grd

x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)

projection(x) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
crs(x)

## make rasters of probability

  x<-raster::rasterize(probsx_mean, x, field="MeanProb", na.rm =TRUE, sp = TRUE)
  projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"

  ## save
writeRaster(x, "ignore/ModelResults/Gridded/09_mean_probs_RB9.tif", format="GTiff", crs="+proj=geocent +ellps=GRS80 +units=m +no_defs", overwrite=TRUE)
  
## get obs and format. makespatial 
obs <- NewDataObsSub %>% as.data.frame() %>%
  dplyr::select(PresAbs, cells, COMID) %>%
  inner_join(probs_sf, by = c("cells", "COMID")) %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=crs(nhd), remove=F) 

## map
library(viridis)

probs_sf

my_breaks <- c(0, 0.25, 0.5, 0.75, 1)

map1 <- ggplot() +
  geom_sf(data = probs_sf, aes(color = MeanProb)) + 
  scale_colour_gradientn(colours=inferno(6), name="Probability of Occurrence",
                         labels = my_breaks) +
  geom_sf(data = subset(obs, PresAbs == 1), shape = 1, size = 2) 


map1

file.name1 <- "ignore/ModelResults/Gridded/09_prob_occs_map_gridded.jpg"
ggsave(map1, filename=file.name1, dpi=300, height=5, width=8)

library(mapview)
library(sf)
library(RColorBrewer)
# webshot::install_phantomjs()

## map
# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, vector.palette = colorRampPalette(c(  "red", "green")) , fgb = FALSE)


m1 <- mapview(probs_sf, zcol = "MeanProb",  legend = TRUE, layer.name = "Probability of Occurrence") +
  mapview(subset(obs, PresAbs == 1), col.regions = "black",cex = 2, layer.name = "Observations")

m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

mapshot(m1, url = paste0(getwd(), "/ignore/ModelResults/Gridded/map_no_clim_gridded.html"),
        file = paste0(getwd(), "/ignore/ModelResults/Gridded/map_no_clim_gridded.png"))
getwd()
