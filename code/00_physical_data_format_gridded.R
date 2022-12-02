### physical data

## data formatting

##workflow
## upload raw data
## get comids to match later
## create rasters
## add hydro
## update remote sensing
## update climate
## 200m grid
## per comid

# packages
library(tidyverse)
library(tidylog)
library(sf)
library(raster)
library(nhdplusTools)
library(readr)
library(mapview)
library(stars)
library(rgdal)
library(rgeos)
library(terra)
library(zoo)
library(terra)


setwd("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")
getwd()

# raw env data ------------------------------------------------------------
path <- "original_model/FullData/"

## Mike's data - env vars with presence absence
data <- read.csv(paste0(path, "200mCells_FullData_PresAbs_Complete_ThinnedCols.csv"))
head(data)

## get coords change crs
orig_grids <- data %>%
  dplyr::select(ID.2, X, Y) %>%
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) %>%
  st_transform(crs=9001)

# Remove climate and remote sensing variables (updated below) -------------

## remove old climate data (also removed remote sensing)

data_red <- data %>%
  dplyr::select(-c(X81pptCr1:FINAL...24..4, PresAbs2005, Presence2))

head(data_red)

## make spatial, crs 4269 to match climate polygons
data_sf<- data_red %>%
  dplyr::select(X,Y, ID.2) %>%
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) 

head(data_sf)
names(data_sf)

st_crs(data_sf)

## upload raster mask - this raster is incorrect, but has the correct crs for extracting
rmask <- raster("ignore/02_mask_raster_network.tif")
crs(rmask)

# Format Climate data -------------------------------------------------

## RB9 boundary


gridsR <- raster(ncol=129, nrow=148, xmn= -2483239, xmx=-2380039, ymn=-4822142 , ymx=-4703742)

projection(gridsR) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
res(gridsR)=c(800,800)

gridsR ## use as template to extract point data
# plot(gridsR)

#### upload annual climate data

## precipitation ###
ppt_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average ppt_30yr.shp") %>%
  rename(ppt_annx = PRISM__) 
class(ppt_ann)
## make spatial and transformCRS
ppt_annSP <- as(ppt_ann, Class = "Spatial")
ppt_annSP <- spTransform(ppt_annSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
ppt_annR <- rasterize(ppt_annSP, gridsR,  'ppt_annx',  na.rm =TRUE, sp = TRUE)
class(ppt_annR)
plot(ppt_annR)

## temperature ###
tmax_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmax_30yr.shp") %>% 
  rename(tmax_annx = PRISM__)
## make spatial and transformCRS
tmax_annSP <- as(tmax_ann, Class = "Spatial")
tmax_annSP <- spTransform(tmax_annSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
tmax_annR <- rasterize(tmax_annSP, gridsR,  'tmax_annx',  na.rm =TRUE, sp = TRUE)

tmin_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmin_30yr.shp") %>% 
  rename(tmin_annx = PRISM__)

## make spatial and transformCRS
tmin_annSP <- as(tmin_ann, Class = "Spatial")
tmin_annSP <- spTransform(tmin_annSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
tmin_annR <- rasterize(tmin_annSP, gridsR,  'tmin_annx',  na.rm =TRUE, sp = TRUE)

## stack all together
annStack <- stack(ppt_annR, tmax_annR, tmin_annR)

## name layers
names(annStack) <- c("ppt_ann", "tmax_ann", "tmin_ann")

## scale to 200m, takes value of original larger cell
annStack200 <- disaggregate(annStack, 4)


## monthly data
ppt_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average ppt_30yr.shp")
## pivot wider by month
ppt_mon <- ppt_mon %>%
  pivot_wider(names_from = "Month", values_from = "pptMonth")

## make spatial and transformCRS
ppt_monSP <- as(ppt_mon, Class = "Spatial")
ppt_monSP <- spTransform(ppt_monSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

#Create list of column names you want to rasterize
fields <- names(ppt_monSP) 
fields
x <- gridsR

## make rasters of each month
for (i in fields){
  x[[i]]<-raster::rasterize(ppt_monSP, x, field=i, na.rm =TRUE, sp = TRUE)
  projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"
  x <- stack(x)
}

x@layers ## check
names(x) <- paste0("pptMon", names(x))
## ppt mon stack
ppt_monR <- x
ppt_monR
# plot(x[[2]]) ## to check

## create raster 
# ppt_monR <- rasterize(ppt_monSP, gridsR,  field = c('pptMonth','Month' ),  na.rm =TRUE, sp = TRUE)


### temp
tmax_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmax_30yr.shp")
## pivot wider by month
tmax_mon <- tmax_mon %>%
  pivot_wider(names_from = "Month", values_from = "tmaxMonth")

## make spatial and transformCRS
tmax_monSP <- as(tmax_mon, Class = "Spatial")
tmax_monSP <- spTransform(tmax_monSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

#Create list of column names you want to rasterize
fields <- names(tmax_monSP) 
fields
x <- gridsR

## make rasters of each env var
for (i in fields){
  x[[i]]<-raster::rasterize(tmax_monSP, x, field=i, na.rm =TRUE, sp = TRUE)
  projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"
  x <- stack(x)
}

x@layers ## check
names(x) <- paste0("tmaxMon", names(x))

## mon stack
tmax_monR <- x
tmax_monR

tmin_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmin_30yr.shp")
## pivot wider by month
tmin_mon <- tmin_mon %>%
  pivot_wider(names_from = "Month", values_from = "tminMonth")

## make spatial and transformCRS
tmin_monSP <- as(tmin_mon, Class = "Spatial")
tmin_monSP <- spTransform(tmin_monSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

##Create list of column names you want to rasterize
fields <- names(tmin_monSP) 
fields
x <- gridsR

## make rasters of each month
for (i in fields){
  x[[i]]<-raster::rasterize(tmin_monSP, x, field=i, na.rm =TRUE, sp = TRUE)
  projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"
  x <- stack(x)
}

x@layers ## check
names(x) <- paste0("tminMon", names(x))
## mon stack
tmin_monR <- x

## stack all together
monStack <- stack(ppt_monR, tmax_monR, tmin_monR)
monStack
## name layers
# names(monStack) <- c("ppt_mon", "tmax_mon", "tmin_mon")

## scale to 200m, takes value of original larger cell
monStack200 <- disaggregate(monStack, 4)
monStack200
## join altogether

climStack <- stack(annStack200, monStack200)
climStack
## save out

writeRaster(climStack, "ignore/00_clim_raster_stack.grd", format="raster", crs="+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", overwrite=TRUE)

## upload
climStack <- stack("ignore/00_clim_raster_stack.grd")

## change crs
projection(climStack) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"

clim_Origx <- cbind(data_red, raster::extract(climStack, orig_grids))
dim(clim_Origx)
names(clim_Origx)[18:56]
head(clim_Origx)
## interpolate the NA values
?na.approx
clim_Origx[18:56] <- na.approx(clim_Origx[18:56])
## 226 missing values changed

head(clim_Origx)

save(clim_Origx, file = "ignore/00_new_clim_orig_env_gridded.Rdata")
load(file = "ignore/00_new_clim_orig_env_gridded.Rdata")

# Remote Sensing data -----------------------------------------------------
## create base raster

sept <- brick("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/TC_2014_RB9/TC_092014_RB9.tif")
sept ## is 30m grids, is extract ok? Mike uses Median (Med) and Variance (Var) within analysis pixel

## aggregate to 180m by median
sept180Med <- aggregate(sept, 6, fun = "median")
names(sept180Med) <- paste0(names(sept180Med), "_Med")

## aggregate to 180m by variance
sept180Var <- aggregate(sept, 6, fun = "var")
names(sept180Var) <- paste0(names(sept180Var), "_Var")

## stack rasters
sept180 <- stack(sept180Med, sept180Var)

## match projection with mask raster
crs(sept180)<-crs(rmask)

## resample to mask layer (200m)
septr <- resample(sept180, rmask, method = "bilinear")
names(septr)

## extract raster values in only 1-3 raster, bind to other data
sept_values <- cbind(clim_Origx, raster::extract(septr[[c(1:3, 7:9)]], orig_grids))
head(sept_values)

## upload wet season data
april <- brick("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/TC_2014_RB9/TC_042014_RB9.tif")
# plot(april)

## aggregate to 180m by median
april180Med <- aggregate(april, 6, fun = "median")
names(april180Med) <- paste0(names(april180Med), "_Med")

## aggregate to 180m by variance
april180Var <- aggregate(april, 6, fun = "var")
names(april180Var) <- paste0(names(april180Var), "_Var")

## stack rasters
april180 <- stack(april180Med, april180Var)

## match projection with mask raster
crs(april180)<-crs(rmask)

## resample to mask layer (200m)
aprilr <- resample(april180, rmask, method = "bilinear")
names(aprilr)

## extract raster values in only 1-3 raster, bind to other data
april_values <- cbind(clim_Orig, raster::extract(aprilr[[c(1:3, 7:9)]], orig_grids))
head(april_values)
## extract raster values in only 1-3 raster, bind to other data
tass_sp <- cbind(sept_values, raster::extract(aprilr[[c(1:3, 7:9)]], orig_grids))
head(tass_sp)

save(tass_sp, file = "ignore/00_tass_cap_climate_original_data.RData")


# Create rasters -----------------------------------------------------------

## upload data
load(file = "ignore/00_tass_cap_climate_original_data.RData")
head(tass_sp)

## format shape file
data_sf <- na.omit(tass_sp) %>%
  dplyr::select(ID.2, X, Y, MRVBF.Mx:TC_042014_RB9.3_Var) #%>%
names(data_sf)
sum(is.na(data_sf))
## make spatial and transformCRS
coordinates(data_sf) <- ~X+Y

projection(data_sf) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"

## create template raster

## dims etc from CurrentGridFeb14.grd

x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)

projection(x) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
crs(x)

#Create list of column names you want to rasterize
fields <- names(data_sf) [2:65]
fields

## make rasters of each env var
for (i in fields){
  x[[i]]<-raster::rasterize(data_sf, x, field=i, na.rm =TRUE, sp = TRUE)
  projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"
  x <- stack(x)
}

x@layers ## check
## define and save layer names
layerNames <- names(x)
save(layerNames, file = "output_data/00_raster_layer_names.RData")

# plot(x[[1]]) ## to check

## save out
writeRaster(x, "ignore/00_raw_new_data_raster.tif", format="GTiff", crs="+proj=geocent +ellps=GRS80 +units=m +no_defs", overwrite=TRUE)

# Add COMIDs - with nhd shapefile -----------------------------------------

## upload nhd shape
nhd <- st_read("/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
## simplify
nhd <- nhd %>%
  st_as_sf %>%
  st_simplify(dTolerance = 0.5, preserveTopology = T)

crs(nhd) 

## convert to points
nhdPoints <- st_cast(nhd, "POINT")
nhdPoints 

# test <- nhdPoints %>%
#   filter(COMID %in% c(22549169, 20348307))

# plot(test)

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, vector.palette = colorRampPalette(c(  "red", "green")) , fgb = FALSE)


m1 <- mapview(nhd, col.regions = "green",cex = 2,layer.name = "nhd reaches (line)") +
  mapview(test, col.regions = "purple",cex = 2,layer.name = "test points")
m1

## looks like same point in the points file is 2 different reaches next to each other.

## upload raster stack made above

x <- stack("ignore/00_raw_new_data_raster.tif")
crs(x) <- crs(rmask)
x
## load names
load( file = "output_data/00_raster_layer_names.RData")

## extract raster values at points
rasterAtPts <- raster::extract(x, nhdPoints, cellnumbers=TRUE)
# rasterAtPts <- na.omit(rasterAtPts)
rasterAtPts
# names(rasterAtPts) <- layerNames
# rasterAtPts
# names(rasterAtPts)
## get coords

coords <- as.data.frame(cbind(nhdPoints, rasterAtPts)) %>% 
  mutate(Longitude = unlist(map(geometry,1)),
         Latitude = unlist(map(geometry,2))) %>%
  dplyr::select(cells, Longitude, Latitude) %>%
  distinct(cells, .keep_all=T)

head(coords)
dim(coords)

rownames(coords) <- seq(1, nrow(coords), 1)

length(unique(coords$cells))
## join with points and remove duplicates

DataComs <- as.data.frame(cbind(nhdPoints, rasterAtPts)) %>% dplyr::select(-geometry) %>% distinct()
head(DataComs)

layerNames
names(DataComs)[7:70] <- layerNames

## get ionly variables and fill NAs
NaApp <- DataComs %>% 
  group_by(COMID) %>%
  dplyr::select(layerNames) %>%
  na.approx() ## this is a quick fix until new data available - some reaches have NAs for all points in reach (all comid)

## get site info
SiteInfo <-  DataComs %>%
  dplyr::select(CLASS, REACHCODE:cells)

## join back together with site info and coords by cell number
DataComsx <- cbind(SiteInfo, NaApp) %>%
  inner_join(coords, by = "cells") 
  
head(DataComsx)

## make spatial and change crs to get data_sf coords

DataComsx <- DataComsx %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=crs(nhdPoints), remove=F) %>%
  st_transform(crs=crs(x)) %>% as.data.frame() %>%
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2))) 
  

head(data_sf)
length(unique(DataComsx$COMID)) ## 2179

names(DataComs)
sum(is.na(DataComsx))

## add in coords from data_sf

# ## upload data
# load(file = "ignore/00_tass_cap_climate_original_data.RData")
# head(tass_sp)
# 
# ## format shape file
# data_sf <- na.omit(tass_sp) %>%
#   dplyr::select(ID.2, X, Y, MRVBF.Mx:TC_042014_RB9.3_Var) #%>%
# names(data_sf)
# 
# DataComs <- cbind(data_sf$X, data_sf$Y, DataComs)

## save out

save(DataComsx, file = "ignore/00_raw_new_data_raster_df_coms.RData")

# Format and join hydro ---------------------------------------------------

delta <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/SD Hydro Vulnerability Assessment - General/Data/RawData/Data_for_SDSU/R/data/Output_09/2022-11-29_predicted_abs_FFMq99_deltaFFMq99_SD_COMIDS_medianDelta_test2q99_test12FFM_allgages.csv")
head(delta)

## change names of columns
delta_long <- delta %>% 
  rename(FlowMetric = metric, MetricValue = abs_FFM_median_cfs, DeltaH = delta_FFM_median_cfs) %>%
  filter(!is.na(FlowMetric)) %>%
  distinct()

unique(delta_long$FlowMetric)

## change names of metrics
delta_long <- delta_long %>%
  # filter(FlowMetric %in% c("fa_mag", "wet_bfl_mag_10", "ds_mag_90", "ds_mag_50")) %>%
  mutate(hydro.endpoint = case_when(FlowMetric == "ds_mag_50" ~ "DS_Mag_50",
                                    FlowMetric == "q99" ~ "Q99",
                                    FlowMetric == "ds_mag_90" ~ "DS_Mag_90",
                                    FlowMetric == "fa_mag" ~ "FA_Mag",
                                    FlowMetric == "peak_10" ~ "Peak_10",
                                    FlowMetric == "peak_2" ~ "Peak_2",
                                    FlowMetric == "peak_5" ~ "Peak_5",
                                    FlowMetric == "sp_mag" ~ "SP_Mag",
                                    FlowMetric == "wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                    FlowMetric == "wet_bfl_mag_50" ~ "Wet_BFL_Mag_50")) 

## values are median, reformat wide
delta_med <- delta_long %>%
  group_by(comid, FlowMetric, hydro.endpoint) %>%
  summarise(MedDelta = MetricValue) %>%
  ungroup() %>%
  rename(COMID = comid) %>%
  dplyr::select(-FlowMetric)  %>%
  pivot_wider(names_from = hydro.endpoint, values_from = MedDelta) %>% dplyr::select(COMID:Wet_BFL_Mag_50)

## join with all data by comid - all RB9
data_hyd_sf <- inner_join(DataComsx, delta_med, by = "COMID") ## 209 reaches don't match

length(unique(DataComsx$COMID)) ## 2179
length(unique(delta_med$COMID)) ## 2116
length(unique(data_hyd_sf$COMID)) ## 2116

head(data_hyd_sf)
dim(data_hyd_sf) ## 16891
## save out

save(data_hyd_sf, file = "ignore/00_RB9_grdded_data.RData")
# load(file = "ignore/03_all_env_data_gridded_comid.RData")
# 
# ## join with all data by comid - observations
# data_hyd_sf_obs <- inner_join(NewDataObsSub, delta_med, by = "COMID") ## 209 reaches don't match
# 
# length(unique(NewDataObsSub$COMID)) ## 322
# length(unique(delta_med$COMID)) ## 2117
# length(unique(data_hyd_sf_obs$COMID)) ## 318
# 
# head(data_hyd_sf_obs)
# 
# ## save out
# save(data_hyd_sf_obs, file = "ignore/03_RB9_grdded_data_observations.RData")
# 
# test <- data_hyd_sf_obs %>%
#   distinct(ID, .keep_all =T)
# 
# length(unique(test$COMID))


# Make raster of all data & hydro -----------------------------------------

## make spatial and transformCRS
coordinates(data_hyd_sf) <- ~X+Y
data_hyd_sf
projection(data_hyd_sf) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"

# data_hyd_sf <- spTransform(data_hyd_sf, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs")) 

## create template raster

## dims etc from CurrentGridFeb14.grd

x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)

projection(x) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
crs(x)

#Create list of column names you want to rasterize
fields <- names(data_hyd_sf) [c(7:70, 76:85)]
fields

## make rasters of each env var
for (i in fields){
  x[[i]]<-raster::rasterize(data_hyd_sf, x, field=i, na.rm =TRUE, sp = TRUE)
  projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"
  x <- stack(x)
}

x@layers ## check
## define and save layer names
layerNames <- names(x)
save(layerNames, file = "output_data/00_final_raster_layer_names.RData")

# plot(x[[1]]) ## to check

## save out
writeRaster(x, "ignore/00_raw_final_data_raster.tif", format="GTiff", crs="+proj=geocent +ellps=GRS80 +units=m +no_defs", overwrite=TRUE)


# Add Comids - not working properly   --------------------------------------------------

## make raster mask
# rmask <- x[[1]]
# crs(rmask) <- "+proj=utm +zone=11 +datum=NAD83"
# crs(x) <-  "+proj=utm +zone=11 +datum=NAD83"
# 
# coms <- raster("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification_PolylineToRaster_v3.tif")
# coms ## 
# names(coms) <- "COMID"
# # coms<-projectRaster(coms, crs = crs(rmask))
# # coms
# 
# 
# ## resample to mask layer
# comsRE <- resample(coms, rmask, method = "ngb")
# comsRE
# 
# comx <- na.omit(as.data.frame(comsRE))
# length(unique(comx$COMID)) ## 2000
# # test <- stack(comsRE, rmask)
# # test <- na.omit(as.data.frame(test, xy=T))
# 
# head(test)
# 
# ## stack with other rasters
# rstack <- stack(x, comsRE)
# plot(rstack)
# 
# ## check and plot
# # rDF <- na.omit(as.data.frame(rstack, xy=T))
# # length(unique(rDF$NHDplus_RB9_csc_Raster))
# # rDFSP <- rDF %>%
# #   st_as_sf(coords=c("x", "y"), crs=26911, remove=F) %>%
# #   filter(!NHDplus_RB9_csc_Raster == 0)
# # 
# # plot(comsRE)
# # plot(rDFSP, add=T)
# 
# ## save raster stack
# writeRaster(rstack, "ignore/00_raw_new_data_raster_Coms_zero.tif", format="GTiff", crs="+proj=geocent +ellps=GRS80 +units=m +no_defs", overwrite=TRUE)
# 
# ## save df
# save(rDF, file="ignore/00_raw_new_data_coms_zero.RData")
# 
# ## for now, get comids from nhdplustools
# 
# rDF <- rDF %>% mutate(ID = 1:nrow(rDF)) 
# 
# orig.sdata.segs <- rDF %>%
#   dplyr::select(x,y, ID) %>%
#   st_as_sf(coords=c("x", "y"), crs=26911, remove=F)%>%
#   st_transform(crs=32611) %>%
#   arrange(ID)
# 
# # use nhdtools to get comids
# data_all_coms <- orig.sdata.segs %>%
#   group_split(ID) %>%
#   set_names(., orig.sdata.segs$ID) %>%
#   map(~discover_nhdplus_id(.x$geometry))
# 
# # flatten into single dataframe instead of list
# data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
#   as.data.frame() %>%
#   rename("COMID"=V1) %>% rownames_to_column(var = "ID") %>%
#   mutate(ID = as.integer(ID))
# data_segs_df
# 
# ## join back to DF
# rDFComs <- full_join(rDF, data_segs_df, by = "ID")
# 
# head(rDFComs)
# 
# 
# save(rDFComs, file="ignore/00_raw_new_data_coms_nhd.RData")
# 
# ## save comids as raster to stack with other env
# 
# rstack <- stack("ignore/00_raw_new_data_raster_Coms_zero.tif")
# names(rstack)
# rstack
# ## layer names
# load(file = "output_data/00_raster_layer_names.RData")
# layerNames
# names(rstack) <- c(layerNames, "COMIDGIS")
# 
# ## env df with comids
# load(file="ignore/00_raw_new_data_coms_nhd.RData")
# head(rDFComs)
# names(rDFComs)
# 
# ## make spatial and transformCRS
# coordinates(rDFComs) <- ~x+y
# 
# projection(rDFComs) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
# ## make raster and add to stack
# 
# x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)
# projection(x) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
# crs(x)
# 
# #Create list of column names you want to rasterize
# fields <- names(rDFComs) [c(1:66)]
# fields
# 
# 
# ## make rasters of each env var
# for (i in fields){
#   x[[i]]<-raster::rasterize(rDFComs, x, field=i, na.rm =TRUE, sp = TRUE)
#   projection(x)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"
#   x <- stack(x)
# }
# 
# 
# layerNames <- names(x)
# layerNames
# x
# 
# save(layerNames, file = "output_data/00_raster_layer_names.RData")
# 
# writeRaster(x, "ignore/00_raw_new_data_raster_Coms_zero.tif", format="GTiff", crs="+proj=geocent +ellps=GRS80 +units=m +no_defs", overwrite=TRUE)
# # 
# 
# # stuff not working -------------------------------------------------------
# ## upload nhd points
# nhdPts <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_points/NHD_reaches_RB9_points.shp")
# dim(nhdPts)
# 
# crs(x) <- "+proj=utm +zone=11 +datum=NAD83"
# ?extract
# test <- raster::extract(x, nhdPts)
# dim(test)
# dim(na.omit(test))
# ## upload raster mask
# rmask <- raster("ignore/02_mask_raster_network.tif")
# crs(rmask)
# 
# 
# ## convert to df
# r_df <- as.data.frame(rmask, xy=T)
# r_df <- na.omit(r_df)
# 
# ## make spatial
# r_dfSP <- r_df %>%
#   st_as_sf(coords=c("x", "y"), crs=4269, remove=F) 
# 
# head(r_dfSP)
# dim(r_dfSP)
# 
# coms <- raster("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHDplus_RB9_Raster_v2/NHDplus_RB9_csc_Raster.tif")
# coms ## no values are zeros, slightly different res and ncells
# plot(coms)
# ## resample to mask layer
# comsRE <- resample(coms, rmask, method = "ngb")
# comsRE
# 
# test <- stack(comsRE, rmask)
# test
# test <- na.omit(test)
# test
# ## make dataframe
# comsDF <- as.data.frame(coms, xy=T) 
# ## change name of comid raster
# names(comsDF)[3] <-"COMID"
# 
# ## some comids are 0, replacewith NA
# comsDF$COMID[comsDF$COMID == 0] <- NA
# comsDF <- na.omit(comsDF)
# length(unique(comsDF$COMID)) ## 1993
# 
# dim(comsDF)
# head(comsDF)
