## add hydro

## workflow
## add hydro median
## add hydro lag year of presence


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

## upload data

# testing data ------------------------------------------------------------

## upload nhd points
# nhdPts <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_points/NHD_reaches_RB9_points.shp")
# nhdPts

## read in dfcreated in 00 - RB9 data
load(file="ignore/00_raw_new_data_coms_nhd.RData")
head(rDFComs)

load(file = "ignore/02_all_data_for_model_gridded.RData") ## data only at pres/abs
head(NewDataObsSub)

## upload comid raster and mask
rmask <- raster("ignore/02_mask_raster_network.tif")
crs(rmask)

coms <- raster("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHDplus_RB9_Raster_v2/NHDplus_RB9_csc_Raster.tif")
coms 

## resample to mask layer
comsRE <- resample(coms, rmask, method = "ngb")
crs(comsRE)
comsRE

comsREDF <- as.data.frame(comsRE, xy=T) %>% rename(comid = NHDplus_RB9_csc_Raster)
head(comsREDF)

## upload original nhd comids
nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHDplus_RB9.shp")
head(nhd)

nhdSC <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
head(nhdSC)

length(unique(rDFComs$COMID)) ## 3841
length(unique(nhdSC$COMID)) ## 2178
length(unique(nhd$COMID)) ## 4522
length(unique(comsREDF$comid)) ## 4261  
length(unique(comsDF$COMID)) ## 3686
length(unique(nhdPts$COMID)) ## 2178
plot(comsRE)

## compare env data with shape
Etest <- nhdPts %>%
  filter(COMID %in% comsDF$COMID) ## 1758



# Upload all data ---------------------------------------------------------

## read in dfcreated in 00 - RB9 data
load(file="ignore/00_raw_new_data_raster_df_coms.RData")
head(DataComs)

load(file = "ignore/02_all_data_for_model_gridded.RData") ## data only at pres/abs
# head(NewDataObsSub)

# Upload and match FFM---------------------------------------------------------------

delta <- read.csv("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/2022-07-28_predicted_abs_FFM_deltaFFM_SD_COMIDS_medianDelta_test5_norunoff.csv")
head(delta)

comids <- unique(delta$comid)

length(unique(delta$comid)) ## 2117
## compare ffm comids with raster data
# Rtest <- DataComs %>%
#   filter(COMID %in% comids)
# 
# length(unique(Rtest$COMID)) ## 2048
# 
# ## which are missing?
# 
# Missing <- DataComs %>%
#   filter(!COMID %in% comids)
# 
# miscoms <- unique(Missing$COMID)
# 
# ## save list of missing comids
# write.csv(miscoms, "output_data/02_missing_comids.csv")
# miscoms
# plot(Missing)
# 
# ## nhd lines
# 
# ## upload nhd shape
# nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
# ## simplify
# nhd <- nhd %>%
#   st_as_sf %>%
#   st_simplify(dTolerance = 0.5, preserveTopology = T)
# 
# 
# ## get line string of missing coms
# 
# nhdMissing <- nhd %>%
#   filter(COMID %in% miscoms)
# 
# nhdMissing
# ## save for Abel
# st_write(nhdMissing, "output_data/03_missing_coms_multiline.shp", append = F)

# Format and join hydro ---------------------------------------------------

## change names of columns
delta_long <- delta %>% 
  rename(FlowMetric = metric, MetricValue = abs_FFM_median_cfs, DeltaH = delta_FFM_median_cfs) %>%
  distinct()

## change names of metrics
delta_long <- delta_long %>%
  # filter(FlowMetric %in% c("fa_mag", "wet_bfl_mag_10", "ds_mag_90", "ds_mag_50")) %>%
  mutate(hydro.endpoint = case_when(FlowMetric == "ds_mag_50" ~ "DS_Mag_50",
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
  summarise(MedDelta = DeltaH) %>%
  ungroup() %>%
  rename(COMID = comid) %>%
  dplyr::select(-FlowMetric)  %>%
  pivot_wider(names_from = hydro.endpoint, values_from = MedDelta) %>% dplyr::select(COMID:Wet_BFL_Mag_50)

## join with all data by comid - all RB9
data_hyd_sf <- inner_join(DataComs, delta_med, by = "COMID") ## 209 reaches don't match

length(unique(DataComs$COMID)) ## 2107
length(unique(delta_med$COMID)) ## 2117
length(unique(data_hyd_sf$COMID)) ## 2048

head(data_hyd_sf)
dim(data_hyd_sf) ## 15531
## save out

save(data_hyd_sf, file = "ignore/03_RB9_grdded_data.RData")
# load(file = "ignore/03_all_env_data_gridded_comid.RData")
# 
# ## join with all data by comid - observations
data_hyd_sf_obs <- inner_join(NewDataObsSub, delta_med, by = "COMID") ## 209 reaches don't match

length(unique(NewDataObsSub$COMID)) ## 322
length(unique(delta_med$COMID)) ## 2117
length(unique(data_hyd_sf_obs$COMID)) ## 318

head(data_hyd_sf_obs)

## save out
save(data_hyd_sf_obs, file = "ignore/03_RB9_grdded_data_observations.RData")


# Plot comids -------------------------------------------------------------

head(data_hyd_sf_obs)
head(NewDataObsSub)

## make delta data spatial

## upload nhd reaches
nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
crs(nhd)
head(nhd)

test <- nhd %>%
  filter(COMID %in% c(22549169, 20348307))

test$geometry

nhd_lines_delta <- nhd %>%
  dplyr::select(CLASS, COMID) %>%
  filter(COMID %in% delta_med$COMID) %>%
  mutate(COMID = as.integer(COMID))

length(unique(nhd_lines_delta$COMID)) ## 2116
names(NewDataObsSub)
names(data_hyd_sf)

## comids matched ffm with all data
data_hyd_sfWG <- nhd %>%
  dplyr::select(CLASS, COMID) %>%
  filter(COMID %in% data_hyd_sf$COMID) %>%
  mutate(COMID = as.integer(COMID))


## observations from all data
NewDataObsSubWG <- NewDataObsSub %>%
  st_transform(crs= crs(nhd_lines_delta)) %>%
  dplyr::select(ID:TC_042014_RB9.3_Var, MRVBF.Mx,pptMonX2, COMID)

## observations with hydro
data_hyd_sf_obsWG <- data_hyd_sf_obs %>%
  st_transform(crs= crs(nhd_lines_delta)) %>%
  dplyr::select(ID:TC_042014_RB9.3_Var, MRVBF.Mx,pptMonX2,DS_Mag_50, COMID)

## map
# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, vector.palette = colorRampPalette(c(  "red", "green")) , fgb = FALSE)


m1 <- mapview(nhd_lines_delta, col.regions = "green",cex = 2,layer.name = "FFM data COMIDs") +
  mapview(data_hyd_sfWG, col.regions = "purple",cex = 2,layer.name = "Matched COMIDs") +
  mapview(NewDataObsSubWG, col.regions = "pink",cex = 3, layer.name = "Observations with all data") +
  mapview(data_hyd_sf_obsWG , col.regions = "red",cex = 3, layer.name = "Observations with FFM")

m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

# mapshot(m1, url = paste0(getwd(), "/ignore/01_full_model_comid_prob_occs_mapview.html"),
#         file = paste0(getwd(), "/ignore/01_full_model_comid_prob_occs_mapview.png"))
getwd()

test <- filter(data_hyd_sf_obsWG, COMID == 	22547077)
test

test2 <- filter(data_hyd_sf_obsWG, LifeStage == "Pseudo")
test2 ## some values missing from env data
length(unique(test2$COMID))

# Add hydro at observations by year -----------------------------------------------
## doesn't work as already median 

head(delta_wide)
head(NewDataObsSub)
str(delta_long)
str(NewDataObsSub)

sort(unique(delta_long$wayr))

delta_wide <- delta_long %>% mutate(wayr = as.character(wayr)) %>%
  dplyr::select(-FlowMetric, - DeltaH) %>%
  pivot_wider(names_from = hydro.endpoint, values_from = MetricValue)

test <- inner_join(NewDataObsSub, delta_long, by = c("COMID" = "comid", "Year" = "wayr"))
test
