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

\
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
load(file="ignore/00_raw_new_data_coms_nhd.RData")
head(rDFComs)

load(file = "ignore/02_all_data_for_model_gridded.RData") ## data only at pres/abs
head(NewDataObsSub)

# Add hydro (median)---------------------------------------------------------------

delta <- read.csv("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/2022-07-28_predicted_abs_FFM_deltaFFM_SD_COMIDS_medianDelta_test5_norunoff.csv")
head(delta)

comids <- unique(delta$comid)

length(unique(delta$comid))
## compare ffm comids with raster data
Rtest <- rDFComs %>%
  filter(COMID %in% comids)

length(unique(Rtest$COMID)) ## 1907

## compare ffm comids with shape
Stest <- nhdSC %>%
  filter(COMID %in% delta$comid)

length(unique(Stest$COMID)) ## 2116

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

## get median delta H per metric (median of years)
delta_med <- delta_long %>%
  group_by(comid, FlowMetric, hydro.endpoint) %>%
  summarise(MedDelta = median(MetricValue)) %>%
  ungroup() %>%
  rename(COMID = comid) %>%
  dplyr::select(-FlowMetric)  %>%
  pivot_wider(names_from = hydro.endpoint, values_from = MedDelta) %>% dplyr::select(COMID:Wet_BFL_Mag_50)

## join with all data by comid - all RB9
data_hyd_sf <- inner_join(comsDF, delta_med, by = "COMID") ## 209 reaches don't match

length(unique(comsDF$COMID)) ## 3686
length(unique(delta_med$COMID)) ## 2117
length(unique(data_hyd_sf$COMID)) ## 1700??? check this!!


## save out

save(data_hyd_sf, file = "ignore/03_RB9_grdded_data.RData")
# load(file = "ignore/03_all_env_data_gridded_comid.RData")

## join with all data by comid - observations
data_hyd_sf_obs <- inner_join(NewDataObsSub, delta_med, by = "COMID") ## 209 reaches don't match
?inner_join
length(unique(NewDataObsSub$COMID)) ## 330
length(unique(delta_med$COMID)) ## 2117
length(unique(data_hyd_sf_obs$COMID)) ## 90

head(data_hyd_sf_obs)

## save out
save(data_hyd_sf_obs, file = "ignore/03_RB9_grdded_data_observations.RData")


# Add hydro at observations by year -----------------------------------------------
## doesn't work as not all years at each comid

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
