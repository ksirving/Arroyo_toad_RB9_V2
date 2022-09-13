### toad observations

### new observation data

# packages
library(tidyverse)
library(tidylog)
library(sf)
library(raster)
library(nhdplusTools)
library(readr)
library(mapview)
library(RStoolbox)


setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")
getwd()

# original data ------------------------------------------------------------
path <- "ignore/FullData/"

## raw data - env vars with presence absence
load(file = "ignore/00_all_env_bio_data_NHD_reach.RData") ## all_data - _abs_hydro
head(all_data)
length(unique(all_data$COMID))
# data <- read.csv(paste0(path, "200mCells_FullData_PresAbs_Complete_ThinnedCols.csv"))
# head(data)
# str(data$ID.2)


# USGS observation data ----------------------------------------------------

bio <- read.csv("input_data/SCCWRP_BUMI_20220816.csv")
head(bio)
dim(bio)
# bio <- bio %>% filter(!Date1 %in% complete.cases(bio$Date1)) ## remove rows with no date

unique(bio$SurveyType) ## includes turtle, remove

sum(bio$SurveyMethod == "Turtle: Visual")

unique(bio$Date1)

## remove rows with no coords
## take only years until 2014
## make dummy ID variable
## remove tuttle from SurveyMethod

bio_sub <- bio %>% 
  filter(!is.na(StartLat), !is.na(StartLong)) %>%
  separate(Date1, into = c("Month", "Day", "Year"), remove = F, sep=" |/") %>% 
  filter(!Year > 14, !SurveyMethod %in% c("Turtle: Visual", "Turtle: Trapping")) %>%
  st_as_sf(coords=c( "StartLong", "StartLat"), crs=4326, remove=F) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("B", ID)) 

dim(bio_sub) ## 2164

head(bio_sub)
range(bio_sub$Year)


## look at sites

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

## plot points
m1 <- mapview(bio_sub, cex=6, col.regions="orange",
              layer.name="Toad Observations USGS")# +
# mapview(bio_sub, cex=6, col.regions="blue",
#         layer.name="Toad Observations Other")  


m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")


# More observations data --------------------------------------------------

bio1 <- read.csv("input_data/arroyo_toad_SDRB9_occurrence_historicaldata_othersources.csv")
# bio2 <- read.csv("input_data/SDRB9_arroyo_toad_CarlsbadFWO_1991_2020.csv")
# dim(bio2)
# bio2
## remove rows with no coords
## take only years between 1990 -  2014
## make dummy ID variable

bio1_sub <- bio1 %>% 
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  separate(EventDate, into = c("Month", "Day", "Year"), remove = F, sep="/") %>% 
  filter(Year %in% 1990:2014) %>%
  # dplyr::select(Latitude:Longitude)
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("H", ID)) 

# get COMIDs --------------------------------------------------------------

# Create dataframe for looking up COMIDS (here use all stations)
data_segs <- bio_sub %>%
  dplyr::select(ID, StartLat, StartLong) %>%
  distinct(ID, StartLat, StartLong) %>% 
  st_as_sf(coords=c("StartLong", "StartLat"), crs=4326, remove=F) %>%
  st_transform(crs=32611) %>%
  arrange(ID)

crs(data_segs)
head(data_segs)
dim(data_segs)

# use nhdtools to get comids
data_all_coms <- data_segs %>%
  group_split(ID) %>%
  set_names(., data_segs$ID) %>%
  map(~discover_nhdplus_id(.x$geometry))

data_all_coms

# flatten into single dataframe instead of list
data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
  as.data.frame() %>%
  rename("COMID"=V1) %>% rownames_to_column(var = "ID") 

head(data_segs_df)

bio_data2 <- full_join(bio_sub, data_segs_df, by = "ID")
object.size(bio_data2)

length(unique(bio_data2$COMID)) ## 331
head(bio_data2)

save(bio_data2, file = "ignore/00_al_bio_data_COMIDs.RData")
load(file = "ignore/00_al_bio_data_COMIDs.RData") # bio_data2

## other obs

# Create dataframe for looking up COMIDS (here use all stations)
data_segs <- bio1_sub %>%
  dplyr::select(ID, Latitude, Longitude) %>%
  distinct(ID, Latitude, Longitude) %>%
  st_as_sf(coords=c("Longitude",  "Latitude"), crs=4269, remove=F) %>%
  st_transform(crs=32611) %>%
  arrange(ID)

head(data_segs)

# # use nhdtools to get comids
data_all_coms <- data_segs %>%
  group_split(ID) %>%
  set_names(., data_segs$ID) %>%
  map(~discover_nhdplus_id(.x$geometry))

data_all_coms

# # flatten into single dataframe instead of list
data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
  as.data.frame() %>%
  rename("COMID"=V1) %>% rownames_to_column(var = "ID")

head(data_segs_df)

bio_data3 <- full_join(bio1_sub, data_segs_df, by = "ID")
object.size(bio_data3)

length(unique(bio_data3$COMID)) ## 19

save(bio_data3, file = "ignore/00_al_bio_data_COMIDs_SDRB9.RData")
load(file = "ignore/00_al_bio_data_COMIDs_SDRB9.RData") # bio_data3

# Format presence absence etc ---------------------------------------------

head(bio_data2)

## format: selected columns, make life stage counts long, 
unique(bio$Count)

bio <- bio_data2 %>%
  dplyr::select(BUMIcount, ID, COMID) %>%
  pivot_longer(BUMIcount, names_to = "LifeStage", values_to = "Count") %>%
  mutate(PresAbs = ifelse(Count < 1, 0, 1))

## make sure if reach has at least one presence then keeps it, if none then absence
bio <- bio %>%
  as.data.frame() %>%
  group_by(COMID, LifeStage) %>% 
  summarise(PresAbs = max(PresAbs)) %>%
  filter(LifeStage == "BUMIcount") %>%
  drop_na() %>%
  mutate(LifeStage = "ALL")

head(bio)
dim(bio) ## 318

unique(bio$PresAbs)

## format other observations
head(bio_data3)

bio1 <- bio_data3 %>%
  as.data.frame() %>%
  dplyr::select(COMID) %>%
  mutate(LifeStage = "ALL", PresAbs = 1) %>%
  distinct()

bio1
bio

# join obs together

bioAll <- bind_rows(bio, bio1) %>% distinct()
bioAll

sum(bioAll$PresAbs ==1)
sum(bioAll$PresAbs ==0)

# plot comids -------------------------------------------------------------

nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHDplus_RB9.shp")
# nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/SD_RB9_boundary.shp")


obs <- DataObs %>%
  dplyr::select(COMID)

nhd_lines_rb9 <- nhd %>%
  dplyr::select(SHAPE_LENG, COMID) %>%
  filter(COMID %in% all_data$COMID) %>% ## subset to orginial data with RB9 extent
  mutate(COMID = as.integer(COMID))



nhd_lines_rb9

dim(nhd_lines_rb9)
str(nhd_lines_rb9)

# nhd_lines_prob <- full_join(nhd_lines_rb9, pred_df, by = "COMID")
nhd_lines_rb9 <- st_zm(nhd_lines_rb9)


## join in obs

nhd_lines_obs <- right_join(nhd_lines_rb9, DataObs, by = "COMID")
nhd_lines_obs <- st_zm(nhd_lines_obs) %>% mutate(nhd_lines_obs, NewObs = as.factor(NewObs)) #%>%
drop_na(Shape_Leng)
dim(nhd_lines_obs)


# library(viridis)
map1 <- ggplot() +
  geom_sf(data = nhd_lines_rb9) +
  geom_sf(data = nhd_lines_obs, aes(colour = NewObs))
# scale_fill_gradientn(colours=rev(magma(6))) ## colours not working

map1

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/01_rb9_map.jpg"
# ggsave(map1, filename=file.name1, dpi=300, height=5, width=6)



# add Mike's observations -------------------------------------------------

## this is a quick fix for now, will need to be fixed

## bio comids
load(file= "ignore/01_all_env_bio_data_NHD_reach_original.RData") ## all_data_obs
head(all_data_obs)
names(all_data_obs)

orig_obs <- all_data_obs %>%
  dplyr::select(COMID, PresAbs200)

head(orig_obs)

unique(orig_obs$COMID)

## join with new obs

DataObs <- full_join(orig_obs, bioAll, by = "COMID") %>% dplyr::select(-LifeStage) %>%
  mutate(NewObs = ifelse(is.na(PresAbs), PresAbs200, PresAbs)) %>% ## use new presences, if NA add old presence/absences 
  dplyr::select(COMID, NewObs)

head(DataObs)

length(unique(DataObs$COMID)) ## 550

