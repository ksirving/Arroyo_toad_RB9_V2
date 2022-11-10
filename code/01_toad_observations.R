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
  mutate(ID = paste0("U", ID)) 

dim(bio_sub) ## 2164
head(bio_sub)

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

head(bio1_sub)


# Observations from Chad --------------------------------------------------

chad1 <- read.csv("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/RawData/Toad_Data/FromChad/MOMsAL.csv.csv")
chad2 <- read.csv("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/RawData/Toad_Data/FromChad/SO_arroyo_toad.csv.csv")
chad3 <- read.csv("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/RawData/Toad_Data/FromChad/USFWS_R9Only.kml.csv")

head(chad1)

chad1_sub <- chad1 %>% 
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  rename(Year = SurYr, SurveyDate = SurDate, Species_1 = CName) %>% 
  filter(Year %in% 1990:2014) %>%
  mutate(Age = NA) %>%
  dplyr::select(SurveyDate, Year, Latitude, Longitude, Species_1, Age) %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F) %>%
  mutate(ID = 1:length(geometry))  %>%
    mutate(ID = paste0("C", ID)) 

chad2_sub <- chad2 %>% 
  filter(!is.na(Lat), !is.na(Long)) %>%
  rename(Latitude = Lat, Longitude = Long) %>%
  separate(SurveyDate, into = c("Month", "Day", "Year"), remove = F, sep="/") %>% 
  mutate(Year = as.integer(Year)) %>%
  filter(Year %in% 1990:2014)  %>%
  dplyr::select(SurveyDate, Year, Latitude, Longitude, Species_1, Age) %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("Ct", ID)) 


chad3_sub <- chad3 %>% ### OCC here might be abundance, not sure though!
  filter(!is.na(Y), !is.na(X)) %>%
  rename(Latitude = Y, Longitude = X) %>%
  rename(Year = YEAR, SurveyDate = DATE_, Species_1 = CNAME) %>% 
  mutate(SurveyDate = as.character(SurveyDate)) %>%
  filter(Year %in% 1990:2014) %>%
  mutate(Age = NA) %>%
  dplyr::select( SurveyDate, Year, Latitude, Longitude, Species_1, Age) %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("Cth", ID)) 

## join altogether

chad_sub <- bind_rows(chad1_sub, chad2_sub, chad3_sub)

## remove other species

unique(chad_sub$Species_1)

chad_sub <- chad_sub %>%
  filter(Species_1 %in% c("arroyo toad", "Arroyo toad", "Arroyo_toad"))

length(unique(chad_sub$ID))

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

save(bio_data2, file = "ignore/01_al_bio_data_COMIDs.RData")
load(file = "ignore/01_al_bio_data_COMIDs.RData") # bio_data2

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

save(bio_data3, file = "ignore/01_al_bio_data_COMIDs_SDRB9.RData")
load(file = "ignore/01_al_bio_data_COMIDs_SDRB9.RData") # bio_data3

## chad obs

# Create dataframe for looking up COMIDS (here use all stations)
data_segs <- chad_sub %>%
  dplyr::select(ID, Latitude, Longitude) %>%
  distinct(ID, Latitude, Longitude) %>% 
  st_as_sf(coords=c("Longitude",  "Latitude"), crs=4326, remove=F) %>%
  st_transform(crs=32611) %>%
  arrange(ID)

crs(data_segs)
head(data_segs)
dim(data_segs)

length(unique(data_segs$ID))

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

bio_data4 <- full_join(chad_sub, data_segs_df, by = "ID")
object.size(bio_data2)

length(unique(bio_data4$COMID)) ## 331
head(bio_data4)

save(bio_data4, file = "ignore/01_chad_bio_data_COMIDs.RData")
load(file = "ignore/01_chad_bio_data_COMIDs.RData") # bio_data4



# add Mike's observations -------------------------------------------------

## add in Mikes observations, remove any that are duplicated

# original pres/abs
inshape="original_model/Current/randomForests/PARTITIONING/DATA3/200mCells_PresAbs2005 PCA PresAbs Pts_MLT.shp" ## presence absence
## format
all_data_obs <- st_read(inshape) %>%
  mutate(Longitude = unlist(map(geometry,1)),
         Latitude = unlist(map(geometry,2))) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("M", ID)) ## p/a

head(all_data_obs)

# Create dataframe for looking up COMIDS (here use all stations)
data_segs <- all_data_obs %>%
  dplyr::select(ID, Latitude,Longitude) %>%
  distinct(ID, Latitude,Longitude) %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=32611, remove=F) %>%
  # st_transform(crs=32611) %>%
  arrange(ID)

# use nhdtools to get comids
data_all_coms <- data_segs %>%
  group_split(ID) %>%
  set_names(., data_segs$ID) %>%
  map(~discover_nhdplus_id(.x$geometry))

# flatten into single dataframe instead of list
data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
  as.data.frame() %>%
  rename("COMID"=V1) %>% rownames_to_column(var = "ID") 


bio_data0 <- full_join(all_data_obs, data_segs_df, by = "ID")
object.size(bio_data0)

length(unique(bio_data0$COMID)) ## 218

save(bio_data0, file = "ignore/01_original_bio_data_COMIDs.RData")
load(file = "ignore/01_original_bio_data_COMIDs.RData")

# Join all presence/absence together (points) ---------------------------------------------

head(bio_data0)

## format: selected columns, make life stage counts long, 
unique(bio$Count)
unique(bio$Year)

head(bio_data2)
## USGS
bio <- bio_data2 %>%
  dplyr::select(BUMIcount, ID, COMID, StartLat, StartLong, Year) %>%
  rename(Latitude = StartLat, Longitude = StartLong) %>%
  pivot_longer(BUMIcount, names_to = "LifeStage", values_to = "Count") %>%
  mutate(PresAbs = ifelse(Count < 1, 0, 1), Year = paste0("20", Year)) %>%
  filter(Year >= 2000) %>%
  distinct()

head(bio)
unique(bio2$LifeStage)

## other
bio1 <- bio_data3 %>%
  dplyr::select(COMID, Latitude, Longitude, ID, Year) %>%
  mutate(LifeStage = "ALL", PresAbs = 1) %>%
  filter(Year >= 2000) %>%
  distinct()

## Chad
bio2 <- bio_data4 %>%
  dplyr::select(COMID, Latitude, Longitude, ID, Age, Year) %>%
  rename(LifeStage = Age) %>%
  mutate(PresAbs = 1, Year = as.character(Year)) %>%
  filter(Year >= 2000) %>%
  distinct()

bio3 <- bio_data0 %>%
  dplyr::select(COMID, Latitude, Longitude, ID, PresAbs200) %>%
  mutate( Year = NA,LifeStage = "ALL" ) %>%
  rename(PresAbs = PresAbs200) %>%
  distinct()
  
## join together

bio_points <- bind_rows(bio, bio1, bio2, bio3) %>% distinct()

dim(bio_points) ## 5130 - unique points
length(unique(bio_points$COMID)) ## 534 - unique comid
# subs <- bio_points %>% filter(Year >= 2000) 
# length(unique(subs$COMID)) ## 485 2000-2014

## map points

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

## plot points
m1 <- mapview(bio_points, cex=6, col.regions="orange",
              layer.name="Toad Observations") #+
  # mapview(subs, cex=6, col.regions="blue",
  #         layer.name="Toad Observations > 2000")


m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

## save out

st_write(bio_points, "ignore/01_toad_obs_points_RB9.shp", append=FALSE)

