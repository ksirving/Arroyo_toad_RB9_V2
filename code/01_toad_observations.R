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


# Join all presence/absence together (points) ---------------------------------------------

head(bio_data2)

## format: selected columns, make life stage counts long, 
unique(bio$Count)
unique(bio$Year)
## USGS
bio <- bio_data2 %>%
  dplyr::select(BUMIcount, ID, COMID, StartLat, StartLong, Year) %>%
  rename(Latitude = StartLat, Longitude = StartLong) %>%
  pivot_longer(BUMIcount, names_to = "LifeStage", values_to = "Count") %>%
  mutate(PresAbs = ifelse(Count < 1, 0, 1), Year = paste0("20", Year)) %>%
  distinct()

## other
bio1 <- bio_data3 %>%
  dplyr::select(COMID, Latitude, Longitude, ID, Year) %>%
  mutate(LifeStage = "ALL", PresAbs = 1) %>%
  distinct()

## Chad
bio2 <- bio_data4 %>%
  dplyr::select(COMID, Latitude, Longitude, ID, Age, Year) %>%
  rename(LifeStage = Age) %>%
  mutate(PresAbs = 1, Year = as.character(Year)) %>%
  distinct()


## join together

bio_points <- bind_rows(bio, bio1, bio2) %>% distinct()

dim(bio_points) ## 4558 - unique points
length(unique(bio_points$COMID)) ## 541 - unique comid
subs <- bio_points %>% filter(Year >= 2000) 
length(unique(subs$COMID)) ## 485 2000-2014

## map points

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

## plot points
m1 <- mapview(bio_points, cex=6, col.regions="orange",
              layer.name="Toad Observations") +
mapview(subs, cex=6, col.regions="blue",
        layer.name="Toad Observations > 2000")


m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

## save out

st_write(bio_points, "ignore/01_toad_obs_points_RB9.shp", append=FALSE)

# Join all presence/absence together (COMID) ------------------------------

## remove pre 2000 if needed**

## make sure if reach has at least one presence then keeps it, if none then absence
bio_coms <- bio %>%
  as.data.frame() %>%
  group_by(COMID, LifeStage) %>% 
  summarise(PresAbs = max(PresAbs)) %>%
  filter(LifeStage == "BUMIcount") %>%
  drop_na() %>%
  mutate(LifeStage = "ALL")

head(bio_coms)
dim(bio_coms) ## 318

unique(bio$PresAbs)

## format other observations
head(bio_data3)

bio1_coms <- bio_data3 %>%
  as.data.frame() %>%
  dplyr::select(COMID) %>%
  mutate(LifeStage = "ALL", PresAbs = 1) %>%
  distinct(COMID, .keep_all = T)
?distinct
bio1_coms
bio

head(bio_data4)
## Chad data
bio2_coms <- bio_data4 %>%
  as.data.frame() %>%
  dplyr::select(COMID, Age) %>%
  rename(LifeStage = Age) %>%
  mutate(PresAbs = 1) %>%
  distinct(COMID, .keep_all = T)

# join obs together

bioAll <- bind_rows(bio_coms, bio1_coms, bio2_coms) %>% distinct()
bioAll

sum(bioAll$PresAbs ==1)
sum(bioAll$PresAbs ==0)

# plot comids -------------------------------------------------------------

nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHDplus_RB9.shp")
# nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/SD_RB9_boundary.shp")

head(nhd)
obs <- bioAll %>%
  dplyr::select(COMID)

nhd_lines_rb9 <- nhd %>%
  dplyr::select(SHAPE_LENG, COMID) %>%
  # filter(COMID %in% all_data$COMID) %>% ## subset to orginial data with RB9 extent
  mutate(COMID = as.integer(COMID))

nhd_lines_rb9

dim(nhd_lines_rb9)
str(nhd_lines_rb9)

# nhd_lines_prob <- full_join(nhd_lines_rb9, pred_df, by = "COMID")
nhd_lines_rb9 <- st_zm(nhd_lines_rb9)

## join in obs

nhd_lines_obs <- right_join(nhd_lines_rb9, bioAll, by = "COMID")
nhd_lines_obs <- st_zm(nhd_lines_obs) %>% mutate(nhd_lines_obs, PresAbs = as.factor(PresAbs)) #%>%
drop_na(Shape_Leng)
dim(nhd_lines_obs)


# library(viridis)
map1 <- ggplot() +
  geom_sf(data = nhd_lines_rb9) +
  geom_sf(data = nhd_lines_obs, aes(colour =PresAbs))
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

