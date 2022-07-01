### new observation data

# packages
library(tidyverse)
library(tidylog)
library(sf)
library(raster)
library(nhdplusTools)
library(readr)
library(mapview)


setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9")
getwd()

# original data ------------------------------------------------------------
path <- "ignore/FullData/"

## raw data - env vars with presence absence
data <- read.csv(paste0(path, "200mCells_FullData_PresAbs_Complete_ThinnedCols.csv"))
head(data)
str(data$ID.2)


# New observation data ----------------------------------------------------

bio <- read.csv("input_data/SCCWRP_BUMI_20220624.csv")
head(bio)
bio <- bio[-c(1562:1565),] ## remove rows with no date


## remove rows with no coords
## take only years until 2014
## make dummy ID variable

bio_sub <- bio %>% 
  filter(!is.na(StartLat), !is.na(StartLong)) %>%
  separate(Date1, into = c("Month", "Day", "Year", "Time"), remove = F, sep=" |/") %>% 
  filter(!Year > 2014) %>%
  st_as_sf(coords=c( "StartLong", "StartLat"), crs=4326, remove=F) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("B", ID)) 

dim(bio_sub)
  
  head(bio_sub)
sum(is.na(bio_sub$StartLat))
sum(is.na(bio_sub$StartLong))


## look at sites

# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)

## plot points
m1 <- mapview(bio_sub, cex=6, col.regions="orange",
              layer.name="Toad Observations") 


m1
m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

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

bio_data2 <- full_join(data_segs, data_segs_df, by = "ID")
object.size(bio_data2)

length(unique(bio_data2$COMID)) ## 214


# Add to predictors  ------------------------------------------------------

head(data)

st_crs(bio_data_sf) == st_crs(data_cds)

## make sure have same coordinates
data_cds <- data %>%
  dplyr::select(ID.2, X, Y, PresAbs2005,  Presence2) %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326, remove=F) 

bio_data_sf <- st_transform(bio_data2, crs=4326) #


## spatial join

test <- st_join(bio_data_sf, data_cds) ## nothing matches
test


