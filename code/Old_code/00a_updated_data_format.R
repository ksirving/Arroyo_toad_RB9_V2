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

save(bio_data2, file = "ignore/00a_al_bio_data_COMIDs.RData")
load(file = "ignore/00a_al_bio_data_COMIDs.RData") # bio_data2

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

save(bio_data3, file = "ignore/00a_al_bio_data_COMIDs_SDRB9.RData")
load(file = "ignore/00a_al_bio_data_COMIDs_SDRB9.RData") # bio_data3

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


# Add to predictors  ------------------------------------------------------

head(DataObs) ## taod observations
names(all_data) ## env data
str(all_data)

# st_crs(bio_data_sf) == st_crs(data_cds)
# 
# ## make sure have same coordinates
# data_cds <- data %>%
#   dplyr::select(ID.2, X, Y, PresAbs2005,  Presence2) %>% 
#   st_as_sf(coords=c("X", "Y"), crs=4326, remove=F) 
# 
# bio_data_sf <- st_transform(bio_data2, crs=4326) #


## join by comid, change NAs in pres/abs to -999

data_bio_env <- full_join(DataObs, all_data, by = "COMID") %>%
  mutate(NewObs = as.numeric(NewObs)) %>%
  mutate(NewObs = replace_na(NewObs, -999))

data_bio_env

length(unique(data_bio_env$COMID)) ## 2296

### 


# Format Climate data -------------------------------------------------

## remove old climate data (also removed remote sensing, i think, until updated)
head(data_bio_env)
names(data_bio_env)

data_red <- data_bio_env %>%
  dplyr::select(-c(X81pptCr1:X81TMnCr13))

head(data_red)

#### upload climate data

## precipitation
ppt_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average ppt_30yr.shp")
ppt_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average ppt_30yr.shp")

## temperature
tmax_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmax_30yr.shp")
tmin_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmin_30yr.shp")

tmax_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmax_30yr.shp")
tmin_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmin_30yr.shp")

## upload original data with comids

load(file= "ignore/00_all_env_bio_data.RData") ## data2
head(data2)

range(data2$FINAL...01..3)

orig_grids <- data2 %>%
  dplyr::select(ID.2, X, Y, COMID) %>%
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) %>%
  st_transform(crs=4269) 

### make sure same crs

head(ppt_mon)
  
st_crs(orig_grids) == st_crs(ppt_ann)

## spatial join - Precip 

pptAnn_join <- ppt_ann %>%
  st_join(orig_grids) %>%
  rename(ppt_ann = PRISM__) %>%
  as.data.frame() %>% na.omit() %>%
  group_by(COMID) %>%
  summarise(pptAnnAv = mean(ppt_ann))

## spread months from rows to columns with month number
pptMon_join <- ppt_mon %>%
  st_join(orig_grids) %>%
  as.data.frame() %>% na.omit() %>%
  group_by(COMID, Month) %>%
  summarise(pptMonAv = mean(pptMonth))  %>%
  mutate(Variable = "pptMonAv") %>%
  unite(col = VariableMon, c(Variable, Month), sep="") %>%
  pivot_wider(names_from = VariableMon, values_from = pptMonAv)


pptData <- full_join(pptAnn_join, pptMon_join, by = "COMID")

### spatial join - Temp

TmaxAnn_join <- tmax_ann %>%
  st_join(orig_grids) %>%
  rename(tmax_ann = PRISM__) %>%
  as.data.frame() %>% na.omit() %>%
  group_by(COMID) %>%
  summarise(TmaxAnn = mean(tmax_ann))

## spread months from rows to columns with month number
TmaxMon_join <- tmax_mon %>%
  st_join(orig_grids) %>%
  # rename(ppt_ann = pptMonth) %>%
  as.data.frame() %>% na.omit() %>%
  group_by(COMID, Month) %>%
  summarise(TmaxMon = mean(tmaxMonth)) %>%
  mutate(Variable = "TmaxMon") %>%
  unite(col = VariableMon, c(Variable, Month), sep="") %>%
  pivot_wider(names_from = VariableMon, values_from = TmaxMon)


TmaxData <- full_join(TmaxAnn_join, TmaxMon_join, by = "COMID")


  TminAnn_join <- tmin_ann %>%
    st_join(orig_grids) %>%
    rename(tmin_ann = PRISM__) %>%
    as.data.frame() %>% na.omit() %>%
    group_by(COMID) %>%
    summarise(TminAnn = mean(tmin_ann))
  
  ## spread months from rows to columns with month number
  TminMon_join <- tmin_mon %>%
    st_join(orig_grids) %>%
    as.data.frame() %>% na.omit() %>%
    group_by(COMID, Month) %>%
    summarise(TminMon = mean(tminMonth)) %>%
    mutate(Variable = "TminMon") %>%
    unite(col = VariableMon, c(Variable, Month), sep="") %>%
    pivot_wider(names_from = VariableMon, values_from = TminMon)
  
  TminData <- full_join(TminAnn_join, TminMon_join, by = "COMID")
  
### join all data together
 
  tempData <- full_join(TminData,  TmaxData, by = c("COMID")) 
  allData <- full_join(pptData, tempData, by = c("COMID"))

  head(allData)  
  length(unique(allData$COMID))
  
  sum(unique(allData$COMID) %in% unique(data_red$COMID))

  head(data_red)
  length(unique(data_red$COMID))
  
  ## join back to main df
  NewData <- full_join(allData, data_red, by = "COMID") #%>% drop_na(PresAbs)
  
  head(NewData)
  dim(NewData)
  length(unique(NewData$COMID))
  
  save(NewData, file = "ignore/00a_new_clim_obs_env_abs_flow.Rdata")
  
load(file = "ignore/00a_new_clim_obs_env_abs_flow.Rdata")
 
 
  str(nhd_lines_obs)
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
  

# Join back to main DF and save -------------------------------------------

head(NewData)
names(NewData)  
dim(NewData)
nhd_lines_rb9$geometry

##  join all env data for region to flow lines
NewDataObs <- NewData %>%
 right_join(nhd_lines_rb9, by = "COMID") %>%
  mutate(COMID = as.integer(COMID)) %>%
  distinct() %>%
    mutate(NewObs = as.numeric(NewObs)) %>%
    mutate(NewObs = replace_na(NewObs, -999)) %>%
    st_as_sf(sf_column_name = "geometry")



length(unique(NewDataObs$COMID)) ## all comids - 1992
length(unique(NewData$COMID))

## join data with only known occurrences - for model build

NewDataObsSub <- NewDataObs %>%
  filter(COMID %in% nhd_lines_obs$COMID) 


length(unique(NewDataObsSub$COMID)) ## occurrence comids 353
length(unique(NewDataObs$COMID)) ## 1992???

# Save out ----------------------------------------------------------------

class(NewDataObsSub)

save(NewDataObsSub, file =  "ignore/00a_data_for_model_abs_flow.RData")
save(NewDataObs, file =  "ignore/00a_data_for_prediction_abs_flow.RData")

# Remote Sensing data -----------------------------------------------------

## upload dry season data

sept <- brick("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/TC_2014_RB9/TC_092014_RB9.tif")
# plot(sept)
sept
# install.packages("RStoolbox")

## convert rb9 comids to points
nhd_points_rb9 <- st_cast(nhd_lines_rb9, "POINT") %>% dplyr::select(-SHAPE_LENG) %>% st_zm()
head(nhd_points_rb9)

## extract raster values
sept_values <- raster::extract(sept, nhd_points_rb9)
# head(sept_values_coord)
sept_values_coord <- cbind(sept_values[,1:3], nhd_points_rb9)

## upload wet season data
april <- brick("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/TC_2014_RB9/TC_042014_RB9.tif")
# plot(april)

## extract raster values
apr_values <- raster::extract(april, nhd_points_rb9)
# head(apr_values)

## join with sept values
tass_sp <- cbind(apr_values[,1:3], sept_values_coord)
head(tass_sp)
dim(tass_sp)

save(tass_sp, file = "ignore/00a_tass_cap_all_sp.RData")

#### summarise by comid - mean and variance

## make longer and group by comid and variable

tass_sp_long <- tass_sp %>%
  pivot_longer(TC_042014_RB9.1:TC_092014_RB9.3, names_to = "Variable", values_to = "Value") %>%
  group_by(COMID, Variable) %>% 
  summarise(MeanVals = mean(Value), VarVals = var(Value)) ## calculations

## make wider rename with mean/var in name

## means
tass_sp_means <- tass_sp_long %>%
  mutate(VarNames = paste0(Variable, "_Var"), Variable = paste0(Variable, "_Mean")) %>%
  pivot_wider(id_cols = "COMID", names_from = Variable, values_from = MeanVals) 

## variance
tass_sp_vars <- tass_sp_long %>%
    mutate(VarNames = paste0(Variable, "_Var"), Variable = paste0(Variable, "_Mean")) %>%
    pivot_wider(id_cols = "COMID", names_from = VarNames, values_from = VarVals) 

## join together

tass_all <- full_join(tass_sp_vars, tass_sp_means, by = "COMID" )

head(tass_all)  

save(tass_all, file="ignore/00a_remote_sensing_data_per_comid.RData")


# Join to main data and save ----------------------------------------------

class(NewDataObsSub)

load(file =  "ignore/00a_data_for_model_abs_flow.RData")
load(file =  "ignore/00a_data_for_prediction_abs_flow.RData")

head(NewDataObs)

## predic data
NewDataObs <- NewDataObs %>% 
  dplyr::select(-contains("FINAL")) %>%
  right_join(tass_all, by = "COMID")

## model data
NewDataObsSub <- NewDataObsSub %>% 
  dplyr::select(-contains("FINAL")) %>%
  left_join(tass_all, by = "COMID") 
  
dim(NewDataObs)
dim(NewDataObsSub)

save(NewDataObsSub, file =  "ignore/00a_data_for_model_abs_flow.RData")
save(NewDataObs, file =  "ignore/00a_data_for_prediction_abs_flow.RData")

