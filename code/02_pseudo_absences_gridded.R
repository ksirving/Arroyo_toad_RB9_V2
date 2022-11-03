### pseudo absences for gridded data

# packages
library(tidyverse)
library(tidylog)
library(sf)
library(raster)
library(nhdplusTools)
library(readr)
library(mapview)

library(sp)
library(rgdal)
library(kernlab)
library(rgl)
library(ks)
library(sm)
library(caret)

getwd()
setwd("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")
## workflow
# upload original and new observations
## snap new obs within 50m buffer
# create pseudo abs with bias surface on gridded data from script 00 
# save gridded data as rasters

# raw env data raster------------------------------------------------------------

## raw data - env vars with presence absence
load(file = "ignore/00_all_env_data_gridded.RData")
head(data_hyd_sf)

envData <- na.omit(data_hyd_sf)
head(envData)

# Physical data -----------------------------------------------------------
getwd()

## read in stack created in 00
xvars <- stack("ignore/00_raw_new_data_raster.grd") ## new rasters @ 200m - change as needed
# xvars <- xvars[[2:97]] ## remove template raster
crs(xvars)
xvars

# Bio data ----------------------------------------------------------------

### add all pres/abs together on previous script, then use that one - snapped. 
## get path for functions
source("original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

# all new pres/abs
# inshape="ignore/01_toad_obs_points_RB9.shp" ## presence absence
# class(sdata)
# orig.sdata<- sdata <- shapefile(inshape) ## p/a

## updated and snapped pres/abs
inshape2 = "output_data/ToadsObs_50m_Final.shp"
orig.sdata<- sdata <- shapefile(inshape2) ## p/a

head(sdata)
crs(sdata)

# NUMBER OF BOOTSTRAP REPLICATES
# b=10001



# Snap occurrence to stream grids -----------------------------------------

## snap points within a 50m buffer
## prep data for snap

## make a raster mask to snap points to
xvars1 <- xvars[[1]]
crs(xvars1) <- "+proj=utm +zone=11 +datum=NAD83"
crs(xvars1)
# ## save out
# writeRaster(xvars1, "ignore/02_mask_raster_network.tif", format="GTiff", overwrite=T)

# snapped in in GIS by Abel, snapped to 50m 

## upload snapped data
bioSnap <- shapefile("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/output_data/ToadsObs_50m_Final.shp")
crs(bioSnap)
## upload raster mask
rmask <- raster("ignore/02_mask_raster_network.tif")
crs(rmask)

## plot grdis and points
# plot(rmask)
# plot(bioSnap, add=T)

## get coordinates of snapped points
bioSnapCds <- bioSnap %>% as.data.frame() %>%
  dplyr::select(ID, COMID, Year:PresAbs, centroid_x, centroid_y) %>%
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs=4326, remove=F) %>%
  # filter(PresAbs == 1) %>%
  st_transform(crs = 9001)

## get cell numbers at pres/abs sites 
cellsPres <- raster::extract(rmask, bioSnapCds, cellnumbers=TRUE)
dim(cellsPres) ## 2797
# head(cellsPres)
# cellsPres
## remove z dimension
bioSnapCds <- st_zm(bioSnapCds)

## change CRS to save, change back to 9001 on upload
bioSnapCds <- bioSnap %>% 
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs=9001, remove=F) %>%
  st_transform(crs = 3310) %>% dplyr::select(PresAbs, centroid_x, centroid_y)

## save out
st_write(bioSnapCds, "output_data/ToadsObs_50m_Final_snapped_coords.shp", append=FALSE) ## check in QGIS


## updated and snapped pres/abs
inshape2 = "output_data/ToadsObs_50m_Final_snapped_coords.shp"
orig.sdata<- shapefile(inshape2) ## p/a
head(orig.sdata)


# KDE Bias Surface --------------------------------------------------------
set.seed(234)

## upload raster for mask (old one as it works for bias)
# xvars <- stack("ignore/00_raw_data_raster.grd")
# xvars <- xvars[[2:97]] ## remove template raster

# develop KDE sampling bias surface
orig.sdata2<-subset(orig.sdata, PresAbs==1)

## format raster
mask <- rmask[[1]]>-1000
crs(mask)
## get cell numbers for bias
bias <- cellsPres[,1]
cells <- unique(sort(bias))
## get coords from bias
kernelXY <- xyFromCell(mask, cells)
samps <- as.numeric(table(bias))

# code to make KDE raster
KDEsur <- sm.density(kernelXY, weights=samps, display="none", ngrid=782, 
                     ylim=c(3600402,3730202), xlim=c(423638,563638), nbins=0)

KDErast=SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
KDErast = SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, 
                                                                 length(KDEsur$estimate))))
KDErast <- raster(KDErast)
KDErast <- resample(KDErast, mask)
KDErast <- KDErast*mask
KDEpts <- rasterToPoints(KDErast) #(Potential) PseudoAbsences are created here 

#Now to integrate Pseudoabsences into Presence Data
a=KDEpts[sample(seq(1:nrow(KDEpts)), size=702, replace=T, prob=KDEpts[,"layer"]),1:2] 

PA.abs<-data.frame(PresAbs=rep(0,nrow(a)))
a.sp<-SpatialPoints(a, proj4string=CRS("+proj=utm +zone=11 +datum=WGS84"))
a.spdf<-SpatialPointsDataFrame(a.sp, PA.abs)
crs(a.spdf)
head(a.spdf)
## format to bind to other data
a.spdf <- a.spdf %>%
  st_as_sf(coords=c("Longitude",  "Latitude"), crs=9001, remove=F) #%>%
  mutate(Longitude = unlist(map(geometry,1)),
         Latitude = unlist(map(geometry,2))) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("P", ID)) %>% 
  mutate(Year = NA, LifeStage = "Pseudo", Count = NA)
head(a.spdf)
crs(a.spdf)

## issue with coords in a.spdf - don't match sdata1 in biodata below so can't extract raster data. 
# Get comids for pseudo obs data ------------------------------------------

orig.sdata.segs <- a.spdf %>%
  st_transform(crs=32611) %>%
  arrange(ID)

orig.sdata.segs

# use nhdtools to get comids
data_all_coms <- orig.sdata.segs %>%
  group_split(ID) %>%
  set_names(., orig.sdata.segs$ID) %>%
  map(~discover_nhdplus_id(.x$geometry))

# flatten into single dataframe instead of list
data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
  as.data.frame() %>%
  rename("COMID"=V1) %>% rownames_to_column(var = "ID") 
data_segs_df

## format to bind with other data. change name of lat/lon to centroidx/y to match snapped data
a.spdf <- full_join(orig.sdata.segs, data_segs_df, by = "ID") %>%
  as.data.frame() %>% 
  mutate(centroid_x = Longitude, centroid_y = Latitude, Join_Count = NA, Id_1 = NA) %>%
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs=9001, remove=F) %>%
  st_transform(crs = 3310)

head(a.spdf)
dim(a.spdf)
class(a.spdf)
crs(a.spdf)
crs(bioSnap)
##  make sdata into sf object
sdata <- bioSnap %>%
  as.data.frame() %>%
  dplyr::select(-coords_x1, -coords_x2) %>%
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs=3310, remove=F) %>%
  mutate(COMID = as.integer(COMID)) 

  dim(sdata)
  crs(sdata)
  head(sdata)
## join with other observations
sdata1<-bind_rows(sdata,a.spdf)
dim(sdata1)
head(sdata1)
sdata1 <- sdata1 %>% dplyr::select(ID, COMID, Year:PresAbs, centroid_x, centroid_y) %>%
## save as spatial object
st_write(sdata1, "ignore/02_pres_abs_pseudo_abs.shp", append = F) ## pseudo_abs
class(sdata)

# Add to env data ---------------------------------------------------------

## presence/absence data
sdata <- st_read("ignore/02_pres_abs_pseudo_abs.shp")

## read in stack created in 00
xvars <- stack("ignore/00_raw_new_data_raster.tif") ## new rasters @ 200m - change as needed
xvars
crs(xvars)
crs(xvars) <- "+proj=utm +zone=11 +datum=NAD83"

## get coordinates of snapped points, change NAs in lifestage so csn drop all nas from DF
bioData <- sdata1 %>% #as.data.frame() %>%
  # dplyr::select(ID, COMID, Year:PresAbs, centroid_x, centroid_y) %>%
  # st_as_sf(coords=c("centroid_x", "centroid_y"), crs=4326, remove=F) %>%
  st_transform(crs = 9001) %>%
  mutate(LifeStage = ifelse(is.na(LifeStage), "Other", LifeStage), 
         Count = ifelse(is.na(Count), 0, Count ))# %>%
  # drop_na()

head(bioData)
crs(bioData)
## remove z dimension
bioData <- st_zm(bioData)

## get raster values at all raster cells
cellsAll <- raster::extract(xvars, 1:ncell(xvars),  df=TRUE)
cellsAll <- na.omit(cellsAll)
## change column names to raster layer names
names(cellsAll) <- c("cells", layerNames)
cellsAll

## get coords at cells
cellsxy <- data.frame(rasterToPoints(xvars))
## add names
names(cellsxy)[3:37] <- layerNames

## join x y with all cells

NewDataObs <- cbind(cellsxy[1:2], cellsAll) ## save as df with coords
dim(NewDataObs) ## 17761
# names(NewDataObs)

## save out
save(NewDataObs, file = "ignore/02_all_data_for_prediction_gridded.RData")

## get raster values at all presence cells
# xvarsRe <- projectRaster(xvars, crs = crs(rmask)) ## reproject raster to match bio data

cellsPres <- raster::extract(xvars, bioData,  cellnumbers = T, df=TRUE)
names(cellsPres)[3:37] <- layerNames
head(cellsPres)
crs(xvarsRe)
dim(xvars)

# test <- filter(cellsAll, cells %in% cellsPres$cells)

## join together biodata and cells with presences

NewDataObsSub <- cbind(bioData, cellsPres)
head(NewDataObsSub)
dim(NewDataObsSub) ## 2797
names(NewDataObsSub)
## upload layer names

load(file = "output_data/00_raster_layer_names.RData") ## layerNames

## change column names to raster layer names
names(NewDataObsSub)[10:44] <- layerNames

## save data for model 
save(NewDataObsSub, file = "ignore/02_all_data_for_model_gridded.RData")

## save data for prediction
# cellsxy <- raster::xyFromCell(xvars[1]) ## might want cell numbers but not working
# head(cellsxy)


