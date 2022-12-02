## testing comid matches and old vs new spatial stuff

## workflow
## 	Check extract with line string - Mike data and new rasters
# Which comids missing - plot
# Use new raster (GIS) to extract

## packages
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

## ## upload raster mask
rmask <- raster("ignore/02_mask_raster_network.tif")
crs(rmask)

# Original data from Mike -------------------------------------------------


path <- "original_model/FullData/"

## Mike's data - env vars with presence absence
data <- read.csv(paste0(path, "200mCells_FullData_PresAbs_Complete_ThinnedCols.csv"))
head(data)

## get coords change crs
orig_grids <- data %>%
  dplyr::select(ID.2, X, Y) 

head(orig_grids)
## make spatial and transformCRS
coordinates(orig_grids) <- ~X+Y
projection(orig_grids) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
## make raster and add to stack

x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)
projection(x) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"

crs(x)

## make raster of id in original grids
  x1 <- raster::rasterize(orig_grids, x, field="ID.2", na.rm =TRUE, sp = TRUE)
  projection(x1)<-"+proj=geocent +ellps=GRS80 +units=m +no_defs"

crs(x1) <- crs(rmask)

# NHD comids --------------------------------------------------------------

nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")

nhd <- nhd %>%
    st_as_sf %>%
    st_simplify(dTolerance = 0.5, preserveTopology = T)

crs(nhd) 
  
nhdLine <- st_cast(nhd, "LINESTRING") #%>% st_transform(crs = 4326)
nhdPoints <- st_cast(nhd, "POINT")
nhdPoints
# Extract values of raster at nhd lines -----------------------------------

# OrigCellsNHD <- raster::extract(x1, nhdLine, cellnumbers=TRUE)
## extract rasters at points
OrigCellsNHDPts <- raster::extract(x1, nhdPoints, cellnumbers=TRUE)

head(OrigCellsNHDPts)

## join together
origCOMIDs <- cbind(nhdPoints, OrigCellsNHDPts)

head(origCOMIDs)
## check no of comids and nas
length(unique(origCOMIDs$COMID)) ## 2179
ind <- which(is.na(origCOMIDs$layer))

origCOMIDs[ind,]

## remove nas from raster layer - this is where no value exists at the nhd locations

origCOMIDsSub <- origCOMIDs %>% drop_na(layer)

length(unique(origCOMIDsSub$COMID)) ## 2107 - only 10 missing!! yay!

plot(origCOMIDsSub)


# extract new raster (GIS) ------------------------------------------------------

## upload 
coms <- raster("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification_PolylineToRaster_v3.tif")
coms ## 
names(coms) <- "COMID"
# coms<-projectRaster(coms, crs = crs(rmask))
# coms


## resample to mask layer
comsRE <- resample(coms, rmask, method = "ngb")

## how many comids?
comx <- na.omit(as.data.frame(comsRE))
length(unique(comx$COMID)) ## 2000

crs(comsRE)

GISCellsNHDPts <- raster::extract(coms, nhdPoints, cellnumbers=TRUE)
GISCellsNHDPts

## join together
GISCOMIDs <- cbind(nhdPoints, GISCellsNHDPts)
GISCOMIDs ## comids don't always match
## remove nas from raster layer - this is where no value exists at the nhd locations

GISCOMIDsSub <- GISCOMIDs %>% drop_na(COMID.1)

length(unique(GISCOMIDsSub$COMID)) ## 2138 hmmm

plot(origCOMIDsSub)


# New data raster (made on 00) --------------------------------------------

n1 <- stack("ignore/00_raw_new_data_raster.tif")
n1
# ## load layer names 
# load(file = "output_data/00_raster_layer_names.RData")
# layerNames

n1 <- n1[[1]] ## get only 1 for test

crs(n1) <- crs(rmask)

NewCellsNHDPts <- raster::extract(n1, nhdPoints, cellnumbers=TRUE)
NewCellsNHDPts

## join together
NewCOMIDs <- cbind(nhdPoints, NewCellsNHDPts)
NewCOMIDs
## remove nas from raster layer - this is where no value exists at the nhd locations

NewCOMIDsSub <- NewCOMIDs %>% drop_na(X00_raw_new_data_raster.1)

length(unique(NewCOMIDsSub$COMID)) ## 2107

plot(origCOMIDsSub)


# rmask extract -----------------------------------------------------------

maskCellsNHDPts <- raster::extract(rmask, nhdPoints, cellnumbers=TRUE)
maskCellsNHDPts

## join together
MaskCOMIDs <- cbind(nhdPoints, maskCellsNHDPts)
MaskCOMIDs
## remove nas from raster layer - this is where no value exists at the nhd locations

MaskCOMIDsSub <- MaskCOMIDs %>% drop_na(X02_mask_raster_network)

length(unique(MaskCOMIDsSub$COMID)) ## 2036

plot(origCOMIDsSub)


# Match with hydro comids -------------------------------------------------

head(NewCOMIDsSub) ## original cells with new data, use this as working from it

## hydro data

delta <- read.csv("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/2022-07-28_predicted_abs_FFM_deltaFFM_SD_COMIDS_medianDelta_test5_norunoff.csv")
head(delta)

length(unique(delta$comid)) ## 2117

comids <- unique(delta$comid)

## compare ffm comids with raster data
Rtest <- NewCOMIDsSub %>%
  filter(COMID %in% comids)

length(unique(Rtest$COMID)) ## 2048

## which are missing?

Missing <- NewCOMIDsSub %>%
  filter(!COMID %in% comids)

miscoms <- unique(Missing$COMID)

## save list of missing comids
write.csv(miscoms, "output_data/00a_missing_comids.csv")

plot(Missing)


# Different comids at points ----------------------------------------------

## read in dfcreated in 00 - RB9 data
load(file = "ignore/03_RB9_grdded_data.RData")
head(data_hyd_sf)

DataComs <- DataComs  %>% distinct()

test1 <- DataComs %>%
  filter(COMID %in% c(22549169, 20348307))

test1

load(file = "ignore/02_all_data_for_model_gridded.RData") ## data only at pres/abs
head(NewDataObsSub)

test2 <- NewDataObsSub %>%
  filter(COMID %in% c(22549169, 20348307))

test2

load(file = "ignore/03_RB9_grdded_data_observations.RData")
head(data_hyd_sf_obs)
dim(data_hyd_sf_obs)
length(unique(data_hyd_sf_obs$COMID)) ## 318

data_hyd_sf_obs2 <- data_hyd_sf_obs %>%
  distinct(ID.1, .keep_all=T)
  
dim(data_hyd_sf_obs2)
length(unique(data_hyd_sf_obs2$COMID)) ## 248
data_hyd_sf_obs2

test3 <- data_hyd_sf_obs %>%
  filter(COMID %in% c(22549169, 20348307))

test3


## upload nhd shape
nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
## simplify
nhd <- nhd %>%
  st_as_sf %>%
  st_simplify(dTolerance = 0.5, preserveTopology = T)

crs(nhd) 

## convert to points
nhdPoints <- st_cast(nhd, "POINT")
nhdPoints 

test <- nhdPoints %>%
  filter(COMID %in% c(22549169, 20348307))

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


# Model tests -------------------------------------------------------------
load(file = "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/original_model/Current/randomForests/PARTITIONING/Cur_10x_Final_Feb2014/Model1/RFClassModel.RData")
mod1 <- model

load(file = "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/original_model/Current/randomForests/PARTITIONING/Cur_10x_Final_Feb2014/Model3/RFClassModel.RData")
mod3 <- model

mod1
mod2
mod3
  