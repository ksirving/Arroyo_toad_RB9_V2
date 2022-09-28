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
library(sf)
library(caret)

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

# Bio data ----------------------------------------------------------------

### add all pres/abs together on previous script, then use that one. currently using original to get the code
## get path for functions
source("original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

# original pres/abs
inshape="ignore/01_toad_obs_points_RB9.shp" ## presence absence
class(sdata)
orig.sdata<- sdata <- shapefile(inshape) ## p/a

head(sdata)

# NUMBER OF BOOTSTRAP REPLICATES
# b=10001

projection(sdata)<-"+proj=utm +zone=11 +datum=WGS84"
projection(orig.sdata)<-"+proj=utm +zone=11 +datum=WGS84"
crs(sdata)



# Physical data -----------------------------------------------------------
getwd()

## read in stack created in 00
xvars <- stack("ignore/00_raw_new_data_raster.grd") ## new raster @ 200m - change as needed
# xvars <- xvars[[2:97]] ## remove template raster
crs(xvars)


# Snap occurrence to stream grids -----------------------------------------

## snap points within a 50m buffer

library('devtools')
install_github("mtalluto/WatershedTools")
library("WatershedTools")
library(geosphere)
?extract

xvars1 <- xvars[[1]]
xvars1
crs(xvars1) <- "+proj=utm +zone=11 +datum=NAD83"
crs(xvars1)
writeRaster(xvars1, "ignore/02_mask_raster_network.tif", format="GTiff", overwrite=T)
?writeRaster
plot(xvars1)
plot(sdata, add=T)

xvarPts <- rasterToPoints(xvars1)
xvarPts <- xvarPts %>%
  as.data.frame() %>%
  st_as_sf(coords=c("x", "y"), crs=9001, remove=F) 

snap_data <- sdata %>%
  as.data.frame() %>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F) %>%
  st_transform(crs = 3310)

st_write(snap_data, "ignore/02_01_toad_obs_points_RB9_3310.shp", append = F)
?st_join
test1 <- st_join(xvarPts, snap_data)
test1
head(snap_data)
dim(snap_data)
dim(xvarPts)
test <- raster::extract(xvars1, sdata)
test
sum(is.na(test))

# KDE Bias Surface --------------------------------------------------------
set.seed(234)

# develop KDE sampling bias surface
orig.sdata2<-subset(orig.sdata, PresAbs==1)

mask <- xvars[[1]]>-1000

bias <- cellFromXY(mask, orig.sdata[-1])

cells <- unique(sort(bias))

kernelXY <- xyFromCell(mask, cells)
samps <- as.numeric(table(bias))

length(samps)

# code to make KDE raster
KDEsur <- sm.density(kernelXY, weights=samps, display="none", ngrid=880, 
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

## format to bind to other data
a.spdf <- a.spdf %>%
  st_as_sf(coords=c("Longitude",  "Latitude"), crs=9001, remove=F) %>%
  mutate(Longitude = unlist(map(geometry,1)),
         Latitude = unlist(map(geometry,2))) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("P", ID)) %>% 
  mutate(Year = NA, LifeStage = "Pseudo", Count = NA)

class(a.spdf)
# Get comids for pseudo obs data ------------------------------------------

orig.sdata.segs <- a.spdf %>%
  # dplyr::select(ID, Latitude, Longitude) %>%
  # distinct(ID, Latitude, Longitude) %>%
  # st_as_sf(coords=c("Longitude",  "Latitude"), crs=9001, remove=F) %>%
  st_transform(crs=32611) %>%
  arrange(ID)

orig.sdata.segs

# use nhdtools to get comids
data_all_coms <- orig.sdata.segs %>%
  group_split(ID) %>%
  set_names(., orig.sdata.segs$ID) %>%
  map(~discover_nhdplus_id(.x$geometry))

data_all_coms

# flatten into single dataframe instead of list
data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
  as.data.frame() %>%
  rename("COMID"=V1) %>% rownames_to_column(var = "ID") #%>%
  # mutate(ID = as.integer(ID))
head(data_segs_df)

a.spdf <- full_join(orig.sdata.segs, data_segs_df, by = "ID")
# object.size(data1)

a.spdf
##  make sdata into sf object
sdata <- sdata %>%
  st_as_sf(coords=c("Longitude",  "Latitude"), crs=9001, remove=F) 
  
## join with other observations
sdata<-bind_rows(sdata,a.spdf)

## save as spatial object
st_write(sdata, "ignore/02_pseudo_abs.shp", append = F) ## pseudo_abs
class(sdata)

# Add to env data ---------------------------------------------------------

?st_write
load(file= "ignore/00_all_env_bio_data.RData")
## make spatial

data_sf<- data2 %>%
  dplyr::select(-PresAbs2005, -Presence2) %>%
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) 

head(data_sf)

## map

# set background basemaps:
# basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
#                   "Esri.NatGeoWorldMap",
#                   "OpenTopoMap", "OpenStreetMap", 
#                   "CartoDB.Positron", "Stamen.TopOSMFeatures")
# 
# mapviewOptions(basemaps=basemapsList, fgb = FALSE)
# 
# m1 <- mapview(data_segs, cex=6, col.regions="orange",
#               layer.name="data points") 
# 
# m1
# m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

## save as shape

st_write(data_sf, "ignore/00_all_env_bio_data.shp", append=FALSE) ## not working??
# 
# check <- st_read("output_data/00_all_env_bio_data.shp")
# plot(check)()
data_sf<- st_read("ignore/00_all_env_bio_data.shp")
data_sf
dim(data_sf)
class(data_sf)

## make long and get mean value per comid

# data_hyd_sf_long <- data_sf %>%
#   pivot_longer(MRVBF.Mx:AvgSlope, names_to = "Variable", values_to = "Value") %>%
#   group_by(COMID, Variable) %>%
#   mutate(Meanvals = mean(Value),
#             Minvals = min(Value),
#             Maxvals = max(Value))
#   
#   save(data_hyd_sf_long, file="ignore/00_all_env_data_scaled.RData")
#   head(data_hyd_sf_long)
#   rm(data_hyd_sf_long)
#   
#   unique(data_hyd_sf_long$Variable)
#   ## elevation = min, catchment area = max, MRVBF = max, VRM1.Mn = min
#   
#   ## format df 
#   data_hyd_sf_longer <- data_hyd_sf_long %>%
#     pivot_longer(Meanvals:Maxvals, names_to= "Stat", values_to = "Values") %>%
#     pivot_wider(names_from = Variable, values_from = Values)
#   
#   head(data_hyd_sf_longer)
#   names(data_hyd_sf_longer)
#   ## take specific stats for each variable
#   meanVars <- data_hyd_sf_longer %>%
#     filter(Stat == "Meanvals") %>%
#     dplyr::select(c(COMID:AvgWaterSt, FINAL...01..3:FINAL...24..4, X81pptCr1:X81TMxCr9))
#   
#   
#   MinVars <- data_hyd_sf_longer %>%
#     filter(Stat == "Minvals") %>%
#     dplyr::select(c(COMID, DEM_10m.Mn, VRM1.Mn:VRM9.Mn))
#   
#   MaxVars <- data_hyd_sf_longer %>%
#     filter(Stat == "Maxvals") %>%
#     dplyr::select(c(COMID, Catchment.A, MRVBF.Mx))
#   
#   ## join toegther
#   
#   all_data <- bind_cols(meanVars, MinVars, MaxVars)
#   
#   head(all_data)
#   names(all_data)
#   
#   all_data <- all_data %>%
#     dplyr::select(-COMID...92, -COMID...99, -Stat) %>%
#     rename(COMID = COMID...1)

# Convert to raster -------------------------------------------------------

## create template raster

## dims etc from CurrentGridFeb14.grd


x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)

projection(x) <- "+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
x

#Create a "rasterBrick" of the template raster
x<-brick(x)

x2 <- x

#Create list of column names you want to rasterize
fields <- names(data_sf) [5:100]
fields
## make rasters of each env var
i
for (i in fields){
  x[[i]]<-rasterize(data_sf, x, field=i)
  projection(x)<-"+proj=utm +zone=11 +datum=WGS84"
  x <- stack(x)
}

x@layers
x2[[i]]

plot(x[[2]]) ## to check
x
## save out
writeRaster(x, "ignore/00_raw_data_raster.grd", format="raster", crs="+proj=utm +zone=11 +datum=WGS84", overwrite=TRUE)
