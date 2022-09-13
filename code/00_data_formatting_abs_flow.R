## data formatting

##workflow
## upload raw data
## get comids to match later
## create rasters
## add hydro

# packages
library(tidyverse)
library(tidylog)
library(sf)
library(raster)
library(nhdplusTools)
library(readr)
library(mapview)


setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")
getwd()

# raw env data ------------------------------------------------------------
path <- "original_model/FullData/"

## raw data - env vars with presence absence
data <- read.csv(paste0(path, "200mCells_FullData_PresAbs_Complete_ThinnedCols.csv"))
head(data)
str(data$ID.2)

# Get COMIDs --------------------------------------------------------------

# Create dataframe for looking up COMIDS (here use all stations)
data_segs <- data %>%
  dplyr::select(ID.2, X, Y) %>%
  distinct(ID.2, X, Y) %>% 
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) %>%
  arrange(ID.2)

head(data_segs)
dim(data_segs)

# use nhdtools to get comids
data_all_coms <- data_segs %>%
  group_split(ID.2) %>%
  set_names(., data_segs$ID.2) %>%
  map(~discover_nhdplus_id(.x$geometry))

data_all_coms

# flatten into single dataframe instead of list
data_segs_df <-data_all_coms %>% flatten_dfc() %>% t() %>%
  as.data.frame() %>%
  rename("COMID"=V1) %>% rownames_to_column(var = "ID.2") %>%
  mutate(ID.2 = as.integer(ID.2))
head(data_segs_df)

data2 <- full_join(data, data_segs_df, by = "ID.2")
object.size(data2)


# Bio data ----------------------------------------------------------------
getwd()

## define and set path
path=paste0(wsdir, "/PARTITIONING/DATA3") 
setwd(path)

## get path for functions
source("Functions.R")

# Training data
inshape="200mCells_PresAbs2005 PCA PresAbs Pts_MLT.shp" ## presence absence

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

# READ SPECIES OBSERVATION DATA 

orig.sdata<- sdata <- shapefile(inshape) ## p/a
# str(sdata@data)
proj4string(sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")
proj4string(orig.sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")
orig.sdata
#sdata <- sdata[1]
# str(sdata)
proj4string(sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")

head(sdata)


# Physical data -----------------------------------------------------------
getwd()

setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")
## read in stack created in 00
xvars <- stack("ignore/00_raw_data_raster.grd") ## original raster @ 200m - change as needed
xvars <- xvars[[2:97]] ## remove template raster

# KDE Bias Surface --------------------------------------------------------
set.seed(234)

# develop KDE sampling bias surface
orig.sdata2<-subset(orig.sdata, PresAbs200==1)

mask <- xvars[[1]]>-1000

bias <- cellFromXY(mask, orig.sdata[,-1])

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

PA.abs<-data.frame(PresAbs200=rep(0,nrow(a)))
head(PA.abs)
a.sp<-SpatialPoints(a, proj4string=CRS("+proj=utm +zone=11 +datum=WGS84"))
a.spdf<-SpatialPointsDataFrame(a.sp, PA.abs)

sdata<-rbind(sdata,a.spdf)

sdata

### get comids for obs data

orig.sdata.segs <- sdata %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(ID = rowname) %>%
  mutate(ID = as.integer(ID)) %>%
  st_as_sf(coords=c("coords.x1", "coords.x2"), crs=32611, remove=F) %>%
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
  rename("COMID"=V1) %>% rownames_to_column(var = "ID") %>%
  mutate(ID = as.integer(ID))
head(data_segs_df)

data1 <- full_join(orig.sdata.segs, data_segs_df, by = "ID")
object.size(data1)


## save as object
save(data2, file= "ignore/00_all_env_bio_data.RData") ## env data
save(data1, file= "ignore/00_pseudo_abs.RData") ## pseudo_abs

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
# data_sf<- st_read("ignore/00_all_env_bio_data.shp")
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


# add hydro ---------------------------------------------------------------

## in RB9 repo

## currently delta ffm but can calculate current if needed

delta <- read.csv("ignore/2022-07-28_predicted_abs_FFM_deltaFFM_SD_COMIDS_medianDelta_test5_norunoff.csv")
head(delta)

## remove duplicates
delta <- delta %>% distinct()

delta_long <- delta %>%
  rename(FlowMetric = metric, AbsFlow = abs_FFM_median_cfs)
  # select(comid region, year, flow_metric, deltah_cur_ref_final, deltaH_watercon_ref_final) %>%
  # pivot_longer(d_ds_dur_ws:d_wet_tim, names_to = "FlowMetric", values_to = "DeltaH")

head(delta_long)

## keep only high R2s
# fa_mag
# wet_bfl_mag_10
# peak_10
# ds_mag_90
# peak_2
# peak_5
# ds_mag_50

unique(delta_long$FlowMetric)
unique(delta$comid)

delta_long <- delta_long %>%
  # filter(FlowMetric %in% c("d_fa_mag", "d_wet_bfl_mag_10", "d_ds_mag_90", "d_ds_mag_50")) %>%
  mutate(hydro.endpoint = case_when(FlowMetric == "ds_mag_50" ~ "DS_Mag_50",
                                    FlowMetric == "ds_mag_90" ~ "DS_Mag_90",
                                    FlowMetric == "fa_mag" ~ "FA_Mag",
                                    FlowMetric == "peak_10" ~ "Peak_10",
                                    FlowMetric == "peak_2" ~ "Peak_2",
                                    FlowMetric == "peak_5" ~ "Peak_5",
                                    FlowMetric == "sp_mag" ~ "SP_Mag",
                                    FlowMetric == "wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                    FlowMetric == "wet_bfl_mag_50" ~ "Wet_BFL_Mag_50")) 
names(delta_long)
## get median delta H 

delta_med <- delta_long %>%
  # group_by(comid, FlowMetric, hydro.endpoint) %>%
  # summarise(MedDelta = median(DeltaH)) %>%
  # ungroup() %>%
  rename(COMID = comid) %>%
  filter(!is.na(FlowMetric)) %>%
  dplyr::select(-c(FlowMetric, wayr, WYT, comid_wy, delta_FFM_median_cfs)) %>%
  pivot_wider(names_from = hydro.endpoint, values_from = AbsFlow)
  

head(delta_med)
dim(delta_med)
length(unique(delta_med$COMID))
dim(data_sf)

names(all_data)
## join spatial sites

data_sub <- data_sf %>%
  dplyr::select(c(ID:Y, COMID))

head(data_sub)

data_hyd_sf <- right_join(data_sf, delta_med, by = "COMID") ## 209 reaches don't match
names(data_hyd_sf)
head(data_sf)

length(unique(data_sub$COMID)) ## 3845
length(unique(delta_med$COMID)) ## 2116
length(unique(data_hyd_sf$COMID)) ## 2116

## scale to nhd reach
data_hyd_sf_long <- data_hyd_sf %>%
    pivot_longer(c(MRVBF.Mx:AvgSlope, DS_Mag_50:Wet_BFL_Mag_50), names_to = "Variable", values_to = "Value") %>%
    group_by(COMID, Variable) %>%
    mutate(Meanvals = mean(Value),
              Minvals = min(Value),
              Maxvals = max(Value))

    # save(data_hyd_sf_long, file="ignore/00_all_env_data_scaled.RData")
    head(data_hyd_sf_long)
    # rm(data_hyd_sf_long)
   
  unique(data_hyd_sf_long$Variable)
  ## elevation = min, catchment area = max, MRVBF = max, VRM1.Mn = min
  ## format df
  data_hyd_sf_longer <- data_hyd_sf_long %>%
    dplyr::select(-Value) %>%
    pivot_longer(Meanvals:Maxvals, names_to= "Stat", values_to = "Values") %>%
    pivot_wider(names_from = Variable, values_from = Values)

  head(data_hyd_sf_longer)
  names(data_hyd_sf_longer)
  str(data_hyd_sf_longer)
  
  ## take specific stats for each variable
  meanVars <- data_hyd_sf_longer %>%
    ungroup() %>%
    filter(Stat == "Meanvals") %>%
    dplyr::select(c(COMID,X81pptCr1:Wet_BFL_Mag_10)) %>%
    distinct()
  

  MinVars <- data_hyd_sf_longer %>%
    ungroup() %>%
    filter(Stat == "Minvals") %>%
    dplyr::select(c(COMID, DEM_10m.Mn, VRM1.Mn:VRM9.Mn))  %>%
    distinct()

  MaxVars <- data_hyd_sf_longer %>%
    ungroup() %>%
    filter(Stat == "Maxvals") %>%
    dplyr::select(c(COMID, Catchment.A, MRVBF.Mx))  %>%
    distinct()

  ## join toegther

  all_data <- bind_cols(meanVars, MinVars[,-1], MaxVars[-1])

  head(all_data)
  names(all_data)
  dim(all_data)
  ## save out
  
  save(all_data, file = "ignore/00_all_env_bio_data_NHD_reach_abs_hydro.RData")

# Raster stuff ------------------------------------------------------------


### extract other physical info

varR <- raster("ignore/00_raw_data_raster.grd")
varR

?extract

test <- raster::extract(varR, data_hyd_sf)

### make into rasters

## dims etc from CurrentGridFeb14.grd

x <- raster(ncol=701, nrow=649, xmn=494438, xmx=495438, ymn=3600402 , ymx=3600402)

projection(x) <- "+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
x

#Create a "rasterBrick" of the template raster
x<-brick(x)

x
names(data_hyd_sf)
#Create list of column names you want to rasterize
fields <- names(data_hyd_sf)[c(6:10)]
fields
## make rasters of each env var
i
for (i in fields){
  x[[i]]<-rasterize(data_hyd_sf, x, field=i)
  projection(x)<-"+proj=utm +zone=11 +datum=WGS84"
  x <- stack(x)
}

x@layers
x2[[i]]

plot(x[[2]]) ## to check
x


## save out
writeRaster(x, "ignore/00_raw_data_raster.grd", format="raster", crs="+proj=utm +zone=11 +datum=WGS84", overwrite=TRUE)


