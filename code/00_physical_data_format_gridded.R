### physical data

## data formatting

##workflow
## upload raw data
## get comids to match later
## create rasters
## add hydro
## update remote sensing
## update climate
## 200m grid
## per comid

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


# setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")
getwd()

# raw env data ------------------------------------------------------------
path <- "original_model/FullData/"

## Mike's data - env vars with presence absence
data <- read.csv(paste0(path, "200mCells_FullData_PresAbs_Complete_ThinnedCols.csv"))
head(data)

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

load(file= "ignore/00_all_env_bio_data.RData") ## data2 - from old code, run and save new
head(data2)

# Remove climate and remote sensing variables (updated below) -------------

## remove old climate data (also removed remote sensing)

data_red <- data2 %>%
  dplyr::select(-c(X81pptCr1:FINAL...24..4))

head(data_red)

## make spatial, crs 4269 to match climate polygons
data_sf<- data_red %>%
  dplyr::select(X,Y, ID.2) %>%
  st_as_sf(coords=c("X", "Y"), crs=4269, remove=F) 

head(data_sf)
names(data_sf)

st_crs(data_sf)

# Convert to raster -------------------------------------------------------

## make raster to extarct climate, hydro and remote sense data
## create template raster

## dims etc from CurrentGridFeb14.grd

x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)

projection(x) <- "+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
x

#Create a "rasterBrick" of the template raster
x<-brick(x)

x2 <- x

#Create list of column names you want to rasterize
fields <- names(data_sf) [5:17]
fields
## make rasters of each env var
i
for (i in fields){
  x[[i]]<-rasterize(data_sf, x, field=i)
  projection(x)<-"+proj=utm +zone=11 +datum=WGS84"
  x <- stack(x)
}

x@layers
# x2[[i]]

plot(x[[2]]) ## to check
x
## save out
writeRaster(x, "ignore/00_raw_original_data_raster.grd", format="raster", crs="+proj=utm +zone=11 +datum=WGS84", overwrite=TRUE)


# Format Climate data -------------------------------------------------

gridsR <- raster(ncol=129, nrow=148, xmn= -2483239, xmx=-2380039, ymn=-4822142 , ymx=-4703742)

projection(gridsR) <- "+proj=geocent +ellps=GRS80 +units=m +no_defs"
res(gridsR)=c(800,800)

gridsR ## use as template to extract point data
plot(gridsR)

#### upload annual climate data

## precipitation ###
ppt_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average ppt_30yr.shp") %>%
  rename(ppt_annx = PRISM__) 
## make spatial and transformCRS
ppt_annSP <- as(ppt_ann, Class = "Spatial")
ppt_annSP <- spTransform(ppt_annSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
ppt_annR <- rasterize(ppt_annSP, gridsR,  'ppt_annx',  na.rm =TRUE, sp = TRUE)

## temperature ###
tmax_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmax_30yr.shp") %>% 
  rename(tmax_annx = PRISM__)
## make spatial and transformCRS
tmax_annSP <- as(tmax_ann, Class = "Spatial")
tmax_annSP <- spTransform(tmax_annSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
tmax_annR <- rasterize(tmax_annSP, gridsR,  'tmax_annx',  na.rm =TRUE, sp = TRUE)

tmin_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmin_30yr.shp") %>% 
  rename(tmin_annx = PRISM__)

## make spatial and transformCRS
tmin_annSP <- as(tmin_ann, Class = "Spatial")
tmin_annSP <- spTransform(tmin_annSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
tmin_annR <- rasterize(tmin_annSP, gridsR,  'tmin_annx',  na.rm =TRUE, sp = TRUE)

## stack all together
annStack <- stack(ppt_annR, tmax_annR, tmin_annR)
## name layers
names(annStack) <- c("ppt_ann", "tmax_ann", "tmin_ann")

## monthly data
ppt_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average ppt_30yr.shp")
## make spatial and transformCRS
ppt_monSP <- as(ppt_mon, Class = "Spatial")
ppt_monSP <- spTransform(ppt_monSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))
## create raster 
ppt_monR <- rasterize(ppt_monSP, gridsR,  field = c('pptMonth','Month' ),  na.rm =TRUE, sp = TRUE)


tmax_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmax_30yr.shp")
## make spatial and transformCRS
tmax_monSP <- as(tmax_mon, Class = "Spatial")
tmax_monSP <- spTransform(tmax_monSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))

## create raster 
tmax_monR <- rasterize(tmax_monSP, gridsR,  field = c('tmaxMonth','Month' ),  na.rm =TRUE, sp = TRUE)


tmin_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmin_30yr.shp")
## make spatial and transformCRS
tmin_monSP <- as(tmin_mon, Class = "Spatial")
tmin_monSP <- spTransform(tmin_monSP, CRS("+proj=geocent +ellps=GRS80 +units=m +no_defs"))
head(tmin_monSP)
## create raster 
tmin_monR <- rasterize(tmin_monSP, gridsR,  field = c('tminMonth','Month' ),  na.rm =TRUE, sp = TRUE)
tmin_monR
plot(tmin_monR)

## stack all together
monStack <- stack(ppt_monR, tmax_monR, tmin_monR)
## name layers
names(monStack) <- c("ppt_mon", "tmax_mon", "tmin_mon")

## join altogether

climStack <- stack(annStack, monStack)
climStack@layers
plot(climStack)
## change CRS 

projection(climStack) <- "+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
crs(climStack)
## save out

writeRaster(climStack, "ignore/00_clim_raster_stack.grd", format="raster", crs="+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", overwrite=TRUE)


## get coords and comids, change crs
orig_grids <- data2 %>%
  dplyr::select(ID.2, X, Y, COMID) %>%
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) %>%
  st_transform(crs=9001)

## extract raster values at grids
clim_Orig <- raster::extract(climStack, orig_grids)
head(clim_Orig) ### all NAs!!!!


save(clim_Orig, file = "ignore/00a_new_clim_orig_env.Rdata")
load(file = "ignore/00a_new_clim_orig_env.Rdata")

# RB9 NHD -----------------------------------------------------------------

nhd <- st_read("/Users/katieirving/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/NHDPlus_V2_Flowline_CA.shp")

head(nhd)

nhd_lines_rb9 <- nhd %>%
  dplyr::select(Shape_Leng, COMID) %>%
  # filter(COMID %in% pred_df$COMID) %>%
  mutate(COMID = as.integer(COMID))

## convert rb9 comids to points
nhd_points_rb9 <- st_cast(nhd_lines_rb9, "POINT") %>% dplyr::select(-Shape_Leng) %>% st_zm()
head(nhd_points_rb9)
## check this!!!!!
# Warning message:
#   In st_cast.sf(nhd_lines_rb9, "POINT") :
#   repeating attributes for all sub-geometries for which they may not be constant

# Remote Sensing data -----------------------------------------------------

sept <- brick("/Users/katieirving/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/TC_2014_RB9/TC_092014_RB9.tif")
sept ## is 30m grids, is extract ok? Mike uses Median (Med) and Variance (Var) within analysis pixel

## extract raster values
sept_values <- raster::extract(sept, orig_grids)
head(sept_values)
sept_values_coord <- cbind(sept_values[,1:3], orig_grids)
head(sept_values_coord)
## upload wet season data
april <- brick("/Users/katieirving/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/TC_2014_RB9/TC_042014_RB9.tif")
# plot(april)

## extract raster values
apr_values <- raster::extract(april, orig_grids)
# head(apr_values)

## join with sept values
tass_sp <- cbind(apr_values[,1:3], sept_values_coord)
head(tass_sp)
dim(tass_sp)

save(tass_sp, file = "ignore/00_tass_cap_all_sp.RData")

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



# Add hydro ---------------------------------------------------------------

delta <- read.csv("/Users/katieirving/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/ignore/2022-07-28_predicted_abs_FFM_deltaFFM_SD_COMIDS_medianDelta_test5_norunoff.csv")
head(delta)

## remove duplicates
delta <- delta %>% distinct()

delta_long <- delta %>%
  # select(comid region, year, flow_metric, deltah_cur_ref_final, deltaH_watercon_ref_final) %>%
  pivot_longer(d_ds_dur_ws:d_wet_tim, names_to = "FlowMetric", values_to = "DeltaH")

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

delta_long <- delta_long %>%
  filter(FlowMetric %in% c("d_fa_mag", "d_wet_bfl_mag_10", "d_ds_mag_90", "d_ds_mag_50")) %>%
  mutate(hydro.endpoint = case_when(FlowMetric == "d_ds_mag_50" ~ "DS_Mag_50",
                                    FlowMetric == "d_ds_mag_90" ~ "DS_Mag_90",
                                    FlowMetric == "d_fa_mag" ~ "FA_Mag",
                                    # FlowMetric == "d_peak_10" ~ "DS_Mag_50",
                                    FlowMetric == "d_peak_2" ~ "Peak_2",
                                    # FlowMetric == "d_peak_5" ~ "DS_Mag_50",
                                    # FlowMetric == "sp_mag" ~ "SP_Mag",
                                    FlowMetric == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                    FlowMetric == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50")) 
names(delta_long)
## get median delta H 

delta_med <- delta_long %>%
  group_by(comid, FlowMetric, hydro.endpoint) %>%
  summarise(MedDelta = median(DeltaH)) %>%
  ungroup() %>%
  rename(COMID = comid) %>%
  dplyr::select(-FlowMetric) %>%
  pivot_wider(names_from = hydro.endpoint, values_from = MedDelta)


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
  pivot_longer(c(MRVBF.Mx:AvgSlope, DS_Mag_50:Wet_BFL_Mag_10), names_to = "Variable", values_to = "Value") %>%
  group_by(COMID, Variable) %>%
  mutate(Meanvals = mean(Value),
         Minvals = min(Value),
         Maxvals = max(Value))

# save(data_hyd_sf_long, file="ignore/00_all_env_data_scaled.RData")
head(data_hyd_sf_long)
rm(data_hyd_sf_long)

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

## save out

save(all_data, file = "ignore/00_all_env_bio_data_NHD_reach.RData")

