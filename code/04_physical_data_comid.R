### data per comid

# packages
library(tidyverse)
library(tidylog)
library(sf)
library(raster)
library(nhdplusTools)
library(readr)
library(mapview)

# devtools::install_github('ceff-tech/ffc_api_client/ffcAPIClient')
# library(ffcAPIClient)
# 
# get_comid_for_lon_lat(-117.3211, 33.38939, online = F)

## upload model and prediction gridded data

load(file =  "ignore/03_RB9_grdded_data.RData") # for predictions
load(file =  "ignore/03_RB9_grdded_data_observations.RData") # for model build

head(data_hyd_sf_obs)
head(data_hyd_sf)

# 22549169, 20348307 have same coord and data values but different hydro value

names(data_hyd_sf_obs)
## separate to bio and env

bioMod <- data_hyd_sf_obs %>%
  as.data.frame() %>%
  dplyr::select(COMID, ID, Year:PresAbs) ## 
  
envMod <- data_hyd_sf_obs %>%
  as.data.frame() %>% 
  dplyr::select(COMID, ID.1, MRVBF.Mx:TC_042014_RB9.3_Var, DS_Mag_50:Wet_BFL_Mag_50) 

envPred <- data_hyd_sf %>% as.data.frame() 


# Physical Data -----------------------------------------------------------

#### summarise by comid - mean, min, max and variance

## make longer and group by comid and variable, just remote sense data first
## median values
RemSense_longMed <- envPred %>%
  dplyr::select(COMID, ends_with("Med")) %>%
  pivot_longer(c(TC_092014_RB9.1_Med:TC_092014_RB9.3_Med, TC_042014_RB9.1_Med:TC_042014_RB9.3_Med), names_to = "Variable", values_to = "Value") %>%
  group_by(COMID, Variable) %>% 
  summarise(Vals = median(Value)) #%>% ## median values

# sum(is.na(RemSense_longMed))
# ind1 <- which(is.na(RemSense_longMed))
# RemSense_longMed[ind1,]

## variance
RemSense_longVar <- envPred %>%
  dplyr::select(COMID, ends_with("Var")) %>%
  pivot_longer(c(TC_092014_RB9.1_Var:TC_092014_RB9.3_Var, TC_042014_RB9.1_Var:TC_042014_RB9.3_Var), names_to = "Variable", values_to = "Value") %>%
  group_by(COMID, Variable) %>% 
  mutate(Vals = var(Value)) %>% ## variance values
  mutate(Vals = ifelse(is.na(Vals), Value, Vals)) %>%
  dplyr::select(COMID, Variable, Vals) %>%
  distinct()
## variance not always calculated - when only 1 comid, in those cases, take original value

# sum(is.na(RemSense_longVar)) ## 6 missing comids
# ind <- which(is.na(RemSense_longVar))
# che <- RemSense_longVar[ind,]

## join together

tass_all <- bind_rows(RemSense_longMed, RemSense_longVar)

## make wider rename with mean/var in name
## means
# tass_sp_means <- RemSense_long %>%
#   mutate(VarNames = paste0(Variable, "_Var"), Variable = paste0(Variable, "_Mean")) %>%
#   pivot_wider(id_cols = "COMID", names_from = Variable, values_from = MeanVals) 
# 
# ## variance
# tass_sp_vars <- RemSense_long %>%
#   mutate(VarNames = paste0(Variable, "_Var"), Variable = paste0(Variable, "_Mean")) %>%
#   pivot_wider(id_cols = "COMID", names_from = VarNames, values_from = VarVals) 

## join together

# tass_all <- full_join(tass_sp_vars, tass_sp_means, by = "COMID" )
# 
head(tass_all)

tass_all <- tass_all %>% 
  pivot_wider(names_from = Variable, values_from = Vals)

save(tass_all, file="ignore/03_remote_sensing_data_per_comid.RData")

## other variables
names(envPred)
## scale to nhd reach 
data_hyd_sf_long <- envPred %>%
  dplyr::select(-c(TC_092014_RB9.1_Med:TC_042014_RB9.3_Var,  DS_Mag_50:Wet_BFL_Mag_50)) %>%
  pivot_longer(c(MRVBF.Mx:tminMonX12), names_to = "Variable", values_to = "Value") %>%
  group_by(COMID, Variable) %>%
  summarise(Meanvals = mean(Value),
         Minvals = min(Value),
         Maxvals = max(Value)) %>%
  distinct()

head(data_hyd_sf_long)
# rm(data_hyd_sf_long)

unique(data_hyd_sf_long$Variable)
## elevation = min, catchment area = max, MRVBF = max, VRM1.Mn = min
## reshape to make sum stats longformat df
data_hyd_sf_longer <- data_hyd_sf_long %>%
  pivot_longer(Meanvals:Maxvals, names_to= "Stat", values_to = "Values") %>%
  pivot_wider(names_from = Variable, values_from = Values)

head(data_hyd_sf_longer)
names(data_hyd_sf_longer)
str(data_hyd_sf_longer)

## take specific stats for each variable
meanVars <- data_hyd_sf_longer %>%
  ungroup() %>%
  filter(Stat == "Meanvals") %>%
  dplyr::select(c(COMID, ppt_ann:pptMonX9, AvgClay:AvgWaterSt)) %>%
  distinct()


MinVars <- data_hyd_sf_longer %>%
  ungroup() %>%
  filter(Stat == "Minvals") %>%
  dplyr::select(c(COMID, DEM_10m.Mn, VRM1.Mn:VRM9.Mn, tmin_ann:tminMonX9))  %>%
  distinct()

MaxVars <- data_hyd_sf_longer %>%
  ungroup() %>%
  filter(Stat == "Maxvals") %>%
  dplyr::select(c(COMID, Catchment.A, MRVBF.Mx, tmax_ann:tmaxMonX9))  %>%
  distinct()

## join dfs
all_data <- bind_cols(meanVars, MinVars[,-1], MaxVars[-1])

## hydro

names(envPred) 

hydro <- envPred %>%
  dplyr::select(COMID, DS_Mag_50:Wet_BFL_Mag_50) %>%
  distinct(COMID, .keep_all = T)

head(hydro)

## join together

hyd_tass <- full_join(hydro, tass_all, by = "COMID")
allData <- full_join(hyd_tass, all_data, by = "COMID")

head(allData)
names(allData)

# allData <- distinct(allData)

## save out

save(allData, file = "ignore/04_all_env_data_NHD_COMID.RData")

# Join presence/absence to env data (COMID) ------------------------------

names(bioMod)
dim(bioMod)

unique(bioMod$LifeStage)
## make sure if reach has at least one presence then keeps it, if none then absence
bio_coms <- bioMod %>%
  as.data.frame() %>%
  group_by(COMID) %>% 
  summarise(PresAbs = max(PresAbs)) %>%
  # rename(COMID = COMID2) %>%
  drop_na() 


head(bio_coms)
dim(bio_coms) ## 318

## how many presense/absences
sum(bio_coms$PresAbs ==1) ## 235
sum(bio_coms$PresAbs ==0) ## 83

head(allData)
## join with env data
NewDataObsComid <- full_join(allData, bio_coms, by = "COMID")

## replace NAs in presabs with -999

NewDataObsComid$PresAbs[is.na(NewDataObsComid$PresAbs)] <- -999


## saveout
save(NewDataObsComid, file = "ignore/04_all_env_bio_data_NHD_COMID.RData")


# plot comids -------------------------------------------------------------

# nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHDplus_RB9.shp")
# # nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/SD_RB9_boundary.shp")
# 
# head(nhd)
# obs <- bioAll %>%
#   dplyr::select(COMID)
# 
# nhd_lines_rb9 <- nhd %>%
#   dplyr::select(SHAPE_LENG, COMID) %>%
#   # filter(COMID %in% all_data$COMID) %>% ## subset to orginial data with RB9 extent
#   mutate(COMID = as.integer(COMID))
# 
# nhd_lines_rb9
# 
# dim(nhd_lines_rb9)
# str(nhd_lines_rb9)
# 
# # nhd_lines_prob <- full_join(nhd_lines_rb9, pred_df, by = "COMID")
# nhd_lines_rb9 <- st_zm(nhd_lines_rb9)
# 
# ## join in obs
# 
# nhd_lines_obs <- right_join(nhd_lines_rb9, bioAll, by = "COMID")
# nhd_lines_obs <- st_zm(nhd_lines_obs) %>% mutate(nhd_lines_obs, PresAbs = as.factor(PresAbs)) #%>%
# drop_na(Shape_Leng)
# dim(nhd_lines_obs)
# 
# 
# # library(viridis)
# map1 <- ggplot() +
#   geom_sf(data = nhd_lines_rb9) +
#   geom_sf(data = nhd_lines_obs, aes(colour =PresAbs))
# # scale_fill_gradientn(colours=rev(magma(6))) ## colours not working
# 
# map1
# 
# file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/01_rb9_map.jpg"
# # ggsave(map1, filename=file.name1, dpi=300, height=5, width=6)
# 
# 
# 
# load(file= "ignore/00_all_env_bio_data.RData")
# ## make spatial
# 
# data_sf<- data2 %>%
#   dplyr::select(-PresAbs2005, -Presence2) %>%
#   st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) 
# 
# head(data_sf)
# 
# ## map
# 
# # set background basemaps:
# # basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
# #                   "Esri.NatGeoWorldMap",
# #                   "OpenTopoMap", "OpenStreetMap", 
# #                   "CartoDB.Positron", "Stamen.TopOSMFeatures")
# # 
# # mapviewOptions(basemaps=basemapsList, fgb = FALSE)
# # 
# # m1 <- mapview(data_segs, cex=6, col.regions="orange",
# #               layer.name="data points") 
# # 
# # m1
# # m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")
# 
# ## save as shape
# 
# st_write(data_sf, "ignore/00_all_env_bio_data.shp", append=FALSE) ## not working??
# # 
# # check <- st_read("output_data/00_all_env_bio_data.shp")
# # plot(check)()
# data_sf<- st_read("ignore/00_all_env_bio_data.shp")
# data_sf
# dim(data_sf)
# class(data_sf)
# 
# ## make long and get mean value per comid
# 
# # data_hyd_sf_long <- data_sf %>%
# #   pivot_longer(MRVBF.Mx:AvgSlope, names_to = "Variable", values_to = "Value") %>%
# #   group_by(COMID, Variable) %>%
# #   mutate(Meanvals = mean(Value),
# #             Minvals = min(Value),
# #             Maxvals = max(Value))
# #   
# #   save(data_hyd_sf_long, file="ignore/00_all_env_data_scaled.RData")
# #   head(data_hyd_sf_long)
# #   rm(data_hyd_sf_long)
# #   
# #   unique(data_hyd_sf_long$Variable)
# #   ## elevation = min, catchment area = max, MRVBF = max, VRM1.Mn = min
# #   
# #   ## format df 
# #   data_hyd_sf_longer <- data_hyd_sf_long %>%
# #     pivot_longer(Meanvals:Maxvals, names_to= "Stat", values_to = "Values") %>%
# #     pivot_wider(names_from = Variable, values_from = Values)
# #   
# #   head(data_hyd_sf_longer)
# #   names(data_hyd_sf_longer)
# #   ## take specific stats for each variable
# #   meanVars <- data_hyd_sf_longer %>%
# #     filter(Stat == "Meanvals") %>%
# #     dplyr::select(c(COMID:AvgWaterSt, FINAL...01..3:FINAL...24..4, X81pptCr1:X81TMxCr9))
# #   
# #   
# #   MinVars <- data_hyd_sf_longer %>%
# #     filter(Stat == "Minvals") %>%
# #     dplyr::select(c(COMID, DEM_10m.Mn, VRM1.Mn:VRM9.Mn))
# #   
# #   MaxVars <- data_hyd_sf_longer %>%
# #     filter(Stat == "Maxvals") %>%
# #     dplyr::select(c(COMID, Catchment.A, MRVBF.Mx))
# #   
# #   ## join toegther
# #   
# #   all_data <- bind_cols(meanVars, MinVars, MaxVars)
# #   
# #   head(all_data)
# #   names(all_data)
# #   
# #   all_data <- all_data %>%
# #     dplyr::select(-COMID...92, -COMID...99, -Stat) %>%
# #     rename(COMID = COMID...1)
# 
# # Convert to raster -------------------------------------------------------
# 
# ## create template raster
# 
# ## dims etc from CurrentGridFeb14.grd
# 
# 
# x <- raster(ncol=701, nrow=649, xmn=423638.013766974, xmx=563838.013766974, ymn=3600402.14370233 , ymx=3730202.14370233)
# 
# projection(x) <- "+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# x
# 
# #Create a "rasterBrick" of the template raster
# x<-brick(x)
# 
# x2 <- x
# 
# #Create list of column names you want to rasterize
# fields <- names(data_sf) [5:100]
# fields
# ## make rasters of each env var
# i
# for (i in fields){
#   x[[i]]<-rasterize(data_sf, x, field=i)
#   projection(x)<-"+proj=utm +zone=11 +datum=WGS84"
#   x <- stack(x)
# }
# 
# x@layers
# x2[[i]]
# 
# plot(x[[2]]) ## to check
# x
# ## save out
# writeRaster(x, "ignore/00_raw_data_raster.grd", format="raster", crs="+proj=utm +zone=11 +datum=WGS84", overwrite=TRUE)
# 
# 
# 


