## precipitation
ppt_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average ppt_30yr.shp") %>%
  rename(ppt_ann = PRISM__) 
## convert to spatial
ppt_annSP <- as(ppt_ann, Class = "Spatial")
ppt_annSP

ppt_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average ppt_30yr.shp")
## convert to spatial
ppt_monSP <- as(ppt_mon, Class = "Spatial")
## temperature
tmax_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmax_30yr.shp") %>% 
  rename(tmax_ann = PRISM__)
## convert to spatial
tmax_annSP <- as(tmax_ann, Class = "Spatial")
tmin_ann <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_annual average tmin_30yr.shp") %>% 
  rename(tmin_ann = PRISM__)
## convert to spatial
tmin_annSP <- as(tmin_ann, Class = "Spatial")
tmax_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmax_30yr.shp")
## convert to spatial
tmax_monSP <- as(tmax_mon, Class = "Spatial")
tmin_mon <- st_read("/Users/katieirving/SCCWRP/PRISM - General/Data/prism_monthly average tmin_30yr.shp")
## convert to spatial
tmin_monSP <- as(tmin_mon, Class = "Spatial")

## create empty raster - match extent, CRS and make 800m grids
gridsR <- raster()
extent(gridsR) <- extent(ppt_annSP)
projection(gridsR) <- "+proj=longlat +datum=NAD83 +no_defs"
res(gridsR) <- c(800,800)

gridsR ## use as template to make rasters from clim data

## create raster for each clim variable
ppt_annR <- rasterize(ppt_annSP, gridsR,  'ppt_ann',  na.rm =TRUE, sp = TRUE)
ppt_annR
ppt_monR <- rasterize(ppt_monSP, gridsR,  field = c('ppt_mon', 'month'),  na.rm =TRUE, sp = TRUE)

tmax_annR <- rasterize(tmax_annSP, gridsR,  'tmax_ann',  na.rm =TRUE, sp = TRUE)
tmax_monR <- rasterize(tmax_monSP, gridsR,  'tmax_mon',  na.rm =TRUE, sp = TRUE)

tmin_annR <- rasterize(tmin_annSP, gridsR,  'tmin_ann',  na.rm =TRUE, sp = TRUE)
tmin_monR <- rasterize(tmin_monSP, gridsR,  'tmin_mon',  na.rm =TRUE, sp = TRUE)

## make brick to extract values
