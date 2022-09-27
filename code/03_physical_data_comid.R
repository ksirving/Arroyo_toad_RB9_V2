### data per comid

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








