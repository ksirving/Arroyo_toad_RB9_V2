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
