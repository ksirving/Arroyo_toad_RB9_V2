### automated gridded model

##packages

library(sp)
library(rgdal)
library(raster)
library(randomForest)
library(kernlab)
library(rgl)
library(ks)
library(sm)
library(sf)
library(tidyverse)
library(mapview)
library(caret)
library(rpart)

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## get path for functions
source("original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

# Upload bio data -------------------------------------------------------------

## updated and snapped pres/abs
bioSnap <- shapefile("/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/output_data/ToadsObs_50m_Final.shp")
head(bioSnap)

orig.sdata <-  bioSnap %>% as.data.frame() %>% ## get coordinates of snapped points
  dplyr::select(ID, COMID, Year:PresAbs, centroid_x, centroid_y) %>%
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs=4326, remove=F) 


# Upload physical data ----------------------------------------------------
## comid env df

load(file = "ignore/04_all_env_data_NHD_COMID.RData") 
head(allData)

## gridded env df with comids
load(file = "ignore/00_RB9_grdded_data.RData") 
head(data_hyd_sf)

## read in env data stack
xvars <- stack("ignore/00_raw_final_data_raster.tif") ## new rasters @ 200m 

## upload and add layer names 
load(file = "output_data/00_final_raster_layer_names.RData")
names(xvars) <- layerNames

## upload old raster mask for crs
rmaskOld <- raster("ignore/02_mask_raster_network.tif")

## use this one for now
mask <- raster("ignore/02_mask_raster_network_new.tif")

## change crs to match mask
crs(xvars) <- crs(rmaskOld)

## format raster - use this one when points are snapped
# mask <- xvars[[65]]
# crs(mask)

## get cell numbers at pres/abs sites 
cellsPres1 <- raster::extract(mask, orig.sdata, cellnumbers=TRUE)
# dim(cellsPres) ## 2797
head(cellsPres1) ## nas here related to final raster, check with new snap
dim(cellsPres1)

## no nas are created with the "new" raster, but are with the "final" raster - fix this issue

# Begin loop here ---------------------------------------------------------

## format loop stuff

models <- paste0("model",seq(1, 10,1))
m=2

for(m in models) {
   
  gridFile <- paste0("ignore/ModelResults/Gridded/", models[m],"/")
  # COMIDFile <- paste0("ignore/ModelResults/COMID/", models[m],"/")


# KDE Bias Surface --------------------------------------------------------

# develop KDE sampling bias surface
orig.sdata2<-subset(orig.sdata, PresAbs==1)

## get cell numbers for bias
# bias <- cellFromXY(mask, orig.sdata[,-1])
bias <- cellsPres1[,1]
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
### check here for pseudo absences - should be ~50/50
## format to bind to other data 
a.spdfx <- a.spdf %>% as.data.frame() %>%
  st_as_sf(coords=c("x",  "y"), crs =32611, remove=F) %>%
  mutate(Longitude = unlist(map(geometry,1)),
         Latitude = unlist(map(geometry,2))) %>%
  mutate(ID = 1:length(geometry)) %>%
  mutate(ID = paste0("P", ID)) %>% 
  mutate(Year = NA, LifeStage = "Pseudo", Count = NA) %>%
  mutate(centroid_x = Longitude, centroid_y = Latitude, Join_Count = NA, Id_1 = NA) %>%
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs=32611) %>%
  st_transform(crs = 4236) %>% 
  mutate(centroid_x = unlist(map(geometry,1)),
         centroid_y = unlist(map(geometry,2)))


##  make sdata into sf object
sdata <- bioSnap %>%
  as.data.frame() %>%
  dplyr::select(-coords_x1, -coords_x2) %>%
  st_as_sf(coords=c("centroid_x", "centroid_y"), crs = 4326, remove=F) %>%
  mutate(COMID = as.integer(COMID)) 

## join with other observations
sdata1<-bind_rows(sdata,a.spdfx)
sdata1 <- sdata1 %>% dplyr::select(ID, COMID, Year:PresAbs, centroid_x, centroid_y) #%>%

save(sdata1, file=paste0(gridFile, "all_presAbs.RData"))

# Add to env data ---------------------------------------------------------

##  change NAs in lifestage so csn drop all nas from DF
bioData <- sdata1 %>% 
  mutate(LifeStage = ifelse(is.na(LifeStage), "Other", LifeStage), 
         Count = ifelse(is.na(Count), 0, Count ))

## remove z dimension
bioData <- st_zm(bioData)

## get raster values at all presence cells
cellsPres <- raster::extract(xvars, bioData,  cellnumbers = T, df=TRUE)
cellsPres ## nas here too, should be fixed with new snap

## join together biodata and cells with presences, remove comid
NewDataObsSub <- cbind(bioData[,-2], cellsPres) 
names(NewDataObsSub)

# dim(NewDataObsSub) # 3579
# length(unique(NewDataObsSub$cells)) ## 1533
# head(NewDataObsSub)

NewDataObsSub <- na.omit(NewDataObsSub)

## join pres/abs with comids from env data - join by grid cell number

COMs <- as.data.frame(data_hyd_sf) %>% dplyr::select(COMID, cells) 

NewDataObsSub <- inner_join(NewDataObsSub, COMs, by = "cells") #%>%
  # distinct(cells, .keep_all=T) 

save(NewDataObsSub, file=paste0(gridFile, "all_presAbs_env_data.RData"))

# Multicolinearality ------------------------------------------------------

## gridded
all_data_obs <- NewDataObsSub %>%
  as.data.frame() %>% 
  dplyr::select(COMID, PresAbs,-c(ppt_ann:tminMonX12, -geometry), MRVBF.Mx:AvgSlope, TC_092014_RB9.1_Med:TC_042014_RB9.3_Var, DS_Mag_50:Wet_BFL_Mag_50) %>%
  drop_na()

cl <- MultiColinear(all_data_obs[,c(3:ncol(all_data_obs))], p=0.05)
xdata <- all_data_obs[,c(3:ncol(all_data_obs))]

for(l in cl) {
  cl.test <- xdata[,-which(names(xdata)==l)]
  print(paste("REMOVE VARIABLE", l, sep=": "))
  MultiColinear(cl.test, p=0.05)
}


# # REMOVE MULTI-COLINEAR VARIABLES
for(l in cl) { all_data_obs <- all_data_obs[,-which(names(all_data_obs)==l)] }


# Random forests ----------------------------------------------------------

ydata <- factor(all_data_obs$PresAbs, levels = c(1,0))
xdata <- all_data_obs[,c(3:ncol(all_data_obs))] 

# class(ydata)
# length(ydata)

# PERCENT OF PRESENCE OBSERVATIONS
( dim(all_data_obs[all_data_obs$PresAbs == 1, ])[1] / dim(all_data_obs)[1] ) * 100 ## 93%

# RUN RANDOM FORESTS MODEL SELECTION FUNCTION
#Also provides variable importance and such
#It is important to look at highest class error - Look at Global OOB errors
#TEST Object is super important; Row number corresponds to paramerers, in output
( rf.model <- rf.modelSel(x=xdata, y=ydata, imp.scale="mir", ntree=b, nodesize=5) ) 

# CREATE NEW XDATA BASED ON SELECTED MODEL AND RUN FINAL RF MODEL (Model 4)
#RF runs differently when you use symbolis languate (using ~ as in an Lin. Model... use the indexing approacy [y=rf.data[,1]...]

sel.vars <- rf.model$PARAMETERS[[1]]# set to use 1 - lowest error rate and has all hydro vars

rf.data <- data.frame(y=ydata, xdata[,sel.vars])	

save(rf.data, file=paste0(gridFile, "data_for_rf_model.RData"))

(rf.final <- randomForest(y=rf.data[,1], x=rf.data[,2:ncol(rf.data)], ntree=b, mtry = 18,nodesize=5,
                          importance=TRUE, norm.votes=TRUE, proximity=TRUE) )

save(rf.final, file=paste0(gridFile, "rf_model.RData"))

# Proximity ------------------------------------------------------------------

### MDS scaling from proximity

# CREATE CMD SCALED DISTANCES OF RF PROXMITIES
distanceMatrix <- dist(1 - rf.final$proximity)

rf.cmd <- cmdscale(distanceMatrix, eig=TRUE, k=4, x.ret=T)  

mdsVarPer <- round(rf.cmd$eig/sum(rf.cmd$eig)*100, 1)

mdsValues <- rf.cmd$points

mdsData <- data.frame(Sample = rownames(mdsValues),
                      X=mdsValues[,1],
                      Y=mdsValues[,2],
                      Status = rf.data$y)

m1 <- ggplot(data = mdsData, aes(x=X, y=Y)) +
  geom_point(aes(color=Status)) +
  theme_bw()+
  xlab(paste("MDS1 - ", mdsVarPer[1], "%", sep=" ")) +
  ylab(paste("MDS2 - ", mdsVarPer[2], "%", sep=" ")) +
  ggtitle(paste("Pres/Abs", "PROXIMITY MATRIX", sep=" - "))


file.name1 <- paste0(gridFile, "no_clim_model_mds_scale_gridded.jpg") 
ggsave(m1, filename=file.name1, dpi=300, height=5, width=6)


# Validation --------------------------------------------------------------

## make y data compatible
rf.data.val <- rf.data %>%
  mutate(y = ifelse(y==1, "Present", "Absent")) %>%
  mutate(y = factor(y, levels = c("Present", "Absent")))

## split into training and testing

setsize <- floor(nrow(rf.data)*0.8)
index <- sample(1:nrow(rf.data), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

trcontrol = trainControl(method='cv', number=10, savePredictions = T,
                         classProbs = TRUE,summaryFunction = twoClassSummary,returnResamp="all")

model = train(y ~ . , data=training, method = "rf", trControl = trcontrol,metric="ROC") 

conMat <- confusionMatrix(predict(model,testing),testing$y)

save(model, file=paste0(gridFile, "validation.RData"))
save(conMat, file=paste0(gridFile, "confusion_matrix.RData"))


pdf(paste0(gridFile, "Trees_no_clim_gridded.pdf"), width=25, height=15)

full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data)
plot(full_tree)
text(full_tree, use.n = T)

dev.off()


# Coeficients and importance ----------------------------------------------

# PLOT VARIABLE IMPORTANCE
pdf(paste0(gridFile, "var_imp_no_clim_gridded.pdf"), width=25, height=15)

varImpPlot(rf.final)

dev.off()

## scaled var imp
p <- as.matrix(rf.final$importance)   
p
ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]]) 

pdf(paste0(gridFile, "var_imp_scaled_no_clim_gridded.pdf"), width=25, height=15)
dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)
dev.off()

## partial dependence plots
pdf(paste0(gridFile, "PartialPlots_6_no_clim_gridded.pdf"), width=8, height=8)
p <- as.matrix(rf.final$importance)    
ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]])  
dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)
plot(rf.final, main="Bootstrap Error Convergence")

for (i in 1:length(names(rf.data[,2:ncol(rf.data)])) ) {
  p <- partialPlot(rf.final, rf.data[,2:ncol(rf.data)], 
                   names(rf.data[,2:ncol(rf.data)])[i], which.class="1", plot=FALSE)   
  p$y <- (p$y - min(p$y)) / (max(p$y) - min(p$y)) 
  plot( y=lowess(y=p$y, x=p$x, f=sp)$y, x=p$x, type="l", ylab="p",  
        xlab=names(rf.data[,2:ncol(rf.data)])[i],
        main=paste("PARTIAL PLOT", names(rf.data[,2:ncol(rf.data)])[i], sep=" - ")) 
}
dev.off() 

} # end loop

# Tuning ------------------------------------------------------------------
## check tuning here, if error can be improved through mtry, rerun models 

##  check error rates

rf.final$err.rate[,1]

oob.error.data <- data.frame(
  Trees = rep(1:nrow(rf.final$err.rate), times = 3),
  Type = rep(c("OOB", "0", "1"), each = nrow(rf.final$err.rate)),
  Error = c(rf.final$err.rate[,"OOB"],
            rf.final$err.rate[, "0"],
            rf.final$err.rate[,"1"])
)


t1 <- ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))
t1
file.name2 <- paste0(gridFile, "oob_error_gridded.jpg") 
ggsave(t1, filename=file.name2, dpi=300, height=5, width=6)

## check split

oob.values <- vector(length=10)
oob.values

for(i in 1:10) {
  temp.model <- randomForest(y~., data=rf.data, mtry = i, ntree=b)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate), 1]
  
}

oob.values
rf.final$mtry
mean(rf.final$err.rate[,1])


# Predict on all rb9 region-------------------------------------------------------

#Index refers to the right column of probabilities - in this model the second column, which is probs of "1"
all_data <- na.omit(NewDataObs)
# all_data <- all_data[,c("COMID", sel.vars)] 
head(all_data)
str(all_data)

pred <- predict(rf.final, all_data,filename="output_data/Current/Model1/SppProbs_no_clim_gridded.img", type="prob",  index=2, 
                na.rm=TRUE, overwrite=TRUE, progress="window")


pred_df <- as.data.frame(predict(rf.final, all_data,filename="output_data/Current/Model1/SppProbs_no_clim_gridded.img", type="prob",  index=2, 
                                 na.rm=TRUE, overwrite=TRUE, progress="window"))
## add comids
pred_df$COMID <- all_data$COMID

head(pred_df)

## get probability, known occs and env data all together 

# obs <- NewDataObsSub 

pred_env <- pred_df %>%
  full_join(all_data, by = "COMID") %>%
  rename(probOcc = 1) %>%
  dplyr::select(probOcc, COMID, NewObs, DS_Mag_50:Wet_BFL_Mag_10)

write.csv(pred_env, "06_hydro_obs_preds_no_clim_model_gridded.csv")

head(pred_env)



# Response curves/relationships ------------------------------------------------

## change in delta

## test change in delta curves with absolute values
## may need to be run with actual values, no delta

## change to absolute
rf.data.ch <- abs(rf.data[, -1]) %>%
  mutate(Y = rf.data$Y)
names(rf.data.ch)
## define ffm 

i=3
m=2

ffm <- names(rf.data.ch)[11:16]
ffm

## observation data with comids
all_data_obs <- na.omit(all_data_obs)
length(all_data_obs$COMID)

dim(rf.data)
names(rf.data)

## get means of all variables

rf.data.mean <- rf.data %>%
  pivot_longer(AvgClay:SP_Mag, names_to = "Variable", values_to="Value") %>%
  group_by(Variable) %>%
  summarise(MeanVal = mean(Value))

## make data long - raw values of ffm

rf.data.long <- rf.data %>%
  pivot_longer(AvgClay:SP_Mag, names_to = "Variable", values_to="Value")

rf.data.long

full_data <- NULL

for(m in 1:length(ffm)) {
  
  metric <- ffm[m]
  metric
  
  metricVals <- rf.data %>%
    pivot_longer(AvgClay:SP_Mag, names_to = "Variable", values_to="Value") %>%
    filter(Variable %in% metric)
  
  head(metricVals)
  
  ## get sequence to predict on
  incr <- seq(min(metricVals$Value), max(metricVals$Value), (max(metricVals$Value)/20))
  incr <- rev(incr)
  
  datax <- data.frame(incr)
  datax
  
  pred_df_incrx <- NULL
  
  ## replace mean of ffm with increment and predi using rf model
  
  for(i in 1: length(incr)) {
    
    
    # data <- rf.data.mean %>%
    #   mutate(MeanVal = ifelse(Variable == metric, incr[i] MeanVal)) %>%
    #   pivot_wider(names_from = Variable, values_from = MeanVal)
    
    ## replace metric values with increment
    data <- rf.data %>%
      select(-paste(metric)) %>%
      mutate(increment = incr[i])
    names(data)
    
    ## change name for prediction
    colnames(data)[17] <- metric
    
    ## predict on increments one at a time
    pred_df_incr <- as.data.frame(predict(rf.final, data, type="prob",  index=2, 
                                          na.rm=TRUE, overwrite=TRUE, progress="window"))
    
    pred_df_incr
    pred_df_incr[,3] <- incr[i]
    pred_df_incr[,4] <- i
    
    pred_df_incrx <- bind_rows(pred_df_incrx, pred_df_incr)
    
    # colnames(pred_df_incr) [i+2] <- paste0("Increment_", i)
    # datax[i,2] <- pred_df_incr[2] ## proability of occurrence
    # datax[,4] <- metric ## metric name
    
    
  }
  
  
  
  pred_df_incrx$COMID <- all_data_obs$COMID
  pred_df_incrx$FFM <- metric
  
  ## change names and combine
  colnames(pred_df_incrx)[1:4] <- c("ProbAbs", "ProbPres", "IncrementValue", "IncrementNumber")
  
  full_data <- bind_rows(full_data, pred_df_incrx)
  
}


head(full_data)
unique(full_data$FFM)

coms <- unique(full_data$COMID)[c(1,4,46,87,245)]

full_data_sub <- full_data %>%
  filter(COMID %in% coms)


## plot

p1 <- ggplot(full_data_sub, aes(y= ProbPres, x = IncrementValue)) +
  geom_path() +
  facet_grid(rows = vars(COMID), cols = vars(FFM), scales = "free_x")

p1

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/06_incremental_preds_sep_coms.jpg"
ggsave(p1, filename=file.name1, dpi=300, height=5, width=6)

# p1 <- ggplot(all_data, aes(y= ProbOcc, x = increment)) +
#   geom_smooth() +
#   facet_wrap(~FFM, scales = "free_x")
# 
# p1
# 
# file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_v2/Figures/02_incremental_preds.jpg"
# ggsave(p1, filename=file.name1, dpi=300, height=5, width=6)


# Plotting probability of occurrence --------------------------------------

nhd <- st_read("/Users/katieirving/SCCWRP/SD Hydro Vulnerability Assessment - General/Data/SpatialData/NHD_reaches_RB9_castreamclassification.shp")
head(pred_df)

head(nhd)

names(all_data)

obs <- all_data_obs %>%
  dplyr::select(COMID, NewObs)

obs

nhd_lines_rb9 <- nhd %>%
  dplyr::select(Shape_Leng, COMID) %>%
  filter(COMID %in% pred_df$COMID) %>%
  mutate(COMID = as.integer(COMID))

head(nhd_lines_prob)
dim(nhd_lines_rb9)
str(nhd_lines_rb9)

nhd_lines_prob <- full_join(nhd_lines_rb9, pred_df, by = "COMID")
nhd_lines_prob <- st_zm(nhd_lines_prob)

# plot(nhd_lines_prob[5])

### round prob values - change later
nhd_lines_prob <- nhd_lines_prob %>%
  rename(probOcc = 3) %>%
  mutate(probRound = round(probOcc, digits=1))

## join in obs

nhd_lines_obs <- right_join(nhd_lines_rb9, obs, by = "COMID")
nhd_lines_obs <- st_zm(nhd_lines_obs) %>% mutate(nhd_lines_obs, NewObs = as.factor(NewObs)) %>%
  drop_na(Shape_Leng) %>% st_centroid()

nhd_lines_obs
nhd_lines_prob

library(viridis)

my_breaks <- c(0, 0.25, 0.5, 0.75, 1)

map1 <- ggplot() +
  geom_sf(data = nhd_lines_prob, aes(color = probOcc)) + 
  scale_colour_gradientn(colours=inferno(6), name="Probability of Occurrence",
                         labels = my_breaks) +
  geom_sf(data = subset(nhd_lines_obs, NewObs == 1), shape = 1, size = 2) 


map1

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/01_full_model_prob_occs_map_no_clim.jpg"
ggsave(map1, filename=file.name1, dpi=300, height=5, width=8)

library(mapview)
library(sf)
library(RColorBrewer)
webshot::install_phantomjs()

## map
# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, vector.palette = colorRampPalette(c(  "red", "green")) , fgb = FALSE)


m1 <- mapview(nhd_lines_prob, zcol = "probOcc",  legend = TRUE, layer.name = "Probability of Occurrence") +
  mapview(subset(nhd_lines_obs, NewObs == 1), col.regions = "black",cex = 2, layer.name = "Observations")

m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")

mapshot(m1, url = paste0(getwd(), "/ignore/map_no_clim_gridded.html"),
        file = paste0(getwd(), "/ignore/map_no_clim_gridded.png"))
getwd()



