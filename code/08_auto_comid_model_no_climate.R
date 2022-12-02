### models on comids

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

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

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
m=1

for(m in models) {
  
  # gridFile <- paste0("ignore/ModelResults/Gridded/", models[m],"/")
  COMIDFile <- paste0("ignore/ModelResults/COMID/", models[m],"/")

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

save(sdata1, file=paste0(COMIDFile, "all_presAbs.RData"))

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

save(NewDataObsSub, file=paste0(COMIDFile, "all_presAbs_env_data.RData"))


# Calculate COMID values --------------------------------------------------

bio_coms <- NewDataObsSub %>%
  as.data.frame() %>%
  group_by(COMID) %>%
  summarise(PresAbs = max(PresAbs)) %>%
  # rename(COMID = COMID2) %>%
  drop_na()

# ## how many presense/absences
sum(bio_coms$PresAbs ==1) ## 235
sum(bio_coms$PresAbs ==0) ## 40

head(allData)
head(bio_coms)
## join with env data
NewDataObsComid <- inner_join(allData, bio_coms, by = "COMID")

## replace NAs in presabs with -999

# NewDataObsComid$PresAbs[is.na(NewDataObsComid$PresAbs)] <- -999

## saveout
save(NewDataObsComid, file=paste0(COMIDFile, "all_presAbs_env_data.RData"))


# Multicolinearality ------------------------------------------------------

## comid
# names(all_data_obs_coms)
all_data_obs <- NewDataObsComid %>%
  as.data.frame() %>% 
  dplyr::select(-c(starts_with(c("ppt", "tmin", "tmax")))) %>%
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

class(ydata)
length(ydata)

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

save(rf.data, file=paste0(COMIDFile, "data_for_rf_model.RData"))

(rf.final <- randomForest(y=rf.data[,1], x=rf.data[,2:ncol(rf.data)], ntree=b, mtry = 18,nodesize=5,
                          importance=TRUE, norm.votes=TRUE, proximity=TRUE) )



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


file.name1 <- paste0(COMIDFile, "no_clim_model_mds_scale_gridded.jpg") 
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

save(model, file=paste0(COMIDFile, "validation.RData"))
save(conMat, file=paste0(COMIDFile, "confusion_matrix.RData"))


pdf(paste0(COMIDFile, "Trees_no_clim_gridded.pdf"), width=25, height=15)

full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data)
plot(full_tree)
text(full_tree, use.n = T)

dev.off()


# Coeficients and importance ----------------------------------------------

# PLOT VARIABLE IMPORTANCE
pdf(paste0(COMIDFile, "var_imp_no_clim_gridded.pdf"), width=25, height=15)

varImpPlot(rf.final)

dev.off()

## scaled var imp
pdf(paste0(COMIDFile, "var_imp_scaled_no_clim_gridded.pdf"), width=25, height=15)

dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)

dev.off()

## partial dependence plots
pdf(paste0(COMIDFile, "PartialPlots_6_no_clim_gridded.pdf"), width=8, height=8)
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
