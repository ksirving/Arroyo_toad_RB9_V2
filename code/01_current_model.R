InstallPkgs <- function(pkgs) {
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss, dep=TRUE, repos="http://cran.us.r-project.org")
  }  
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss, dep=TRUE, repos="http://cran.us.r-project.org")
  }
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
  }
}
InstallPkgs(c("sp","rgdal","raster","randomForest","kernlab","rgl","ks"))
#Normally need to install these packages,"rgl","ks"
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

setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")

wsdir="original_model/Current/randomForests"

## define and set path
path=paste0(wsdir, "/PARTITIONING/DATA3") 
setwd(path)

## get path for functions
source("Functions.R")

# Bio data ----------------------------------------------------------------
getwd()

## define and set path
path=paste0(wsdir, "/PARTITIONING/DATA3") 
setwd(path)



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
xvars <- stack("ignore/00_raw_data_raster.grd")
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


# Join to scaled data -----------------------------------------------------

## this is a quick fix for now, will need to be fixed

## bio comids
load(file= "ignore/00_pseudo_abs.RData") ## data1
head(data1)

## make spatial
data1 <- data1 %>%
  st_as_sf(coords=c("X", "Y"), crs=32611, remove=F) %>%
  dplyr::select(-ID)

## make sure if reach has at least one presence then keeps it, if none then absence
data1a <- data1 %>%
  as.data.frame() %>%
  group_by(COMID) %>% 
  summarise(PresAbs200 = max(PresAbs200))
  
data1a


## join presence/absence to env df

load(file = "ignore/00_all_env_bio_data_NHD_reach.RData") ## all_data
head(all_data)
names(all_data2)

## join all data with all data with pres/absence
all_data2 <- full_join(all_data, data1a, by = "COMID") %>% ### doesn't match all comids - fix!
  mutate(PresAbs200 = as.numeric(PresAbs200)) %>%
  mutate(PresAbs200 = replace_na(PresAbs200, -999))


all_data_obs  <- right_join(all_data, data1a, by = "COMID")
sum(is.na(all_data_obs))

all_data_obs <- na.omit(all_data_obs) ### fix this issue!!! 
names(all_data_obs)

save(all_data_obs, file = "ignore/01_all_env_bio_data_NHD_reach_original.RData")

dim(all_data_obs)

# Join observations and mulitcolinearlity -------------------------------------------------------

# FETCH RASTER VALUES FOR X VARABLES AND JOIN TO SPECIES OBSERVATION DATA (Assigns raster values to points)
# sum(is.na(all_data_obs))
# 
# which(is.na(all_data_obs))
# 
# all_data_obs[1672:1702,100:105]

cl <- MultiColinear(all_data_obs[,2:101], p=0.05)
xdata <- all_data_obs[,2:101]

for(l in cl) {
cl.test <- xdata[,-which(names(xdata)==l)]
print(paste("REMOVE VARIABLE", l, sep=": "))
MultiColinear(cl.test, p=0.05)
}

l
# # REMOVE MULTI-COLINEAR VARIABLES
for(l in cl) { all_data_obs <- all_data_obs[,-which(names(all_data_obs)==l)] }

### not sure how this works, check and edit

names(all_data_obs)

# Random forest model -----------------------------------------------------

ydata <- as.factor(all_data_obs$PresAbs200)
xdata <- all_data_obs[,2:59] ## do not include template layer

class(ydata)
length(ydata)

# PERCENT OF PRESENCE OBSERVATIONS
( dim(all_data_obs[all_data_obs$PresAbs200 == 1, ])[1] / dim(all_data_obs)[1] ) * 100 ## 48%

# RUN RANDOM FORESTS MODEL SELECTION FUNCTION
#Also provides variable importance and such
#It is important to look at highest class error - Look at Global OOB errors
#TEST Object is super important; Row number corresponds to paramerers, in output
( rf.model <- rf.modelSel(x=xdata, y=ydata, imp.scale="mir", ntree=b, nodesize=5) ) 

# CREATE NEW XDATA BASED ON SELECTED MODEL AND RUN FINAL RF MODEL (Model 4)
#RF runs differently when you use symbolis languate (using ~ as in an Lin. Model... use the indexing approacy [y=rf.data[,1]...]

sel.vars <- rf.model$PARAMETERS[[3]]# set to use 3 - lowest error rate and has all hydro vars

rf.data <- data.frame(y=ydata, xdata[,sel.vars])	

(rf.final <- randomForest(y=rf.data[,1], x=rf.data[,2:ncol(rf.data)], ntree=b, nodesize=5,
                          importance=TRUE, norm.votes=TRUE, proximity=TRUE) )
names(rf.data )

## split into training and testing

setsize <- floor(nrow(rf.data)*0.8)
index <- sample(1:nrow(rf.data), size = setsize)
training <- rf.data[index,]
testing <- rf.data[-index,]

(rf.train <- randomForest(y=training[,1], x=training[,2:ncol(training)], ntree=b, nodesize=5,
                          importance=TRUE, norm.votes=TRUE, proximity=TRUE) )

plot(rf.train)

dim(testing)

result <- as.data.frame( predict(rf.train, testing[, 2:16], type = "response"))
result$y <- testing$y

plot(result)


### visualise trees

library(rpart)

pdf("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9/Figures/01_Trees.pdf", width=25, height=15)

full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data)
plot(full_tree)
text(full_tree, use.n = T)

dev.off()

##  check error rates

rf.final$err.rate

oob.error.data <- data.frame(
  Trees = rep(1:nrow(rf.final$err.rate), times = 3),
  Type = rep(c("OOB", "0", "1"), each = nrow(rf.final$err.rate)),
  Error = c(rf.final$err.rate[,"OOB"],
            rf.final$err.rate[, "0"],
            rf.final$err.rate[,"1"])
)
oob.error.data 

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

## check split

oob.values <- vector(length=10)

for(i in 1:10) {
  temp.model <- randomForest(y~., data=rf.data, mtry = i, ntree=b)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate), 1]
  
}

oob.values

### MDS scaling from proximity

# CREATE CMD SCALED DISTANCES OF RF PROXMITIES
distanceMatrix <- dist(1 - rf.final$proximity)

rf.cmd <- cmdscale(distanceMatrix, eig=TRUE, k=4, x.ret=T)  
# pa.col <- cRamp(ydata)  
# rf.cmd <- data.frame(rf.cmd$points)

mdsVarPer <- round(rf.cmd$eig/sum(rf.cmd$eig)*100, 1)
mdsVarPer

mdsValues <- rf.cmd$points

mdsData <- data.frame(Sample = rownames(mdsValues),
                      X=mdsValues[,1],
                      Y=mdsValues[,2],
                      Status = rf.data$y)

mdsData

m1 <- ggplot(data = mdsData, aes(x=X, y=Y)) +
  geom_point(aes(color=Status)) +
  theme_bw()+
  xlab(paste("MDS1 - ", mdsVarPer[1], "%", sep=" ")) +
  ylab(paste("MDS2 - ", mdsVarPer[2], "%", sep=" ")) +
  ggtitle(paste("Pres/Abs", "PROXIMITY MATRIX", sep=" - "))

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9/Figures/01_full_model_mds_scale.jpg"
ggsave(m1, filename=file.name1, dpi=300, height=5, width=6)


# Coeficients and importance ----------------------------------------------

rf.final$importance[,1]
varImpPlot(rf.final)

rf.final$oob.times
rf.final$test

# PLOT BOOTSTRAP ERROR CONVERGENCE
plot(rf.final, main="Bootstrap Error Convergence")

# PLOT VARIABLE IMPORTANCE
p <- as.matrix(rf.final$importance)    
ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]])  
dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)

# PLOT PROTOTYPES (MULTIVARIATE CLASS CENTERS)					  
rf.p <- classCenter(xdata, ydata, rf.final$prox)
v1=4; v2=1
plot(xdata[,v1], xdata[,v2], pch=21, xlab=names(xdata)[v1], ylab=names(xdata)[v2],
     bg=c("black", "grey94")[as.numeric(factor(ydata))], main="Data with Prototypes")
points(rf.p[,v1], rf.p[,v2], pch=21, cex=2.5, bg=c("blue", "red"))					  



# Predictions and model performance ----------------------------------------------------------------

library(ROCR)

all_data_obs <- all_data_obs %>% mutate(PresAbs200 = as.factor(PresAbs200))

#Index refers to the right column of probabilities - in this model the second column, which is probs of "1"

pred1 <- predict(rf.final, all_data_obs,filename="output_data/Current/Model1/SppProbs.img", type="prob",  index=2, 
        na.rm=TRUE, overwrite=TRUE, progress="window")

pred1

pred2 = prediction(pred1[,2], all_data_obs$PresAbs200)
pred2

# # 1. Area under curve
# performance(perf, "auc")
# 
# 2. True Positive and Negative Rate
perf = performance(pred2, "tpr","fpr")
# 3. Plot the ROC curve
plot(perf,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray") ## model is too good!!!

perf <- performance(pred2,"tpr","fpr")
perf@x.name
plot(perf)
perf@x.values

# precision/recall curve (x-axis: recall, y-axis: precision)
perf <- performance(pred2, "prec", "rec")
perf
plot(perf)

# sensitivity/specificity curve (x-axis: specificity,
# y-axis: sensitivity)
perf <- performance(pred2, "sens", "spec")
perf
plot(perf)

rf.final$confusion
rf.final$pred[order(model_rf$pred$rowIndex),2]

#output

# Predict on all rb9 region-------------------------------------------------------

#Index refers to the right column of probabilities - in this model the second column, which is probs of "1"
all_data <- na.omit(all_data) ## fix this issue!!!

pred <- predict(rf.final, all_data,filename="output_data/Current/Model1/SppProbs.img", type="prob",  index=2, 
                              na.rm=TRUE, overwrite=TRUE, progress="window")


pred_df <- as.data.frame(predict(rf.final, all_data,filename="output_data/Current/Model1/SppProbs.img", type="prob",  index=2, 
                              na.rm=TRUE, overwrite=TRUE, progress="window"))
## add comids
pred_df$COMID <- all_data$COMID



# partial dependence plots ------------------------------------------------

path
sp=0.6
pdf("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9/Figures/PartialPlots_6.pdf", width=8, height=8)
p <- as.matrix(rf.final$importance)    
ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]])  
dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)
plot(rf.final, main="Bootstrap Error Convergence")

names(rf.data)
for (i in 1:length(names(rf.data[,2:ncol(rf.data)])) ) {
  p <- partialPlot(rf.final, rf.data[,2:ncol(rf.data)], 
                   names(rf.data[,2:ncol(rf.data)])[i], which.class="1", plot=FALSE)   
  p$y <- (p$y - min(p$y)) / (max(p$y) - min(p$y)) 
  plot( y=lowess(y=p$y, x=p$x, f=sp)$y, x=p$x, type="l", ylab="p",  
        xlab=names(rf.data[,2:ncol(rf.data)])[i],
        main=paste("PARTIAL PLOT", names(rf.data[,2:ncol(rf.data)])[i], sep=" - ")) 
}
dev.off() 

## partial plots only hydro
# 
# pdf("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9/Figures/PartialPlots_hydro_fa_mag.pdf", width=8, height=8)
# 
# hyd_vars <- rf.data[,44:47]
# hyd_vars
# i= 45
# names(rf.data)
# 
# for (i in 44:length(names(hyd_vars)) ) {
#   p <- partialPlot(rf.final, rf.data[,2:ncol(rf.data)], 
#                    names(rf.data[,2:ncol(rf.data)])[i], which.class="1", plot=FALSE)   
#   p$y <- (p$y - min(p$y)) / (max(p$y) - min(p$y)) 
#   plot( y=lowess(y=p$y, x=p$x, f=sp)$y, x=p$x, type="l", ylab="p",  
#         xlab=names(rf.data[,2:ncol(rf.data)])[i],
#         main=paste("PARTIAL PLOT", names(rf.data[,2:ncol(rf.data)])[i], sep=" - ")) 
# }
# dev.off() 



# Plotting probability of occurrence --------------------------------------

nhd <- st_read("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/NHDPlus_V2_Flowline_CA.shp")
head(pred)

head(nhd)

names(all_data)

obs <- all_data_obs %>%
  dplyr::select(COMID, PresAbs200)

nhd_lines_rb9 <- nhd %>%
  dplyr::select(Shape_Leng, COMID) %>%
  filter(COMID %in% pred_df$COMID) %>%
  mutate(COMID = as.integer(COMID))

head(nhd_lines_prob)
dim(nhd_lines_rb9)
str(nhd_lines_rb9)

nhd_lines_prob <- full_join(nhd_lines_rb9, pred_df, by = "COMID")
nhd_lines_prob <- st_zm(nhd_lines_prob)

plot(nhd_lines_prob[5])

### round prob values - change later
nhd_lines_prob <- nhd_lines_prob %>%
  rename(probOcc = 4) %>%
  mutate(probRound = round(probOcc, digits=1))

## join in obs

nhd_lines_prob <- full_join(nhd_lines_prob, obs, by = "COMID")

nhd_lines_prob

library(viridis)
map1 <- ggplot() +
  geom_sf(data = nhd_lines_prob, aes(color = probOcc)) + 
  geom_sf(data = subset(data1, PresAbs200 == 1)) +
  scale_fill_gradientn(colours=rev(magma(6))) ## colours not working

map1

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9/Figures/01_full_model_prob_occs_map.jpg"
ggsave(map1, filename=file.name1, dpi=300, height=5, width=6)

library(mapview)
library(sf)
library(RColorBrewer)
## map
# set background basemaps:
basemapsList <- c("Esri.WorldTopoMap", "Esri.WorldImagery",
                  "Esri.NatGeoWorldMap",
                  "OpenTopoMap", "OpenStreetMap", 
                  "CartoDB.Positron", "Stamen.TopOSMFeatures")

mapviewOptions(basemaps=basemapsList, fgb = FALSE)


m1 <- mapview(nhd_lines_prob, zcol = "probRound", col.regions=brewer.pal(11, "YlGn"))

m1@map %>% leaflet::addMeasure(primaryLengthUnit = "meters")


##################################################
# PROXIMITY MULTIDIMENSIONAL SCALING (MDS) PLOTS #
##################################################

######### FUNCTION TO CREATE COLOR VECTOR #########
#AWESOME!!! SAVE THIS FOREVER!!! - Orders color ramp with data#
cRamp <- function(x,d=c("blue", "red")){
  crange <- function(x)(as.numeric(x)-min(as.numeric(x)))/diff(range(as.numeric(x)))
  cols <- colorRamp(d)(crange(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  


###################################################

# CREATE CMD SCALED DISTANCES OF RF PROXMITIES
rf.cmd <- cmdscale(1 - rf.final$proximity, eig=TRUE, k=4)  
pa.col <- cRamp(ydata)  
rf.cmd <- data.frame(rf.cmd$points)

# 2D PLOT OF FIRST TWO MDS DIMENSIONS	
plot(rf.cmd[,1:2], ylab="DIM 1", xlab="DIM 2", pch=16, col=pa.col,
     main= paste("Pres/Abs", "PROXIMITY MATRIX MDS d=2", sep=" - "))
legend("bottomright", pch=c(16,16), col=c("blue", "red"), 
       legend=c("Absent","Present")) 

# 3D PLOT OF FIRST THREE MDS DIMENSIONS - Interactive using Windows RGL Driver (Library RGL)	   
plot3d(rf.cmd[,1],rf.cmd[,2],rf.cmd[,3], col=pa.col,
       pch=18, size=1.25, type="s", xlab="MDS dim 1", ylab="MDS dim 2", 
       zlab="MDS dim 3")	

# SAVE R IMAGE 
save.image( paste(path, "Cur_10x_Final_Feb2014/Model1/RFClassModel.RData", sep="/") ) 