############################
##Load Necessary Packages###
############################
library(sp)
library(rgdal)
library(raster)
library(randomForest)
library(kernlab)
library(rgl)
library(ks)
library(sm)

### added by Katie 
library(sf)

#Import averaged rasters
#havg<-raster("D:/AT DistributionModeling/Historic/randomForests/PARTITIONING/Hist_10x_Final/hAvg.img")
cavg<-raster("ignore/Current/randomForests/PARTITIONING/Cur_10x_Final_Feb2014/cbAvg.img")

wsdir="ignore/Current/randomForests"
#wsdir="E:/MLT DATA/AT DistributionModeling/Historic/randomForests"
path=paste(wsdir, "PARTITIONING", sep="/")
#  load( paste(path, "RFModel_Hist.RData", sep="/"))
setwd(path)
path
getwd()
############################################
##Import Points and Create Pseudoabsences###
############################################

sf::sf_use_s2(FALSE) ## switch off spherical geom

inshape="200mCells_PresAbs2005 PCA PresAbs Pts_MLT"



## data not loading: changed DATA2 to DATA3, removed "path"
orig.sdata<-sdata <- readOGR(dsn=paste("DATA3", sep="/"), layer=inshape)

## another option if needed to update code
# orig.sdata<-sdata <- st_read("DATA3/200mCells_PresAbs2005 PCA PresAbs Pts_MLT.shp")

str(sdata@data)
## warning message
proj4string(sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")
proj4string(orig.sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")

## the below is repeated code?
#sdata <- sdata[1]
# str(sdata)
# proj4string(sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")


## KDE Bias Surface
##################################################
## KDE Bias Surface
##################################################
# develop KDE sampling bias surface
orig.sdata2<-subset(orig.sdata, PresAbs200==1)
mask <- cavg[[1]]>-1000
bias <- cellFromXY(mask, orig.sdata[,-1])
cells <- unique(sort(bias))
kernelXY <- xyFromCell(mask, cells)
samps <- as.numeric(table(bias))

# code to make KDE raster
KDEsur <- sm.density(kernelXY, weights=samps, display="none", ngrid=812, 
                     ylim=c(3600402,3730202), xlim=c(423638,563638), nbins=0)
KDErast=SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
KDErast = SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, 
                                                    length(KDEsur$estimate))))
KDErast <- raster(KDErast)
KDErast <- resample(KDErast, mask)
KDErast <- KDErast*mask
KDEpts <- rasterToPoints(KDErast, spatial=F) #(Potential) PseudoAbsences are created here 

#Now to integrate Pseudoabsences into Presence Data
a=KDEpts[sample(seq(1:nrow(KDEpts)), size=702, replace=T, prob=KDEpts[,"layer"]),1:2]

PA.abs<-data.frame(PresAbs200=rep(0,nrow(a)))
head(PA.abs)
a.sp<-SpatialPoints(a, proj4string=CRS("+proj=utm +zone=11 +datum=WGS84"))
a.spdf<-SpatialPointsDataFrame(a.sp, PA.abs)
sdata<-rbind(sdata,a.spdf)
plot(sdata)

#############################
#############################
##############################
training.pred<-data.frame(extract(cavg, sdata))
training.pred<-cbind(sdata@data[,1], training.pred)
colnames(training.pred)<-c("PresAbs", "Pred")
training.pred$ID<-seq(1,nrow(training.pred),1)
training.pred<-training.pred[,c(3,1,2)]
training.pred$Class<-training.pred$Pred
training.pred$Class[training.pred$Class>=0.6357364]<-1
training.pred$Class[training.pred$Class<0.6357364]<-0
head(training.pred)#Make sure it looks right

#Try out PresenceAbsence functions
library(PresenceAbsence)
# install.packages("PresenceAbsence")

str(training.pred)

## change to numeric
training.pred$PresAbs <- as.numeric(training.pred$PresAbs)

auc(training.pred)
auc.roc.plot(training.pred)
calibration.plot(training.pred)
cmx1<-cmx(training.pred)
error.threshold.plot(training.pred[1:880,], opt.methods=c(seq(1,11,1)))
Kappa(cmx1)
predicted.prevalence(training.pred)
presence.absence.accuracy(training.pred)
presence.absence.hist(training.pred)
presence.absence.summary(training.pred)
roc.plot.calculate(training.pred)

##############################
###Explore effects of different thresholds
#############################
#With MaxPCC and Sens=Spec (same) from PlotROC function
predPA<-cavg
predPA[predPA>=0.49228]<-1
predPA[predPA<0.49228]<-0

#Calculate Minimum Area Occupied
freq(predPA) #Provides data to calculate Min. Occupied Area
freq(predPA)[2,2]/(freq(predPA)[2,2]+freq(predPA)[1,2])



##############################
###Validation from Jeff Evans
#############################
path=paste(wsdir, "PARTITIONING", sep="/")
#  load( paste(path, "RFModel_Hist.RData", sep="/"))
setwd(path)
#InstallPkgs(c("randomForest","caTools","psych"))
library(randomForest)
library(caTools)
library(psych)
# install.packages("psych")

source(paste(paste( "DATA3", sep="/"), "Functions.R", sep="/"))


#ConfusionMatrix
confu<-table(training.pred$PresAbs,training.pred$Class)

Kappa <- function(x) {
    n <- sum(x)
      ni <- apply(x, 1, sum)
        nj <- apply(x, 2, sum)
          m <- min(length(ni), length(nj))
        p0 <- sum(diag(x[1:m, 1:m]))/n
      pc <- sum(ni[1:m] * nj[1:m])/n^2
    return( (p0 - pc)/(1 - pc) )
	}

# Kappa 
Kappa(confu)
#0.8874842

# weighted kappa is 
#  (probability of observed matches - probability of expected matches)/
#  (1 - probability of expected matches)
# 
#  Unweighted Kappa just considers the matches on the main diagonal whereas weighted kappa  
#     considers off diagonal elements as well. 
# Weighted Kappa
library(psych) #Used for cohen.kappa
 wts <- matrix(c(0,1,1,0),ncol=2)
    
 cohen.kappa(confu[,1:2], wts)
  wts <- matrix(c(0.9,0.1,0.1,0.9),ncol=2)
    cohen.kappa(confu[,1:2], wts)

##########################################################################
# CHUNK 2
# Backprediction
##########################################################################
##########################################################################
##########################################################################
  
######### do not have rf.final in r environ - need code from RF_CLASS_MODEL?
    
# Predict dependent variable and class probabilities 
#Relies on Jeff Evans' functions
rf.pred <- predict(rf.final, rf.data[,2:ncol(rf.data)], type="response")
rf.prob <- as.data.frame(predict(rf.final, rf.data[,2:ncol(rf.data)], type="prob"))
ObsPred <- data.frame(cbind(Observed=as.numeric(as.character(ydata)), 
                      PRED=as.numeric(as.character(rf.pred)), Prob1=rf.prob[,2], 
					  Prob0=rf.prob[,1]) )
V0=rf.final$votes[,1]					  
V1=rf.final$votes[,2]
					  
op <- (ObsPred$Observed == ObsPred$PRED)
( pcc <- (length(op[op == "TRUE"]) / length(op))*100 )
			
# Calculate variety of validation statistics using backprediction
PlotROC(ObsPred[,2], ObsPred[,1])#, summaryFile="Models_Oct14/AUC-Results_backpred.txt")


#### USE THIS CODE!!!!
PlotROC(training.pred[,3], training.pred[,2], summaryFile="Cur_10x_Final_Feb2014/AUC-Results_backpred_prespseudo.txt")

PlotROC(training.pred[,3], training.pred[,2], cutoff="manual", cutoffValue=0.492, summaryFile="Cur_10x_Final_Feb2014/AUC-Results_backpred_prespseudo_LowSens.txt")
0.492280772
############################################################
#The below calculates the statistics for ONLY true pres/abs
PlotROC(training.pred[1:880,3], training.pred[1:880,2], summaryFile="Cur_10x_Final_Feb2014/AUC-Results_backpred_TruepresAbs.txt")

confu2<-table(training.pred$PresAbs[1:880],training.pred$Class[1:880])

Kappa <- function(x) {
    n <- sum(x)
      ni <- apply(x, 1, sum)
        nj <- apply(x, 2, sum)
          m <- min(length(ni), length(nj))
        p0 <- sum(diag(x[1:m, 1:m]))/n
      pc <- sum(ni[1:m] * nj[1:m])/n^2
    return( (p0 - pc)/(1 - pc) )
	}
Kappa(confu2)
#0.8610655


#save.image( paste(path, "Cur_10x_Final_Feb2014/Model10/CurrentAssessment.RData", sep="/") ) 


############################
###Testing wityh ROCR
############################
library(ROCR)
    pred1 <- prediction(training.pred[1:880,3], training.pred[1:880,2])
        pred2 <- prediction(training.pred[,3], training.pred[,2])

    # Create the ROC performance object. 
    perf <- performance(pred1, "tpr", "fpr")
perf

(performance(pred, "auc"))

    auc <- performance(pred, "auc")@y.values[[1]]
plot(performance(pred,"sens","spec"))

#To plot Sensitivity/Specifity
x<-(performance(pred1,"sens"))
y<-(performance(pred1,"spec"))
y
plot(x)
plot(y, add=T)
plot(

#Find Sens=Spec value
x1<-(rbind(x@x.values, x@y.values, y@y.values))
x1<-(cbind(x@x.values[[1]], x@y.values[[1]], y@y.values[[1]]))
x1

[787,] 0.492280772 1.000000000 1.00000000