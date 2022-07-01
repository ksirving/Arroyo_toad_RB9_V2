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

#Import averaged rasters
#havg<-raster("D:/AT DistributionModeling/Historic/randomForests/PARTITIONING/Hist_10x_Final/hAvg.img")
#cavg<-raster("D:/AT DistributionModeling/Current/randomForests/PARTITIONING/Cur_10x_Final/cAvg.img")
havg<-raster("D:/AT DistributionModeling/Historic/randomForests/PARTITIONING/Hist_10x_Final_1981/hbAvg.img")

wsdir="D:/AT DistributionModeling/Historic/randomForests"
#wsdir="E:/MLT DATA/AT DistributionModeling/Historic/randomForests"
path=paste(wsdir, "PARTITIONING", sep="/")
#  load( paste(path, "RFModel_Hist.RData", sep="/"))
setwd(path)
############################################
##Import Points and Create Pseudoabsences###
############################################
inshape="200mCells_PresPts"
sdata <- readOGR(dsn=paste(path, "DATA2_1981", sep="/"), layer=inshape)
  str(sdata@data)

orig.sdata<-sdata<-sdata[3]

proj4string(sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")
proj4string(orig.sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")

#sdata <- sdata[1]
str(sdata)
proj4string(sdata)<-CRS("+proj=utm +zone=11 +datum=WGS84")

## KDE Bias Surface
# develop KDE sampling bias surface
mask <- havg[[1]]>-1000
bias <- cellFromXY(mask, sdata[,-1])
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

#Now to integrate Pseuooabsences into Presence Data
a=KDEpts[sample(seq(1:nrow(KDEpts)), size=1037, replace=T, prob=KDEpts[,"layer"]),1:2]

Presence<-data.frame(Presence=rep(0,nrow(a)))
a.sp<-SpatialPoints(a, proj4string=CRS("+proj=utm +zone=11 +datum=WGS84"))
sdata@data$Presence<-as.numeric(sdata@data$Presence) #Need to make sure Presence Colum in species data is numeric, matching with the created pseudo-absence data
a.spdf<-SpatialPointsDataFrame(a.sp, Presence)
sdata<-rbind(sdata[1],a.spdf)
plot(sdata)
#############################
#############################
##############################
training.pred<-data.frame(extract(havg, sdata))
training.pred<-cbind(sdata@data[,1], training.pred)
colnames(training.pred)<-c("PresAbs", "Pred")
training.pred$ID<-seq(1,nrow(training.pred),1)
training.pred<-training.pred[,c(3,1,2)]
training.pred$Class<-training.pred$Pred
#training.pred$Class[training.pred$Class>=0.659524]<-1
#training.pred$Class[training.pred$Class<0.659524]<-0
head(training.pred)#Make sure it looks right

#Try out PresenceAbsence functions
library(PresenceAbsence)
auc(training.pred)
auc.roc.plot(training.pred)
calibration.plot(training.pred)
cmx1<-cmx(training.pred)
error.threshold.plot(training.pred, opt.methods=c(seq(1,11,1)))
Kappa(cmx1)
predicted.prevalence(training.pred)
presence.absence.accuracy(training.pred)
presence.absence.hist(training.pred)
presence.absence.summary(training.pred)
#roc.plot.calculate(training.pred)

##############################
###Explore effects of different thresholds
#############################
#With MaxPCC and Sens=Spec (same) from PlotROC function
predPA<-havg
predPA[predPA>=0.680961904 ]<-1
predPA[predPA<0.680961904]<-0

#Calculate Minimum Area Occupied
freq(predPA) #Provides data to calculate Min. Occupied Area
freq(predPA)[2,2]/(freq(predPA)[2,2]+freq(predPA)[1,2])


#With Lowest Value for Sens=1 (Spec = 0.801002893)
predPA<-havg
predPA[predPA>=0.43512684]<-1
predPA[predPA<0.43512684]<-0

#Calculate Minimum Area Occupied
freq(predPA) #Provides data to calculate Min. Occupied Area
freq(predPA)[2,2]/(freq(predPA)[2,2]+freq(predPA)[1,2])

##############################
###Validation from Jeff Evans
#############################
#path=paste(wsdir, "PARTITIONING", sep="/")
#  load( paste(path, "RFModel_Hist.RData", sep="/"))
#setwd(path)
#InstallPkgs(c("randomForest","caTools","psych"))
library(randomForest)
library(caTools)
library(psych)

source(paste(paste(path, "DATA2_1981", sep="/"), "Functions.R", sep="/"))


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

confu<-table(training.pred2$PresAbs,training.pred2$Class)

	
# Kappa 
Kappa(confu)
#0.9479

# weighted kappa is 
#  (probability of observed matches - probability of expected matches)/
#  (1 - probability of expected matches)
# 
#  Unweighted Kappa just considers the matches on the main diagonal whereas weighted kappa  
#     considers off diagonal elements as well. 


##########################################################################
# CHUNK 2
# Backprediction
##########################################################################
##########################################################################
##########################################################################
  
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

#USE THIS CODE
PlotROC(training.pred[,3], training.pred[,2], summaryFile="Hist_10x_Final_1981/AUC-Results_backpred_prespseudo.txt")

training.pred2<-training.pred
training.pred2$Class[training.pred2$Class>=0.43512684]<-1
training.pred2$Class[training.pred2$Class<0.43512684]<-0
head(training.pred2)#Make sure it looks right
#PlotROC(training.pred2[,4], training.pred2[,2], summaryFile="Hist_10x_Final_1981/AUC-Results_backpred_prespseudoMinTPR1.txt")
PlotROC(training.pred2[,3], training.pred2[,2], cutoff="Manual", cutoffValue=0.43512684,summaryFile="Hist_10x_Final_1981/AUC-Results_backpred_prespseudoMinTPR1.txt")



#save.image( paste(path, "Hist_10x_Final_1981/Historic1981Assessment.RData", sep="/") ) 






############################
###Testing wityh ROCR
############################
pred1 <- prediction(training.pred[,3], training.pred[,2])

    # Create the ROC performance object. 
perf <- performance(pred, "tpr", "fpr")
perf

(performance(pred1, "sens", "spec"))

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
x1<-(cbind(x@x.values[[1]], x@y.values[[1]], y@y.values[[1]]))
x1
t(x1)



#0.692110789 
(0.9353905497+0.9035679846)
#0.689961004 
(0.9363548698+0.9035679846)
#0.687601240 
(0.9363548698+0.9026036644)
#0.681011899 
(0.9373191900+0.9026036644)
#THIS IS IT 0.680961904 
(0.9382835101+0.9026036644)
#0.678162184 
(0.9382835101+0.9016393443)
#0.677712229 
(0.9382835101+0.9006750241)
#0.677272273 
(0.9392478303+0.9006750241)
