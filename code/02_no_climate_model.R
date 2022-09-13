## model with no climate

## random forest - full model

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

setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2")

## upload data

load(file =  "ignore/00a_data_for_model_abs_flow.RData") # NewDataObsSub - for model build
load(file =  "ignore/00a_data_for_prediction_abs_flow.RData") #NewDataObs - for prediction

NewDataObs <- NewDataObs %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

dim(NewDataObs)

## get path for functions
source("original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

# Join observations and mulitcolinearlity -------------------------------------------------------
names(NewDataObsSub)
## remove climate vars
## remove nas and rearrange for easy indexing
## remove ds 50 and change neg ffm values to zero
all_data_obs <- NewDataObsSub %>%
  select(-DS_Mag_50) %>%
  # pivot_longer(DS_Mag_90:Wet_BFL_Mag_10, names_to = "ffm", values_to = "value") %>%
  # mutate(value = ifelse(value < 0, 0, value)) %>%
  # pivot_wider(names_from = ffm, values_from = value) %>%
  select(COMID, NewObs, SHAPE_LENG,-c(pptAnnAv:TmaxMon12), AvgClay:MRVBF.Mx, TC_042014_RB9.1_Var:TC_092014_RB9.3_Mean, geometry) %>%
  drop_na()



sum(is.na(all_data_obs))
dim(all_data_obs)
names(all_data_obs)


cl <- MultiColinear(all_data_obs[,c(4:35)], p=0.05)
xdata <- all_data_obs[,c(4:35)]
xdata

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

set.seed(234)
# NUMBER OF BOOTSTRAP REPLICATES
b=10001

ydata <- factor(all_data_obs$NewObs, levels = c(1,0))
xdata <- all_data_obs[,c(4:34)] ## do not include template layer

class(ydata)
length(ydata)

# PERCENT OF PRESENCE OBSERVATIONS
( dim(all_data_obs[all_data_obs$NewObs == 1, ])[1] / dim(all_data_obs)[1] ) * 100 ## 55%

# RUN RANDOM FORESTS MODEL SELECTION FUNCTION
#Also provides variable importance and such
#It is important to look at highest class error - Look at Global OOB errors
#TEST Object is super important; Row number corresponds to paramerers, in output
( rf.model <- rf.modelSel(x=xdata, y=ydata, imp.scale="mir", ntree=b, nodesize=5) ) 

# CREATE NEW XDATA BASED ON SELECTED MODEL AND RUN FINAL RF MODEL (Model 4)
#RF runs differently when you use symbolis languate (using ~ as in an Lin. Model... use the indexing approacy [y=rf.data[,1]...]

sel.vars <- rf.model$PARAMETERS[[3]]# set to use 3 - lowest error rate and has all hydro vars

rf.data <- data.frame(y=ydata, xdata[,sel.vars])	

(rf.final <- randomForest(y=rf.data[,1], x=rf.data[,2:ncol(rf.data)], ntree=b, mtry = 2,nodesize=5,
                          importance=TRUE, norm.votes=TRUE, proximity=TRUE) )
# names(rf.data)
# rf.data

#OOB estimate of  error rate: 19.14%


# Tuning ------------------------------------------------------------------

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
oob.values

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
m1

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9/Figures/01_full_model_mds_scale_no_clim.jpg"
ggsave(m1, filename=file.name1, dpi=300, height=5, width=6)


# Validation --------------------------------------------------------------

## make y data compatible
rf.data.val <- rf.data %>%
  mutate(y = ifelse(y==1, "Present", "Absent")) %>%
  mutate(y = factor(y, levels = c("Present", "Absent")))

str(rf.data.val)

## split into training and testing

setsize <- floor(nrow(rf.data)*0.8)
index <- sample(1:nrow(rf.data), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]
training


trcontrol = trainControl(method='cv', number=10, savePredictions = T,
                         classProbs = TRUE,summaryFunction = twoClassSummary,returnResamp="all")

model = train(y ~ . , data=training, method = "rf", trControl = trcontrol,metric="ROC") 

model$resample

model

# mtry  ROC        Sens       Spec     
# 2    0.8343590  0.7866667  0.7615385
# 9    0.8264103  0.7600000  0.7461538
# 16    0.8264103  0.7666667  0.7538462
# 
# ROC was used to select the optimal model using the largest value.
# The final value used for the model was mtry = 2.


confusionMatrix(predict(model,testing),testing$y)

# Confusion Matrix and Statistics
# 
# Reference
# Prediction Present Absent
# Present      35      8
# Absent        5     22
# 
# Accuracy : 0.8143          
# 95% CI : (0.7034, 0.8972)
# No Information Rate : 0.5714          
# P-Value [Acc > NIR] : 1.539e-05       
# 
# Kappa : 0.616           
# 
# Sensitivity : 0.875           
# Specificity : 0.7333          
# Pos Pred Value : 0.814           
# Neg Pred Value : 0.8148          
# Prevalence : 0.5714          
# Detection Rate : 0.5             
# Detection Prevalence : 0.6143          
# 
# 'Positive' Class : Present  
### visualise trees

library(rpart)

pdf("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/01_Trees.pdf", width=25, height=15)

full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data)
plot(full_tree)
text(full_tree, use.n = T)

dev.off()


# Coeficients and importance ----------------------------------------------

rf.final$importance[,1]

pdf("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/var_imp_no_clim.pdf", width=20, height=10)

varImpPlot(rf.final)

dev.off()

vp1

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


# Predict on all rb9 region-------------------------------------------------------

#Index refers to the right column of probabilities - in this model the second column, which is probs of "1"
all_data <- na.omit(NewDataObs)
# all_data <- all_data[,c("COMID", sel.vars)] 
head(all_data)
str(all_data)

pred <- predict(rf.final, all_data,filename="output_data/Current/Model1/SppProbs_no_clim.img", type="prob",  index=2, 
                na.rm=TRUE, overwrite=TRUE, progress="window")


pred_df <- as.data.frame(predict(rf.final, all_data,filename="output_data/Current/Model1/SppProbs_no_clim.img", type="prob",  index=2, 
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

write.csv(pred_env, "02_hydro_obs_preds_no_clim_model.csv")

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

file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/02_incremental_preds_sep_coms.jpg"
ggsave(p1, filename=file.name1, dpi=300, height=5, width=6)

# p1 <- ggplot(all_data, aes(y= ProbOcc, x = increment)) +
#   geom_smooth() +
#   facet_wrap(~FFM, scales = "free_x")
# 
# p1
# 
# file.name1 <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_v2/Figures/02_incremental_preds.jpg"
# ggsave(p1, filename=file.name1, dpi=300, height=5, width=6)

## partial dependence plots
sp=0.6
pdf("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Arroyo_toad_RB9_V2/Figures/PartialPlots_6_no_clim.pdf", width=8, height=8)
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



# Plotting probability of occurrence --------------------------------------

nhd <- st_read("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/Cannabis_Eflows/ignore/NHDPlus_V2_Flowline_CA.shp")
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

mapshot(m1, url = paste0(getwd(), "/ignore/map_no_clim.html"),
        file = paste0(getwd(), "/ignore/map_no_clim.png"))
getwd()



# Checking hydro data -----------------------------------------------------

## the hydro curves aren't polynomial over zero delta, why?

head(pred_env)

## check predictions, prob of occurece
pred_ev <- pred_env %>%
  # dplyr::select(probOcc, DS_Mag_50:Wet_BFL_Mag_10) %>%
  pivot_longer(DS_Mag_50:Wet_BFL_Mag_10, names_to = "FFM", values_to = "Value")

head(pred_ev)

ggplot(data = pred_ev, aes(y=probOcc, x=Value)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~FFM, scales="free_x")

?geom_smooth

plot(pred_ev$DS_Mag_50, pred_ev$probOcc)

## check known occurrences

orgi_dat_wide <- NewDataObsSub %>%
  dplyr::select(NewObs, DS_Mag_50:Wet_BFL_Mag_10)

orgi_dat <- NewDataObsSub %>%
  dplyr::select(NewObs, DS_Mag_50:Wet_BFL_Mag_10) %>%
  pivot_longer(DS_Mag_50:Wet_BFL_Mag_10, names_to = "FFM", values_to = "Value")

ggplot(data = subset(pred_ev, !NewObs == -999), aes(y=NewObs, x=Value)) +
  geom_point() +
  geom_smooth( method = loess, se = F) +
  facet_wrap(~FFM, scales="free_x")



