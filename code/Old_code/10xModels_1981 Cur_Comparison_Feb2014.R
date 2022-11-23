## final rasters current to historic comparison 

library(sp)
library(rgdal)
library(raster)

havg<-raster("ignore/Historic/randomForests/PARTITIONING/Hist_10x_Final_1981/hbAvg.img")
cavg<-raster("ignore/Current/randomForests/PARTITIONING/Cur_10x_Final_Feb2014/cbAvg.img")

havg
C.predPA<-cavg
C.predPA[C.predPA>=0.492280772]<-1
C.predPA[C.predPA<0.492280772]<-0

freq(C.predPA)
     # value  count
# [1,]     0  41443
# [2,]     1   4862
# [3,]    NA 408644

H.predPA<-havg
H.predPA[H.predPA>=0.43512684]<-1
H.predPA[H.predPA<0.43512684]<-0

freq(H.predPA)
     # value  count
# [1,]     0  39650
# [2,]     1   6655
# [3,]    NA 408644

S.predPA<-stack(H.predPA, C.predPA)

freq(S.predPA)

plot(S.predPA)


#Create Transition Raster
trans <- S.predPA[[1]]-S.predPA[[2]]
freq(trans)
freq(trans)
     # value  count
# [1,]    -1   1467
# [2,]     0  41578
# [3,]     1   3260
# [4,]    NA 408644


##############
#-1 =  1467/46305 = 0.03168124  OF TOTAL STREAM SECTIONS
# 1 =  3260/46305 = 0.07040276
# -1 and 1 = 0.102084
#If all 3260 pixels transformed to habitat, that would be an increase in habitat of 67.02303% (3260/4864)
pal<-heat.colors(3)
pal<-terrain.colors(3, alpha = 1)
pal<-cm.colors(3, alpha = 1)
pal<-bpy.colors(3, alpha = 1)

## object n not found
topo.colors(n, alpha = 1)
cm.colors(n, alpha = 1)

pal<-colorRampPalette(c("red", "black", "green"))( 3 ) ## (n)

plot(trans, col=pal, colNA="gray")

writeRaster(trans, "D:/AT DistributionModeling/TransitionRast_Feb2014.img")

plot(trans, col=pal, colNA="gray")
freq(trans)
help(plot)
