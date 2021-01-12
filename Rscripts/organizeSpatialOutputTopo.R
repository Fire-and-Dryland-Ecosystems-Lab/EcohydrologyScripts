#title: "organizeSpatialOutputTopo"
#author: "ErinJHanan"

#################################################################################################
# Creates a data fram with separate columns for over and understory and adds vegetation type
# Full basin after initializing soils for each veg type

library(plyr)

setwd("...")

patch = readin_rhessys_output("yourOutputPrefix", c=1)
overstory = subset(patch$cdg, patch$cdg$stratumID < 1000000)
understory = subset(patch$cdg, patch$cdg$stratumID >= 1000000)


#--------

# read in, add basin, zone, stratum. Replace names with the names for your maps
hill_ID=scan(file="../../auxdata.subset/hillMap.asc", skip=6, na.strings="*")
patch_ID=scan(file="../../auxdata.subset/patchMap.asc", skip=6, na.strings="*")
zone_ID = patch_ID
stratum_ID = patch_ID
landcover=scan(file="../../auxdata.subset/landcoverMap.asc", skip=6, na.strings="*")
elevation=scan(file="../../auxdata.subset/dem30m.asc", skip=6, na.strings="*")
slope=scan(file="../../auxdata.subset/slope.asc", skip=6, na.strings="*")
aspect=scan(file="../../auxdata.subset/aspect.asc", skip=6, na.strings="*")

x=rep(1,length(hill_ID)) # uses the number of cases
basin_ID=x

#combine them
tmp = as.data.frame(cbind(basin_ID,hill_ID,zone_ID,patch_ID,stratum_ID,landcover,elevation,slope,aspect))
tmp2 = subset(tmp, is.na(tmp$patch_ID)==F)
tmp = aggregate(tmp2, by=list(tmp2$aspect,tmp2$slope,tmp2$elevation,tmp2$landcover, tmp2$stratum_ID, tmp2$patch_ID, tmp2$zone_ID, tmp2$hill_ID, tmp2$basin_ID), mean)

#reorder them, remove subset_key
tmp=tmp[,c(10:18)]
locations = tmp

# create data frame with outputs
stratum_ID = patch$cdg$stratumID
patch_ID = patch$cdg$patchID
LAI = patch$cdg$proj_lai
layer = rep(1:2, length(stratum_ID))

data = as.data.frame(cbind(patch_ID, stratum_ID,LAI,layer))
over.sub = subset(data, data$layer==1)
under.sub = subset(data, data$layer==2)
names(over.sub)[names(over.sub)=="LAI"] <- "overstoryLAI"
names(under.sub)[names(under.sub)=="LAI"] <- "understoryLAI"
data = merge(over.sub,under.sub,by="patch_ID")
data=data[, c(1,2,3,5,6)]

# merge outputs with locations
results=merge(locations,data,by="patch_ID")
results=results[, c(2,3,4,1,5,6,7,8,9,11,12,13)]
#names(results)[names(results)=="stratum_ID.x"] <- "Stratum_ID"
names(results)[names(results)=="stratum_ID.y"] <- "understory_stratum_ID"

#-------

under_patches=subset(results, results$landcover==42 | results$landcover==41)

stats.results <- ddply(results, c("landcover"), summarise,
                       N    = length(understoryLAI),
                       mean = mean(understoryLAI),
                       sd   = sd(understoryLAI),
                       se   = sd / sqrt(N))
stats.results


pine=subset(results, results$landcover==42)
plot(pine$patch_ID, pine$overstoryLAI, col="green", ylim=c(0,6))
points(pine$patch_ID, pine$understoryLAI, col="blue")

deciduous=subset(results, results$landcover==41)
plot(deciduous$patch_ID, deciduous$overstoryLAI, col="green", ylim=c(0,6))
points(deciduous$patch_ID, deciduous$understoryLAI, col="blue")

plot(pine$elevation, pine$understoryLAI, col="green")
plot(pine$slope, pine$understoryLAI, col="green")
plot(pine$aspect, pine$understoryLAI, col="green")

good=subset(deciduous, deciduous$understoryLAI>1)
good=subset(good, good$overstoryLAI>2) 
summary(good$overstoryLAI)
