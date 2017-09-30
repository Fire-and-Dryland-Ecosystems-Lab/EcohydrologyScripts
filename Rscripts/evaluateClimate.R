#####################################################################################
# Just Seeders NoFire, read in data and plot vegetation vs. precip --- this is for runs with 2 canopy layers


sdglNF=read.table("NoFire_grow_stratum.daily", header=T)
sdglNF=mkdate(sdglNF)
s1NF = subset(sdglNF, sdglNF$stratumID == 39393)
s2NF = subset(sdglNF, sdglNF$stratumID != 39393)
bdglNF=read.table("NoFire_grow_basin.daily", header=T)
bdglNF=mkdate(bdglNF)
bdlNF=read.table("NoFire_basin.daily", header=T)
bdlNF=mkdate(bdlNF)
zdglNF=read.table("NoFire_grow_zone.daily", header=T)
zdglNF=mkdate(zdglNF)
pdglNF=read.table("NoFire_grow_zone.daily", header=T)
pdglNF=mkdate(pdglNF)

par(mar=c(5,5,5,5))
plot(s1NF$date, s1NF$proj_lai, type="l", col="darkgreen", yaxt="n", main="May No Fire Seeders", ylim=c(0,7), xlab="Date", ylab="")
lines(s2NF$date, s2NF$proj_lai, col="green3")
axis(side=2, at=c(0,1,2,3,4,5))
mtext("LAI", side=2, line=2.75, at=0.25, col="darkgreen")
mtext("Shrubs", side=2, line=2, at=2.5, cex=0.75, col="darkgreen")
mtext("Herbs", side=2, line=2, at=0.25, cex=0.75, col="green3")
par(new=T)
plot(zdglNF$date, zdglNF$rain, type="l", yaxt="n", xaxt="n", ylim=c(500,0), col="dodgerblue2", xlab="", ylab="")
axis(side=4, at=c(0, 100, 200))
mtext("Rainfall (mm)", side = 4, line=2.5, at=125, col="dodgerblue2")

plot(s1NF$date, s1NF$root_depth, type="l")

#############################################################################################
# plot relationship between precip and transpiration
plot(bdlNF$precip, bdlNF$trans)

# Aggregate precip by year and subset for 2007-2021
precip.wy = aggregate(bdlNF$precip, by=list(bdlNF$wy), sum)
precip.wy <- subset(precip.wy, Group.1 == 2007:2021, select=c(Group.1, x))
plot(precip.wy$Group.1,precip.wy$x, type = "h")

# Aggregate transpiration by year and subset for 2007-2021
trans.wy = aggregate(bdlNF$trans, by=list(bdlNF$wy), sum)
trans.wy <- subset(trans.wy, Group.1 == 2007:2021, select=c(Group.1, x))
plot(trans.wy$Group.1,trans.wy$x, type = "h")

# Aggregate streamflow by year and subset for 2007-2021
streamflow.wy = aggregate(bdlNF$streamflow, by=list(bdlNF$wy), sum)
streamflow.wy <- subset(streamflow.wy, Group.1 == 2007:2021, select=c(Group.1, x))
# plot(streamflow.wy$Group.1,streamflow.wy$x, type = "h")

# Aggregate evaporation by year and subset for 2007-2021
evap.wy = aggregate(bdlNF$evap, by=list(bdlNF$wy), sum)
evap.wy <- subset(evap.wy, Group.1 == 2007:2021, select=c(Group.1, x))
# plot(evap.wy$Group.1,evap.wy$x, type = "h")


# Merge by water year
tmp = merge(precip.wy, trans.wy, by=c("Group.1"))
tmp2 = merge(streamflow.wy, evap.wy, by=c("Group.1"))
tmp3 = merge(tmp, tmp2, by=c("Group.1"))

names(tmp3)[1]<-"wy"
names(tmp3)[2]<-"precip"
names(tmp3)[3]<-"trans"
names(tmp3)[4]<-"streamflow"
names(tmp3)[5]<-"evap"

tmp3$et=tmp3$evap+tmp3$trans
climate=tmp3


#calculate ET
et.wy = trans.wy+evap.wy

# Plot relationships
names = c("2007","2008","2009","2010","2011","2012","2013","2014","2015","2016",
          "2017", "2018", "2019", "2020", "2021")
plot(climate$precip, climate$trans, xlim=c(-10,1450), ylim=c(400,550))
# text(climate$precip, climate$trans, labels=names, cex=0.6, pos=4, col="red") 
# spread.labels(climate$precip, climate$trans, labels=names, offset = F, ony= TRUE, 
#          between=FALSE, linecol=par("fg"), srt=0)

plot(climate$precip, climate$streamflow, xlim=c(-10,1300))
text(climate$precip, climate$streamflow, labels=names, cex=0.6, pos=4, col="red") 



