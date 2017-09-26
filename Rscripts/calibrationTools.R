# from Tague Lab
#Read in observed

library(chron)
setwd("...")

obs=read.table("obs2.txt", head=T)
# sometimes the date gets messed up if so, run this:
dt = seq.dates(from="10/1/1928", length=length(obs$mm))
obs$date=as.Date(dt)

# to get columns for year month day
# obs$year <- format(obs[,1], "%Y") # year
# obs$month <- format(obs[,1], "%m") # month
# obs$day <- format(obs[,1], "%d") # day

source("computecalstatsExtended.R")

# Run static calibration, for static trying to capture totals, so look at percent error
head(stats3d)

bst.stat.err = subset(stats3d, stats3d$perr.total > -1 & stats3d$perr.total < 1)

summary(stats3d)
summary(stats3d[,8:10])
bad = subset(stats3d, stats3d$perr.total > 50)

save.image("StatCal.RData")
ls()

plot(stats3d$m, stats3d$nse)
plot(stats3d$m, stats3d$nse, ylim=c(0,1))
plot(stats3d$m, stats3d$perr.total, ylim=c(-1,1))
plot(stats3d$k, stats3d$nse, ylim=c(0,1))
plot(stats3d$m, stats3d$nse)
plot(result3d$date, log(result3d$V1), type="l")
lines(result3d$date, log(result3d$mm), col="red")

###################################################################################################
#choose best 2 paramenter sets, run model 10 years

setwd("/work.flamenco/Projects/JohnsonCreek/out/spinup/cal")

bdg85=read.table("p85_grow_basin.daily", head=T)
bdg162=read.table("p162_grow_basin.daily", head=T)
bd162=read.table("p162_basin.daily", head=T)
bdg85=mkdate(bdg85)
bdg162=mkdate(bdg162)
bd162=mkdate(bd162)

plot(bdg85$date, bdg85$lai, type="l")

plot(bdg162$date, bdg162$lai, type="l")

plot(bd162$date, bd162$streamflow, type="l")
lines(bd162$date, bd162$precip, col="blue")

###################################################################################################
# Run dynamic calibration

bst = subset(stats3d, stats3d$nse > .6 & stats3d$nselog > .5 & stats3d$perr.total > -10 & stats3d$perr.total < 10)
bstnse = subset(stats3d, stats3d$nse > .6)
ts.plot(result3d[,bstnse$row+1])


plot(stats3d$m, stats3d$nse)
plot(stats3d$m, stats3d$nse, ylim=c(0,1))
plot(stats3d$m, stats3d$perr.total)
plot(stats3d$m, stats3d$perr.total, ylim=c(-1,1))

plot(stats3d$k, stats3d$nse, ylim=c(0,1))
plot(stats3d$k, stats3d$perr.total)

plot(stats3d$po, stats3d$nse)
plot(stats3d$po, stats3d$nse, ylim=c(0,1))
plot(stats3d$po, stats3d$perr.total)
plot(stats3d$po, stats3d$perr.total, ylim=c(-5,5))

plot(stats3d$pa, stats3d$nse)
plot(stats3d$pa, stats3d$nse, ylim=c(0,1))
plot(stats3d$pa, stats3d$perr.total)
plot(stats3d$pa, stats3d$perr.total, ylim=c(-5,5))

plot(stats3d$gw1, stats3d$nse)
plot(stats3d$gw1, stats3d$nse, ylim=c(0,1))
plot(stats3d$gw1, stats3d$perr.total)
plot(stats3d$gw1, stats3d$perr.total, ylim=c(-5,5))

plot(stats3d$gw2, stats3d$nse)
plot(stats3d$gw2, stats3d$nse, ylim=c(0,1))
plot(stats3d$gw2, stats3d$perr.total)
plot(stats3d$gw2, stats3d$perr.total, ylim=c(-5,5))

plot(result3d$date, (result3d$V1.29), type="l", ylim=c(0,15)) #171
#lines(result3d$date, (result3d$V2.44), col="blue") #171
lines(result3d$date, (result3d$mm), col="red") #observed

plot(cumsum(result3d$V1.29), type="l", ylim=c(0,3000)) #171
#lines(result3d$date, (result3d$V2.44), col="blue") #171
lines(cumsum(result3d$mm), col="red") #observed


plot(result3d$date, (result3d$V1.58), type="l", ylim=c(0,15)) #171
#lines(result3d$date, (result3d$V2.44), col="blue") #171
lines(result3d$date, (result3d$mm), col="red") #observed

plot(cumsum(result3d$V1.58), type="l", ylim=c(0,3000)) #171
#lines(result3d$date, (result3d$V2.44), col="blue") #171
lines(cumsum(result3d$mm), col="red") #observed
###################################################################################################
# plot outputs from last cal run
setwd("/work.flamenco/Projects/JohnsonCreek/out/cal/sims")
bdsen3=read.table("sen3_basin.daily", head=T)
bdsen3=mkdate(bdsen3)
plot(bdsen3$date, bdsen3$snowpack/100, type="l", ylim=c(0,15))
lines(result3d$date, (result3d$V22.2), type="l", col="green") # best = sim 72
lines(result3d$date, (result3d$mm), col="red") #observed

plot(bdsen3$date, log(bdsen3$snowpack), type="l")
lines(result3d$date, log(result3d$V22.3), type="l", col="green") # simulation 25 in cal file 1
lines(result3d$date, log(result3d$mm), col="red") #observed

###################################################################################################
# test runs varying snow vs. rain and timing of melt
setwd("/work.flamenco/Projects/JohnsonCreek/out/test/snow")
bd1.1=read.table("1.1_basin.daily", head=T)
bd2.1=read.table("2.1_basin.daily", head=T)
bd3.1=read.table("3.1_basin.daily", head=T)
bd4.1=read.table("4.1_basin.daily", head=T)
bd5.1=read.table("5.1_basin.daily", head=T)

bd1.2=read.table("1.2_basin.daily", head=T)
bd2.2=read.table("2.2_basin.daily", head=T)
bd3.2=read.table("3.2_basin.daily", head=T)
bd4.2=read.table("4.2_basin.daily", head=T)
bd5.2=read.table("5.2_basin.daily", head=T)

bd1.3=read.table("1.3_basin.daily", head=T)
bd2.3=read.table("2.3_basin.daily", head=T)
bd3.3=read.table("3.3_basin.daily", head=T)
bd4.3=read.table("4.3_basin.daily", head=T)
bd5.3=read.table("5.3_basin.daily", head=T)

bd1.1=mkdate(bd1.1)
bd2.1=mkdate(bd2.1)
bd3.1=mkdate(bd3.1)
bd4.1=mkdate(bd4.1)
bd5.1=mkdate(bd5.1)

bd1.2=mkdate(bd1.2)
bd2.2=mkdate(bd2.2)
bd3.2=mkdate(bd3.2)
bd4.2=mkdate(bd4.2)
bd5.2=mkdate(bd5.2)

bd1.3=mkdate(bd1.3)
bd2.3=mkdate(bd2.3)
bd3.3=mkdate(bd3.3)
bd4.3=mkdate(bd4.3)
bd5.3=mkdate(bd5.3)

bdg1.1=read.table("1.1_grow_basin.daily", head=T)
bdg2.1=read.table("2.1_grow_basin.daily", head=T)
bdg3.1=read.table("3.1_grow_basin.daily", head=T)
bdg4.1=read.table("4.1_grow_basin.daily", head=T)
bdg5.1=read.table("5.1_grow_basin.daily", head=T)

bdg1.2=read.table("1.2_grow_basin.daily", head=T)
bdg2.2=read.table("2.2_grow_basin.daily", head=T)
bdg3.2=read.table("3.2_grow_basin.daily", head=T)
bdg4.2=read.table("4.2_grow_basin.daily", head=T)
bdg5.2=read.table("5.2_grow_basin.daily", head=T)

bdg1.3=read.table("1.3_grow_basin.daily", head=T)
bdg2.3=read.table("2.3_grow_basin.daily", head=T)
bdg3.3=read.table("3.3_grow_basin.daily", head=T)
bdg4.3=read.table("4.3_grow_basin.daily", head=T)
bdg5.3=read.table("5.3_grow_basin.daily", head=T)

bdg1.1=mkdate(bdg1.1)
bdg2.1=mkdate(bdg2.1)
bdg3.1=mkdate(bdg3.1)
bdg4.1=mkdate(bdg4.1)
bdg5.1=mkdate(bdg5.1)

bdg1.2=mkdate(bdg1.2)
bdg2.2=mkdate(bdg2.2)
bdg3.2=mkdate(bdg3.2)
bdg4.2=mkdate(bdg4.2)
bdg5.2=mkdate(bdg5.2)

bdg1.3=mkdate(bdg1.3)
bdg2.3=mkdate(bdg2.3)
bdg3.3=mkdate(bdg3.3)
bdg4.3=mkdate(bdg4.3)
bdg5.3=mkdate(bdg5.3)

plot(bd1.1$date, bd1.1$snowpack/100, type="l", ylim=c(0,15))
lines(bd1.1$date, (bd1.1$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd1.1$year < 2005], result3d$mm, col="red") #observed

plot(bd2.1$date, bd2.1$snowpack/100, type="l", ylim=c(0,15))
lines(bd2.1$date, (bd2.1$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd2.1$year < 2005], result3d$mm, col="red") #observed

plot(bd3.1$date, bd3.1$snowpack/100, type="l", ylim=c(0,15))
lines(bd3.1$date, (bd3.1$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd3.1$year < 2005], result3d$mm, col="red") #observed

plot(bd4.1$date, bd4.1$snowpack/100, type="l", ylim=c(0,15))
lines(bd4.1$date, (bd4.1$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd4.1$year < 2005], result3d$mm, col="red") #observed

plot(bd5.1$date, bd5.1$snowpack/100, type="l", ylim=c(0,15))
lines(bd5.1$date, (bd5.1$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd5.1$year < 2005], result3d$mm, col="red") #observed

plot(bd1.2$date, bd1.2$snowpack/100, type="l", ylim=c(0,15))
lines(bd1.2$date, (bd1.2$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd1.2$year < 2005], result3d$mm, col="red") #observed

plot(bd2.2$date, bd2.2$snowpack/100, type="l", ylim=c(0,15))
lines(bd2.2$date, (bd2.2$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd2.2$year < 2005], result3d$mm, col="red") #observed

plot(bd3.2$date, bd3.2$snowpack/100, type="l", ylim=c(0,15))
lines(bd3.2$date, (bd3.2$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd3.2$year < 2005], result3d$mm, col="red") #observed

plot(bd4.2$date, bd4.2$snowpack/100, type="l", ylim=c(0,15))
lines(bd4.2$date, (bd4.2$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd4.2$year < 2005], result3d$mm, col="red") #observed

plot(bd5.2$date, bd5.2$snowpack/100, type="l", ylim=c(0,15))
lines(bd5.2$date, (bd5.2$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd5.2$year < 2005], result3d$mm, col="red") #observed

plot(bd1.3$date, bd1.3$snowpack/100, type="l", ylim=c(0,15))
lines(bd1.3$date, (bd1.3$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd1.3$year < 2005], result3d$mm, col="red") #observed

plot(bd2.3$date, bd2.3$snowpack/100, type="l", ylim=c(0,15))
lines(bd2.3$date, (bd2.3$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd2.3$year < 2005], result3d$mm, col="red") #observed

plot(bd3.3$date, bd3.3$snowpack/100, type="l", ylim=c(0,15))
lines(bd3.3$date, (bd3.3$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd3.3$year < 2005], result3d$mm, col="red") #observed

plot(bd4.3$date, bd4.3$snowpack/100, type="l", ylim=c(0,15))
lines(bd4.3$date, (bd4.3$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd4.3$year < 2005], result3d$mm, col="red") #observed

plot(bd5.3$date, bd5.3$snowpack/100, type="l", ylim=c(0,15))
lines(bd5.3$date, (bd5.3$streamflow), type="l", col="blue") # modeled streamflow
lines(result3d$date[bd5.3$year < 2005], result3d$mm, col="red") #observed

plot(bdg5.1$date, bdg5.1$lai, type="l", col="black", ylim=c(0,2.5))
lines(bdg5.2$date, bdg5.2$lai, type="l", col="blue")
lines(bdg5.3$date, bdg5.3$lai, type="l", col="red")

plot(bd5.1$date, bd5.1$evap+bd5.1$trans, type="l", col="black")
lines(bd5.2$date, bd5.1$evap+bd5.1$trans, type="l", col="blue")
lines(bd5.3$date, bd5.1$evap+bd5.1$trans, type="l", col="red")

bd1.1$unfilled_cap = bd1.1$sat_def-bd1.1$rz_storage-bd1.1$unsat_stor
plot(bd1.1$date, bd1.1$unfilled_cap, type="l", col="black")


plot(bd5.1$date, bd5.1$precip, type="l", col="dodgerblue")
lines(bd5.1$date, bd5.1$streamflow, type="l", col="darkblue")
lines(result3d$date[bd5.3$year < 2005], result3d$mm, col="red") #observed

plot(bd5.2$date, bd5.2$precip, type="l", col="dodgerblue")
lines(bd5.2$date, bd5.2$streamflow, type="l", col="darkblue")
lines(result3d$date[bd5.3$year < 2005], result3d$mm, col="red") #observed

plot(bd5.3$date, bd5.3$precip, type="l", col="dodgerblue")
lines(bd5.3$date, bd5.3$streamflow, type="l", col="darkblue")
lines(result3d$date[bd5.3$year < 2005], result3d$mm, col="red") #observed

# calculate NSE individually for these runs

tmp = merge(bd4.2, obs, by=c("date"))
nse4.2 = nse(tmp$streamflow, tmp$mm)
nse4.2
log4.2 = lognse(tmp$streamflow, tmp$mm)
log4.2
err4.2 = mper.err(tmp$streamflow, tmp$mm)
err4.2

tmp = merge(bd3.2, obs, by=c("date"))
nse3.2 = nse(tmp$streamflow, tmp$mm)
nse3.2
log3.2 = lognse(tmp$streamflow, tmp$mm)
log3.2
err3.2 = mper.err(tmp$streamflow, tmp$mm)
err3.2

tmp = merge(bd2.2, obs, by=c("date"))
nse2.2 = nse(tmp$streamflow, tmp$mm)
nse2.2
log2.2 = lognse(tmp$streamflow, tmp$mm)
log2.2
err2.2 = mper.err(tmp$streamflow, tmp$mm)
err2.2

tmp = merge(bd5.2, obs, by=c("date"))
nse5.2 = nse(tmp$streamflow, tmp$mm)
nse5.2
log5.2 = lognse(tmp$streamflow, tmp$mm)
log5.2
err5.2 = mper.err(tmp$streamflow, tmp$mm)
err5.2

tmp = merge(bd5.3, obs, by=c("date"))
nse5.3 = nse(tmp$streamflow, tmp$mm)
nse5.3
log5.3 = lognse(tmp$streamflow, tmp$mm)
log5.3
err5.3 = mper.err(tmp$streamflow, tmp$mm)
err5.3


###########################################################
# Test different values for GW1

source("computecalstatsGW.Rscpt")

stats3d
head(result3d)

plot(result3d$date, (result3d$V1.29), type="l", col="blue", ylim=c(0,15)) #Original
lines(result3d$date, (result3d$V1.1), col="blue") #.12
lines(result3d$date, (result3d$V1.2), col="blue") #.012
lines(result3d$date, (result3d$V1.1), col="green") #.12 .. low snow T
lines(result3d$date, (result3d$V1.2), col="green") #.012 .. low snow T
lines(result3d$date, (result3d$V1.2), col="green") #.012 .. low snow T
lines(result3d$date, (result3d$mm), col="red") #observed
