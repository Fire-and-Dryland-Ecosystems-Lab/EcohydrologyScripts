nbalance = function(rhout)
 
 rhout = res2

 beg.bdg = subset(rhout$bdg, rhout$bdg$day==1 & rhout$bdg$month==1)
 end.bdg = subset(rhout$bdg, rhout$bdg$day==31 & rhout$bdg$month==12)
nb = nrow(beg.bdg)
ne = nrow(end.bdg)
if (beg.bdg$wy[1] > end.bdg$wy[1])
	end.bdg = end.bdg[2:ne,]
nb = nrow(beg.bdg)
ne = nrow(end.bdg)
if (end.bdg$wy[1] > beg.bdg$wy[1])
	beg.bdg = beg.bdg[2:ne,]
nb = nrow(beg.bdg)
ne = nrow(end.bdg)
if (beg.bdg$wy[nb] > end.bdg$wy[ne])
	beg.bdg = beg.bdg[1:ne,]
nb = nrow(beg.bdg)
ne = nrow(end.bdg)
if (beg.bdg$wy[nb] < end.bdg$wy[ne])
	end.bdg = end.bdg[1:nb,]
nb = nrow(beg.bdg)
ne = nrow(end.bdg)

 delta.bdg = end.bdg-beg.bdg
 delta.bdg$wy = beg.bdg$wy
rhout$bdg.wy = subset(rhout$bdg.wy, rhout$bdg.wy$wy %in% delta.bdg$wy)
 
 years = unique(rhout$bdg.wy$year)
 nyears = length(years)
 nbalance = as.data.frame(matrix(nrow=nyears), ncol=1)
colnames(nbalance) = "year"
nbalance$year = years

nbalance$streamN =  rhout$bdg.wy$streamflow_NO3+rhout$bdg.wy$streamflow_DON + rhout$bdg.wy$streamflow_NH4
nbalance$atmN = rhout$bdg.wy$N_dep + rhout$bdg.wy$nfix
nbalance$dsoil = delta.bdg$soiln*1000.0 + delta.bdg$sminn + delta.bdg$litrn*1000.0 + delta.bdg$nitrate + delta.bdg$DON*1000.0 + delta.bdg$surfaceN
nbalance$dplant = delta.bdg$plantn*1000.0
nbalance$dgw = (delta.bdg$gwNO3 + delta.bdg$gwNH4+ delta.bdg$gwDON)*1000.0
nbalance$den = rhout$bdg.wy$denitrif

nbalance$orgN = (delta.bdg$plantn + delta.bdg$soiln + delta.bdg$litrn) * 1000.0
nbalance$Nout = nbalance$streamN + nbalance$den
nbalance$dstore = nbalance$dplant + nbalance$dgw + nbalance$dsoil
nbalance$net = nbalance$atmN - nbalance$Nout - nbalance$dstore

tmp = subset(rhout$bdg, rhout$bdg.wy==rhout$bdg.wy$wy[1])
plot(cumsum(tmp$streamflow_NO3 + tmp$streamflow_DON + tmp$streamflow_NH4 + tmp$denitrif))
lines(cumsum(tmp$ndep*1000), col="red")


nbalance.day = tmp[,c("year","month","day")]
tmp2 = as.data.frame(apply(tmp[,c("plantn","litrn","soiln","sminn","nitrate",
                                  "DON","surfaceN","gwNO3","gwNH4","gwDON")], 2, diff) )

nbalance.day$atmN = tmp$N_dep*1000 + tmp$nfix
nbalance.day$plantn = c(0, tmp2$plantn*1000.0)
nbalance.day$osoiln = c(0, tmp2$litrn*1000.0 + tmp2$soiln*1000.0)
nbalance.day$msoiln = c(0, tmp2$sminn + tmp2$nitrate + tmp2$DON*1000.0) 
nbalance.day$surfn = c(0, tmp2$surfaceN)
nbalance.day$gwn = c(0, (tmp2$gwNO3 + tmp2$gwNH4 + tmp2$gwDON)*1000)
nbalance.day$gwout = (tmp$gwNO3out + tmp$gwNH4out + tmp$gwDONout)
nbalance.day$Nout= tmp$streamflow_NO3 + tmp$streamflow_DON + tmp$streamflow_NH4 + tmp$denitrif + nbalance.day$gwout


nbalance.day$store = (nbalance.day$plantn + nbalance.day$osoiln + nbalance.day$msoiln + nbalance.day$gwn + nbalance.day$surfn)
nbalance.day$net = nbalance.day$atmN - nbalance.day$Nout - nbalance.day$store
nbalance.day = nbalance.day[2:365,]





