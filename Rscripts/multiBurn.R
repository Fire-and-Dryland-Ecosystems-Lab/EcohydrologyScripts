#title: "multiBurn"
#authors: "ErinJHanan and NaomiTague"

# creates box plots for all simulations (for each year) -- shows inter-simulation variability over time

nday = 365*20+5
nsims = 15
start.yrs = seq(from=2007, to=2021)
nyrs = 20
lai.all = as.data.frame(matrix(nrow=nday, ncol=length(start.yrs)))

for (i in 1:length(start.yrs)) {

curr.year = start.yrs[i]
iname = sprintf("bdgml%d",start.yrs[i])
curr = get(iname)

begin.year = curr.year - 5
end.year = begin.year+20-1

tmp = subset(curr, curr$year %in% c(begin.year:end.year))

lai.all[,i] = tmp$lai
}


lai.all$date = seq.dates(from="1/1/2000", length=length(lai.all$V1))
lai.all$month = as.integer(months(lai.all$date))
lai.all$year = as.integer(as.character(years(lai.all$date)))

tmp = aggregate(lai.all, by=list(lai.all$month,lai.all$year), mean)
lai.all.myr = tmp
lai.all.myr$date = 15
#lai.all.myr$date = as.Date(paste(lai.all.myr$year, lai.all.myr$month, lai.all.myr$day, sep=","))
lai.all.myr$date = lai.all.myr$year*12+lai.all.myr$month

tmp = stack(lai.all.myr[,3:(nsims+2)])
tmp$date = rep(lai.all.myr$date, times=nsims)

boxplot(tmp$values~tmp$date)


tmp = aggregate(lai.all, by=list(lai.all$year), mean)
lai.all.yr = tmp
tmp = stack(lai.all.yr[,2:(nsims+1)])
tmp$date = rep(lai.all.yr$year, times=nsims)
tmp$yr = tmp$date-2000-5
boxplot(tmp$values~tmp$yr)




