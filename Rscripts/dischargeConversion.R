# from Tague Lab

setwd("...")

# conversion from cfs to mm

obs = read.table("streamflow.jc.1928_2015.txt", header=T)
#obs = flow
area.mile2 = 218 
mil_to_foot = 27878400
ft_to_mm = 304.8

obs$mm = obs$cfs/(area.mile2*mil_to_foot)*ft_to_mm*60*60*24

obs$date = seq.dates(from="10/01/1928",length=length(obs$mm))

obs$year = as.numeric(as.character(years(obs$date)))
obs$month = as.numeric(months(obs$date))
obs$day = as.numeric(days(obs$date))
obs$wateryear = ifelse((obs$month >= 10), obs$year+1, obs$year) 

# output for caliba
write.table(obs[,c(1,6)], file="obs2.txt", row.names=F, col.names=F)

