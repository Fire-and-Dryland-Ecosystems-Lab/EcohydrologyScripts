# clim_compare
# 
# get snowtel data, compare to netcdf clim and aggregated clim
# Added compare with Nikkis generated clim

# install.packages("snotelr")

library(RHESSysIOinR)
library(snotelr)
library(tidyverse)
library(ncdf4)
library(data.table)
library(cowplot)
options(scipen = 999)
# ============================== Inputs ==============================

# basestation clim (aggregated from netcdf)
input_base_clim = "clim/ward"

#snotel
 sno_meta = snotel_info()
find_site = sno_meta[which(grepl("ward creek", sno_meta$site_name)),]

# netcdf
input_nc_pcp = "clim/crop_agg_met_gridmet_pr_1979-01-01_2024-09-30.nc"
input_nc_tmin = "clim/crop_agg_met_gridmet_tmmn_1979-01-01_2024-09-30.nc"
input_nc_tmax = "clim/crop_agg_met_gridmet_tmmx_1979-01-01_2024-09-30.nc"

# daymet
input_nc_pcp_day = "clim/agg_met_daymet_prcp_1980-01-01_2023-12-31.nc"
input_nc_tmin_day = "clim/agg_met_daymet_tmin_1980-01-01_2023-12-31.nc"
input_nc_tmax_day = "clim/agg_met_daymet_tmax_1980-01-01_2023-12-31.nc"

# daymet
input_daymet = "clim/daymet/wepp_cli_edit.txt"

# ============================== Data ingest ==============================
good_vars = function(df) {
  if (any(!c("date", "rain", "tmin", "tmax","source") %in% names(df))) {
    stop("wrong col names")
  }
  df = df[,names(df) %in% c("date", "rain", "tmin", "tmax","source")]
  if (!isa(df$date, "Date")) {
    df$date = as.Date(df$date)
  }
  if (max(df$rain, na.rm = T) < 1) { #assume in meters
    df$rain = df$rain * 1000
  }
  if (mean(df$tmin, na.rm = T) > 200) {
    df$tmin =  df$tmin - 273.15
  }
  if (mean(df$tmax, na.rm = T) > 200) {
    df$tmax =  df$tmax - 273.15
  }
  return(df[,c("date", "rain", "tmin", "tmax","source")])
}

# ------------------------------ basestation/aggregated from netcdf ------------------------------ #
agg_clim = read_clim(input_base_clim)
agg_clim$source = "basestation_gridmet"

# ------------------------------ snotel ------------------------------ #
snotel = snotel_download(find_site$site_id, internal = T)
names(snotel)[names(snotel) == "precipitation"] = "rain"
names(snotel)[names(snotel) == "temperature_min"] = "tmin"
names(snotel)[names(snotel) == "temperature_max"] = "tmax"
snotel$source = "snotel"

# ------------------------------ netcdf ------------------------------ #
pr_edit_nc = nc_open(input_nc_pcp)
tmin_edit_nc = nc_open(input_nc_tmin)
tmax_edit_nc = nc_open(input_nc_tmax)

dat = ncvar_get(tmax_edit_nc, attributes(tmax_edit_nc$var)$names[1])
# dat[2:4,2:3,]
tmaxavg = apply(dat, 3, mean)
# tmax_df = as.data.frame((t(dat)))
dat = ncvar_get(tmin_edit_nc, attributes(tmin_edit_nc$var)$names[1])
tminavg = apply(dat, 3, mean)
# tmin_df = as.data.frame((t(dat)))
dat = ncvar_get(pr_edit_nc, attributes(pr_edit_nc$var)$names[1])
# rain_df = as.data.frame((t(dat)))
pravg = apply(dat, 3, mean)

# nc_clim = data.frame(rain = (rain_df$V1+rain_df$V2)/2,
#                      tmin = (tmin_df$V1+tmin_df$V2)/2,
#                      tmax = (tmax_df$V1+tmax_df$V2)/2)
nc_clim = data.frame(rain = pravg,
                     tmin = tminavg,
                     tmax = tmaxavg)
nc_clim$date = NA

days = ncvar_get(pr_edit_nc, "time")
dates =  as.Date("1900-01-01") + days
nc_clim$date = dates
nc_clim$source = "netcdf_gridmet"

nc_close(pr_edit_nc)
nc_close(tmin_edit_nc)
nc_close(tmax_edit_nc)

# ======================== DAYMET ================================
pr_edit_ncday = nc_open(input_nc_pcp_day)
tmin_edit_ncday = nc_open(input_nc_tmin_day)
tmax_edit_ncday = nc_open(input_nc_tmax_day)

dat = ncvar_get(tmax_edit_nc, attributes(tmax_edit_nc$var)$names[4])
tmaxavg = apply(dat, 3, mean)
# tmax_df = as.data.frame((t(dat)))
dat = ncvar_get(tmin_edit_nc, attributes(tmin_edit_nc$var)$names[4])
tminavg = apply(dat, 3, mean)
# tmin_df = as.data.frame((t(dat)))
dat = ncvar_get(pr_edit_nc, attributes(pr_edit_nc$var)$names[4])
# rain_df = as.data.frame((t(dat)))
pravg = apply(dat, 3, mean)

# nc_clim = data.frame(rain = (rain_df$V1+rain_df$V2)/2,
#                      tmin = (tmin_df$V1+tmin_df$V2)/2,
#                      tmax = (tmax_df$V1+tmax_df$V2)/2)
nc_clim_day = data.frame(rain = pravg,
                     tmin = tminavg,
                     tmax = tmaxavg)
                     nc_clim_day$date = NA

days = ncvar_get(pr_edit_nc, "time")
dates =  as.Date("1950-01-01") + days
nc_clim_day$date = dates
nc_clim_day$source = "netcdf_daymet"

nc_close(pr_edit_ncday)
nc_close(tmin_edit_ncday)
nc_close(tmax_edit_ncday)

# ====================================================================

writeclim = F
if (writeclim) {
  nc_clim = good_vars(nc_clim)
  outname = "clim/ward_netcdfgridmet_agg"
  startdate = "1979 01 01 01"
  
  write(paste0(startdate,"\n", paste0(format(nc_clim$tmin,nsmall = 1),collapse = "\n")) , file = paste0(outname,".tmin" ))
  write(paste0(startdate,"\n", paste0(format(nc_clim$tmax,nsmall = 1),collapse = "\n")) , file = paste0(outname,".tmax" ))
  write(paste0(startdate,"\n", paste0(nc_clim$rain/1000,collapse = "\n")) , file = paste0(outname,".rain" ))
  
  base_source = "clim/ward.base"
  base_in = readLines(base_source)
  base_out = gsub(gsub(".base","", base_source), outname, base_in)
  writeLines(base_out,con = paste0(outname,".base"))
}

# ------------------------------ DAYMET ------------------------------ #
daymet = read.table(file = input_daymet, header = T)
daymet$date =  as.Date(paste(daymet$year, daymet$month, daymet$day, sep = "-"))
names(daymet)[names(daymet) == "prcp_mm"] = "rain"
daymet$source = "daymet"

writeclim = F
if (writeclim) {
  mindate = min(daymet$date)
  date_head = "1990 01 01 01"
  clim_files = c("rain", "tmin", "tmax")
  fname = paste0("clim/ward_daymet", ".", clim_files)
  lapply(fname, write, x = date_head)
  data.table::fwrite(as.data.frame(daymet$prcp_mm/1000), file = fname[1], append = T)
  data.table::fwrite(as.data.frame(daymet$tmin), file = fname[2], append = T)
  data.table::fwrite(as.data.frame(daymet$tmax), file = fname[3], append = T)
  base_in = readLines("clim/ward.base")
  base_out = gsub("ward", "ward_daymet", base_in)
  writeLines(base_out,con = "clim/ward_daymet.base")
}

# ============================== Combine ==============================

agg_clim = good_vars(agg_clim)
snotel = good_vars(snotel)
nc_clim = good_vars(nc_clim)
daymet = good_vars(daymet)

climate = rbind(agg_clim, snotel, nc_clim, daymet)

climate = rbind(nc_clim, nc_clim_day)

climate$date = as.POSIXlt(climate$date)
climate$year = climate$date$year + 1900
climate$month = climate$date$mon + 1
climate$day = climate$date$mday
climate$wy = ifelse(climate$month >= 10, climate$year + 1, climate$year)
climate$yd = lubridate::yday(climate$date)
climate$wyd = get_waterYearDay(climate$date)
climate$yearmon = zoo::as.yearmon(climate$date, "%Y %m")

# ============================== Plots ==============================
pr = climate %>% group_by(source, year) %>% summarise(Precipitation = sum(rain)) %>%
  ggplot() +
  geom_line(aes(x = year, y = Precipitation, color = source) ) +
  labs(title = "Cumulative Annual Precipitation", x = "Year" ) +
  theme(legend.title=element_blank())

tmin = climate %>% group_by(source, year) %>% summarise(`Minimum Temperature` = mean(tmin)) %>%
  ggplot() +
  geom_line(aes(x = year, y = `Minimum Temperature`, color = source) ) +
  labs(title = "Average Daily Minimum Temperature", x = "Year" ) +
  theme(legend.title=element_blank())

tmax = climate %>% group_by(source, year) %>% summarise(`Maximum Temperature` = mean(tmax)) %>%
  ggplot() +
  geom_line(aes(x = year, y = `Maximum Temperature`, color = source) ) +
  labs(title = "Average Daily Maximum Temperature", x = "Year" ) +
  theme(legend.title=element_blank())

grid = plot_grid(pr, tmin, tmax, ncol=1)

saveplots = F
if (saveplots) {
  ggsave("plots/clim_compare.jpg", grid, width = 10, height = 12)
}

# ============================== Added gridmet-based ml data from Nikki ==============================
# gridmet-based data from nikki
input_ml_gridmet = "../data/20230427_will_results.csv"
gridmet_ml = read.csv(input_ml_gridmet)
unique(gridmet_ml$name)

gridmet_ml_ward = gridmet_ml %>% 
  filter(name == "WARD C AT HWY 89 NR TAHOE PINES CA", quantile == 50) %>%
  select(!c("Num"))
gridmet_ml_ward$source = "Gridmet_ml"
gridmet_ml_ward[,c("X","Y","name","quantile")] = NULL
gridmet_ml_ward = pivot_wider(gridmet_ml_ward,names_from = "variable",values_from = "prediction")
names(gridmet_ml_ward)[names(gridmet_ml_ward) == "pr"] = "rain"
gridmet_ml_ward$tmin = (gridmet_ml_ward$tmin - 32) * (5/9)
gridmet_ml_ward$tmax = (gridmet_ml_ward$tmax - 32) * (5/9)
summary(gridmet_ml_ward)
summary(climate)

climate = nc_clim
climate = add_dates(climate)

climate_mn = climate %>% group_by(source, year, month) %>% summarise(tmax = mean(tmax), tmin = mean(tmin), rain = mean(rain))

climate_mn = climate %>% group_by(source, year, month) %>% summarise(tmax = mean(tmax), tmin = mean(tmin), rain = sum(rain))

climate_mn = rbind(climate_mn, gridmet_ml_ward)
climate_mn$year_month = zoo::as.yearmon(paste0(climate_mn$year,"-",climate_mn$month))

climate_yr = climate %>% group_by(source, year) %>% summarise(tmax = mean(tmax), tmin = mean(tmin), rain = sum(rain))

climate_yr %>% group_by(source) %>% summarise(Precipitation = mean(rain))


pr2 = climate_mn %>% group_by(source, month) %>% summarise(Precipitation = mean(rain)) %>%
  ggplot() +
  geom_col(aes(x = month, y = Precipitation, fill = source),position="dodge" ) +
  labs(title = "Mean Cumulative Monthly Precipitation", x = "Year" ) +
  theme(legend.title=element_blank())

tmin = climate_mn %>% group_by(source, year) %>% filter(year%in%c(1980,1990,2000,2010,2020)) %>% summarise(`Minimum Temperature` = mean(tmin)) %>%
  ggplot() +
  geom_col(aes(x = year, y = `Minimum Temperature`, fill = source),position="dodge" ) +
  labs(title = "Average Daily Minimum Temperature", x = "Year" ) +
  theme(legend.title=element_blank())

tmax = climate_mn %>% group_by(source, year) %>% filter(year%in%c(1980,1990,2000,2010,2020)) %>% summarise(`Maximum Temperature` = mean(tmax)) %>%
  ggplot() +
  geom_col(aes(x = year, y = `Maximum Temperature`, fill = source),position="dodge" ) +
  labs(title = "Average Daily Maximum Temperature", x = "Year" ) +
  theme(legend.title=element_blank())

cols_multi = plot_grid(pr, tmin, tmax, ncol=1)

pr2 = climate_mn %>% group_by(source, year, month) %>% filter(year %in% c(1980,1990,2000,2010,2020)) %>% summarise(Precipitation = sum(rain)) %>%
  ggplot() +
  geom_line(aes(x = month, y = Precipitation, color = source), linewidth=1) +
  labs(title = "Cumulative Annual Precipitation", x = "Year Month" ) +
  facet_wrap(vars(year)) +
  theme(legend.title=element_blank())

tmin2 = climate_mn %>% group_by(source, year, month) %>% filter(year%in%c(1980,1990,2000,2010,2020)) %>% summarise(`Minimum Temperature` = mean(tmin)) %>%
  ggplot() +
  geom_line(aes(x = month, y = `Minimum Temperature`, color = source), linewidth=1) +
  labs(title = "Average Daily Minimum Temperature", x = "Year" ) +
  facet_wrap(vars(year)) +
  theme(legend.title=element_blank())

tmax2 = climate_mn %>% group_by(source, year, month) %>% filter(year%in%c(1980,1990,2000,2010,2020)) %>% summarise(`Maximum Temperature` = mean(tmax)) %>%
  ggplot() +
  geom_line(aes(x = month, y = `Maximum Temperature`, color = source), linewidth=1 ) +
  labs(title = "Average Daily Maximum Temperature", x = "Year" ) +
  facet_wrap(vars(year)) +
  theme(legend.title=element_blank())

#grid2 = plot_grid(pr2, tmin2, tmax2, ncol=1)
saveplots = F
if (saveplots) {
  ggsave("plots/clim_annual.jpg", cols_multi, width = 10, height = 12)
  ggsave("plots/clim_pr_mn.jpg", pr2, width = 10, height = 12)
  ggsave("plots/clim_tmin_mn.jpg", tmin2, width = 10, height = 12)
  ggsave("plots/clim_tmax_mn.jpg", tmax2, width = 10, height = 12)
}


pr = nc_clim %>% group_by(source, year) %>% summarise(Precipitation = sum(rain)) %>%
  ggplot() +
  geom_line(aes(x = year, y = Precipitation, color = source) ) +
  labs(title = "Cumulative Annual Precipitation", x = "Year" ) +
  theme(legend.title=element_blank())