# Get basin outlines from USGS gauge IDs
# for potential CARB sites

library(dataRetrieval)
library(waterData)
library(tidyverse)
library(nhdplusTools)
library(sf)
library(openxlsx)
library(terra)

# Ward Creek
# 10336676
site_ID = "10336676"

studysites = read.xlsx(xlsxFile = "~/../Box/CARB/SiteSelection/2022_OurProject/StudySitesAllProjects.xlsx", sheet = "All")

# ----- Get Streamflow -----
get_streams = F
if (get_streams) {
  site_ID = "10336676"
  sdate <- "1950-01-01" 
  edate   <- "2022-09-30" 
  code  <- "00060" # for streamflow - cubic ft/sec
  stat = "00003" # daily mean
  
  site_info <- readNWISsite(site_ID)
  # Download data
  q_cfs <- importDVs(site_ID, code = code, stat = stat,  
                     sdate = sdate, edate = edate)
  # mi2 to km2 
  basin_km2 = site_info$drain_area_va * 2.58998811
  # convert unit  first cubic feet /s to cubic m /s then to mm / day
  # convert cubic feet to cubic meter val*0.0283168
  observation = q_cfs %>% select(c(dates,val))
  observation = observation %>% mutate(mm=val*0.0283168*(1000*24*3600/basin_km2))
  observation$dates = as.Date(observation$dates)# select the time series or not
  # just to check
  na.count = which(is.na(observation$val))
  obs_file_out = ""
  write.table(observation,obs_file_out,sep="\t",col.names=F,row.names=F) ## val is cfs/s and  mm is mm/day
}

# ----- Get DEM & Basin -----
# https://rdrr.io/github/USGS-R/nhdplusTools/f/vignettes/plot_nhdplus.Rmd

# nwissite <- list(featureSource = "nwissite", featureID = paste0("USGS-",site_ID))
# basin <- get_nldi_basin(nwissite)
# plot(basin)

studysites = read.xlsx(xlsxFile = "~/../Box/CARB/SiteSelection/2022_OurProject/StudySitesAllProjects.xlsx", sheet = "All")
basin_list = list()
cents = list()
for (i in 1:nrow(studysites)) {
  nwissite <- list(featureSource = "nwissite", featureID = paste0("USGS-",studysites$STAID[i]))
  basin_list[[i]] <- get_nldi_basin(nwissite)
  #st_write(obj = basin_list[[i]], dsn = paste0("spatial_source/basin_outlines/", stringr::str_extract(studysites$STANAME[i], '\\w*') ),driver = "ESRI Shapefile"  )
  cents[[i]] = sf::st_centroid(basin_list[[i]])
  cents[[i]] = st_transform(cents[[i]], "EPSG:3310")
}

names(cents) = studysites$STANAME

centdf = data.table::rbindlist(cents)
centdf$name = studysites$STANAME
centdf$X = sapply(centdf$geometry, function(X) {X[1]} )
centdf$Y = sapply(centdf$geometry, function(X) {X[2]} )

write.csv(centdf, file = "../data/basin_centroids_epsg3310.csv")

# v <- vect(basin_list[[1]])
# crs(v) <- "EPSG:32611"
# #crs(v) <- "EPSG:26911"
# 
# r <- rast(v, nrow = 1000, ncol = 1000 )
# 
# #r <- rast(v, resolution = 30, units = "meter")
# plot(r)
# z <- rasterize(v, r)
# plot(z)
# cellSize(z)


