# Get Climate Data
#
# This script shows how to download gridmet or daymet using wget scripts and your input basin to spatially subset the climate data
# Can get either only precip, tmin, tmax, or the full set of clim inputs

library(terra)
library(ncdf4)
library(utils)

# --------------------------------------------------------------------
# ------------------------------ INPUTS ------------------------------
# --------------------------------------------------------------------
plot = T
basin = rast("preprocessing/whitebox/basin.tif")

# get basin extent in lat lon, with buffer
basin_proj = project(basin, "EPSG:4326")
basin_extent = ext(basin_proj)
# this buffer should work for gridmet which has ~4km grid, will be big for daymet but should not be an issue
basin_extent_buffer = extend(basin_extent, 0.03)
basin_extent_buffer = round(basin_extent_buffer, digits = 4)

if (plot) {
  plot(basin_extent_buffer)
  plot(basin_extent,add=T)
  plot(basin_proj,add=T)
}

north = as.numeric(ymax(basin_extent_buffer))
west = as.numeric(xmin(basin_extent_buffer))
east = as.numeric(xmax(basin_extent_buffer))
south = as.numeric(ymin(basin_extent_buffer))

# --------------------------------------------------------------------
# ------------------------------ DAYMET ------------------------------
# --------------------------------------------------------------------

# -- DAYMET DOWNLOAD HAS A LARGE DELAY WHEN GETTING MANY YEARS, MAY NEED TO USE A LARGE VALUE FOR TIMEOUT OR TRY MANUALLY

# If you just want to do it all online
# source: https://thredds.daac.ornl.gov/thredds/ncss/grid/daymet-v4-agg/na.ncml/dataset.html

# min datetime 1980-01-01T12:00:00Z max 2023-12-31T12:00:00Z
# vars = dayl prcp srad (daylight average incident shortwave radiation) swe tmax tmin vp
# example
# https://thredds.daac.ornl.gov/thredds/ncss/daymet-v4-agg/na.ncml?var=lat&var=lon&var=tmax&north=46.10139&west=-118.18491&east=-117.87436&south=45.89437&disableProjSubset=on&horizStride=1&time_start=1980-01-01T12%3A00%3A00Z&time_end=1980-12-31T12%3A00%3A00Z&timeStride=1&accept=netcdf

# this is the aggregated dataset for the daymet data
dataset = "/thredds/ncss/daymet-v4-agg/na.ncml"
# rhessys expects climate inputs in separate files, so will generate urls/files for each var separately
vars = c("tmax", "tmin","prcp")

startdate = "1980-01-01"
enddate = "2023-12-31"
# enddate = "1980-12-31"

urltext = paste0("https://thredds.daac.ornl.gov",dataset,"?var=lat&var=lon&var=",vars,
"&north=",north,"&west=",west,"&east=",east,"&south=",south,"&disableProjSubset=on&horizStride=1&time_start=",
startdate,"T12%3A00%3A00Z&time_end=",enddate,"T12%3A00%3A00Z&timeStride=1&accept=netcdf")

destnames = paste0("clim/agg_met_daymet_",vars,"_",startdate,"_",enddate, ".nc")

# --------------- DOWNLOAD -----------------
# download is slow and has big delay at start for many years
# On windows, have had both success and failure with both libcurl and curl. - most recently curl worked
# downmethod = "libcurl"
downmethod = "curl"
options(timeout = max(5000, getOption("timeout")))

checkdown1 = download.file(url = urltext[1],destfile = destnames[1], method = downmethod, cacheOK = F)
checkdown2 = download.file(url = urltext[2],destfile = destnames[2], method = downmethod, cacheOK = F)
checkdown3 = download.file(url = urltext[3],destfile = destnames[3], method = downmethod, cacheOK = F)

# --------------------------------------------------------------------
# ------------------------------ GridMET -----------------------------
# --------------------------------------------------------------------

# source: http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
# then select netcdf subset

# all vars: c("pr", "rmax", "rmin", "sph", "srad", "th", "tmmn", "tmmx", "vpd","vs")
# all varnames c("precipitation_amount", "daily_maximum_relative_humidity", "daily_minimum_relative_humidity", "daily_mean_specific_humidity", "daily_mean_shortwave_radiation_at_surface", "daily_mean_wind_direction", "daily_minimum_temperature", "daily_maximum_temperature", "daily_mean_vapor_pressure_deficit", "daily_mean_wind_speed")

vars = c("pr", "rmax", "rmin", "sph", "srad", "th", "tmmn", "tmmx", "vpd","vs")
varnames = c("precipitation_amount", "daily_maximum_relative_humidity", "daily_minimum_relative_humidity", "daily_mean_specific_humidity", "daily_mean_shortwave_radiation_at_surface", "daily_mean_wind_direction", "daily_minimum_temperature", "daily_maximum_temperature", "daily_mean_vapor_pressure_deficit", "daily_mean_wind_speed")

datasets = paste0("/thredds/ncss/agg_met_",vars,"_1979_CurrentYear_CONUS.nc")

# dates- starts 1979-1-1, ends current date 
startdate = "1979-01-01"
enddate = "2024-09-30"

urltext = paste0("http://thredds.northwestknowledge.net:8080",datasets,"?var=",varnames,
"&north=",north,"&west=",west,"&east=",east,"&south=",south,"&disableProjSubset=on&horizStride=1&time_start=",
startdate,"T12%3A00%3A00Z&time_end=",enddate,"T00%3A00%3A00Z&timeStride=1&accept=netcdf")

destnames = paste0("clim/agg_met_gridmet_",vars,"_",startdate,"_",enddate, ".nc")

# --------------- DOWNLOAD -----------------
# download is slow and has big delay at start for many years
# On windows, have had both success and failure with both libcurl and curl. - most recently curl worked
# downmethod = "libcurl"
downmethod = "curl"
options(timeout = max(1800, getOption("timeout")))

checkdownlist = list()
for (i in seq_along(urltext)) {
  checkdownlist[i] = download.file(url = urltext[i],destfile = destnames[i], method = downmethod, cacheOK = F)
}

