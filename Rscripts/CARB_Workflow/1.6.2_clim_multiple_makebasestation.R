# 1.6.2 Separate basestation setup

library(ncdf4)
library(terra)

plots = T

# ------------------------------ INPUTS ------------------------------
input_nc_pcp = "clim/crop_agg_met_pr_1979_CurrentYear_CONUS.nc"
input_nc_tmin = "clim/crop_agg_met_tmmn_1979_CurrentYear_CONUS.nc"
input_nc_tmax = "clim/crop_agg_met_tmmx_1979_CurrentYear_CONUS.nc"

# USE SAME BASIN RASTER AS PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
basin = rast("preprocessing/spatial90m/basin.tif")
DEM = rast("preprocessing/spatial90m/dem.tif")

# location for maps to be used with base station creation
map_dest = "clim/netcdfmaps"

# where to output a new zone map based on Netcdf grid
zone_dest = "preprocessing/spatial90m/"

# specify climate files individualy or via pattern, USE PROCESSED CLIMATE INPUTS FROM PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
clim_files = list.files(path = "clim",pattern = "crop_agg_met_", full.names = T)

# CHECK THE EDITS TO THE BASE STATION FILE AT THE END, FILE NAMES MAY NEED TO BE CORRECTED

# ------------------------------ OUTPUT NETCDF AS TIF FOR USE IN CREATING GRID ------------------------------
# we know extent should be since we just cropped using it
basin_vect = as.polygons(basin)
basin_vect_unproj = project(basin_vect, "+proj=longlat +datum=WGS84 +no_defs ")

clim_tmp = rast(clim_files[1])
# tmp = nc_open(clim_files[1])
# tmp$var$precipitation_amount$dim[[1]]$vals
# tmp$var$precipitation_amount$dim[[2]]$vals

ext(clim_tmp) = ext(basin_vect_unproj)
# output geotif raster using first day of data from ncdf
# writeRaster(clim_tmp[[1]], "preprocessing/spatial_source/GRIDSOURCE_crop_agg_met_pr_1979_CurrentYear_CONUS.tif", overwrite=T)
nc_grid = clim_tmp[[1]]

# ------------------------------ GENERATE GRID MAP IN PROJECTED AND UNPROJECTED CRS  ------------------------------
#nc_grid = rast("preprocessing/spatial_source/GRIDSOURCE_crop_agg_met_pr_1979_CurrentYear_CONUS.tif")
id_grid = rast(nc_grid)
values(id_grid) = seq_along(values(nc_grid))
names(id_grid) = "ID"
id_grid_proj = project(id_grid, basin) # back to original projection

# ------------------------------ GET NETCDF DATA ------------------------------

# ------------------------------ netcdf ------------------------------ #
pr_edit_nc = nc_open(input_nc_pcp)
tmin_edit_nc = nc_open(input_nc_tmin)
tmax_edit_nc = nc_open(input_nc_tmax)

dat = ncvar_get(tmax_edit_nc, attributes(tmax_edit_nc$var)$names[1])
tmax_df = as.data.frame((t(dat)))
dat = ncvar_get(tmin_edit_nc, attributes(tmin_edit_nc$var)$names[1])
tmin_df = as.data.frame((t(dat)))
dat = ncvar_get(pr_edit_nc, attributes(pr_edit_nc$var)$names[1])
rain_df = as.data.frame((t(dat)))

# CORRECT CLIM DATA
if (max(df$rain, na.rm = T) < 1) { #assume in meters
  df$rain = df$rain * 1000
}
if (mean(df$tmin, na.rm = T) > 200) {
  df$tmin =  df$tmin - 273.15
}
if (mean(df$tmax, na.rm = T) > 200) {
  df$tmax =  df$tmax - 273.15
}


# ------------------------------ WRITE FILES ------------------------------
# something like:
zoneIDS = unique(zone_grid)

for (i in zoneIDS) {
  
  outname = "clim/ward_netcdfgridmet_agg"
  startdate = "1979 01 01 01"
  
  # WRITE CLIMATE FILES
  write(paste0(startdate,"\n", paste0(format(nc_clim$tmin,nsmall = 1),collapse = "\n")) , file = paste0(outname,".tmin" ))
  write(paste0(startdate,"\n", paste0(format(nc_clim$tmax,nsmall = 1),collapse = "\n")) , file = paste0(outname,".tmax" ))
  write(paste0(startdate,"\n", paste0(nc_clim$rain/1000,collapse = "\n")) , file = paste0(outname,".rain" ))
  
  # WRITE BASESTATION
  base_source = "clim/ward.base"
  base_in = readLines(base_source)
  base_out = gsub(gsub(".base","", base_source), outname, base_in)
  writeLines(base_out,con = paste0(outname,".base"))
}

