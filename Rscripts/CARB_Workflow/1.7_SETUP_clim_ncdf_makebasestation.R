# clim_ncdf_makebasestation
# 
# Create a RHESSys basestation file based on Gridmet netcdf inputs
# 
# Requires
#   - createbaseinfo_netcdf.c compiled binary
# 
# Inputs:
#   - Processed/edited/subset aggregated (through time), netcdf climate data, via: #http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
#     (easies to use/download spatial subset to begin with, via the NetcdfSubset download option)
#     Has already been processed and modified via clim_ncdf_processgridmet
#   - gridment elevation metadata netcdf (from same source)
#   - Basin and DEM raster maps

library(ncdf4)
library(terra)
library(rhutils)

plots = T

# ------------------------------ INPUTS ------------------------------
# uncompiled source of the createbaseinfo_netcdf.c code
# ncbase_src = "../../Util/createbaseinfo_netcdf.c"

# USE SAME BASIN RASTER AS PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
basin = rast("preprocessing/whitebox/basin.tif")
DEM = rast("preprocessing/whitebox/dem.tif")

# location for maps to be used with base station creation
# map_dest = "clim/netcdfmaps"

# where to output a new zone map based on Netcdf grid
# zone_dest = "preprocessing/whitebox/"
zone_dest = "preprocessing/whitebox/clim_grid.tif"

# BASESTATION DEST
dest = file.path("clim/netcdf.base")

# specify climate files individualy or via pattern, USE PROCESSED CLIMATE INPUTS FROM PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
clim_files = list.files(path = "clim",pattern = "crop_agg_met_", full.names = T)
# clim_files = list.files(path = "clim",pattern = "agg_met_daymet", full.names = T)
# clim_files = list.files(path = "clim",pattern = "agg_met_daymet.*test", full.names = T)
# backup since sometimes the clim netcdf can not work to create the grid
# gridmet_metadata = "../data/gridmet/metdata_elevationdata.nc"
# CHECK THE EDITS TO THE BASE STATION FILE AT THE END, FILE NAMES MAY NEED TO BE CORRECTED

# ------------------------------ END INPUTS ------------------------------
# for later use
basin_vect = as.polygons(basin)
basin_vect_latlon = project(basin_vect, "EPSG:4326")
# ------------------------------ GET GRID MAP FROM NETCDF ------------------------------
# for comparison, exploring original data
# tmp = nc_open(clim_files[1])
# latitude <- ncvar_get(tmp, "lat")
# longitude <- ncvar_get(tmp, "lon")

# there are slight differences between getting the grid from the metadata elevation map and the 
# netcdf already clipped using NCO. Using the clipped data from NCO as a default for now, but 
# could easily be better the other way, code should be easy to toggle

# Get CRS from netcdf
nctmp = nc_open(clim_files[1])
globalatts = ncatt_get(nctmp, varid = 0)
# lcc_atts = ncatt_get(nctmp, varid = "lambert_conformal_conic")
nc_close(nctmp)

# ---------- Grid from clipped climate netcdf files ----------
grid_from_clim = T
if (grid_from_clim) {
  # read in netcdf as raster, get only 1 day of precip
  nc_grid_lonlat = rast(clim_files[1])[[1]]
  # set values to sequential, chnage name
  values(nc_grid_lonlat) = seq_along(values(nc_grid_lonlat))
  names(nc_grid_lonlat) = "ID"
  # project to utm, mask to get zone/basestation map
  nc_grid = project(nc_grid_lonlat, basin, method="near")
  id_grid = mask(nc_grid, basin)
}

grid_from_elevmeta = F
if (grid_from_elevmeta) {
  gridmeta_proj = project(rast(gridmet_metadata), crs(basin))
  # ----- method using gridmet elevation metadata -----
  # clip/crop the grid
  nc_grid  = crop(gridmeta_proj, basin, extend = T, touches = T, snap = "out")
  # fix name, set values to sequential
  names(nc_grid) = "ID"
  values(nc_grid) = seq_along(values(nc_grid))
  nc_grid_lonlat = project(nc_grid, "EPSG:4326")
  # resample the grid to basin res (from 4km to 90m or whatever)
  basin_extended = extend(basin, nc_grid)
  nc_grid_resample = resample(nc_grid, basin_extended, method = "near")
  # mask then crop to original basin extent
  nc_grid_mask = mask(nc_grid_resample, basin_extended)
  id_grid = crop(nc_grid_mask, basin)
}

if (plots) {
  par(mfrow = c(2, 1))
  # layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
  # projected grid
  plot(nc_grid_lonlat)
  plot(basin_vect_latlon, add=T)
  # masked map (zone basestation map)
  plot(id_grid)
  plot(basin_vect, add=T)
  dev.copy2pdf(file = file.path("clim/", "nc_clim_grid.pdf"), width = 6, height = 10)
}

# --------------- output the zone grid for use in preprocessing ---------------
# zone = mask(x = id_grid_proj, mask = basin)
writeRaster(id_grid, zone_dest, overwrite=T)

# ------------------------------ GENERATE INPUTS FOR NETCDF BASE STATION CREATION  ------------------------------
# Maps need to be based on the netcdf latlong grid resolution + extent

# demproj = project(DEM, "EPSG:4326")
demproj = project(DEM, crs(nc_grid_lonlat)) # can just project based on read in data
dem_nc = resample(demproj, nc_grid_lonlat,method="average")

plot(demproj)
plot(dem_nc)

# ------------------------------ GET VAR NAMES  ------------------------------
varnames = rep(NA, length(clim_files))
for (i in seq_along(clim_files)) {
  tmp = nc_open(clim_files[i])
  vartmp = names(tmp$var)
  if (length(vartmp) > 1) {
    vartmp <- vartmp[!vartmp %in% c("lat", "lon", "lambert_conformal_conic")]
  }
  varnames[i] = vartmp
  nc_close(tmp)
}

# =========================== WRITE BASE STATION ==========================
# ID,lat,lon etc inputs are maps, matrix, single value, or values - where length(values) == length(values(map)) == length(c(t(matrix)))

write_basestation(
  dest = dest, year_start_index = 1900, search_dist = 0.04, 
  netcdf_var_x = "lon", netcdf_var_y = "lat", 
  precip_multiplier = 0.001, rhum_multiplier = 0.01, temp_unit = "K", 
  IDs = nc_grid_lonlat, dem = dem_nc, lat = crds(nc_grid_lonlat, df=T)$y, lon = crds(nc_grid_lonlat, df=T)$x, 
  x = NULL, y = NULL, lai = 3.0, screen_height = 10.0,
  netcdf_tmax_filename = clim_files[grepl("tmmx|tmax", clim_files)], 
  netcdf_var_tmax = varnames[grepl("tmmx|tmax", clim_files)], 
  netcdf_tmin_filename = clim_files[grepl("tmmn|tmin", clim_files)], 
  netcdf_var_tmin = varnames[grepl("tmmn|tmin", clim_files)], 
  netcdf_rain_filename = clim_files[grepl("pr|prcp", clim_files)],
  netcdf_var_rain = varnames[grepl("pr|prcp", clim_files)]
)


