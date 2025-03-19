# clim_ndf_processgridmet
# 
# This workflow takes aggregated (through time) GridMet netcdf climate data as input, and 
# subsets and processes the netcdf data to be usable in RHESSys
# 
# Requires
#   - NCO (installed on linux with apt with 'sudo apt install NCO')
# 
# Inputs:
#   - Aggregated (through time), netcdf climate data, via: #http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
#     (easies to use/download spatial subset to begin with, via the NetcdfSubset download option)
#   - Basin raster

library(ncdf4)
library(terra)

plot = T
basin = rast("preprocessing/whitebox/basin.tif")
# ext(basin)

# clim_files = list.files(path = "../data/gridmet/",pattern = "agg_met_", full.names = T)
# clim_files = list.files(path = "clim/",pattern = "agg_met_gridmet", full.names = T)
clim_files = list.files(path = "clim/",pattern = "agg_met_daymet", full.names = T)

# destination folder for cropped + edited ncdf clim files
clim_dest = "clim"

# ------------------------------ CROP + PROCESS NETCDF CLIMATE DATA USING BASIN EXTENT AND NCO ------------------------------
# FOR NETCDF TO WORK WITH RHESSYS:
#   Must use 'time' as temporal dimension, below we change every instance of 'day' to 'time' to be safe
#   Must have variable as floating point/double data type, NOT INTEGER

basin_vect = as.polygons(basin)
# writeVector(basin_vect, "preprocessing/spatial_source/basin_vect.shp", overwrite=T)
# projection for netcdf data is: +proj=longlat +datum=WGS84 +no_defs 
# basin_vect_latlon = project(basin_vect, "+proj=longlat +datum=WGS84 +no_defs ")
basin_vect_latlon = project(basin_vect, "EPSG:4326")
# writeVector(basin_vect_latlon, "preprocessing/spatial_source/basin_vect_latlon.shp", overwrite=T)

basin_extent = ext(basin_vect_latlon)
# cell size is 4 km, buffer should be at least 3 km at most latitudes 
# so should cover the centroids of netcdf cells not completely inside basin
basin_extent_buffer = extend(basin_extent, 0.03)

# this should be in lon lat wgs84 to match ncdf - see xy axes
if (plot) {
  plot(basin_extent_buffer)
  plot(basin_extent,add=T)
  plot(basin_vect_latlon,add=T)
}

basin_extent = basin_extent_buffer

# iterate through each climate file
for (infile in clim_files) {
  # destination for cropped file
  outfile = paste0(file.path(clim_dest, "crop_"), basename(infile))
  # get the variable name for when modifying data type
  nctmp = nc_open(infile)
  varname = names(nctmp$var)
  globalatts = ncatt_get(nctmp, varid = 0)
  nc_close(nctmp)
  


  # ---------- DAYMET ----------
  if (!is.null(globalatts$source) && grepl("daymet|Daymet", globalatts$source)) {
    # for daymet
    if (length(varname) > 1) {
      varname <- varname[!varname %in% c("lat", "lon", "lambert_conformal_conic")]
    }

    # ncatt_get(nc = nctmp, varid ="prcp")
    # projvar <- nctmp$var$lambert_conformal_conic
    # proj <- ncatt_get(nc = nctmp, varid = "lambert_conformal_conic")

    crop_cmd <- paste0(
      "ncks -O -d x,", basin_extent[1], ",", basin_extent[2],
      " -d y,", basin_extent[3], ",", basin_extent[4], " ", infile, " ", outfile, "\n",
      # "ncatted -O -a coordinates,,m,c,'time lat lon' ", outfile, "\n",
      # "ncap2 -O -s '", varname, "=double(", varname, ")' ", outfile, " ", outfile
    )
    # ---------- GRIDMET ----------
  } else {
    crop_cmd <- paste0(
      "ncks -O -d lon,", basin_extent[1], ",", basin_extent[2],
      " -d lat,", basin_extent[3], ",", basin_extent[4], " ", infile, " ", outfile, "\n",
      "ncrename -d day,time -v day,time -O ", outfile, " ", outfile, "\n",
      "ncatted -O -a coordinates,,m,c,'time lat lon' ", outfile, "\n",
      "ncap2 -O -s '", varname, "=double(", varname, ")' ", outfile, " ", outfile
    )
  }


  test = F
  if (test) {
    gridnc = nc_open("clim/crop_agg_met_gridmet_pr_1979-01-01_2024-09-30.nc")
    
    ncatt_get(nc = gridnc, varid ="precipitation_amount")
    
    nc_close(gridnc)

    # ----

    # nctmp = nc_open("clim/agg_met_daymet_tmax_1980-01-01_2023-12-31.nc")
    # nctmp = nc_open("clim/agg_met_daymet_tmin_1980-01-01_2023-12-31.nc")
    # nctmp = nc_open("clim/agg_met_daymet_prcp_1980-01-01_2023-12-31.nc")
    filesin = c("clim/agg_met_daymet_tmax_1980-01-01_2023-12-31.nc", "clim/agg_met_daymet_tmin_1980-01-01_2023-12-31.nc", "clim/agg_met_daymet_prcp_1980-01-01_2023-12-31.nc")

    for (i in seq_along(filesin)) {
      nctmp = nc_open(filesin[i])

      lat = ncdim_def(name = "lat", units = "km", vals = ncvar_get(nctmp,varid = "y"), longname = "y coordinate")
      lon = ncdim_def(name = "lon", units = "km", vals = ncvar_get(nctmp,varid = "x"), longname = "x coordinate")
      time = ncdim_def(name = "time", units = "days since 1950-01-01 00:00:00", vals = ncvar_get(nctmp,varid = "time"), unlim = T, calendar = "standard", longname = "24-hour day based on local time")
      
      var = nctmp$var[!names(nctmp$var) %in% c("lambert_conformal_conic","lat","lon")]
      vardef = ncvar_def(name = var[[1]]$name, longname = var[[1]]$longname, units = var[[1]]$units, dim = list(lat,lon,time), prec = "float" )
      newnc = nc_create(filename = paste0("clim/agg_met_daymet_",var[[1]]$name,"_test.nc"),vars = vardef, force_v4=T)
      ncvar_put(nc = newnc, varid = var[[1]]$name, vals = ncvar_get(nctmp,varid = var[[1]]$name))
  
      nc_close(nctmp)
      nc_close(newnc)
    }


    test = nc_open("clim/agg_met_daymet_prcp_test.nc")
    ncvar_get(test,"prcp")
    nc_close(test)

  }

  system()
  
  # ASSUMES EITHER UNIX OR IF ON WINDOWS USING WSL
  if (.Platform$OS.type == "windows") {
    # tmp = noquote(paste("bash -c \"", crop_cmd, "\"", sep = ""))
    tmp = noquote(paste("wsl ", crop_cmd, sep = ""))
  } else {
    tmp = crop_cmd
  }
  
  system(tmp)
}


  # TO CHECK THAT LAT AND LON IS -180 TO 180 (IT SEEMS TO BE - this is different from MACA data)
  # dims = unlist(strsplit(ncatt_get(nctmp, "precipitation_amount", "dimensions")$value, " "))
  # lonvals = nctmp$var$precipitation_amount$dim[[which(dims=="lon")]]$vals
  # debug = F
  # if (debug){
  #   tmp = rast(infile)[[1]]
  #   plot(tmp)
  #   tmp2 = rast("../../CARB/data/pr_2023.nc")[[1]]
  #   plot(tmp2)
  # }
