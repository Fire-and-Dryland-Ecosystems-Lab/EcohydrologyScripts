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
basin = rast("preprocessing/spatial90m/basin.tif")

# source netcdfs (these are default file names from THREADS aggregate download)
# pr = "../data/gridmet/agg_met_pr_1979_CurrentYear_CONUS.nc"
# tmin = "../data/gridmet/agg_met_tmmn_1979_CurrentYear_CONUS.nc"
# tmax = "../data/gridmet/agg_met_tmmx_1979_CurrentYear_CONUS.nc"
# rmax = "../data/gridmet/agg_met_rmax_1979_CurrentYear_CONUS.nc"
# rmin = "../data/gridmet/agg_met_rmin_1979_CurrentYear_CONUS.nc"
# vpd = "../data/gridmet/agg_met_vpd_1979_CurrentYear_CONUS.nc"
# sph = "../data/gridmet/agg_met_sph_1979_CurrentYear_CONUS.nc"
# srad = "../data/gridmet/agg_met_srad_1979_CurrentYear_CONUS.nc"
# wdir = "../data/gridmet/agg_met_th_1979_CurrentYear_CONUS.nc"
# wspd = "../data/gridmet/agg_met_vs_1979_CurrentYear_CONUS.nc"
# clim_files = c(pr, tmin, tmax, rmax, rmin, vpd, sph, srad, wdir, wspd)

# OR
clim_files = list.files(path = "../data/gridmet/",pattern = "agg_met_", full.names = T)

# destination folder for cropped + edited ncdf clim files
clim_dest = "clim"

# ------------------------------ CROP + PROCESS NETCDF CLIMATE DATA USING BASIN EXTENT AND NCO ------------------------------
# FOR NETCDF TO WORK WITH RHESSYS:
#   Must use 'time' as temporal dimension, below we change every instance of 'day' to 'time' to be safe
#   Must have variable as floating point/double data type, NOT INTEGER

basin_vect = as.polygons(basin)
# writeVector(basin_vect, "preprocessing/spatial_source/basin_vect.shp", overwrite=T)
# projection for netcdf data is: +proj=longlat +datum=WGS84 +no_defs 
basin_vect_unproj = project(basin_vect, "+proj=longlat +datum=WGS84 +no_defs ")
# writeVector(basin_vect_unproj, "preprocessing/spatial_source/basin_vect_unproj.shp", overwrite=T)

# this should be in lon lat wgs84 to match ncdf - see xy axes
if (plot) {plot(basin_vect_unproj)}

basin_extent = ext(basin_vect_unproj)

# iterate through each climate file
for (infile in clim_files) {
  # get the variable name for when modifying data type
  nctmp = nc_open(infile)
  
  # TO CHECK THAT LAT AND LON IS -180 TO 180 (IT SEEMS TO BE - this is different from MACA data)
  # dims = unlist(strsplit(ncatt_get(nctmp, "precipitation_amount", "dimensions")$value, " "))
  # lonvals = nctmp$var$precipitation_amount$dim[[which(dims=="lon")]]$vals
  
  varname = names(nctmp$var)
  nc_close(nctmp)
  outfile = paste0(file.path(clim_dest, "crop_"), basename(infile))
  crop_cmd = paste0("ncks -O -d lon,",basin_extent[1],",",basin_extent[2]," -d lat,",basin_extent[3],",",basin_extent[4]," ",infile," ",outfile,"\n",
                    "ncrename -d day,time -v day,time -O ", outfile, " ", outfile, "\n",
                    "ncatted -O -a coordinates,,m,c,'time lat lon' ", outfile, "\n",
                    "ncap2 -O -s '",varname,"=double(",varname,")' ", outfile," ", outfile)
  tmp = noquote(paste("bash -c \"", crop_cmd, "\"", sep = ""))
  system(tmp)
}


