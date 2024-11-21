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

plots = T

# ------------------------------ INPUTS ------------------------------
# uncompiled source of the createbaseinfo_netcdf.c code
ncbase_src = "../exmaples/createbaseinfo_netcdf.c"

# USE SAME BASIN RASTER AS PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
basin = rast("preprocessing/whitebox/basin.tif")
DEM = rast("preprocessing/whitebox/dem.tif")

# location for maps to be used with base station creation
map_dest = "clim/netcdfmaps"

# where to output a new zone map based on Netcdf grid
zone_dest = "preprocessing/whitebox/"

# specify climate files individualy or via pattern, USE PROCESSED CLIMATE INPUTS FROM PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
clim_files = list.files(path = "clim",pattern = "crop_agg_met_", full.names = T)
# backup since sometimes the clim netcdf can not work to create the grid
gridmet_metadata = "../data/gridmet/metdata_elevationdata.nc"
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
  par(mfrow = c(1, 2))
  # projected grid
  plot(nc_grid_lonlat)
  plot(basin_vect_latlon, add=T)
  # masked map (zone basestation map)
  plot(id_grid)
  plot(basin_vect, add=T)
}


# --------------- output the zone grid for use in preprocessing ---------------
# zone = mask(x = id_grid_proj, mask = basin)
writeRaster(id_grid, file.path(zone_dest, "zone.tif"), overwrite=T)

# ------------------------------ GENERATE INPUTS FOR NETCDF BASE STATION CREATION  ------------------------------
# Maps need to be based on the netcdf latlong grid resolution + extent
demproj = project(DEM, "EPSG:4326")
dem_nc = resample(demproj, nc_grid_lonlat,method="average")
# LAI
lai_nc = dem_nc
values(lai_nc) = 3.0
names(lai_nc) = "lai"

if (!file.exists(map_dest)) {
  dir.create(map_dest)
}

writeRaster(dem_nc, filename = file.path(map_dest,"dem.asc"), filetype="AAIGrid", gdal = c("FORCE_CELLSIZE=TRUE"), NAflag=-9999, overwrite=T)
writeRaster(lai_nc, filename = file.path(map_dest,"lai.asc"), filetype="AAIGrid", gdal = c("FORCE_CELLSIZE=TRUE"), NAflag=-9999, overwrite=T)
writeRaster(nc_grid_lonlat, filename = file.path(map_dest,"cellid.asc"), filetype="AAIGrid", datatype =  "INT4S", gdal = c("FORCE_CELLSIZE=TRUE"), NAflag=-9999, overwrite=T)
file.remove(list.files(path = map_dest, pattern = ".prj|.aux.xml", full.names = T))

# format: ID ID Y X X Y

xyid_loc = cbind(1:nrow(crds(nc_grid_lonlat)), 1:nrow(crds(nc_grid_lonlat)), crds(nc_grid_lonlat)[,c("y","x")], crds(nc_grid_lonlat))
write.table(xyid_loc, file.path(map_dest,"xyid_loc.txt"), row.names = F, col.names = F)

# ------------------------------ RUN C BIN TO GET BASE STATION  ------------------------------
# COMPILE
system(noquote(paste("bash -c \"", paste0("gcc ", ncbase_src, " -o ", file.path(map_dest,"create_netcdfbase")), "\"", sep = "")))

# RUN -- YOU MAY NEED TO EDIT THIS, command should look like: ./create_netcdfbase cellid.asc lai.asc dem.asc xyid_loc.txt clim/ netcdf.base
cmd = paste0(file.path(map_dest,"create_netcdfbase"), " ", file.path(map_dest,"cellid.asc"), " ", file.path(map_dest,"lai.asc"), " ",
             file.path(map_dest,"dem.asc"), " ", file.path(map_dest,"xyid_loc.txt"), " ", dirname(clim_files[1]), "/ ", file.path(map_dest,"netcdf.base"))
system(noquote(paste("bash -c \"",cmd ,"\"", sep = "")))

# ------------------------------ EDIT/FIX BASE STATION  ------------------------------

varnames = rep(NA, length(clim_files))
for (i in seq_along(clim_files)) {
  tmp = nc_open(clim_files[i])
  varnames[i] = names(tmp$var)
    
}

ncbase = read.table(file.path(map_dest,"netcdf.base"))

ncbase[ncbase[,2]=="year_start_index", 1] = "1900"

ncbase[ncbase[,2]=="netcdf_tmax_filename", 1] = clim_files[grepl("tmmx", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_tmax", 1] = varnames[grepl("tmmx", clim_files)]

ncbase[ncbase[,2]=="netcdf_tmin_filename", 1] = clim_files[grepl("tmmn", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_tmin", 1] = varnames[grepl("tmmn", clim_files)]

ncbase[ncbase[,2]=="netcdf_rain_filename", 1] = clim_files[grepl("pr", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_rain", 1] = varnames[grepl("pr", clim_files)]

ncbase[ncbase[,2]=="netcdf_huss_filename", 1] = clim_files[grepl("daily_mean_specific_humidity", varnames)]
ncbase[ncbase[,2]=="netcdf_var_huss", 1] = varnames[grepl("daily_mean_specific_humidity", varnames)]

ncbase[ncbase[,2]=="netcdf_rmax_filename", 1] = clim_files[grepl("rmax", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_rmax", 1] = varnames[grepl("rmax", clim_files)]

ncbase[ncbase[,2]=="netcdf_rmin_filename", 1] = clim_files[grepl("rmin", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_rmin", 1] = varnames[grepl("rmin", clim_files)]

ncbase[ncbase[,2]=="netcdf_rsds_filename", 1] = clim_files[grepl("daily_mean_shortwave_radiation_at_surface", varnames)]
ncbase[ncbase[,2]=="netcdf_var_rsds", 1] = varnames[grepl("daily_mean_shortwave_radiation_at_surface", varnames)]

ncbase[ncbase[,2]=="netcdf_was_filename", 1] = clim_files[grepl("daily_mean_wind_speed", varnames)]
ncbase[ncbase[,2]=="netcdf_var_was", 1] = varnames[grepl("daily_mean_wind_speed", varnames)]

write.table(ncbase, file.path(map_dest,"netcdf.base"), row.names = F, col.names = F, quote = F)


# ------------------------------ OLD GRASS VERSION ------------------------------
grass = F
if (grass) {
  # GRASS COMMANDS
  # 
  # Create location with tif version of cropped netcdf, with extent fixed by vector basin file
  # r.proj location=Wardproj mapset=PERMANENT input=basin
  # g.region n=39.162039 s=39.111069 w=-120.247294 e=-120.156854 rows=1 cols=2 -p
  # g.region res=00:00:02.7 -a -p
  # r.proj --overwrite --verbose location=Wardproj mapset=PERMANENT input=basin
  # g.region n=39.162750 s=39.111000 w=-120.247500 e=-120.156750 rows=1 cols=2 -p
  
  # r.mapcalc expression=xloc=x()                                                   
  # r.mapcalc expression=yloc=y()                                                   
  # r.mapcalc expression=xmap=col()                                                 
  # r.mapcalc expression=ymap=row() 
  
  # r.mapcalc expression=xyid=(ymap-1)*12+xmap
  # r.mapcalc expression=screen_height=2.0
  # r.proj -g location=Wardproj mapset=PERMANENT input=dem30f
  # g.region n=39:10:13.361221N s=39:03:07.027948N w=120:15:27.127436W e=120:08:08.47742W rows=427 cols=337
  # g.region res=00:00:02.7 -a -p
  # r.proj --overwrite --verbose location=Wardproj mapset=PERMANENT input=dem30f
  # g.region rast=xyid res=00:00:02.7
  # r.mapcalc expression=dem2_7sec=dem30f
  # r.mapcalc expression=lai2_7sec=3.0
  # r.mapcalc expression=basin2_7sec=basin
  # g.region rast=xyid res=00:02:30
  # r.resamp.stats input=basin output=test method=sum --overwrite
  # r.mapcalc expression=basinmask=if(test >= 1, 1)
  # g.remove type=raster name=test@PERMANENT -f
  # r.mask raster=basinmask maskcats=1
  # r.mapcalc expression=xyid_msk = xyid
  # r.resamp.stats input=dem2_7sec output=dem2m30s method=average --overwrite
  # r.resamp.stats input=lai2_7sec output=lai2m30s method=average --overwrite
  # r.resamp.stats input=screen_height output=screenheight2m30s method=average --overwrite
  
  # set base="C:/Users/burke/Documents/Carb/Ward/clim/"
  # 
  # r.out.ascii input=xyid_msk output="%base%/cellid_msk.asc"
  # r.out.ascii input=dem2m30s output="%base%/dem_msk.asc"
  # r.out.ascii input=lai2m30s output="%base%/lai_msk.asc"
  # r.out.ascii input=screenheight2m30s output="%base%/screen_msk.asc"
  # r.out.ascii input=xloc output="%base%/xloc.asc"
  # r.out.ascii input=yloc output="%base%/yloc.asc"
  
  # FIX HEADERS
  # r.out.gdal input=xyid_msk output="%base%/cellid.asc" format=AAIGrid type=Int16 createopt=FORCE_CELLSIZE=TRUE nodata=-9999 --overwrite
  # r.out.gdal input=dem2m30s output="%base%/dem.asc" format=AAIGrid type=Float64 createopt=FORCE_CELLSIZE=TRUE nodata=-9999 --overwrite
  # r.out.gdal input=lai2m30s output="%base%/lai.asc" format=AAIGrid type=Float64 createopt=FORCE_CELLSIZE=TRUE nodata=-9999 --overwrite
  # r.out.gdal input=screenheight2m30s output="%base%/screen.asc" format=AAIGrid type=Float64 createopt=FORCE_CELLSIZE=TRUE nodata=-9999 --overwrite
  # r.out.gdal input=xloc output="%base%/xloc.asc" format=AAIGrid type=Float64 createopt=FORCE_CELLSIZE=TRUE nodata=-9999 --overwrite
  # r.out.gdal input=yloc output="%base%/yloc.asc" format=AAIGrid type=Float64 createopt=FORCE_CELLSIZE=TRUE nodata=-9999 --overwrite
  
  # file.remove(list.files(path = "clim/", pattern = ".prj|.aux.xml", full.names = T))
  # 
  # xloc = read.table("clim/xloc.asc", skip = 6)
  # yloc = read.table("clim/yloc.asc", skip = 6)
  # # COMBINE - by rows, eg ids are seq along each row
  # # format: ID ID Y X X Y
  # xyid_loc = data.frame(seq_along(as.vector(unlist(xloc))), seq_along(as.vector(unlist(xloc))),as.vector(unlist(yloc)), 
  #                       as.vector(unlist(xloc)), as.vector(unlist(xloc)), as.vector(unlist(yloc)), fix.empty.names = F)
  # write.table(xyid_loc, "clim/xyid_loc.txt", row.names = F, col.names = F)
  # 
  # # ./create_netcdfbase cellid.asc lai.asc dem.asc xyid_loc.txt clim/ netcdf.base
  # 
  # # CHECK VARS
  # pr_nc = nc_open("clim/crop_agg_met_pr_1979_CurrentYear_CONUS.nc")
  # tmin_nc = nc_open("clim/crop_agg_met_tmmn_1979_CurrentYear_CONUS.nc")
  # tmax_nc = nc_open("clim/crop_agg_met_tmmx_1979_CurrentYear_CONUS.nc")
}




