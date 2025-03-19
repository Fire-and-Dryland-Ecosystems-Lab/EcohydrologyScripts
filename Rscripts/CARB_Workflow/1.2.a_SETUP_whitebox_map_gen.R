# 1.2.4_whitebox_TEST.R
# 
# https://whiteboxr.gishub.org/articles/demo.html

# install.packages("whitebox")
# whitebox::install_whitebox()

# for soil: ?
# appendTextureclass {soilassessment}

library(whitebox)
library(terra)
library(rhutils)
library(leaflet)
library(leafem)
library(sf)

wbt_version()
wbt_init()

# ==================== Inputs ====================
# source DEM
source_dem = "preprocessing/spatial_source/ned_dem_basinclip.tif"
# source gauge location shapefile, snap dist in map unit (m)
source_gauge = "preprocessing/spatial_source/gauge_loc.shp"
gauge_snap_dist = 90

res = 90

# stream_threshold=100
stream_threshold=150
# stream_threshold=170


# for final and temp output
output_dir = "preprocessing/whitebox"
tmp_dir = file.path(output_dir, "wb_tmp")

plots = T
writeplots = T
testing = F

# check + add folders
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
if (!dir.exists(tmp_dir)) {
  dir.create(tmp_dir)
}

# ==================== TRIM, RESAMPLE, Fill, d8 direction, accumulation, Streams, Basin, basin outlet  ====================
maps_out = wbox_dem2streams_gauge_basin(source_dem, src_gauge, res, stream_threshold = stream_threshold, gauge_snap_dist = gauge_snap_dist, output_dir, plots = T, writeplots = T, overwrite = T)

# ==================== Subbasins ====================
subbasin = wbox_subbasins(dem_brch = maps_out[["DEM"]], streams = maps_out[["streams"]], output_dir = output_dir, stream_threshold = stream_threshold, tmp_dir = tmp_dir, plots = T, writeplots = T, overwrite = T) 

sub_sizes = wbox_subbasins_vis(subbasin, maps_out[["streams"]])

# ==================== Aggregate outlier/small subbasins ====================
# IF THERE ARE SUBBASINS THAT ARE TOO SMALL OR STRANGE SHAPES AND NEED TO BE LUMPED INTO OTHERS
# USE LEAFLET MAP OF SUBBASINS OVERLAID WITH STREAMS, HOVER FOR SUBBASIN ID
# THEN RECLASS HILLSLOPES AS NEEDED, AND CHECK CHANGES AFTERWARDS
# NORMAL PLOT INCLUDES COUNT OF SUBBASINS, MIN AND MAX SIZES/CELL COUNTS
reclasshill = F
if (reclasshill) {
  s = rast(subbasin)
  rcl = as.matrix(data.frame(is = c(31,32, 25), becomes = c(4,4,11)))
  s_new <- classify(s, rcl)
  writeRaster(s_new, filename = file.path(tmp_dir, "subbasins_new.tif"),overwrite = T)

  subsizes = wbox_subbasins_vis(file.path(tmp_dir, "subbasins_new.tif"), maps_out[["streams"]])
  
  wbox_subbasins_plot(file.path(tmp_dir, "subbasins_new.tif"), streams, output_dir, stream_threshold)
  

  # ================= IF IT IS GOOD ======================
  file.rename(from = file.path(output_dir, "subbasins.tif"), file.path(output_dir, "subbasins_original.tif"))
  file.rename(from = file.path(tmp_dir, "subbasins_new.tif"), file.path(output_dir, "subbasins.tif"))
}

# ==================== Slope and aspect ====================
maps_out2 = wbox_slope_aspect_horizons(dem_brch = maps_out[["DEM"]], plots = T, writeplots = T, overwrite = T)

# ==================== Soils ====================
soil_texture = polaris2texture(
  basin = maps_out[["basin"]], 
  sand = "preprocessing/spatial_source/POLARISOut/mean/sand/0_5/lat3738_lon-120-119.tif",
  clay = "preprocessing/spatial_source/POLARISOut/mean/clay/0_5/lat3738_lon-120-119.tif",
  plot_out = file.path(output_dir, "SoilTextures_baseline.pdf")
)

reclasssoils = F
if (reclasssoils) {
  # soil_texture
  rcl = as.matrix(data.frame(is = c(11,8), becomes = c(12,9)))
  soil_texture_new <- classify(soil_texture, rcl)

  soil_texture_plot(soil_texture_new, file.path(output_dir, "SoilTextures_reclass.pdf"))
  
  # ================= IF IT IS GOOD ======================
  writeRaster(soil_texture_new, file.path(output_dir, "soils.tif"),overwrite=T)
}


# copy other maps
# patches, soils, rules
# writeRaster(trim(rast("preprocessing/spatial90m/patches.tif")), file.path(output_dir,"patches.tif"), overwrite=T)
# writeRaster(trim(rast("preprocessing/spatial90m/rules_90m.tif")), file.path(output_dir,"rules_LPC_90m.tif"), overwrite=T)



# ==================== TROUBLESHOOTING -- Check Maps ====================
if (testing) {
  template = "preprocessing/template/carb_msr.template"
  map_dir = "preprocessing/whitebox/"
  streams = "streams.tif"
  temp_read = template_read(template = template)
  maps_in = c(unique(temp_read[[5]][,2]), streams)
  # map_paths = file.path(map_dir, maps_in)
  
  # ==================== Check exists ====================
  file_paths = vector(mode = "character")
  for (i in maps_in) { #simple check
    file = list.files(path = map_dir, pattern = paste("^",i,"$",sep = ""),full.names = TRUE)
    file_paths = c(file_paths, file)
  }
  length(maps_in) == length(file_paths)
  
  # ==================== try read ====================
  read_stack = try(terra::rast(x = file_paths))
  names(read_stack) = maps_in
  plot(read_stack)
  
  # ==================== Check projection + source driver ====================
  # Check projections (read_stack will error if proj is different, but arguments might be different) -----
  p = vector(mode = "character",length = length(read_stack[1]))
  d = p
  for (i in 1:length(read_stack[1])) {
    p[i] = terra::crs(read_stack[[i]])
    text = sf::gdal_utils(util = "info", source = file_paths[i], quiet = T)
    pattern <- "Driver: (.+?)(?=\n|$)"
    d[i] <- sub("Driver: ","", unlist(regmatches(text, gregexpr(pattern, text, perl = TRUE))))
  }
  cat("Unique map projections (including args): ",length(unique(p)))
  
  # ==================== NaN and NA handling ====================
  # Handle NaNs - set to NA
  cat("Setting NaNs to NA.\n")
  terra::values(read_stack)[is.nan(terra::values(read_stack))] = NA
  
  cat("Trimming NAs.\n")
  read_stack = terra::trim(read_stack) #get rid of extra background

  read_stack = mask(read_stack,read_stack$basin.tif)
  
  read_stack$rules_LPC_90m.tif[!is.na(read_stack$basin.tif) & is.na(read_stack$rules_LPC_90m.tif) ]
  plot(read_stack$basin.tif)
  plot(!is.na(read_stack$basin.tif) & is.na(read_stack$rules_LPC_90m.tif), col = c("transparent","black"), add=T)
  plot(read_stack$streams.tif, add=T, col = "red")
  
  plot(trim(rast("preprocessing/spatial90m/rules_LPC_90m.tif")), col = "black")
  plot(trim(rast("preprocessing/spatial90m/basin.tif")), add=T)
  plot(trim(rast("preprocessing/spatial90m/dem.tif")), add=T)
  x = trim(rast(dem_source_path))
  plot(x)
  plot(x == min(values(x),na.rm = T), add = T, col =  c("transparent","black"))
  
  # Check for missing data (within world map mask) - no fix, just an error since I think this will break things if left unchecked
  cat("Checking for missing data within bounds of world map.\n")
  if (!is.null("map_info")) {
    wrld_vals = !is.na(terra::values(read_stack[["basin.tif"]]))
    NAs_in_wrld = lapply(as.data.frame(terra::values(read_stack)), function(X) {sum(is.na( X[wrld_vals]))})
    NAs_in_wrld[[which(grepl("streams",names(NAs_in_wrld))) ]] = NULL
    if (any(NAs_in_wrld > 0) ) {
      cat("One or more maps have NAs within the bounds of the world map, see maps and counts of NAs below:\n")
      print(NAs_in_wrld[NAs_in_wrld > 0])
      
    }
  }
  
  
  map_df = as.data.frame(read_stack)
}


# ==================== Testing ====================
testing = F
if (testing) {
  ### TEST ###
  demf = rast("preprocessing/spatial90m/dem_fill_grass.tif")
  par(mfrow = c(1,1))
  plot(trim(demf), main = "Filled DEM (GRASS)")
  
  whor_grass = rast("preprocessing/spatial90m/west_horizon_grasstest.tif")
  ehor_grass = rast("preprocessing/spatial90m/east_horizon_grasstest.tif")
  
  wbt_horizon_angle(dem = "preprocessing/spatial90m/dem_fill_grass.tif", output = file.path(wb_tmp, "horizon_west_grasstest.tif"), azimuth = 270, max_dist = 100000)
  wbt_horizon_angle(dem = "preprocessing/spatial90m/dem_fill_grass.tif", output = file.path(wb_tmp, "horizon_east_grasstest.tif"), azimuth = 90, max_dist = 100000)
  whor_test = rast(file.path(wb_tmp, "horizon_west_grasstest.tif"))
  ehor_test = rast(file.path(wb_tmp, "horizon_east_grasstest.tif"))
  
  whor_grass_deg = (whor_grass) * (180/pi)
  ehor_grass_deg = (ehor_grass) * (180/pi)
  
  par(mfrow = c(2, 2))
  plot(trim(whor_grass_deg), main = "GRASS West Horizon", breaks = seq(-20, 50, length.out = 11))
  plot(trim(ehor_grass_deg), main = "GRASS East Horizon", breaks = seq(-20, 50, length.out = 11))
  plot(trim(whor_test), main = "WhiteBox West Horizon", breaks = seq(-20, 50, length.out = 11))
  plot(trim(ehor_test), main = "WhiteBox East Horizon", breaks = seq(-20, 50, length.out = 11))
}
