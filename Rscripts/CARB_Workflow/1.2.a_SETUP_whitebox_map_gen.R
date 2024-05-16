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

wbt_version()
wbt_init()

# ==================== Inputs ====================
# source DEM
dem_source_path = "preprocessing/spatial_source/ned_dem_basinclip.tif"
# source gauge location shapefile, snap dist in map unit (m)
gauge_source = "preprocessing/spatial_source/gauge_loc.shp"
gauge_snap_dist = 90

stream_threshold=100

# for final and temp output
output_dir = "preprocessing/whitebox"
wb_tmp = file.path(output_dir, "wb_tmp")

plots = T
writeplots = T
testing = F

# check + add folders
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}if (!dir.exists(wb_tmp)) {
  dir.create(wb_tmp)
}


# ==================== starting DEM ====================
# copy and trim source DEM -- could clip manually here too
dem_source = rast(dem_source_path)
writeRaster(trim(dem_source), file.path(wb_tmp,"dem_trim.tif"), overwrite=T)

# ==================== CHANGE RESOLUTION HERE ====================
dem_resample = resample(dem_source, rast(ext(dem_source),resolution = 90))
# resample seems to work just as well
# dem_agg = aggregate(dem_source,fact = 3)

# if (plots) {
#   par(mfrow = c(3,2))
#   plot(dem_source, main = "Source DEM")
#   hist(dem_source)
#   plot(dem_resample, main = "Resample DEM")
#   hist(dem_resample)
#   plot(dem_agg, main = "Agg DEM")
#   hist(dem_agg)
#   # if (writeplots) {
#   #   dev.copy2pdf(file = file.path(output_dir, "plot0_DEMres_change.pdf"), width = 8, height = 6)
#   # }
#   summary(dem_source)
#   summary(dem_resample)
#   summary(dem_agg)
# }
writeRaster(trim(dem_resample), file.path(wb_tmp,"dem_trim.tif"), overwrite=T)

# ==================== Fill, d8 direction, accumulation ====================
## Breach depressions to ensure continuous flow
wbt_breach_depressions(dem = file.path(wb_tmp,"dem_trim.tif"), output = file.path(wb_tmp, "dem_brch.tif"))
# wbt_breach_depressions_least_cost(dem = file.path(wb_tmp,"dem_trim.tif"), output = file.path(wb_tmp, "dem_brch.tif"), dist = 50)
# wbt_fill_depressions(dem = file.path(wb_tmp,"dem_trim.tif"), output = file.path(wb_tmp, "dem_brch.tif"))

## Generate d8 flow pointer (note: other flow directions are available)
wbt_d8_pointer(dem = file.path(wb_tmp, "dem_brch.tif"), output = file.path(wb_tmp, "dem_brch_ptr_d8.tif"))
# this is fractional? d8, more equivilent to previous grass method
# wbt_fd8_pointer(dem = file.path(wb_tmp, "dem_brch.tif"), output = file.path(wb_tmp, "dem_brch_ptr_fd8.tif"))
## Generate d8 flow accumulation in units of cells (note: other flow directions are available)
wbt_d8_flow_accumulation(input = file.path(wb_tmp, "dem_brch.tif"), output = file.path(wb_tmp, "dem_brch_accum_d8.tif"), out_type = "cells")
# wbt_fd8_flow_accumulation(dem = file.path(wb_tmp, "dem_brch.tif"), output = file.path(wb_tmp, "dem_brch_accum_fd8.tif"), out_type = "cells")

# ==================== Streams ====================
## Generate streams with a stream initiation threshold of based on threshold param above
wbt_extract_streams(flow_accum = file.path(wb_tmp, "dem_brch_accum_d8.tif"), output = file.path(wb_tmp, "streams.tif"), threshold = stream_threshold)

# ==================== Basin, basin outlet ====================
wbt_jenson_snap_pour_points(pour_pts = gauge_source, 
                            streams = file.path(wb_tmp, "streams.tif"), 
                            output = file.path(output_dir, "gauge_loc_snap.shp"), 
                            snap_dist = gauge_snap_dist)

wbt_watershed(d8_pntr = file.path(wb_tmp, "dem_brch_ptr_d8.tif"), 
              pour_pts = file.path(output_dir, "gauge_loc_snap.shp"), 
              output = file.path(output_dir, "basin.tif"))

if (plots) {
  par(mfrow=c(1,1), mar=c(3,3,3,7))
  plot(rast(file.path(output_dir, "basin.tif")))
  par(mfrow=c(1,1), mar=c(3,3,3,7), new=TRUE)
  plot(rast(file.path(wb_tmp, "streams.tif")), add=T, col="black")
  par(mfrow=c(1,1), mar=c(3,3,3,7), new=TRUE)
  plot(vect(file.path(output_dir, "gauge_loc_snap.shp")), add=T,col="red")
}

# wbt_basins(d8_pntr = file.path(wb_tmp, "dem_brch_ptr_d8.tif"), output = file.path(output_dir, "basin.tif"))


# ==================== Crop and Trim by basin ====================
dem_brch_mask = mask(rast(file.path(wb_tmp, "dem_brch.tif")),rast(file.path(output_dir, "basin.tif")))
writeRaster(dem_brch_mask,filename = file.path(output_dir, "dem.tif"), overwrite = T)

streams_mask = mask(rast(file.path(wb_tmp, "streams.tif")),rast(file.path(output_dir, "basin.tif")))
writeRaster(streams_mask,filename = file.path(output_dir, "streams.tif"), overwrite = T)


# ==================== Subbasins/hillslopes ====================
# HILLSLOPES DOESNT WORK SINCE IT EXCLUDES THE STREAM PIXELS
wbt_subbasins(d8_pntr = file.path(wb_tmp, "dem_brch_ptr_d8.tif"), streams = file.path(wb_tmp, "streams.tif"), output = file.path(wb_tmp, "subbasins.tif"))

# wbt_hillslopes(d8_pntr = file.path(wb_tmp, "dem_brch_ptr_d8.tif"), streams = file.path(wb_tmp, "streams.tif"), output = file.path(wb_tmp, "hillslopes.tif"))
# hills = rast(file.path(wb_tmp, "hillslopes.tif"))
# hills = mask(hills,rast(file.path(output_dir, "basin.tif")))
# hillssize = summary(as.factor(values(hills)))
# hillssize[hillssize <10]
# plot(hills)

subbasins = rast(file.path(wb_tmp, "subbasins.tif"))
subbasins = mask(subbasins,rast(file.path(output_dir, "basin.tif")))
# values(subbasins) = as.factor(values(subbasins))

colors <- rainbow(length(seq(1,27)))
plot(subbasins, main = "Subbasins", col=colors)
# plot(rast( file.path(wb_tmp, "streams.tif")), main = "Streams",add=T,col="black")

# check if there are tiny subbasins
tmp = summary(as.factor(values(subbasins)))
tmp[tmp<20]
length(tmp[tmp<20])


reclasshill = F
if (reclasshill) {
  # reclass
  # 33 into 32
  # 34 into 36

  subset = subbasins
  subset[subset != 34] = NA
  tarsub = subbasins
  tarsub[tarsub != 36] = NA
  
  plot(subbasins, main = "Subbasins", col=colors)
  plot(subset, col = "black",add=T)
  plot(streams_mask, col="grey",add=T)
  plot(tarsub, col="brown",add=T)
  
  tmp = rast(file.path(wb_tmp, "subbasins.tif"))

  tmp[tmp == 33] = 32
  tmp[tmp == 34] = 36
  
  tmp2  = summary(as.factor(values(tmp)))
  tmp2[tmp2<20]
  
  writeRaster(x = tmp, filename = file.path(output_dir, "subbasins.tif"), overwrite = T)
}



# ==================== Slope and aspect ====================
wbt_slope(dem = file.path(output_dir, "dem.tif"), output = file.path(output_dir, "slope.tif"), units = 'degrees')
# aspect is standard (0 == 360 == NORTH )
wbt_aspect(dem = file.path(output_dir, "dem.tif"), output = file.path(output_dir, "aspect.tif"))

if (plots) {
  par(mfrow = c(2, 2))
  plot(rast(file.path(wb_tmp, "dem_brch.tif")), main = "Filled DEM")
  plot(rast(file.path(wb_tmp, "dem_brch_accum_fd8.tif")), main = "Flow Accumulation")
  plot(rast(file.path(output_dir, "slope.tif")), main = "Slope")
  plot(rast(file.path(output_dir, "aspect.tif")), main = "Aspect")
  if (writeplots) {
    dev.copy2pdf(file = file.path(output_dir, "plots1_dem_acc_slope_aspect.pdf"), width = 8, height = 6)
  }
}

# ==================== Horizons ====================
wbt_horizon_angle(dem = file.path(output_dir, "dem.tif"), output = file.path(wb_tmp, "horizon_east.tif"), azimuth = 270, max_dist = 100000)
wbt_horizon_angle(dem = file.path(output_dir, "dem.tif"), output = file.path(wb_tmp, "horizon_west.tif"), azimuth = 90, max_dist = 100000)
# sin of radian horizon
horizon_east_sin = sin(rast(file.path(wb_tmp, "horizon_east.tif")))
horizon_west_sin = sin(rast(file.path(wb_tmp, "horizon_west.tif")))
writeRaster(horizon_east_sin, file.path(output_dir, "e_horizon.tif"))
writeRaster(horizon_west_sin, file.path(output_dir, "w_horizon.tif"))



if (plots) {
  par(mfrow = c(2, 2))
  plot(rast(file.path(output_dir, "e_horizon.tif")), main = "East Horizon")
  plot(rast(file.path(output_dir, "w_horizon.tif")), main = "West Horizon")
  plot(rast( file.path(output_dir, "streams.tif")), main = "Streams")
  plot(rast(file.path(output_dir, "subbasins.tif")), main = "Subbasins")
  if (writeplots) {
    dev.copy2pdf(file = file.path(output_dir, "plots2_horizons_streamfs_subbasins.pdf"), width = 8, height = 6)
  }
}


# ==================== Soils ====================
# echo "import soils layers, originally from R-polaris workflow"
# r.import input=C:\Users\burke\Documents\CARB\Pitman\preprocessing\spatial_source\POLARISOut\mean\clay\0_5\lat3738_lon-120-119.tif output=clay_0_5
# r.import input=C:\Users\burke\Documents\CARB\Pitman\preprocessing\spatial_source\POLARISOut\mean\sand\0_5\lat3738_lon-120-119.tif output=sand_0_5
# r.soils.texture sand=sand_0_5@PERMANENT clay=clay_0_5@PERMANENT scheme=C:\Users\burke\Documents\CARB\data\USDA.dat output=soil_texture
# echo "import the soil texture"
# r.proj input=soil_texture location=soil mapset=rhessys output=soil_texture method=nearest --v --o
# r.out.gdal in=soil_texture output="%base%/soil_texture.tif" format=GTiff --o

# copy other maps
# patches, soils, rules
writeRaster(trim(rast("preprocessing/spatial90m/patches.tif")), file.path(output_dir,"patches.tif"), overwrite=T)
writeRaster(trim(rast("preprocessing/spatial90m/rules_LPC_90m.tif")), file.path(output_dir,"rules_LPC_90m.tif"), overwrite=T)
writeRaster(trim(rast("preprocessing/spatial90m/soils.tif")), file.path(output_dir,"soils.tif"), overwrite=T)
writeRaster(trim(rast("preprocessing/spatial90m/basin.tif")), file.path(output_dir,"basin.tif"), overwrite=T)



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
