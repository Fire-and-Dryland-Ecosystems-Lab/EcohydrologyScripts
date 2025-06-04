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
source_dem = "preprocessing/spatial_source/DEM_Sooke.tif"
# source gauge location shapefile, snap dist in map unit (m)
source_gauge = "preprocessing/rithet_outlet.shp"
gauge_snap_dist = 90

res = 90

# 90 res
stream_threshold=75
# stream_threshold=100
# stream_threshold=150
# stream_threshold=175
# 30 res
# stream_threshold=250

# for final and temp output
output_dir = "preprocessing/whitebox90"
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
maps_out = wbox_dem2streams_gauge_basin(source_dem, source_gauge, res, stream_threshold = stream_threshold, gauge_snap_dist = gauge_snap_dist, output_dir, plots = T, writeplots = T, overwrite = T)

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
  rcl = as.matrix(data.frame(is = c(15), becomes = c(19)))
  s_new <- classify(s, rcl)
  writeRaster(s_new, filename = file.path(tmp_dir, "subbasins_new.tif"),overwrite = T)

  subsizes = wbox_subbasins_vis(file.path(tmp_dir, "subbasins_new.tif"), maps_out[["streams"]])
  
  wbox_subbasins_plot(file.path(tmp_dir, "subbasins_new.tif"),  maps_out[["streams"]], output_dir, stream_threshold)
  

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

# ==================== PATCH MAP ====================

p_map = rast(file.path(output_dir, "dem.tif"))
names(p_map) = "patches"
# resample to get unique patches 1 per cell
values(p_map)[!is.nan(values(p_map))] = seq_along(values(p_map)[!is.nan(values(p_map))])
# writeRaster(p_map, "preprocessing/spatial90m/patches.tif", overwrite=T)
writeRaster(p_map, file.path(output_dir, "patches.tif"), overwrite=T)

