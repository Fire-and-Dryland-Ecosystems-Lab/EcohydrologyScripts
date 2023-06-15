# preprocessing of spatial inputs

library(RHESSysPreprocessing)
library(terra)


# map generation/modification
dem_map = rast("preprocessing/spatial90m/dem.tif")
#dem_map = rast("preprocessing/spatial180m/dem.tif")

mask_map = rast("preprocessing/spatial90m/basin.tif")
mask_vect = as.polygons(mask_map)
writeVector(mask_vect, "preprocessing/spatial_source/basin_vect90m.shp")

p_map = dem_map
names(p_map) = "patches"
# resample to get unique patches 1 per cell
values(p_map)[!is.nan(values(p_map))] = seq_along(values(p_map)[!is.nan(values(p_map))])
writeRaster(p_map, "preprocessing/spatial90m/patches.tif", overwrite=T)
#writeRaster(p_map, "preprocessing/spatial180m/patches.tif", overwrite=T)

# soils corrections/reclass
# soils_map = rast("preprocessing/spatial180m/soil_texture.tif")
soils_map = rast("preprocessing/spatial90m/soil_texture.tif")

summary(as.factor(values(soils_map)))

#
# 1 clay		
# 2 silty-clay	
# 3 silty-clay-loam
# 4 sandy-clay	
# 5 sandy-clay-loam
# 6 clay-loam	
# 7 silt
# 8 silt-loam	
# 9 loam		
# 10 sand		
# 11 loamy-sand	
# 12 sandy-loam	

m = matrix(c(0, 12,
             8, 9,
             11, 12), 
           ncol = 2, byrow = T)

new_soils_map = classify(soils_map, m)
summary(as.factor(values(new_soils_map)))

# writeRaster(new_soils_map, "preprocessing/spatial90m/soils.tif", overwrite=T)
writeRaster(new_soils_map, "preprocessing/spatial90m/soils.tif", overwrite=T)

#soils_map2 = rast("preprocessing/spatial180m/soils.tif")
