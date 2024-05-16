# preprocessing of spatial inputs

library(RHESSysPreprocessing)
library(terra)
# terra::aggregate()
# terra::resample()

# map generation/modification
dem_map = rast("preprocessing/whitebox/dem.tif")
#dem_map = rast("preprocessing/spatial180m/dem.tif")

mask_map = dem_map
mask_vect = as.polygons(mask_map)
writeVector(mask_vect, "preprocessing/spatial_source/basin_vect90m.shp")

p_map = dem_map
names(p_map) = "patches"
# resample to get unique patches 1 per cell
values(p_map)[!is.nan(values(p_map))] = seq_along(values(p_map)[!is.nan(values(p_map))])
writeRaster(p_map, "preprocessing/whitebox/patches.tif", overwrite=T)
#writeRaster(p_map, "preprocessing/spatial180m/patches.tif", overwrite=T)






# soils corrections/reclass
# soils_map = rast("preprocessing/spatial180m/soil_texture.tif")
# soils_map = rast("preprocessing/spatial_source/soil_texture.tif")
# 
# clay = rast("preprocessing/spatial_source/POLARISOut/mean/clay/0_5/lat3738_lon-121-120.tif")
# sand = rast("preprocessing/spatial_source/POLARISOut/mean/sand/0_5/lat3738_lon-121-120.tif")
# 
# claymask = crop(project(clay, dem_map), dem_map)
# claymask = mask(claymask, dem_map)
# sandmask = crop(project(sand, dem_map), dem_map)
# sandmask = mask(sandmask, dem_map)
# writeRaster(claymask, "preprocessing/spatial_source/POLARISOut/claymask.tif", overwrite=T)
# writeRaster(sandmask, "preprocessing/spatial_source/POLARISOut/sandmask.tif", overwrite=T)

soils_map = rast("preprocessing/spatial90m/soil_texture_masked.tif")
soils_map = mask(soils_map, rast("preprocessing/whitebox/patches.tif"))
summary(as.factor(values(soils_map)))
plot(soils_map)
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

m = matrix(c(12, 9), 
           ncol = 2, byrow = T)

new_soils_map = classify(soils_map, m)
summary(as.factor(values(new_soils_map)))

plot(new_soils_map)

# writeRaster(new_soils_map, "preprocessing/spatial90m/soils.tif", overwrite=T)
writeRaster(new_soils_map, "preprocessing/whitebox/soils.tif", overwrite=T)

#soils_map2 = rast("preprocessing/spatial180m/soils.tif")
