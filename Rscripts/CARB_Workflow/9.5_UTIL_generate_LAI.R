
library(terra)
library(stars)
library(stringr)
source("fun_targetdriven_LAI_NDVI.R")

# -------------------- Inputs --------------------
LPC_ver = T
nlcd_ver = F
basin_path = "preprocessing/spatial90m/basin.tif"
#landsat
landsat_dir = "preprocessing/landsat/"
landsat_pattern = "T1_B\\d\\.TIF$" # at end file name, gets only the bands
# maybe check https://github.com/logan-berner/LandsatTS in future

#veg cover
nlcd_path = "~/Projects/data/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img"
lpc_veg_path = "preprocessing/spatial90m/LPC_veg_cover.tif"

# input basin
basin = rast(basin_path)
basin = trim(basin)
basinnoprj = project(basin, "epsg:4326")
writebasin = F
if (writebasin){
  basin.poly = as.polygons(basin)
  writeVector(basin.poly,paste0(gsub(".tif","", basin_path),".shp"), overwrite=T)
}

# -------------------- NDVI input/parsing --------------------
ndvi_clean = getndvi(basin_path = basin_path, landsat_dir = landsat_dir, landsat_pattern = landsat_pattern, plots = T)

# -------------------- Veg cover - NLCD or other --------------------
if (nlcd_ver) {
  # get veg cover type
  nlcd = rast(nlcd_path)
  nlcd_prj = project(nlcd, ndvi_clean)
  nlcd_crop = crop(nlcd_prj, ndvi_clean)
  nlcd_mask = mask(nlcd_crop, ndvi_clean)
  nlcd_mask = droplevels(nlcd_mask)

  summary(nlcd_mask)
  cats(nlcd_mask)
  levels(nlcd_mask)

  # Do some lumping to deal with outliers
  # mixed forest, deciduous,woody wetland -> evergreen forest. emergent herb -> herb
  veg_cover_clean = subst(nlcd_mask,
                          c("Woody Wetlands", "Developed, Open Space", "Developed, Low Intensity",
                            "Developed, Medium Intensity", "Emergent Herbaceous Wetlands", "Barren Land", "Open Water"),
                          c("Evergreen Forest", "Evergreen Forest", "Evergreen Forest", "Evergreen Forest", "Herbaceous", "Herbaceous", "Herbaceous"))

  veg_cover_clean = droplevels(veg_cover_clean)
  names(veg_cover_clean) = "veg_cover"
  names(levels(veg_cover_clean)[[1]])[2] = "veg_cover"
  summary(veg_cover_clean)
  plot(veg_cover_clean)
}

if (LPC_ver) {
  lpc_veg = trim(rast(lpc_veg_path))
  levels(lpc_veg) = data.frame(value = c(1,3,4,5), veg_cover = c("Evergreen Forest", "Herbaceous","Other", "Shrub/Scrub" ))
  plot(lpc_veg)
  veg_cover_clean = lpc_veg
}

# -------------------- NDVI inf and back notes --------------------
# From Hanan et al 2018
# NDVI∞ is the maximum NDVI observed in each region and
# NDVIback is the background NDVI (in the absence of vegetation) for each region.
# k is a parameter corresponding to the extinction of solar radiation through a canopy
# Table 2:
# Vegetation        k       NDVI∞       NDVIback
# Johnson Creek
#     Pine          0.40    0.98        0.07
#     Grass         0.59    0.90        0.07
#     Shrub         0.55    0.90        0.07
#     Deciduous     0.54    0.88        0.07
# Rattlesnake Canyon
#     Shrub         0.55    0.92        0.11

# LAI = -1/k * ln((ndvi.inf - ndvi) / (ndvi.inf - ndvi.back))

## ========== k/epc.ext_coef/beer's law from param defs
## 0.400000	   epc.ext_coef ## pine 1
## 0.550000	   epc.ext_coef ## shrub 5
## 0.371000	   epc.ext_coef chaparrel 5
## 0.48, 0.5   grass
## 0.4, 0.5          evergreen
## 0.54         deciduous
## 0.5          conifer p301


# -------------------- NDVI inf --------------------

# Determine NDVIinf for each veg type (i.e., largest NDVI)
veg_NDVIinf = levels(veg_cover_clean)[[1]]
veg_NDVIinf$min = sapply(veg_NDVIinf$value, function(X){min(ndvi_clean[veg_cover_clean == X], na.rm = T)})
veg_NDVIinf$max = sapply(veg_NDVIinf$value, function(X){max(ndvi_clean[veg_cover_clean == X], na.rm = T)})
# handle outliers, use 2% instead?
veg_NDVIinf$min2pct = sapply(veg_NDVIinf$value, function(X){quantile(ndvi_clean[veg_cover_clean == X], probs = 0.02, na.rm=T)})
veg_NDVIinf$max2pct = sapply(veg_NDVIinf$value, function(X){quantile(ndvi_clean[veg_cover_clean == X], probs = 0.98, na.rm=T)})
# take mean of values in 2 percentile?
veg_NDVIinf$min2pct_avg = mapply(function(X,Y){mean(ndvi_clean[veg_cover_clean == X & ndvi_clean < Y], na.rm=T)}, veg_NDVIinf$value, veg_NDVIinf$min2pct)
veg_NDVIinf$max2pct_avg = mapply(function(X,Y){mean(ndvi_clean[veg_cover_clean == X & ndvi_clean > Y], na.rm=T)}, veg_NDVIinf$value, veg_NDVIinf$max2pct)

# these values are taken from comment section above
veg_NDVIinf$k[veg_NDVIinf$veg_cover == "Evergreen Forest"] = 0.5
veg_NDVIinf$k[veg_NDVIinf$veg_cover == "Shrub/Scrub"] = 0.55
veg_NDVIinf$k[veg_NDVIinf$veg_cover == "Herbaceous"] = 0.48
if (LPC_ver) {
  veg_NDVIinf$k[veg_NDVIinf$veg_cover  == "Other"] = 0.8
}

# starting from the nlcd map, reclass the values from the veg cover to the ndvi inf
NDVIinf = subst(as.numeric(veg_cover_clean), veg_NDVIinf$value, veg_NDVIinf$max, raw=T)

# k based on the k included for each veg cover
k = subst(as.numeric(veg_cover_clean), veg_NDVIinf$value, veg_NDVIinf$k, raw=T)

# find the NDVIBACK for whole image and brightest value
# Determine NDVIback for image (i.e., average of 20 brightest pixels)
NDVIback = mean(tail(sort(values(ndvi_clean),decreasing=TRUE),20) )

laivis = function(lai) {
  par(mfrow = c(1,2))
  plot(lai)
  hist(lai)
  summary(lai)
}

lai = ndvi2lai(ndvi = ndvi_clean, ndvi_inf = NDVIinf, ndvi_back = NDVIback)
laivis(lai)

NDVIinf2 = subst(as.numeric(veg_cover_clean), veg_NDVIinf$value, veg_NDVIinf$max2pct_avg, raw=T)
lai2 = ndvi2lai(ndvi = ndvi_clean, ndvi_inf = NDVIinf2, ndvi_back = NDVIback)
laivis(lai2)

NDVIinf3 = subst(as.numeric(veg_cover_clean), veg_NDVIinf$value, veg_NDVIinf$max2pct, raw=T)
lai3 = ndvi2lai(ndvi = ndvi_clean, ndvi_inf = NDVIinf3, ndvi_back = NDVIback)
laivis(lai3)


options("scipen"=999, "digits"=3)

# CHOOSE A LAI
# writeRaster(lai,"preprocessing/spatial90m/lai.asc", overwrite=T)
writeRaster(lai,"preprocessing/spatial90m/lai.tif", overwrite=T)
