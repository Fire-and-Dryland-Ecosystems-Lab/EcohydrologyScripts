# plot site locations 
library(xlsx)

sites =  read.xlsx("../Documents/Possible Watersheds.xlsx", sheetIndex = 2)

library(sf)
library(terra)
CA = st_read("../data/ca-state-boundary/CA_State_TIGER2016.shp")

pts = st_as_sf(x = sites, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326)
#pts_proj = sf_project(st_crs(pts), st_crs(CA), pts)

ca_proj = st_transform(CA, crs= st_crs(pts))

deficit = rast("../data/def_1989_2019.tif")

crs(deficit)
crs(pts)

deficit_prj = project(deficit, crs(pts))

plot(deficit_prj,main=NULL, axes = F, legend=F)
plot(st_geometry(pts), pch =16, axes=T, add=T)
p = recordPlot()


plot(st_geometry(ca_proj), main=NULL)
plot(st_geometry(pts), pch =16, axes=T, add=T)

plot(st_geometry(ca_proj), main=NULL, add=T)
plot(pts["STANAME"], pch =16, axes=T)
#plot(pts["STANAME"], pch =16)



patches = trim(rast("preprocessing/spatial90m/patches.tif"))
plot(patches)

psub = patches
values(psub)[values(psub) != 1] = NA
psub = trim(psub)

writeRaster(psub,"preprocessing/spatial90m/ward_subset_vegid1.tif")

terra::project(psub, "EPSG:4269")


