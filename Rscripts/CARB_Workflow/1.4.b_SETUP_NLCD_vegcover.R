# ------------------------------ ALTERNATIVE = GET NLCD ------------------------------
run = F
if (run) {
  # devtools::install_github("ropensci/FedData")
  library(FedData)
  
  maskPolygon <- polygon_from_extent(raster::extent(ext(mask_map)[1:4]), proj4string='+proj=utm +datum=NAD83 +zone=11')
  
  NLCD <- get_nlcd(template=maskPolygon, label='VEPIIN', force.redo = T, extraction.dir = paste0("preprocessing/spatial_source/NLCD/"))
  
  NLCD_proj = project(as(NLCD, "SpatRaster"), mask_map, method = "near")
  NLCD_crop = crop(NLCD_proj, mask_map)
  NLCD_mask = mask(NLCD_crop, mask_map)
  
  plot(NLCD_mask, type = "classes")
  hist(NLCD_mask)
  nlcd_colors()
  
  summary(as.factor(values(NLCD_mask)))
  
  NLCD_mask_reclass = terra::classify(NLCD_mask, data.frame(c(11, 21, 22, 31, 42, 52, 71, 90, 95), c(31, 31, 31, 31, 7, 50, 50, 7, 50) ))
  
  writeRaster(NLCD_mask_reclass, "preprocessing/spatial90m/NLCD_veg_cover.tif", overwrite = T)
  
  # PLOT for saving
  par(mfrow = c(1,2))
  plot(LPC_veg_cover_map, type = "classes", main ="LPC Veg Cover")
  plot(NLCD_mask_reclass, type = "classes", main = "NLCD Veg Cover")
  dev.off()
}