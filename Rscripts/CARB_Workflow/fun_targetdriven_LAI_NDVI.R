# LAI and target driven functions

# -------------------- Functions --------------------
getndvi = function(basin_path, landsat_dir, landsat_pattern, plots =  F) {
  basin = rast(basin_path)
  basin = trim(basin)
  ## read in landsat data
  landsat_files = list.files(path = landsat_dir, recursive = TRUE, pattern = landsat_pattern, full.names = TRUE)
  band3 = rast(landsat_files[str_detect(landsat_files,".*B3.*")])
  band4 = rast(landsat_files[str_detect(landsat_files,".*B4.*")])
  #band5 = rast(landsat_files[str_detect(landsat_files,".*B5.*")])
  # band3 = crop(band3, basin)
  # band4 = crop(band4, basin)
  
  # get the NDVI and band 8 including the time of when LAI is maximum
  # (NIR - R) / (NIR + R) band4 - band1 / band4 + band1
  # https://www.usgs.gov/core-science-systems/nli/landsat/landsat-normalized-difference-vegetation-index?qt-science_support_page_related_con=0#qt-science_support_page_related_con
  # In Landsat 4-7, NDVI = (Band 4 – Band 3) / (Band 4 + Band 3).
  # In Landsat 8-9, NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)
  
  ndvi = (band4 - band3)/(band4 + band3)
  ndvi_prj = project(ndvi, basin)
  ndvi_clean  = mask(ndvi_prj, basin)
  
  if (plots) {
    hist(ndvi_clean)
    plot(ndvi_clean)
  }
  return(ndvi_clean)
}

ndvi2lai = function(ndvi, ndvi_inf, ndvi_back) {
  tmp = (ndvi_inf - ndvi) / (ndvi_inf - ndvi_back)
  tmp[tmp <= 0] = 0
  lai = (-1/k) * log(tmp)
  pct99 = quantile(values(lai)[!is.infinite(values(lai)) & !is.nan(values(lai))], 0.99, na.rm=T)
  lai[lai>pct99] = pct99
  lai[lai<0] = 0
  names(lai) = "LAI"
  return(lai)
}


generate_spinup_targetfile = function(output_file, basin_map, hill_map, zone_map, patch_map, stratum_map, lai_map, correct_msr_patchIDs = F, correct_unique_strataIDs= F) {
  
  read_maps = try(rast(c(basin_map, hill_map, zone_map, patch_map, stratum_map, lai_map)), silent = T)
  if (class(read_maps) == "try-error") {
    maplist = lapply(c(basin_map, hill_map, zone_map, patch_map, stratum_map, lai_map), function(X) {trim(rast(X))})
    read_maps = rast((maplist))
  }
  names(read_maps) = c("basin_ID","hill_ID","zone_ID","patch_ID","stratum_ID","LAI")
  # fix patch and stratum IDs
  if (correct_msr_patchIDs)
    values(read_maps$patch_ID) = values(read_maps$patch_ID) * 100 + 1
  if (correct_unique_strataIDs)
    values(read_maps$stratum_ID) = values(read_maps$patch_ID) * 10 + 1
  
  plot(read_maps)
  
  ID_df = as.data.frame(read_maps)
  
  # summary(ID_df)
  # cat("Any nans in basin: " ,any(is.nan(ID_df$basin_ID)),"\n")
  # cat("Any nans in lai: " ,any(is.nan(ID_df$LAI)),"\n")
  
  ID_df_agg = aggregate(ID_df$LAI, by=list(ID_df$patch_ID, ID_df$zone_ID, ID_df$hill_ID, ID_df$basin_ID), mean)
  
  if (nrow(ID_df) > nrow(ID_df_agg)) {
    ID_df_out = ID_df_agg
  } else {
    ID_df_out = ID_df
  }
  
  # summary(ID_df)
  # summary(ID_df_agg)
  # nrow(ID_df)
  # nrow(ID_df_agg)
  
  # any(is.na(ID_df$LAI))
  # any(is.infinite(ID_df$LAI))
  # any(is.null(ID_df$LAI))
  # any(ID_df$LAI < 0)
  
  header = sprintf("%d num_stratum\n%d num_targets\nLAI", nrow(ID_df_out), 1)
  write(header, file=output_file)
  write.table(ID_df_out, file = output_file, append=T, quote=F, row.names=F, col.names = T)
  cat("Wrote target file:",output_file,"\n") 
}

