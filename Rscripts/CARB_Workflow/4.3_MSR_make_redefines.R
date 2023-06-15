# make redefines for treatments
library(rhutils)
library(RHESSysIOinR)

redef_maps = F
if (redef_maps) {
  # Build redefines for low and high area treated scenarios
  library(terra)
  patch = rast("preprocessing/BC90m_tiff/patch.tif")
  slope = rast("preprocessing/BC90m_tiff/slope.tif")
  # slope > 50% rise, or 26.6 degrees (map is in degrees)
  # going to use 40% or 22 deg
  max_thin = patch[slope < 22]
  max_thin_fam = max_thin[!is.na(max_thin)]
  max_thin_shrub = as.numeric(c(paste0(max_thin_fam, "012"), paste0(max_thin_fam, "022"),  paste0(max_thin_fam, "032")) )
  max_thin_conifer = as.numeric(paste0(max_thin_fam, c("021")))
  max_thin_map = slope < 22
  plot(max_thin_map)
  patch_vec = values(patch, mat = F, na.rm = T)
  min_thin = sample(patch_vec, size = length(patch_vec) * 0.2)
  #save(min_thin,file = "r_obj/min_thin_sample.rdata")
  load("r_obj/min_thin_sample.rdata")
  length(min_thin) / length(patch_vec)
  plot(patch %in% min_thin)
  min_thin_map = patch %in% min_thin
  min_thin_map = mask(min_thin_map, patch)
  plot(min_thin_map)
  min_thin_shrub = as.numeric(c(paste0(min_thin, "012"), paste0(min_thin, "022"), paste0(min_thin, "032")))
  min_thin_conifer = as.numeric(paste0(min_thin, "021") )
  writeRaster(max_thin_map, "preprocessing/BC90m_tiff/max_thin.tif")
  writeRaster(min_thin_map, "preprocessing/BC90m_tiff/min_thin.tif", overwrite=T)
}



# ------------------------------ THIN REDEFINES ------------------------------
# Thin2 rxfire - 100% thin_harvest shrub (0% thin_remain), 0% thin_harvest tree (100% thin_remain), 100% remove litter
worldfile = "worldfiles/Ward_msr90m_thin2.world"
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin2_redef1_remain.world" ,std_thin = "1", veg_parm_ID = 21) # remain
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin2_redef2_harvest.world" ,std_thin = "1", veg_parm_ID = 25) # harvest
litter_vars = c("litter_cs.litr1c", "litter_ns.litr1n",	"litter_cs.litr2c", "litter_cs.litr3c", "litter_cs.litr4c", "cs.cwdc", "ns.cwdn")
build_redefine2(worldfile, out_file = "worldfiles/redefs/Ward_msr90m_thin2_redef3_litter.world" , vars = litter_vars, values = 0)

# thin3 - 90% thin_harvest shrub + tree, 100% thin_remain
worldfile = "worldfiles/Ward_msr90m_thin3.world"
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin3_redef1_harvest.world" ,std_thin = "0.9", veg_parm_ID = c(21,25)) # harvest
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin3_redef2_remain.world" ,std_thin = "1", veg_parm_ID = c(21,25)) # remain

# thin4 - 50% thin_harvest shrub + tree, 100% thin_remain
worldfile = "worldfiles/Ward_msr90m_thin4.world"
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin4_redef1_harvest.world" ,std_thin = "0.5", veg_parm_ID = c(21,25)) # harvest
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin4_redef2_remain.world" ,std_thin = "1", veg_parm_ID = c(21,25)) # remain

# thin5 - 0% thin_harvest shrub + tree, 100% thin_remain
worldfile = "worldfiles/Ward_msr90m_thin5.world"
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin5_redef1_remain.world" ,std_thin = "1", veg_parm_ID = c(21,25)) # remain


# thin6 - 20% thin_harvest shrub, 90% harvest tree, 100% thin_remain
worldfile = "worldfiles/Ward_msr90m_thin6.world"
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin6_redef1_harvest.world" ,std_thin = "0.2", veg_parm_ID = c(25)) # harvest
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin6_redef2_harvest.world" ,std_thin = "0.9", veg_parm_ID = c(21)) # harvest
build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin6_redef3_remain.world" ,std_thin = "1", veg_parm_ID = c(21,25)) # remain

# worldfile = "worldfiles/Ward_msr90m_thin2.world"
# 
# world = read_world(worldfile = worldfile)
# 
# familyIDs = unique(world[world$vars == "family_ID", "values"])
# 
# thin_strata_list = lapply(familyIDs, function(X) {
#   strataIDs = unique(world[world$level == "canopy_strata" & world$ID %in% paste0(X, c("01", "02", "03" ,"04" ,"05" ,"06", "07"),"1") , "ID"])
#   veg_parms = world[world$level == "canopy_strata" & world$ID %in% paste0(X, c("01", "02", "03" ,"04" ,"05" ,"06", "07"),"1") & world$vars == "veg_parm_ID", "values"]
#   thin_strata = strataIDs[duplicated(veg_parms, fromLast = T)]
#   return(thin_strata)
# })
# 
# thin_strata = unlist(thin_strata_list)
# 
# build_redefine(worldfile, "worldfiles/redefs/Ward_msr90m_thin1_redef.world" ,
#                std_thin = "1", strataID = thin_strata)





# ------------------------------  ------------------------------