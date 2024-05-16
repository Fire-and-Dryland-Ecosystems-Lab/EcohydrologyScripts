# 1.11_SETUP_basin_stats
# 
# Compile and write file of stats about the basin
# Basin-scale & hillslope scale
# - Vegparmid pct areas
# - slope(s)
# - elevations(s)

library(RHESSysPreprocessing)
library(tidyverse)
library(knitr)

########## INPUTS ##########

worldfile = "preprocessing/preprocess_out/BigCreek_msr90m_nc.world"
outfile = "preprocessing/basin_hillslope_stats.txt"

#   vegid_tree = 1,
#   vegid_shrub = 5,
#   vegid_herb = 3,
#   vegid_other = 4,

###############

tmp = basin_hillslope_stats(worldfile = worldfile, outfile = outfile, LPC_vegnames = T)

# come back to thi

# basin_hillslope_stats = function(worldfile, outfile = NULL, LPC_vegnames = T) {
#   
#   cat("Reading in worldfile ...\n")
#   world = read_world(worldfile = worldfile,hill_col = T, patch_col = T, zone_col = T)
#   cat("Completed reading worldfile.\n")
#   
#   # Basin-scale - Vegparmid pct areas
#   pjoin = merge(world[world$vars == "area" & world$level == "patch", c("values","patch_ID", "hillslope_ID")], 
#                 world[world$vars == "veg_parm_ID",c("values","patch_ID")], by = "patch_ID", all = T)
#   veg_areas = aggregate(as.numeric(pjoin$values.x), by = list(pjoin$values.y), FUN = sum)
#   names(veg_areas) = c("veg_parm_ID","area")
#   totalarea = sum(as.numeric(world[world$vars == "area" & world$level == "zone","values"]))
#   veg_areas$percent_area = veg_areas$area/totalarea
#   veg_areas$veg_names = veg_areas$veg_parm_ID
#   
#   if (LPC_vegnames) {
#     veg_areas$veg_names = case_match(veg_areas$veg_names, "1" ~ "tree", "5" ~ "shrub","3" ~ "herb","4" ~ "bare")
#   } else {
#     veg_areas$veg_names = paste0("vegID_",veg_areas$veg_names)
#   }
#   
#   # hillslope veg areas
#   hill_veg_areas = aggregate(as.numeric(pjoin$values.x), by = list(pjoin$hillslope_ID, pjoin$values.y), FUN = sum)
#   names(hill_veg_areas) = c("hillslope_ID", "veg_parm_ID","area")
#   hill_veg_areas = hill_veg_areas[order(hill_veg_areas$hillslope_ID,hill_veg_areas$veg_parm_ID),]
#   
#   hill_totalarea = aggregate(as.numeric(world[world$vars == "area" & world$level == "zone", c("values" )]), 
#                              by = list(world[world$vars == "area" & world$level == "zone", c("hillslope_ID" )]),
#                              FUN = sum)
#   names(hill_totalarea) = c("hillslope_ID", "hill_totalarea")
#   
#   hill_veg_areas = merge(hill_veg_areas, hill_totalarea, by = "hillslope_ID")
#   hill_veg_areas$percent_area = hill_veg_areas$area/hill_veg_areas$hill_totalarea
#   hill_veg_areas$veg_names = hill_veg_areas$veg_parm_ID
#   
#   if (LPC_vegnames) {
#     hill_veg_areas$veg_names = case_match(hill_veg_areas$veg_names, "1" ~ "tree", "5" ~ "shrub","3" ~ "herb","4" ~ "bare")
#   } else {
#     hill_veg_areas$veg_names = paste0("vegID_",hill_veg_areas$veg_names)
#   }
#   
#   # hill_veg_areas %>% group_by(veg_names) %>% 
#   #   summarize(mean = mean(percent_area, na.rm = TRUE),
#   #             median = median(percent_area, na.rm = TRUE),
#   #             min = min(percent_area, na.rm = TRUE),
#   #             max = max(percent_area, na.rm = TRUE))
#   
#   # basin - slope
#   zonevars = cbind(world[world$vars == "slope" & world$level == "zone", ],
#                    area = world[world$vars == "area" & world$level == "zone", "values"])
#   names(zonevars)[names(zonevars) =="values"] = "slope"
#   zonevars$slope = as.numeric(zonevars$slope)
#   zonevars$area = as.numeric(zonevars$area)
#   basin_slope = sum(zonevars$slope * (zonevars$area/sum(zonevars$area)) )
#   
#   # hillslope - slope
#   # zonevars = merge(zonevars, hill_totalarea, by = "hillslope_ID")
#   zonevars$slope_mult = zonevars$slope * zonevars$area
#   hillslope_slope = aggregate(zonevars$slope_mult, by = list(zonevars$hillslope_ID), FUN = sum)
#   names(hillslope_slope)[1] = "hillslope_ID"
#   hillslope_slope = merge(hillslope_slope, hill_totalarea)
#   names(hillslope_slope)[names(hillslope_slope) == "hill_totalarea"] = "hillslope_area"
#   hillslope_slope$slope = hillslope_slope$x/hillslope_slope$hillslope_area
#   hillslope_slope$x = NULL
#   
#   # basin - elevations
#   basin_elevation = as.numeric(world[world$var == "z" & world$level == "basin","values"])
#   
#   # hillslope - elevation
#   hillslope_elevation = world[world$vars == "z" & world$level == "hillslope", c("hillslope_ID","values")]
#   rownames(hillslope_elevation) = NULL
#   names(hillslope_elevation)[2] = "elevation"
#   
#   hillslope_df = merge(hillslope_slope, hillslope_elevation)
#   hill_veg_wider = pivot_wider(hill_veg_areas[,c("hillslope_ID", "veg_names","percent_area")], names_from = "veg_names", values_from = "percent_area") 
#   hillslope_df = merge(hillslope_df, hill_veg_wider)
#   
#   basin_df = data.frame(
#     area_sq_m = sum(hill_totalarea$hill_totalarea), area_sq_km = sum(hill_totalarea$hill_totalarea)/1000000, 
#     elevation = basin_elevation, slope = basin_slope,
#     num_hillslopes = nrow(hillslope_df), mean_hill_area_sq_km = mean(hillslope_df$hillslope_area)/1000000, 
#     num_zones = length(unique(world$zone_ID)), mean_zone_area = mean(as.numeric(world[world$vars == "area" & world$level == "zone","values"])), 
#     num_patches = length(unique(world$patch_ID)), mean_patch_area = mean(as.numeric(world[world$vars == "area" & world$level == "patch","values"])),
#     pivot_wider(veg_areas[,c("veg_names","percent_area")], names_from = "veg_names", values_from = "percent_area") 
#   )
#   
#   basin_df_tbl = basin_df
#   hillslope_df_tbl = hillslope_df
#   
#   if (LPC_vegnames) {
#     basin_df_tbl$tree = paste0(format(basin_df_tbl$tree * 100, digits = 1), "%")
#     basin_df_tbl$herb = paste0(format(basin_df_tbl$herb * 100, digits = 1), "%")
#     basin_df_tbl$bare = paste0(format(basin_df_tbl$bare * 100, digits = 1), "%")
#     basin_df_tbl$shrub = paste0(format(basin_df_tbl$shrub * 100, digits = 1), "%")
#     hillslope_df_tbl$tree = paste0(format(hillslope_df_tbl$tree * 100, digits = 1), "%")
#     hillslope_df_tbl$herb = paste0(format(hillslope_df_tbl$herb * 100, digits = 1), "%")
#     hillslope_df_tbl$bare = paste0(format(hillslope_df_tbl$bare * 100, digits = 1), "%")
#     hillslope_df_tbl$shrub = paste0(format(hillslope_df_tbl$shrub * 100, digits = 1), "%")
#     dgts = 1
#   } else  {
#     dgts = 2
#   }
#   
#   kbl_basin = kable(basin_df_tbl, caption = "Basin stats", digits = dgts)
#   kbl_hill = kable(hillslope_df_tbl, caption = "Hillslope stats", digits = dgts)
#   print(kbl_basin)
#   print(kbl_hill)
#   
#   # output
#   if (!is.null(outfile)) {
#     writeLines(capture.output(c(kbl_basin, kbl_hill)), outfile)
#   }
#   
#   return(list(basin_stats = basin_df, hillslopes_stats = hillslope_df))
# }