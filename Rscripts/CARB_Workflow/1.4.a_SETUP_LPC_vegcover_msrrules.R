# pctcover2MSRrules
# 
# Uses an input basin boundary map and set of lifeform percent cover maps (tree, shrub, herb, other)
# using these inputs it creates a set of MSR rules and a corresponding map
library(rhutils)
library(terra)
library(data.table)
source("../R/fun_MSR_thinning.R")
source("R/0_global_vars.R")

# ------------------------------ INPUTS ------------------------------
lpc_data_path = "../data/Landcover/2019/" 

# to get the grid
mask_map_path = "preprocessing/whitebox/basin.tif"

lpc_grid = vect("../data/Landcover/CONUS_C2_ARD_grid/conus_c2_ard_grid.shp")
mask_rast = rast(mask_map_path)
mask_rast_prj = project(mask_rast,lpc_grid)
grid_crop = crop(lpc_grid,mask_rast_prj)
grid_crop_vals = values(grid_crop)

lpc_pattern = paste0("h",sprintf("%0*d", 3, grid_crop_vals$h),"v", sprintf("%0*d", 3, grid_crop_vals$v))
# lpc_pattern = "h003v008"

round_value = 10

output_rules_map = "preprocessing/whitebox/rules_LPC_90m.tif"
output_rules_file = "preprocessing/rules/LPC_90m.rules"
output_std_veg_cover = "preprocessing/whitebox/LPC_veg_cover.tif"

# --------------------------------------------------------------------

# If your basin covers multiple lpc tiles you may need to edit this function
lpc_map_list = lpc_extract(lpc_data_path = lpc_data_path, lpc_pattern = lpc_pattern, mask_map_path = mask_map_path)

lpc_df = lpc_as_df(lpc_map_list)

lpc_bin_df = lpc_bin(lpc_df = lpc_df, round_value = round_value)

outlist = make_std_rules(lpc_bin_df)
lpc_rules_df = outlist[[1]] # cover data and rule ID for each patch
unique_rules = outlist[[2]] # unique cover combinations/rules for the basin

#baseline rules
rules_map_std = make_rules_map(lpc_rules_df, mask_map_path, lpc_map_list)
writeRaster(rules_map_std, output_rules_map, overwrite = T)

# LPC2std_veg_cover
LPC_std_vegcover = LPC2std_veg_cover(lpc_df, mask_map_path, lpc_map_list, vegid_tree = 1, vegid_shrub = 5, vegid_herb = 3, vegid_other = 4)
writeRaster(LPC_std_vegcover, output_std_veg_cover, overwrite = T)

write_rules_file(
  unique_rules,
  output_rules_file,
  vegid_tree = 1,
  vegid_shrub = 5,
  vegid_herb = 3,
  vegid_other = 4,
  strata_ct = NULL
)

