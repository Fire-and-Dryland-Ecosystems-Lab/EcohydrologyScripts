# setup thinning rules
library(rhutils)
library(terra)
library(data.table)
source("R/0_global_vars.R")
source("../r/fun_MSR_thinning.R")


# ------------------------------ INPUTS ------------------------------
lpc_data_path = "../data/Landcover/2019/" 

# to get the grid
mask_map_path = "preprocessing/whitebox/basin.tif"

lpc_grid = vect("../data/Landcover/CONUS_C2_ARD_grid/conus_c2_ard_grid.shp")
mask_rast = rast(mask_map_path)
mask_rast_prj = project(mask_rast,crs(lpc_grid))
grid_crop = crop(lpc_grid,mask_rast_prj)
grid_crop_vals = values(grid_crop)

lpc_pattern = paste0("h",sprintf("%0*d", 3, grid_crop_vals$h),"v", sprintf("%0*d", 3, grid_crop_vals$v))
# lpc_pattern = "h003v008"

round_value = 10

output_rules_map = "preprocessing/whitebox/rules_LPC_90m.tif"
output_rules_file = "preprocessing/rules/LPC_90m.rules"
output_std_veg_cover = "preprocessing/whitebox/LPC_veg_cover.tif"

# --------------------------------------------------------------------

lpc_map_list = lpc_extract(lpc_data_path = lpc_data_path, lpc_pattern = lpc_pattern, mask_map_path = mask_map_path)

lpc_df = lpc_as_df(lpc_map_list)

lpc_bin_df = lpc_bin(lpc_df = lpc_df, round_value = round_value)

outlist = make_std_rules(lpc_bin_df)
lpc_rules_df = outlist[[1]]
unique_rules = outlist[[2]]

rules_map_std = make_rules_map(lpc_rules_df, mask_map_path, lpc_map_list)

# These steps shouldve already happened in previous step

# writeRaster(rules_map_std, output_rules_map, overwrite = T)
# 
# write_rules_file(
#   unique_rules,
#   output_rules_file,
#   vegid_tree = 1,
#   vegid_shrub = 5,
#   vegid_herb = 3,
#   vegid_other = 4,
#   strata_ct = NULL
# )

# ideally turn this into a loop that pulls directly from the excel doc for thin values
# treatments = readxl::read_excel("../data/FuelTreatments/CARB_RHESSys_FuelTreatments.xlsx", sheet = 1)

# 1 is NO TREATMENT

# 2 rx fire
tree_thin_pct = 0.05
shrub_thin_pct = 0.6
output_rules_file_thin2 = "preprocessing/rules/LPC_90m_thin2.rules"
thin2_rules = make_thin_rule(unique_rules = unique_rules, tree_thin_pct = tree_thin_pct, shrub_thin_pct = shrub_thin_pct)
write_rules_file(thin2_rules, output_rules_file_thin2, thin_vegid_mod = 20)

# 3 heavy harvest
tree_thin_pct = 0.4
shrub_thin_pct = 0.9
output_rules_file_thin3 = "preprocessing/rules/LPC_90m_thin3.rules"
thin3_rules = make_thin_rule(unique_rules = unique_rules, tree_thin_pct = tree_thin_pct, shrub_thin_pct = shrub_thin_pct)
write_rules_file(thin3_rules, output_rules_file_thin3, thin_vegid_mod = 20)

# 4 moderate thinning
tree_thin_pct = 0.1
shrub_thin_pct = 0.75
output_rules_file_thin4 = "preprocessing/rules/LPC_90m_thin4.rules"
thin4_rules = make_thin_rule(unique_rules = unique_rules, tree_thin_pct = tree_thin_pct, shrub_thin_pct = shrub_thin_pct)
write_rules_file(thin4_rules, output_rules_file_thin4, thin_vegid_mod = 20)

# 5 mastication
tree_thin_pct = 0.1
shrub_thin_pct = 0.9
output_rules_file_thin5 = "preprocessing/rules/LPC_90m_thin5.rules"
thin5_rules = make_thin_rule(unique_rules = unique_rules, tree_thin_pct = tree_thin_pct, shrub_thin_pct = shrub_thin_pct)
write_rules_file(thin5_rules, output_rules_file_thin5, thin_vegid_mod = 20)

# 6 clearcut
tree_thin_pct = 0.8
shrub_thin_pct = 0.8
output_rules_file_thin6 = "preprocessing/rules/LPC_90m_thin6.rules"
thin6_rules = make_thin_rule(unique_rules = unique_rules, tree_thin_pct = tree_thin_pct, shrub_thin_pct = shrub_thin_pct)
write_rules_file(thin6_rules, output_rules_file_thin6, thin_vegid_mod = 20)

