# Generate target file - target driven spinup
library(terra)
source("fun_targetdriven_LAI_NDVI.R")

basin_map = "preprocessing/spatial90m/basin.tif"
hill_map = "preprocessing/spatial90m/subbasins.tif"
zone_map = "preprocessing/spatial90m/patches.tif"
patch_map = "preprocessing/spatial90m/patches.tif"
stratum_map = patch_map
lai_map = "preprocessing/spatial90m/lai.tif"

output_file = "defs/spinup_LAI_targets.txt"
# output_file = "defs/STD_spinup_LAI_targets.txt"

correct_msr_patchIDs = T
# correct_msr_patchIDs = F
correct_unique_strataIDs = T

generate_spinup_targetfile(output_file, basin_map, hill_map, zone_map, patch_map, stratum_map, lai_map, correct_msr_patchIDs, correct_unique_strataIDs)


# tmp = check_worldfile("preprocessing/preprocess_out/Pitman_STD_90m_2strata.world")
#
# tmp = read_world("preprocessing/preprocess_out/Pitman_STD_90m_2strata.world")
# pids = unique(tmp$values[tmp$vars == "patch_ID"])
# pids = as.numeric(pids)
# min(pids)
# max(pids)
# length(pids)
