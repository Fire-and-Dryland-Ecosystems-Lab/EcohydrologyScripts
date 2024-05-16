# SETUP_thin_worlds
library(RHESSysIOinR)
library(rhutils)
library(tidyverse)
library(data.table)
source("../r/fun_shift_msr_world_vals.R")

# -----------------------------------
# basin_name = "BigCreek"
basin_name = site
thinrules = c(2,3,4,5)
dest = "preprocessing/preprocess_out/"
# -----------------------------------

#map_dir = "preprocessing/spatial180m/"
overwrite = TRUE
streams = "streams.tif"
unique_strata_ID = TRUE
seq_patch_IDs = F

# ---------------- MSR -----------------
convert_aspect = F
asprules = paste0("preprocessing/rules/LPC_90m_thin",thinrules,".rules")
map_dir = "preprocessing/whitebox/"

# ----------------- netcdf -----------------
msrncdf = T
if (msrncdf) {
  config = "msr90m_nc"
  template = "preprocessing/templates/carb_msr_nc.template"
  name = file.path(dest,paste0(basin_name,"_", config, "_thin",thinrules))
}

# source_world = "worldfiles/Ward_msr90m_basinsoilspin400y_reset.world"
# dest_worlds = paste0("preprocessing/preprocess_out/Ward_msr90m_thin",thinrules,".world")
# output_worlds = paste0("worldfiles/Ward_msr90m_thin",thinrules,"_reset.world")

for (i in seq_along(thinrules)) {

  RHESSysPreprocess(template = template,
                    name = name[i],
                    map_dir = map_dir,
                    streams = streams,
                    overwrite = overwrite,
                    asprules = asprules[i],
                    unique_strata_ID = unique_strata_ID,
                    seq_patch_IDs = seq_patch_IDs)
  
  # file.copy(from = paste0(name,".world"), to = file.path("worldfiles",basename(paste0(name,".world"))), overwrite = T)
  # shift_msr_world_vals(source_world, dest_worlds[i], output_worlds[i])
  
}



# This script transposes values between two MSR worldfiles, with the input worldfile containing 
# spun-up values you want to move to (a copy of) the second. This transposition is based on veg paramater ID,
# and assumes the input world only has 1 of each veg parameter, and then transposes each set of patch and strata
# variables based on veg parm, with the output/destination worldfile having 1 or more of each veg parm. Can 
# account for variable patch family compositions within the world.
# 
# Ex: 
# input world patches: ID 1, 2, 3          | veg parms 11, 22, 33
# output world patches: ID 1, 2, 3, 4, 5   | veg parms 11, 11, 22, 22, 33
# input world patch and strata vars from ID 1 get transposed to output world patches 1 & 2
# 

# ==================== START I/O ====================
# world_in_path = "worldfiles/Ward_msr90m_soilvegspun.world"
# world_out_path = "worldfiles/Ward_msr90m_thin2.world"
# world_output_filename = "worldfiles/Ward_msr90m_thin2_init.world"

# ==================== END I/O ====================



# shift_msr_world_vals(source_world, "worldfiles/Ward_msr90m_thin2.world", "worldfiles/Ward_msr90m_thin2_init.world")
# shift_msr_world_vals(source_world, "worldfiles/Ward_msr90m_thin3.world", "worldfiles/Ward_msr90m_thin3_init.world")
# shift_msr_world_vals(source_world, "worldfiles/Ward_msr90m_thin4.world", "worldfiles/Ward_msr90m_thin4_init.world")
# shift_msr_world_vals(source_world, "worldfiles/Ward_msr90m_thin5.world", "worldfiles/Ward_msr90m_thin5_init.world")
# shift_msr_world_vals(source_world, "worldfiles/Ward_msr90m_thin6.world", "worldfiles/Ward_msr90m_thin6_init.world")


