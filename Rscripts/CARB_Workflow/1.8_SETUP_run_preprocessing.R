# Multiscale routing preprocess setup

library(RHESSysPreprocessing)
library(rhutils)
# source("R/0_global_vars.R")

# -----------------------------------
name = "preprocessing/preprocess_out/WallaWallaNCdaymet"
map_dir = "preprocessing/whitebox/"
template = "preprocessing/template/walla_NCdaymet.template"
check_template(template)
# -----------------------------------
overwrite = TRUE
streams = "streams.tif"
unique_strata_ID = TRUE
seq_patch_IDs = F
convert_aspect = F

# ----------------- run -----------------
RHESSysPreprocess(template = template,
                  name = name,
                  map_dir = map_dir,
                  streams = streams,
                  overwrite = overwrite,
                  unique_strata_ID = unique_strata_ID,
                  seq_patch_IDs = seq_patch_IDs,
                  convert_aspect = convert_aspect)

file.copy(from = paste0(name,".world"), to = file.path("worldfiles",basename(paste0(name,".world"))), overwrite = T)
file.copy(from = paste0(name,".flow"), to = file.path("flowtables",basename(paste0(name,".flow"))), overwrite = T)
















 # ---------- DEBUG -----------
roads = NULL
impervious = NULL
roofs = NULL
header = FALSE
meta = FALSE
output_patch_map = FALSE
parallel = TRUE
make_stream = 4
flownet_name = name
road_width = NULL
skip_hillslope_check = TRUE




create_flownet(flownet_name = name,
               template = template,
               map_dir = map_dir,
               asprules = asprules,
               streams = streams,
               overwrite = overwrite)


