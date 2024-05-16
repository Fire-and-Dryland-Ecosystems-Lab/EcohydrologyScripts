# Multiscale routing preprocess setup

library(RHESSysPreprocessing)
source("R/0_global_vars.R")

# -----------------------------------
# basin_name = "BigCreek"
basin_name = site
dest = "preprocessing/preprocess_out/"
# -----------------------------------

#map_dir = "preprocessing/spatial180m/"
overwrite = TRUE
streams = "streams.tif"
unique_strata_ID = TRUE
seq_patch_IDs = F
# asprules = NULL

# ---------------- MSR -----------------
convert_aspect = F
asprules = "preprocessing/rules/LPC_90m.rules"
map_dir = "preprocessing/whitebox/"

# ---------------- std -----------------
msrstd = F
if (msrstd) {
  config = "msr90m"
  template = "preprocessing/template/carb_msr.template"
  name = file.path(dest,paste0(basin_name,"_",config))
}

# ----------------- netcdf -----------------
msrncdf = T
if (msrncdf) {
  config = "msr90m_nc"
  template = "preprocessing/templates/carb_msr_nc.template"
  name = file.path(dest,paste0(basin_name,"_", config))
}

# ----------------- run -----------------
RHESSysPreprocess(template = template,
                  name = name,
                  map_dir = map_dir,
                  streams = streams,
                  overwrite = overwrite,
                  asprules = asprules,
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


