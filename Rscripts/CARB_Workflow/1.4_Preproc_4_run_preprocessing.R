# Multiscale routing preprocess setup

library(RHESSysPreprocessing)

type = "raster"
typepars = "preprocessing/spatial90m/"
#typepars = "preprocessing/spatial180m/"
overwrite = TRUE
streams = "streams.tif"
unique_strata_ID = TRUE
seq_patch_IDs = F
asprules = NULL

# ---------- STD -----------
#name = "preprocessing/preprocess_out/Ward180m"
# name = "preprocessing/preprocess_out/Ward90m"
# template = "preprocessing/template/carb_std.template"
name = "preprocessing/preprocess_out/Ward90m_nc"
template = "preprocessing/template/carb_std_nc.template"

RHESSysPreprocess(template = template,
                  name = name,
                  type = type,
                  typepars = typepars,
                  streams = streams,
                  overwrite = overwrite,
                  asprules = asprules,
                  unique_strata_ID = unique_strata_ID,
                  seq_patch_IDs = seq_patch_IDs)

file.copy(from = paste0(name,".world"), to = file.path("worldfiles",basename(paste0(name,".world"))), overwrite = T)
file.copy(from = paste0(name,".flow"), to = file.path("flowtables",basename(paste0(name,".flow"))), overwrite = T)

# ---------- MSR -----------
# name = "preprocessing/preprocess_out/Ward_msr90m"
# template = "preprocessing/template/carb_msr.template"

# ----- netcdf -----
name = "preprocessing/preprocess_out/Ward_msr90m_nc"
template = "preprocessing/template/carb_msr_nc.template"
asprules = "preprocessing/rules/LPC_90m.rules"


RHESSysPreprocess(template = template,
                  name = name,
                  type = type,
                  typepars = typepars,
                  streams = streams,
                  overwrite = overwrite,
                  asprules = asprules,
                  unique_strata_ID = unique_strata_ID,
                  seq_patch_IDs = seq_patch_IDs)

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
               type = type,
               typepars = typepars,
               asprules = asprules,
               streams = streams,
               overwrite = overwrite)


