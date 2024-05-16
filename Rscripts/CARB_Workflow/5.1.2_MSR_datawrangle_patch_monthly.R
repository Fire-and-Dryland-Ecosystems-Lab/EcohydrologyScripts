# -------------------- data wrangle --------------------
library(RHESSysIOinR)
library(rhutils)
source("R/0_global_vars.R")
source("../R/fun_data_aggregation.R")

# -------------------------------------------------------------
# -------------------- Input Data --------------------
# -------------------------------------------------------------
# out_dir = "output/rh_out_2023-08-25--07-32-12/"
out_dir = "output/rh_out_2024-02-08--11-02-30/"
level = "patch"
timestep = "monthly"
source_name = paste0(level,"_",timestep)
dest_dir = paste0("r_obj/",source_name)
if (!file.exists(dest_dir)) {dir.create(dest_dir)}

# -------------------------------------------------------
# -------------------- CSV to .rdata --------------------
# -------------------------------------------------------
files_in = list.files(path = out_dir, pattern = paste0(".*_",level,".*\\.csv$"), full.names = T)
files_out = file.path(dest_dir, paste0(scenario_df$name,"_",source_name,".rda"))

tmpsummary = list()
for (i in seq_along(files_in)) {
  DT_patch_monthly = fread(files_in[i])
  tmpsummary[[i]] = summary(DT_patch_monthly) # check the data
  world = read_world(scenario_df$worldfiles[i])
  DT_patch_monthly = patch_add_world_vars(DT_patch_monthly, world)
  
  cat(paste0("Saving ",paste0(scenario_df$name[i],"_",source_name,".rda")," to ",dest_dir,"\n"))
  save(DT_patch_monthly, file =  files_out[i])
  gc()
}

# --------------------------------------------------------------
# -------------------- Patch to Hillslope +Agg to single files  --------------------
# --------------------------------------------------------------
patch2hill_agg(patch_rdata = files_out, output_rdata = file.path(dest_dir, paste0("allruns_hill_monthly.rda")), site = site)

# # --------------------------------------------------------------
# # -------------------- Agg to annual  --------------------
# # --------------------------------------------------------------
# 
# keepvars = c("day", "month", "year", "basinID", "hillID", "zoneID", "patchID","familyID", "stratumID", "family_area", "area_prop", "veg_parm_ID", "area", "wy","yr_mn","Scenario")
# 
# aggvars = names(patchfamily_mn_DT)[!names(patchfamily_mn_DT) %in%  keepvars]
# patchfamily_wateryr_DT = patchfamily_mn_DT[, lapply(.SD, mean), by=c("Scenario", "wy", "basinID", "hillID", "zoneID","familyID"), .SDcols = aggvars]
# 
# aggvars = names(hill_mn_DT)[!names(hill_mn_DT) %in%  keepvars]
# hill_wateryr_DT = hill_mn_DT[, lapply(.SD, mean), by=c("Scenario", "wy", "basinID", "hillID"), .SDcols = aggvars]
# 
# aggvars = names(hillvegID_mn_DT)[!names(hillvegID_mn_DT) %in%  keepvars]
# hillvegID_wateryr_DT = hillvegID_mn_DT[, lapply(.SD, mean), by=c("Scenario", "wy", "basinID", "hillID"), .SDcols = aggvars]
# 
# save(patchfamily_wateryr_DT, file = file.path(dest_dir, paste0("allruns_patchfamily_wateryear.rda")))
# save(hill_wateryr_DT, file = file.path(dest_dir, paste0("allruns_hill_wateryear.rda")))
# save(hillvegID_wateryr_DT, file = file.path(dest_dir, paste0("allruns_hill_wateryear_vegID.rda")))
# 
# 







