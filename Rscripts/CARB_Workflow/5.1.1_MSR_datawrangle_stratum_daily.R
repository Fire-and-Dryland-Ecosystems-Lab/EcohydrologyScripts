# data wrangle
library(RHESSysIOinR)
library(rhutils)
source("R/0_global_vars.R")

# -------------------------------------------------------------
# -------------------- Ingest Stratum Data --------------------
# -------------------------------------------------------------
out_dir = "output/rh_out_2024-04-04--01-22-38/"
strata_dest_dir = "r_obj/strata_daily"

# worldfiles = c("worldfiles/Ward_msr90m_soilvegspun.world", 
#                "worldfiles/Ward_msr90m_thin2_init.world",
#                "worldfiles/Ward_msr90m_thin3_init.world",
#                "worldfiles/Ward_msr90m_thin4_init.world",
#                "worldfiles/Ward_msr90m_thin5_init.world",
#                "worldfiles/Ward_msr90m_thin6_init.world")


# -------------------------------------------------------
# -------------------- CSV to .rdata --------------------
# -------------------------------------------------------
if (!file.exists(strata_dest_dir)) {
  dir.create(strata_dest_dir)
}
strata_files_in = list.files(path = out_dir, pattern = ".*_stratum.*\\.csv$", full.names = T)
run_names = gsub("_stratum", "", gsub(".csv", "", basename(strata_files_in)))

for (i in seq_along(strata_files_in)) {
  stratum_daily_DT = fread(strata_files_in[i])
  save(stratum_daily_DT, file = file.path(strata_dest_dir, paste0(run_names[i],".rda") ))
  gc()
}

# ----------------------------------------------------------
# -------------------- Daily to Monthly --------------------
# ----------------------------------------------------------
strata_daily_R = list.files(path = strata_dest_dir, pattern = ".*_daily", full.names = T)
run_names = gsub("_daily.rda", "", basename(strata_daily_R))

for (i in seq_along(strata_daily_R)) {
  load(strata_daily_R[i])
  aggvars = names(stratum_daily_DT)
  aggvars = aggvars[aggvars!=c("day","month","year","basinID","hillID","zoneID","patchID","stratumID")]
  
  stratum_mn_DT = stratum_daily_DT[,lapply(.SD, mean) , 
                                   by = c("basinID", "hillID", "zoneID", "patchID", "stratumID", "year","month"), .SDcols = aggvars]
  stratum_mn_DT$run = run_names[i]
  stratum_mn_DT = RHESSysIOinR::add_dates(stratum_mn_DT)
  stratum_mn_DT = stratum_mn_DT[order(yr_mn)]
  save(stratum_mn_DT, file = file.path(strata_dest_dir, paste0(run_names[i],"_monthly.rda") ))
  gc()
}
  
# --------------------------------------------------------------
# -------------------- Stratum to Hillslope --------------------
# --------------------------------------------------------------
# strata_monthly_R = list.files(path = strata_dest_dir, pattern = ".*_monthly.rda", full.names = T)
# run_names = gsub("_monthly.rda", "", basename(strata_monthly_R))
out_dir = "output/rh_out_2024-03-26--23-44-37/"
strata_dest_dir = "r_obj/strata_monthly"

stratum_files_in = list.files(path = out_dir, pattern = ".*_stratum.*\\.csv$", full.names = T)

if (length(stratum_files_in) == 0) {
  stratum_files_in = list.files(path = out_dir, pattern = ".*\\.csv$", full.names = T)
}
if (length(stratum_files_in) == 0) {
  cat("No stratum files found at:",out_dir)
  return(NULL)
}

 dt_list = list()
 dt_listvegid = list()
 
# for (i in seq_along(strata_monthly_R)) {
for (i in seq_along(stratum_files_in)) {
  runname = gsub("_stratum", "", gsub(".csv", "", basename(stratum_files_in[i])))
  stratum_mn_DT = fread(stratum_files_in[i] )
  climate = scenario_df$climate[scenario_df$name == runname]
  runs = scenario_df$runs[scenario_df$name == runname]
  treat_type = scenario_df$runs_fullnames[scenario_df$name == runname]
  scenario = paste0(site,"_msr_", runs,"_", climate)
  stratum_mn_DT$Scenario = scenario
  stratum_mn_DT$Treat_type = treat_type
  stratum_mn_DT$Climate = climate
  stratum_mn_DT$Site = "Ward"
  stratum_mn_DT = add_dates(stratum_mn_DT)

  # load(strata_monthly_R[i]) #stratum_mn_DT
  # get areas and join by strata id or patch
  world = read_world(scenario_df$worldfiles[scenario_df$name == runname])
  patchareas = world[world$level == "patch" & world$vars == "area",c("values","ID")]
  names(patchareas) = c("area", "patchID")
  patchareas$patchID = as.numeric(patchareas$patchID)
  patchareas$area = as.numeric(patchareas$area)
  length(unique(patchareas$patchID)) == length(unique(stratum_mn_DT$patchID))
  all(patchareas$patchID %in% stratum_mn_DT$patchID)
  stratum_mn_DT = stratum_mn_DT[patchareas, on = .(patchID)]
  # get family ID
  stratum_mn_DT$familyID = floor(stratum_mn_DT$patchID/100)
  family_area = stratum_mn_DT[,sum(area), by = c("familyID", "year","month")]
  names(family_area)[4] = "family_area"
  stratum_mn_DT = stratum_mn_DT[family_area, on = .(familyID, year, month)]
  stratum_mn_DT$area_prop = stratum_mn_DT$area/stratum_mn_DT$family_area
  # get vegparm per strata
  vegparm = world[world$level == "canopy_strata" & world$vars == "veg_parm_ID",c("values","ID")]
  names(vegparm) = c("veg_parm_ID", "stratumID")
  vegparm$veg_parm_ID = as.numeric(vegparm$veg_parm_ID)
  vegparm$stratumID = as.numeric(vegparm$stratumID)
  stratum_mn_DT = stratum_mn_DT[vegparm, on = .(stratumID)]
  
  # # aggregate to patch family using area weights
  # aggvars = c("cs.totalc","cs.cpool","cs.leafc","cs.live_stemc","cs.dead_stemc" , "cs.live_crootc" ,"cs.dead_crootc", "cs.frootc","totalc","stemc","rootc","leafc")
  vars = c("month", "year", "basinID", "hillID", "zoneID", "patchID", "stratumID", "Scenario", 
              "Treat_type", "Climate", "Site", "wy", "yr_mn", "area", "familyID", "family_area", "area_prop", "veg_parm_ID")
  # aggvars = aggvars[aggvars %in% names(stratum_mn_DT)]
  aggvars = names(stratum_mn_DT)[!names(stratum_mn_DT) %in% vars]
  family_tmp = stratum_mn_DT
  family_tmp = family_tmp[,(aggvars) := .SD*area_prop, .SDcols = aggvars]
  family_mn_DT = family_tmp[,lapply(.SD, sum) ,by = c("Site", "Scenario","Treat_type","Climate","basinID", "hillID", "zoneID", "familyID", "year","month", "wy","yr_mn"), .SDcols = c(aggvars,"area")]
  # save(family_mn_DT, file = file.path(strata_dest_dir, paste0(run_names[i],"_monthly_patchfamily.rda") ))
  # aggregate to hillslope
  hill_mn_DT = family_mn_DT[,lapply(.SD, mean) ,by = c("Site", "Scenario","Treat_type","Climate","basinID", "hillID", "year","month", "wy","yr_mn"), .SDcols = c(aggvars,"area")]
  # save(hill_mn_DT, file = file.path(strata_dest_dir, paste0(run_names[i],"_monthly_hill.rda") ))
  
  # aggregate by vegparm
  hillvegID = stratum_mn_DT[,lapply(.SD, mean) ,by = c("Site", "Scenario","Treat_type","Climate","basinID", "hillID","veg_parm_ID", "year","month", "wy","yr_mn"), .SDcols = c(aggvars,"area")]
  dt_listvegid[[i]] = hillvegID
  dt_list[[i]] = hill_mn_DT
}


 hillvegID_mn_DT = rbindlist(dt_listvegid)
 hill_mn_DT = rbindlist(dt_list)
 
save(hillvegID_mn_DT, file = file.path(strata_dest_dir, "hillvegID_mn.rda" ))
save(hill_mn_DT, file = file.path(strata_dest_dir, "hill_mn.rda" ))
