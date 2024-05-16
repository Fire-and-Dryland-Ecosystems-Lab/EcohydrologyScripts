
library(RHESSysIOinR)
library(rhutils)

# WARD CREEK MSR RHESSYS RUN

runs = c("baseline", "thin2", "thin3", "thin4", "thin5", "thin6")
worldfiles = c("worldfiles/Ward_msr90m_soilvegspun.world",
               "worldfiles/Ward_msr90m_thin2_init.world",
               "worldfiles/Ward_msr90m_thin3_init.world",
               "worldfiles/Ward_msr90m_thin4_init.world",
               "worldfiles/Ward_msr90m_thin5_init.world",
               "worldfiles/Ward_msr90m_thin6_init.world")

flowtables = c("flowtables/Ward_msr90m.flow",
               "flowtables/Ward_msr90m_thin2.flow",
               "flowtables/Ward_msr90m_thin3.flow",
               "flowtables/Ward_msr90m_thin4.flow",
               "flowtables/Ward_msr90m_thin5.flow",
               "flowtables/Ward_msr90m_thin6.flow")

i=2

# -------------------- Project/Run Name --------------------
name = paste0("Ward_msr_test_", runs[i])
# -------------------- Input RHESSys --------------------
clim = "clim/netcdf"
clim = "clim/ward_netcdfgridmet_agg"

#dates = read_clim(clim,dates_out = T)
# clim_repeat(clim, "clim/ward_1000y", 1000, "years")

#dates = c("1979 1 1 1", "2020 9 30 24")
dates = c("2000 1 1 1", "2001 1 1 1")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.4",
  tec_file = paste0("tecfiles/ward",name,".tec"),
  world_file = worldfiles[i],
  world_hdr_prefix = paste0("hdr_",name),
  flowtable = flowtables[i],
  start = dates[1],
  end = dates[2],
  output_folder = "output/",
  output_prefix = name,
  commandline_options = "-g -vmort_off -msr"
)

# -------------------- Input Headers --------------------
input_hdr = IOin_hdr(
  basin = "defs/basin.def",
  hillslope = "defs/hill.def",
  zone = "defs/zone_Ward.def",
  soil = c("defs/soil_sandy-loam_Ward.def", "defs/soil_loam_Ward.def"),
  landuse = "defs/lu_Ward.def",
  stratum = c("defs/stratum_evergreen_Ward.def", "defs/stratum_evergreen_Ward_thin.def", 
              "defs/stratum_shrub_Ward.def", "defs/stratum_shrub_Ward_thin.def", 
              "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  basestations = paste0(clim, ".base")
)


# --------------------  Def File Parameters --------------------

# load("r_obj/cal_defpars.rdata")

#input_def_pars = IOin_def_pars_simple(input_def_pars)
input_def_pars = NULL

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "2000 10 1 1", end = dates[2], output_state = F)
if (i >= 2) {
  source_redefs = list.files(path = "worldfiles/redefs",pattern = runs[i],full.names = T)
  ct = 0
  if (any(grepl("litter", source_redefs))) {
    ct = ct + 1
    redef_date = tec_copy_redef(input_redef = source_redefs[grepl("litter", source_redefs)],
                                redef_date = paste0("2005 3 1 ",ct), worldfile = input_rhessys$world_file)
    input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world"))
  }
  # two harvests possible
  if (any(grepl("harvest", source_redefs))) {
    ct = ct + 1
    redef_date = tec_copy_redef(input_redef = source_redefs[grepl("harvest", source_redefs)][1],
                                redef_date = paste0("2005 3 1 ",ct), worldfile = input_rhessys$world_file)
    input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world_thin_harvest"))
  }
  if (sum(grepl("harvest", source_redefs))== 2) {
    ct = ct + 1
    redef_date = tec_copy_redef(input_redef = source_redefs[grepl("harvest", source_redefs)][2],
                                redef_date = paste0("2005 3 1 ",ct), worldfile = input_rhessys$world_file)
    input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world_thin_harvest"))
  }
  if (any(grepl("remain", source_redefs))) {
    ct = ct + 1
    redef_date = tec_copy_redef(input_redef = source_redefs[grepl("remain", source_redefs)],
                                redef_date = paste0("2005 3 1 ",ct), worldfile = input_rhessys$world_file)
    input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world_thin_remain"))
  }
  input_tec_data = input_tec_data[with(input_tec_data, order(year, month, day, hour)), ]
}

# -------------------- Output filter --------------------
source("../R/output_aliases.R")

outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = paste0(name,"_basin"),
  spatial_level = "basin",
  spatial_ID = 1,
  variables = "patch.streamflow"
)
# outfilter2 = build_output_filter(
#   timestep = "daily",
#   output_format = "csv",
#   output_path = "output",
#   output_filename = paste0(name,"_hillslope"),
#   spatial_level = "hill",
#   spatial_ID = 1,
#   variables = c("patch.streamflow", "patch.rz_storage", "patch.unsat_storage", "patch.evaporation", "patch.evaporation_surf","patch.transpiration_unsat_zone", "patch.transpiration_sat_zone")
# )
# outfilter3 = build_output_filter(
#   timestep = "daily",
#   output_format = "csv",
#   output_path = "output",
#   output_filename = paste0(name,"_stratum"),
#   spatial_level = "stratum",
#   spatial_ID = 1,
#   variables = c("stratum.cs.totalc", "stratum.cs.cpool", "stratum.cs.leafc", "stratum.cs.live_stemc", "stratum.cs.dead_stemc", "stratum.cs.live_crootc", "stratum.cs.dead_crootc", "stratum.cs.frootc")
# )
output_filter = IOin_output_filters(outfilter, file_name = paste0("./output/filters/filter",name,".yml") )

# -------------------- Run --------------------

run_rhessys_single(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr,
  def_pars = input_def_pars,
  tec_data = input_tec_data,
  output_filter = output_filter,
  return_cmd = F,
  write_run_metadata = F
)

out_dir = collect_output()

# -------------------- WATER BALANCE --------------------
library(tidyverse)
library(cowplot)

DT = get_basin_daily(out_dir)
bdwatbal = watbal_basin_of(bd = DT)

theme_set(theme_cowplot())

water_balance_daily = bdwatbal %>%
  ggplot() + aes(x = date, y = watbal) + geom_line() + ggtitle("Daily Water Balance") + labs(caption = "For checking data only")

water_balance_annual = bdwatbal %>% group_by(year) %>% summarize(watbal = sum(watbal)) %>%
  ggplot() + aes(x = year, y = watbal) + geom_line() + ggtitle("Annual Water Balance") + labs(caption = "For checking data only")

water_balance_annual

# lapply(cmds, system)
# 
#system(cmds[1])
# 
# 
# 
# out_dir = collect_output()
# beepr::beep(3)

#plotpdf_allvars(out_dir, "msr_run")

