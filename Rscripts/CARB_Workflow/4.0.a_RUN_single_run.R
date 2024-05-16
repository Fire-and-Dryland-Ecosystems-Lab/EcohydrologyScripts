
library(RHESSysIOinR)
library(rhutils)
source("R/0_global_vars.R")
source("../R/output_aliases.R")
options(scipen = 999)

resetveg = F
if (resetveg) {
  # world_path = "worldfiles/Pitman_MSR_soilspin_200y.world"
  world = read_world(worldfile = world_path)
  veg_vars =
    c("cs.cpool", "cs.leafc", "cs.dead_leafc", "cs.leafc_store", "cs.leafc_transfer",
      "cs.live_stemc", "cs.livestemc_store", "cs.livestemc_transfer", "cs.dead_stemc",
      "cs.deadstemc_store", "cs.deadstemc_transfer", "cs.live_crootc", "cs.livecrootc_store",
      "cs.livecrootc_transfer", "cs.dead_crootc", "cs.deadcrootc_store", "cs.deadcrootc_transfer",
      "cs.frootc", "cs.frootc_store", "cs.frootc_transfer", "cs.cwdc", 
      "epv.prev_leafcalloc", "ns.npool", "ns.leafn", "ns.dead_leafn",
      "ns.leafn_store", "ns.leafn_transfer", "ns.live_stemn", "ns.livestemn_store",
      "ns.livestemn_transfer", "ns.dead_stemn", "ns.deadstemn_store", "ns.deadstemn_transfer",
      "ns.live_crootn", "ns.livecrootn_store", "ns.livecrootn_transfer", "ns.dead_crootn",
      "ns.deadcrootn_store", "ns.deadcrootn_transfer", "ns.frootn", "ns.frootn_store",
      "ns.frootn_transfer", "ns.cwdn", "ns.retransn")
  world$values[world$vars %in% veg_vars] = "0.0"
  
  # insert rows for spinup object ID
  world$ind = 1:nrow(world)
  tmp = world[world$vars == "veg_parm_ID",]
  tmp$vars = "spinup_object_ID"
  tmp$values = 1
  tmp$ind = tmp$ind + 0.5
  world = rbind(world, tmp)
  world = world[order(world$ind),]
  world$ind = NULL
  
  write_world(world, gsub(".world", "_reset.world", world_path))
}

resetsnow = F
if (resetsnow) {
  world_path = scenario_df$worldfiles[scenario_df$runs == "baseline"][1]
  world = read_world(worldfile = world_path)
  snow_vars = c("snowpack.water_equivalent_depth", "snowpack.water_depth", "snowpack.T", 
                "snowpack.surface_age", "snowsnowpack.energy_deficit","snow_stored")
  world$values[world$vars %in% snow_vars] = "0.0"
  write_world(world, gsub(".world", "_snowreset.world", world_path))
  snowworld = gsub(".world", "_snowreset.world", world_path)
}

# tmp = read_world("worldfiles/Ward_msr90m_baseline_init.world")
# snow = as.numeric(tmp[tmp$vars == "snowpack.water_equivalent_depth","values"])


# MSR RHESSYS VEG GROW
# -------------------- Project/Run Name --------------------
script = "test_single_run"
name = paste0(site,"_",script)

# -------------------- Input RHESSys --------------------
#clim = "clim/netcdf"
clim = "clim/ward_netcdfgridmet_agg"
 
# dates = read_clim(clim,dates_out = T)
# clim_repeat(clim, "clim/ward_1000y", 1000, "years")

#dates = c("1979 1 1 1", "1999 9 30 24")
# dates = c("1979 1 1 1", "2019 9 30 24")
dates = c("1979 1 1 1", "1999 9 30 24")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.4",
  tec_file = paste0("tecfiles/",site,".tec"),
  # world_file = scenario_df$worldfiles[scenario_df$runs == "baseline"][1],
  world_file = snowworld,
  world_hdr_prefix = paste0("hdr_",name),
  flowtable = scenario_df$flowtables[scenario_df$runs == "baseline"][1],
  start = dates[1],
  end = dates[2],
  output_folder = "output/",
  output_prefix = name,
  commandline_options = "-g -vmort_off -climrepeat -msr"
  # -vegspinup ./defs/spinup_LAI_targets.txt 
)

# -------------------- Input Headers --------------------
# input_hdr = IOin_hdr(
#   basin = "defs/basin.def",
#   hillslope = "defs/hill.def",
#   zone = "defs/zone.def",
#   soil = c("defs/soil_sandy-loam.def", "defs/soil_loam.def"),
#   landuse = "defs/lu.def",
#   stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
#   basestations = paste0(clim, ".base")
# )

input_hdr = IOin_hdr(
  basin = "defs/basin.def", hillslope = "defs/hill.def", zone = paste0("defs/zone_",site,".def"),
  soil = c(paste0("defs/soil_sandy-loam_",site,".def"), paste0("defs/soil_loam_",site,".def")),
  landuse = paste0("defs/lu_",site,".def"),
  stratum = c("defs/stratum_evergreen.def", 
              paste0("defs/stratum_shrub.def"), 
              "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  # spinup = "defs/spinup.def",
  basestations = paste0(clim, ".base")
)
# names(input_hdr)[names(input_hdr) == "spinup"] = "spinup_def"

# --------------------  Def File Parameters --------------------
pars_list = list(
  # SHRUB
  # list(input_hdr$stratum_def[2], "epc.livewood_turnover", 0.2),
  list(input_hdr$stratum_def[2], "epc.leaf_turnover", 0.26),
  # list(input_hdr$stratum_def[2], "epc.froot_turnover", 0.45), # c(0.27, 0.6)
  # list(input_hdr$stratum_def[2], "epc.branch_turnover", 0.03) # 0.2 p301
  
  # EVERGREEN
  # list(input_hdr$stratum_def[1], "epc.froot_turnover", 0.3),
  # list(input_hdr$stratum_def[1], "epc.livewood_turnover", 0.01), # was 0.01
  # list(input_hdr$stratum_def[1], "epc.branch_turnover", 0.03), # 0.01 p301
  list(input_hdr$stratum_def[1], "epc.leaf_turnover", 0.22),
  list(input_hdr$stratum_def[1], "epc.froot_turnover", 0.26),
  list(input_hdr$stratum_def[1], "epc.alloc_stemc_leafc", c(0.2)),
  list(input_hdr$stratum_def[1], "epc.alloc_frootc_leafc", c(1.4)),
  list(input_hdr$stratum_def[1], "epc.alloc_crootc_stemc", c(0.5)),
  list(input_hdr$stratum_def[1], "epc.allocation_flag", c("waring")),
  # ----- soils -----
  list(input_hdr$soil_def[1], "m", 0.22 ), 
  list(input_hdr$soil_def[1], "m_v", c(0.22)), 
  list(input_hdr$soil_def[1], "Ksat_0", 264),
  list(input_hdr$soil_def[1], "Ksat_0_v", 264) 
)
input_def_pars = IOin_def_pars_simple(pars_list)
# input_def_pars = NULL

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = F)

# -------------------- Output filter --------------------
# output_carb from aliases
outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = paste0(name,"_basin"),
  spatial_level = "basin",
  spatial_ID = 1,
  variables = output_carb
)
output_filter = IOin_output_filters(outfilter, file_name = "./output/filters/filter.yml")

# -------------------- Run --------------------
run_rhessys_single(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr,
  def_pars = input_def_pars,
  tec_data = input_tec_data,
  output_filter = output_filter,
  return_cmd = F,
  write_run_metadata = T
)
out_dir = collect_output()


plotpdf_allvars(out_dir, name)

df = get_basin_daily(out_dir)
df = watbal_basin_of(df)
summary(df$watbal)
plot(df$date,df$watbal,type="l")

# statefile = worldstate(input_rhessys$world_file)
# name = world_name(input_rhessys, "_soilspin")
# file.rename(statefile, name)


# -------------------- Reset --------------------
# 
# world_path = name
# world = read_world(worldfile = world_path)
# 
# veg_vars =
#   c(
#     "cs.cpool", "cs.leafc", "cs.dead_leafc", "cs.leafc_store", "cs.leafc_transfer", "cs.live_stemc", "cs.livestemc_store", "cs.livestemc_transfer", "cs.dead_stemc",
#     "cs.deadstemc_store", "cs.deadstemc_transfer", "cs.live_crootc", "cs.livecrootc_store", "cs.livecrootc_transfer", "cs.dead_crootc", "cs.deadcrootc_store", 
#     "cs.deadcrootc_transfer", "cs.frootc", "cs.frootc_store", "cs.frootc_transfer", "cs.cwdc","epv.prev_leafcalloc", "ns.npool", "ns.leafn", "ns.dead_leafn", "ns.leafn_store", 
#     "ns.leafn_transfer", "ns.live_stemn", "ns.livestemn_store", "ns.livestemn_transfer", "ns.dead_stemn", "ns.deadstemn_store", "ns.deadstemn_transfer", 
#     "ns.live_crootn", "ns.livecrootn_store", "ns.livecrootn_transfer", "ns.dead_crootn", "ns.deadcrootn_store", "ns.deadcrootn_transfer", "ns.frootn", 
#     "ns.frootn_store", "ns.frootn_transfer", "ns.cwdn", "ns.retransn"
#   )
# 
# world$values[world$vars %in% veg_vars] = "0.0"
# 
# write_world(world, gsub(".world", "_reset.world", world_path))




