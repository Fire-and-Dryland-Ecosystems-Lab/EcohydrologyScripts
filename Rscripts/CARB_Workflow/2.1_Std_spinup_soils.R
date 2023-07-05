
library(RHESSysIOinR)
library(rhutils)

world_name = function(input_rhessys, addstr = NULL) {
  endyr = regmatches(input_rhessys$end_date, regexpr("^\\d{4}\\>" ,input_rhessys$end_date))
  startyr = regmatches(input_rhessys$start_date, regexpr("^\\d{4}\\>" ,input_rhessys$start_date))
  yrlen = as.numeric(endyr) - as.numeric(startyr)
  name = gsub(".world", paste0("_", yrlen,"yr", addstr,".world"), input_rhessys$world_file)
}

# WARD CREEK STANDARD RHESSYS SOIL SPINUP

# -------------------- Project/Run Name --------------------
name = "Ward_std_soilspinup"
# -------------------- Input RHESSys --------------------
clim = "clim/netcdf"
 
# dates = read_clim(clim,dates_out = T)
# clim_repeat(clim, "clim/ward_1000y", 1000, "years")

#dates = c("1979 1 1 1", "2020 9 30 24")
dates = c("1979 1 1 1", "2979 9 30 24")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.4netcdf",
  tec_file = "tecfiles/ward.tec",
  world_file = "worldfiles/Ward90m_nc.world",
  world_hdr_prefix = name,
  flowtable = "flowtables/Ward90m_nc.flow",
  start = dates[1],
  end = dates[2],
  output_folder = "output/",
  output_prefix = name,
  commandline_options = "-g -netcdfgrid -climrepeat"
)

# -------------------- Input Headers --------------------
input_hdr = IOin_hdr(
  basin = "defs/basin.def",
  hillslope = "defs/hill.def",
  zone = "defs/zone.def",
  soil = c("defs/soil_sandy-loam.def", "defs/soil_loam.def"),
  landuse = "defs/lu.def",
  # stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  basestations = paste0(clim, ".base")
)

# --------------------  Def File Parameters --------------------
pars_list = list(
  # EVERGREEN
  list("defs/stratum_evergreen.def", "epc.allocation_flag", "combined"),
  list("defs/stratum_evergreen.def", "epc.froot_turnover", 0.3),
  list("defs/stratum_evergreen.def", "epc.livewood_turnover", 0.013), # was 0.01
  list("defs/stratum_evergreen.def", "epc.alloc_livewoodc_woodc", 0.41),
  list("defs/stratum_evergreen.def", "epc.alloc_stemc_leafc", 0.8),
  # SHRUB
  list("defs/stratum_shrub.def", "epc.allocation_flag", "combined"),
  list("defs/stratum_shrub.def", "epc.livewood_turnover", 0.27),
  list("defs/stratum_shrub.def", "epc.leaf_turnover", 0.33),
  list("defs/stratum_shrub.def", "epc.froot_turnover", 0.45), # c(0.27, 0.6)
  list("defs/stratum_shrub.def", "epc.alloc_livewoodc_woodc", 0.95)
)
input_def_pars = IOin_def_pars_simple(pars_list)

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = T)

# -------------------- Output filter --------------------
source("../R/output_aliases.R")

outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = "soilspin_basin",
  spatial_level = "basin",
  spatial_ID = 1,
  variables = output_vars_minimal
)
output_filter = IOin_output_filters(outfilter, file_name = "./output/filters/filter.yml")

# -------------------- Run --------------------
start = Sys.time()
run_rhessys_single(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr,
  def_pars = input_def_pars,
  tec_data = input_tec_data,
  output_filter = output_filter,
  return_cmd = F,
  write_run_metadata = T
)
end = Sys.time()
end - start

out_dir = collect_output()
# collect_params(dir = "./")
# out_dir = collect_csvs(dir = "output", dir_base = "spinsoils")

plotpdf_allvars(out_dir, "spinsoils")

statefile = worldstate(input_rhessys$world_file)
name = world_name(input_rhessys, "_soilspin")

file.rename(statefile, name)


# -------------------- Reset --------------------

world_path = name
world = read_world(worldfile = world_path)

veg_vars =
  c(
    "cs.cpool", "cs.leafc", "cs.dead_leafc", "cs.leafc_store", "cs.leafc_transfer", "cs.live_stemc", "cs.livestemc_store", "cs.livestemc_transfer", "cs.dead_stemc",
    "cs.deadstemc_store", "cs.deadstemc_transfer", "cs.live_crootc", "cs.livecrootc_store", "cs.livecrootc_transfer", "cs.dead_crootc", "cs.deadcrootc_store", 
    "cs.deadcrootc_transfer", "cs.frootc", "cs.frootc_store", "cs.frootc_transfer", "cs.cwdc","epv.prev_leafcalloc", "ns.npool", "ns.leafn", "ns.dead_leafn", "ns.leafn_store", 
    "ns.leafn_transfer", "ns.live_stemn", "ns.livestemn_store", "ns.livestemn_transfer", "ns.dead_stemn", "ns.deadstemn_store", "ns.deadstemn_transfer", 
    "ns.live_crootn", "ns.livecrootn_store", "ns.livecrootn_transfer", "ns.dead_crootn", "ns.deadcrootn_store", "ns.deadcrootn_transfer", "ns.frootn", 
    "ns.frootn_store", "ns.frootn_transfer", "ns.cwdn", "ns.retransn"
  )

world$values[world$vars %in% veg_vars] = "0.0"

write_world(world, gsub(".world", "_reset.world", world_path))




