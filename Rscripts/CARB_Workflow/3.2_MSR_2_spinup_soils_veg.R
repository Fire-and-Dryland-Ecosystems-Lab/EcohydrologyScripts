
library(RHESSysIOinR)
library(rhutils)

# WARD CREEK MSR RHESSYS SOIL SPINUP + VEG GROW

resetveg = F
if (resetveg) {
  world_path = "worldfiles/Ward_msr90m_spun.world"
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
}


# -------------------- Project/Run Name --------------------
name = "Ward_msr_spingrow"
# -------------------- Input RHESSys --------------------
#clim = "clim/netcdf"
clim = "clim/ward"
 
# dates = read_clim(clim,dates_out = T)
# clim_repeat(clim, "clim/ward_1000y", 1000, "years")

#dates = c("1979 1 1 1", "1999 9 30 24")
dates = c("1979 1 1 1", "2079 9 30 24")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.4",
  tec_file = "tecfiles/ward.tec",
  world_file = "worldfiles/Ward_msr90m_soilspin.world",
  world_hdr_prefix = name,
  flowtable = "flowtables/Ward_msr90m.flow",
  start = dates[1],
  end = dates[2],
  output_folder = "output/",
  output_prefix = name,
  commandline_options = "-g -vmort_off -climrepeat -msr"
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

load("r_obj/cal_defpars.rdata")
input_def_pars = IOin_def_pars_simple(input_def_pars)

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = T)

# -------------------- Output filter --------------------
source("../R/output_aliases.R")

outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = "msrgrow_basin",
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
beepr::beep(3)

plotpdf_allvars(out_dir, "msr_grow")

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




