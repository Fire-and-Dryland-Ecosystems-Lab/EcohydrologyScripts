
library(RHESSysIOinR)
library(rhutils)

# WARD CREEK STANDARD RHESSYS TEST RUNS

# -------------------- Project/Run Name --------------------
name = "Ward_std_nc"
# -------------------- Input RHESSys --------------------
clim = "clim/netcdf_edited"

#dates = read_clim(clim,dates_out = T)
# clim_repeat(clim, "clim/ward_500y", 500, "years")

# dates = c("1979 10 1 01", "2016 9 30 24")
dates = c("1979 9 1 1", "2000 9 30 24")

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
  commandline_options = "-g -netcdfgrid"
)
# 
# -------------------- Input Headers --------------------
input_hdr = IOin_hdr(
  basin = "defs/basin.def",
  hillslope = "defs/hill.def",
  zone = "defs/zone.def",
  soil = c("defs/soil_sandy-loam.def", "defs/soil_loam.def"),
  landuse = "defs/lu.def",
  stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  basestations = paste0(clim, ".base")
)

# --------------------  Def File Parameters --------------------
# input_def_pars = IOin_def_pars_simple(input_def_pars)
input_def_pars = NULL

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = F)

# -------------------- Output filter --------------------
source("~/CARB/R/output_aliases.R")

outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = "testrun_basin",
  spatial_level = "basin",
  spatial_ID = 1,
  variables = c(output_vars_minimal)
)

output_filter = IOin_output_filters(outfilter, file_name = "./output/filters/filter.yml")

# -------------------- Run --------------------

# file.copy("~/Repos/RHESSys-develop/rhessys/RHESSys7.4","../bin/rhessys7.4netcdf", overwrite = T)

start = Sys.time()
cmd = run_rhessys_single(
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


cmd_out = system(cmd, intern = T)

writeLines(cmd_out, "output/cmdout.txt")



out_dir = collect_output()
# collect_params(dir = "./")
# out_dir = collect_csvs(dir = "output", dir_base = "spinsoils")

out_dir = "~/CARB/Ward/output/rh_out_2022-08-28--14-14-40/"
plotpdf_allvars(out_dir, "spinsoils")

out_dir = "~/Repos/Multiscale_Methods/output/thincompare2022-08-16_10.02.14/"


check_worldfile("worldfiles/Ward90m_nc.world")

world = read_world("worldfiles/Ward90m_nc.world")

patchsizes = sum(as.numeric(world$values[world$vars == "area" & world$level=="patch"]))
zonesizes = sum(as.numeric(world$values[world$vars == "area" & world$level=="zone"]))

# patch areas do not sum to zone area 1081668.518000 for zone 1
izone = which(world$vars == "zone_ID")
index_max = c(izone[2:length(izone)] - 1, length(world$vars))
world$zone = NA
world$zone[izone[1]:length(world$vars)] = unname(unlist(mapply(rep, world$values[izone], (index_max - izone) + 1 )))
world$zone[world$level == "hillslope"] = NA

ihill = which(world$vars == "hillslope_ID")
ihill_max = c(ihill[2:length(ihill)] - 1, length(world$vars))
world$hillslope = NA
world$hillslope[ihill[1]:length(world$vars)] = unname(unlist(mapply(rep, world$values[ihill], (ihill_max - ihill) + 1 )))

zone1size = sum(as.numeric(world$values[world$vars == "area" & world$level=="zone" & world$zone == 1 & world$hillslope == 18]))

pzonesize = sum(as.numeric(world$values[world$vars == "area" & world$level=="patch" & world$zone == 1 & world$hillslope == 18]))
