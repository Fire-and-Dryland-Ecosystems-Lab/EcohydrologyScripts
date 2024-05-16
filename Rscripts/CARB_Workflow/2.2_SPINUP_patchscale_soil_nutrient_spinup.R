# MSR_1_soil_spin
library(RHESSysIOinR)
library(rhutils)
library(data.table)
source("../R/fun_spinup.R")

# Uses a worldfile to generate single patch worlds based on each veg parameter, then runs each world and reincorporates the
# resulting soil nutrient values into the original worldfile

# ---------- Input world ----------
world_path = "worldfiles/Ward_msr90m.world"

# read world
world = as.data.table(read_world(world_path))

# add cols
world = world_add_level_i(world)
world = world_add_patch_vegparmIDcol(world)

veg = unique(world$vegparm[!is.na(world$vegparm)])
output_world_paths = paste0("worldfiles/Soil_spin_worlds/Ward_msr90m_vegID",veg,".world")
output_flow_paths = paste0("flowtables/Soil_spin_flowtables/Ward_msr90m_vegID",veg,".flow")

# ---------- Make 1 Patch Worlds ----------
# selecting patches based on single veg parm IDs, since world is 1 strata
p1_veg = rbindlist(lapply(veg, function(x) world[min(which(world$vegparm == x)), ]))
# make new worlds
newworlds = lapply(p1_veg$unique_ID, FUN = extract_world, world = world)

mapply(write_world,newworlds,output_world_paths)

# ---------- Make 1 Patch Flow Tables ----------

# ASSUMES 1 HILLSLOPE
make_1pflow = function(world) {
  hill = unique(world$ID[world$level == "hillslope"])
  patch = unique(world$pID[!is.na(world$pID)])
  zone = unique(world$ID[world$level == "zone"])
  xcen = world$values[world$level == "patch" & world$vars == "x"]
  ycen = world$values[world$level == "patch" & world$vars == "y"]
  zcen = world$values[world$level == "patch" & world$vars == "z"]
  area = world$values[world$level == "patch" & world$vars == "area"]
  flow_out = paste0(length(hill),"\n",
                    hill, "\t", length(patch),"\n",
                    patch, "\t", zone, "\t", hill, "\t", xcen, "\t", ycen, "\t", zcen, "\t", area, "\t", area, "\t", 1, "\t", area, "\t", 0)
  return(flow_out)
}

newflows = lapply(newworlds, make_1pflow)
mapply(writeLines, newflows, output_flow_paths)

# -------------------- RUN FOR SOIL NUTRIENT SPINUP --------------------

run_cmds = list()

for (i in seq_along(veg)) {
  # -------------------- Project/Run Name --------------------
  name = paste0("Ward_msr90m_vegID",veg[i])
  # -------------------- Input RHESSys --------------------
  clim = "clim/ward"
  
  #dates = c("1979 1 1 1", "2020 9 30 24")
  dates = c("1979 1 1 1", "2750 9 30 24")
  
  input_rhessys = IOin_rhessys_input(
    version = "../bin/rhessys7.4",
    tec_file = "tecfiles/ward.tec",
    world_file = output_world_paths[i],
    world_hdr_prefix = paste0("hdr_",name),
    flowtable = output_flow_paths[i],
    start = dates[1],
    end = dates[2],
    output_folder = "output/",
    output_prefix = name,
    commandline_options = "-g -vmort_off -climrepeat"
  )
  
  # -------------------- Input Headers --------------------
  # input_hdr = IOin_hdr(
  #   basin = "defs/basin.def",
  #   hillslope = "defs/hill.def",
  #   zone = "defs/zone.def",
  #   soil = c("defs/soil_sandy-loam.def", "defs/soil_loam.def"),
  #   landuse = "defs/lu.def",
  #   # stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  #   stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  #   basestations = paste0(clim, ".base")
  # )
  
  input_hdr = IOin_hdr(
    basin = "defs/basin.def",
    hillslope = "defs/hill.def",
    zone = "defs/zone_Ward.def",
    soil = c("defs/soil_sandy-loam_Ward.def", "defs/soil_loam_Ward.def"),
    landuse = "defs/lu_Ward.def",
    stratum = c("defs/stratum_evergreen_Ward.def", "defs/stratum_shrub_Ward.def", 
                "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
    basestations = paste0(clim, ".base")
  )
  
  # --------------------  Def File Parameters --------------------
  # pars_list = list(
  #   # EVERGREEN
  #   list("defs/stratum_evergreen.def", "epc.allocation_flag", "combined"),
  #   list("defs/stratum_evergreen.def", "epc.froot_turnover", 0.3),
  #   list("defs/stratum_evergreen.def", "epc.livewood_turnover", 0.013), # was 0.01
  #   list("defs/stratum_evergreen.def", "epc.alloc_livewoodc_woodc", 0.41),
#   list("defs/stratum_evergreen.def", "epc.alloc_stemc_leafc", 0.8),
  #   # SHRUB
  #   list("defs/stratum_shrub.def", "epc.allocation_flag", "combined"),
  #   list("defs/stratum_shrub.def", "epc.livewood_turnover", 0.27),
  #   list("defs/stratum_shrub.def", "epc.leaf_turnover", 0.33),
  #   list("defs/stratum_shrub.def", "epc.froot_turnover", 0.45), # c(0.27, 0.6)
  #   list("defs/stratum_shrub.def", "epc.alloc_livewoodc_woodc", 0.95)
  # )
  # input_def_pars = IOin_def_pars_simple(pars_list)
  input_def_pars = NULL
  
  # -------------------- Make Tec File --------------------
  input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = T)
  
  # -------------------- Output filter --------------------
  source("../R/output_aliases.R")
  
  outfilter = build_output_filter(
    timestep = "daily",
    output_format = "csv",
    output_path = "output",
    output_filename = paste0("soilspin_basin_veg",veg[i]),
    spatial_level = "basin",
    spatial_ID = 1,
    variables = output_vars_minimal
  )
  output_filter = IOin_output_filters(outfilter, file_name = paste0("./output/filters/filter_veg",veg[i],".yml"))
  
  run_cmds[[i]] = run_rhessys_single(
    input_rhessys = input_rhessys,
    hdr_files = input_hdr,
    def_pars = input_def_pars,
    tec_data = input_tec_data,
    output_filter = output_filter,
    return_cmd = T,
    write_run_metadata = T
  )
  
}

# -------------------- MANUAL PARALLEL RUNS --------------------
library(parallel)

n_cores = detectCores() - 1
# cl = makeCluster(length(run_cmds))
cl = makePSOCKcluster(names = length(run_cmds), port = floor(runif(1,11000,11999)), outfile = "output/clusteroutfile")
clusterExport(cl = cl, varlist = c("run_cmds"), envir = environment())
parLapply(cl = cl, X = seq_along(veg), fun = function(X, Y) { system(Y[[X]])}, Y = run_cmds)
stopCluster(cl)

out_dir = collect_output()

# -------------------- Run --------------------

# run_rhessys_single(
#   input_rhessys = input_rhessys,
#   hdr_files = input_hdr,
#   def_pars = input_def_pars,
#   tec_data = input_tec_data,
#   output_filter = output_filter,
#   return_cmd = F,
#   write_run_metadata = T
# )


# 
# out_dir = collect_output()
# collect_params(dir = "./")
# out_dir = collect_csvs(dir = "output", dir_base = "spinsoils")

plotpdf_allvars(out_dir, "spinsoils")

# statefile = worldstate(input_rhessys$world_file)
# name = world_name(input_rhessys, "_soilspin")
# 
# file.rename(statefile, name)


# -------------------- New World with Spun Soil Nutrients --------------------

# choosing the Y2750 worlds
world_path_spun = list.files("worldfiles/Soil_spin_worlds/","Y2750", full.names = T)

# world_path_spun = c("worldfiles/Soil_spin_worlds/Ward_msr90m_vegID1.world.Y2979M9D30H23.state",
#                "worldfiles/Soil_spin_worlds/Ward_msr90m_vegID3.world.Y2979M9D30H23.state",
#                "worldfiles/Soil_spin_worlds/Ward_msr90m_vegID4.world.Y2979M9D30H23.state",
#                "worldfiles/Soil_spin_worlds/Ward_msr90m_vegID5.world.Y2979M9D30H23.state")

worlds_spun = lapply(world_path_spun, FUN = function(X){as.data.table(read_world(X))})
worlds_spun = lapply(worlds_spun, world_add_level_i)
worlds_spun = lapply(worlds_spun, world_add_patch_vegparmIDcol)

soil_vars = c("soil_cs.soil1c",  "soil_cs.soil2c", "soil_cs.soil3c", 
              "soil_cs.soil4c", "soil_ns.sminn", "soil_ns.nitrate"
              # "litter.rain_stored", "litter_cs.litr1c", "litter_ns.litr1n", 
              # "litter_cs.litr2c", "litter_cs.litr3c", "litter_cs.litr4c", 
              # "rootzone.depth", "snow_stored", "rain_stored", 
              # "epv.wstress_days", "epv.max_fparabs", "epv.min_vwc", 
              # "gw.storage", "gw.NO3"
              )

world_dest = as.data.table(read_world(world_path))
world_dest = world_add_level_i(world_dest)
world_dest = world_add_patch_vegparmIDcol(world_dest)

vegids = as.numeric(gsub("vegID","", regmatches(world_path_spun, regexpr("vegID\\d+", world_path_spun))))

# ew bad
for (vid in vegids) {
  for (v in soil_vars) {
    world_dest[vars %in% soil_vars & vegparm == vid & vars == v, "values"] = 
      worlds_spun[[which(vegids == vid)]][vars %in% soil_vars & vegparm == vid & vars == v, "values"]
  }
}

write_world(world_dest, gsub(".world", "_patchsoilspin.world", world_path))

