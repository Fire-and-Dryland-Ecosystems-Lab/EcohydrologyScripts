
library(RHESSysIOinR)
library(rhutils)
source("R/0_global_vars.R")
# WARD CREEK MSR RHESSYS RUN

cmds = vector(length = length(scenario_df$name), mode = "character")
filelist_allruns = list()

for (i in seq_along(scenario_df$name)) {
  # -------------------- Project/Run Name --------------------
  # name = paste0("Ward_msr", runs[i])
  name = scenario_df$name[i]
  # -------------------- Input RHESSys --------------------
  #clim = "clim/netcdf"
  clim = "clim/ward_netcdfgridmet_agg"
  #dates = read_clim(clim,dates_out = T)
  # clim_repeat(clim, "clim/ward_1000y", 1000, "years")
  
  input_rhessys = IOin_rhessys_input(
    version = "../bin/rhessys7.4",
    tec_file = paste0("tecfiles/ward",name,".tec"),
    world_file = scenario_df$worldfiles[i],
    world_hdr_prefix = paste0("hdr_",name),
    flowtable = scenario_df$flowtables[i],
    start = scenario_df$start[i],
    end = scenario_df$end[i],
    output_folder = "output/",
    output_prefix = name,
    commandline_options = "-g -vmort_off"
    # -msr
  )
  
  # -------------------- Input Headers --------------------
  input_hdr = IOin_hdr(
    basin = "defs/basin.def", hillslope = "defs/hill.def", zone = "defs/zone_Ward.def",
    soil = c("defs/soil_sandy-loam_Ward.def", "defs/soil_loam_Ward.def"), landuse = "defs/lu_Ward.def",
    stratum = c("defs/stratum_evergreen_Ward.def", "defs/stratum_evergreen_Ward_thin.def", "defs/stratum_shrub_Ward.def", 
                "defs/stratum_shrub_Ward_thin.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
    basestations = paste0(clim, ".base")
  )
  
  # --------------------  Def File Parameters --------------------
  # load("r_obj/cal_defpars.rdata")
  #input_def_pars = IOin_def_pars_simple(input_def_pars)
  input_def_pars = NULL
  
  # -------------------- Make Tec File --------------------
  input_tec_data = IOin_tec_std(start = scenario_df$start_output[i], end = scenario_df$end[i], output_state = F)
  input_tec_data[3,] = input_tec_data[2,]
  input_tec_data[3, "hour"] =  as.numeric(input_tec_data[2, "hour"]) + 1
  input_tec_data[3, "name"] = "print_monthly_on"
  treat_date = scenario_df$treat_date[i]
  
  if (scenario_df$runs[i] != "baseline") {
    source_redefs = list.files(path = "worldfiles/redefs",pattern = scenario_df$runs[i],full.names = T)
    ct = 0
    if (any(grepl("litter", source_redefs))) {
      ct = ct + 1
      redef_date = tec_copy_redef(input_redef = source_redefs[grepl("litter", source_redefs)],
                                  redef_date = paste0(treat_date," ",ct), worldfile = input_rhessys$world_file)
      input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world_multiplier"))
    }
    # two harvests possible
    if (any(grepl("harvest", source_redefs))) {
      ct = ct + 1
      redef_date = tec_copy_redef(input_redef = source_redefs[grepl("harvest", source_redefs)][1],
                                  redef_date = paste0(treat_date," ",ct), worldfile = input_rhessys$world_file)
      input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world_thin_harvest"))
    }
    if (sum(grepl("harvest", source_redefs))== 2) {
      ct = ct + 1
      redef_date = tec_copy_redef(input_redef = source_redefs[grepl("harvest", source_redefs)][2],
                                  redef_date = paste0(treat_date," ",ct), worldfile = input_rhessys$world_file)
      input_tec_data = rbind(input_tec_data,c(redef_date, "redefine_world_thin_harvest"))
    }
    if (any(grepl("remain", source_redefs))) {
      ct = ct + 1
      redef_date = tec_copy_redef(input_redef = source_redefs[grepl("remain", source_redefs)],
                                  redef_date = paste0(treat_date," ",ct), worldfile = input_rhessys$world_file)
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
    # variables = c("patch.streamflow")
    variables = output_carb
  )
  outfilter2 = build_output_filter(
    timestep = "monthly",
    output_format = "csv",
    output_path = "output",
    output_filename = paste0(name,"_patch"),
    spatial_level = "patch",
    spatial_ID = 1,
    variables = c("patch.streamflow", "patch.et", "patch.trans", "patch.PET","patch.psn" ,"patch.theta", "patch.rz_storage", "patch.unsat_storage", "patch.sat_deficit", "patch.pcp","patch.rain_thru","patch.snow_thru" )
  )
  outfilter3 = build_output_filter(
    timestep = "monthly",
    output_format = "csv",
    output_path = "output",
    output_filename = paste0(name,"_stratum"),
    spatial_level = "stratum",
    spatial_ID = 1,
    variables = c("stratum.totalc", "stratum.stemc", "stratum.rootc","stratum.leafc")
  )
  output_filter = IOin_output_filters(outfilter, outfilter2, outfilter3, file_name = paste0("./output/filters/filter",name,".yml") )
  
  # -------------------- Run --------------------
  
  cmds[i] = run_rhessys_single(
    input_rhessys = input_rhessys,
    hdr_files = input_hdr,
    def_pars = input_def_pars,
    tec_data = input_tec_data,
    output_filter = output_filter,
    return_cmd = T,
    write_run_metadata = F
  )
  filelist_allruns[[i]] = list_rh_input_files(input_rhessys, input_hdr, input_def_pars, input_tec_data, output_filter)
  
}

# -------------------- HPC --------------------
name = "Ward_treat"
rhout_str = rhout_write_for_hpc(cmds, name)

filelist = unlist(lapply(filelist_allruns,"[[","all_files"))
filelist = c(filelist[!duplicated(filelist)], rhout_str)

# build copy/rsync command
# src = "~/Projects/CARB/Ward/"
dest = "/data/gpfs/assoc/firelab/william.burke/Ward/"
usr = "wburke@pronghorn.rc.unr.edu:"

copy_cmd = paste0("rsync -azPR ", paste0(filelist,collapse = " ") , " ", usr,dest)
dput(copy_cmd)


# -------------------- MANUAL PARALLEL RUNS --------------------
manualpara = F
if (manualpara) {
  library(parallel)
  start = Sys.time()
  cl = parallel::makeCluster(5)
  parallel::clusterExport(cl = cl, varlist = c("cmds"), envir = environment())
  parallel::parLapply(cl = cl, X = seq_along(cmds), fun = function(X, cmds) { system(cmds[X])}, cmds = cmds)
  parallel::stopCluster(cl)
  out_dir = collect_output(alert = F)
  end = Sys.time()
  end - start
}





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

# lapply(cmds, system)
# 
#system(cmds[1])
# 
# 
# 
# out_dir = collect_output()
# beepr::beep(3)

#plotpdf_allvars(out_dir, "msr_run")

