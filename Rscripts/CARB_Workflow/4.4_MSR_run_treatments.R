
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

cmds = vector(length = length(runs), mode = "character")

for (i in seq_along(runs)) {
  # -------------------- Project/Run Name --------------------
  name = paste0("Ward_msr", runs[i])
  # -------------------- Input RHESSys --------------------
  #clim = "clim/netcdf"
  clim = "clim/ward_netcdfgridmet_agg"
  
  #dates = read_clim(clim,dates_out = T)
  # clim_repeat(clim, "clim/ward_1000y", 1000, "years")
  
  #dates = c("1979 1 1 1", "2020 9 30 24")
  dates = c("2000 1 1 1", "2020 9 30 24")
  
  input_rhessys = IOin_rhessys_input(
    version = "../../../repos/rhessys-develop/rhessys/rhessys7.4",
    tec_file = paste0("tecfiles/ward",name,".tec"),
    world_file = worldfiles[i],
    world_hdr_prefix = name,
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
    variables = c("stratum.cs.net_psn", "patch.lai", "patch.totalc", "patch.evaporation", "patch.streamflow")
  )
  output_filter = IOin_output_filters(outfilter, file_name = paste0("./output/filters/filter",name,".yml") )
  
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

}


# -------------------- MANUAL PARALLEL RUNS --------------------
manualpara = F
if (manualpara) {
  library(parallel)
  cl = parallel::makeCluster(11)
  parallel::clusterExport(cl = cl, varlist = c("cmds"), envir = environment())
  parallel::parLapply(cl = cl, X = 1:6, fun = function(X, cmds) { system(cmds[X])}, cmds = cmds)
  # stop the cluster
  parallel::stopCluster(cl)
  
  out_dir = collect_output()
  
}


# lapply(cmds, system)
# 
#system(cmds[1])
# 
# 
# 
# out_dir = collect_output()
# beepr::beep(3)

#plotpdf_allvars(out_dir, "msr_run")

