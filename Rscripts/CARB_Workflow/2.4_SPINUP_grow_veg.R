
library(RHESSysIOinR)
library(rhutils)
source("r/0_global_vars.R")

# WARD CREEK MSR RHESSYS SOIL SPINUP + VEG GROW

resetveg = F
if (resetveg) {
  world_path = "worldfiles/Ward_msr90m_basinsoilspin400y.world"
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

cmds = vector(length = length(spinup_df$runs), mode = "character")

for (i in seq_along(spinup_df$runs)) {
  
  # -------------------- Project/Run Name --------------------
  name = paste0("Ward_msr90m_",spinup_df$runs[i])
  # name = "Ward_msr_spingrow"
  # name = "Ward_msr_targetspinup"
  
  # -------------------- Input RHESSys --------------------
  #clim = "clim/netcdf"
  clim = "clim/ward"
  
  # dates = read_clim(clim,dates_out = T)
  # clim_repeat(clim, "clim/ward_1000y", 1000, "years")
  
  dates = c("1979 1 1 1", "2029 9 30 24")
  # dates = c("1979 1 1 1", "2179 9 30 24")
  
  input_rhessys = IOin_rhessys_input(
    version = "../bin/rhessys7.4",
    tec_file = "tecfiles/ward.tec",
    world_file = spinup_df$reset_worlds[i],
    world_hdr_prefix = paste0("hdr_",name),
    # flowtable = "flowtables/Ward_msr90m.flow",
    flowtable = spinup_df$flowtables[i],
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
    zone = "defs/zone_Ward.def",
    soil = c("defs/soil_sandy-loam_Ward.def", "defs/soil_loam_Ward.def"),
    landuse = "defs/lu_Ward.def",
    stratum = c("defs/stratum_evergreen_Ward.def", "defs/stratum_evergreen_Ward_thin.def", "defs/stratum_shrub_Ward.def", 
                "defs/stratum_shrub_Ward_thin.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
    # spinup = "defs/spinup.def",
    basestations = paste0(clim, ".base")
  )
  # names(input_hdr)[names(input_hdr) == "spinup"] = "spinup_def"
  
  # --------------------  Def File Parameters --------------------
  
  # load("r_obj/cal_defpars.rdata")
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
  # input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = T)
  input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = T)
  

  
  # outputdates = seq(1980, 2170, 10)
  # for (i in outputdates) {
  #   input_tec_data = add_tec(input_tec_data, paste(i, "8 1 1"), tecname = "print_daily_on")
  #   input_tec_data = add_tec(input_tec_data, paste(i, "8 1 2"), tecname = "print_daily_growth_on")
  #   input_tec_data = add_tec(input_tec_data, paste(i, "8 2 1"), tecname = "print_daily_off")
  #   input_tec_data = add_tec(input_tec_data, paste(i, "8 2 2"), tecname = "print_daily_growth_off")
  # }
  
  
  # -------------------- Output filter --------------------
  source("../R/output_aliases.R")
  
  outfilter = build_output_filter(
    timestep = "daily",
    output_format = "csv",
    output_path = "output",
    output_filename = paste0(name,"_basin"),
    spatial_level = "basin",
    spatial_ID = 1,
    variables = c(output_vars_minimal,"stratum.cs.leafc", "stratum.cs.live_stemc")
  )
  outfilter2 = build_output_filter(
    timestep = "daily",
    output_format = "csv",
    output_path = "output",
    output_filename = paste0(name,"_patch"),
    spatial_level = "patch",
    spatial_ID = 1,
    variables = c("patch.lai")
  )
  output_filter = IOin_output_filters(outfilter,outfilter2, file_name = paste0("./output/filters/filter",name,".yml"))
  
  RHESSysIOinR::write_output_filter(output_filter = output_filter)
  
  
  # -------------------- Run --------------------
  cmds[i] = run_rhessys_single(
    input_rhessys = input_rhessys,
    hdr_files = input_hdr,
    def_pars = input_def_pars,
    tec_data = input_tec_data,
    output_filter = output_filter,
    return_cmd = T,
    write_run_metadata = T
  )
  # out_dir = collect_output()
  
}


manualpara = F
if (manualpara) {
  library(parallel)
  start = Sys.time()
  cl = parallel::makeCluster(6)
  parallel::clusterExport(cl = cl, varlist = c("cmds"), envir = environment())
  parallel::parLapply(cl = cl, X = seq_along(cmds), fun = function(X, cmds) { system(cmds[X])}, cmds = cmds)
  parallel::stopCluster(cl)
  out_dir = collect_output(alert = F)
  end = Sys.time()
  end - start
}


plotpdf_allvars(out_dir, "msr_spin_grow_veg")

statefile = sapply(spinup_df$reset_worlds, worldstate, USE.NAMES = F)
# name = world_name(input_rhessys, "_soilspin")

file.rename(statefile, paste0("worldfiles/Ward_msr90m_",spinup_df$runs,"_init.world"))


