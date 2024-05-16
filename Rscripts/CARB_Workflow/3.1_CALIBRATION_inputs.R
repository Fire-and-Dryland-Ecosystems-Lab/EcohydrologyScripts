library(RHESSysIOinR)
library(rhutils)
source("R/0_global_vars.R")
# WARD CREEK MSR RHESSYS CALIBRATION INPUTS

# scenario_df[1,]

# -------------------- Project/Run Name --------------------
# name = "Ward_msr_cal_vegsens"
name = "Ward_msr_cal"
# -------------------- Input RHESSys --------------------
#clim = "clim/netcdf"
# clim = "clim/ward"
clim = "clim/ward_netcdfgridmet_agg"

#dates = read_clim(clim, dates_out = T)
#dates = c("1990 1 1 1", "2000 9 30 24")
dates = c("1989 1 1 1", "2020 9 30 24")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.4",
  tec_file = "tecfiles/ward.tec",
  world_file = scenario_df[1,]$worldfiles,
  world_hdr_prefix = paste0("hdr_",name),
  flowtable = scenario_df[1,]$flowtables,
  start = dates[1],
  end = dates[2],
  output_folder = "output/",
  output_prefix = name,
  commandline_options = "-g -vmort_off"
  #-climrepeat -msr
)

# -------------------- Input Headers --------------------
input_hdr = IOin_hdr(
  basin = "defs/basin.def", hillslope = "defs/hill.def", zone = "defs/zone_Ward.def",
  soil = c("defs/soil_sandy-loam_Ward.def", "defs/soil_loam_Ward.def"), landuse = "defs/lu_Ward.def",
  stratum = c("defs/stratum_evergreen_Ward.def", "defs/stratum_shrub_Ward.def",  "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  basestations = paste0(clim, ".base")
)

# --------------------  Def File Parameters --------------------
old = F
if (old) {
  # pars_list = list(
  # EVERGREEN
  # list(input_hdr$stratum_def[1], "epc.froot_turnover", 0.26),
  # list(input_hdr$stratum_def[1], "epc.livewood_turnover", 0.015), # was 0.01
  # list(input_hdr$stratum_def[1], "epc.alloc_livewoodc_woodc", 0.55),
  # list(input_hdr$stratum_def[1], "epc.alloc_stemc_leafc", c(0.8,0.2, 1.62, 0.8,0.2, 1.62,0.8,0.2, 1.62)),
  # list(input_hdr$stratum_def[1], "epc.alloc_frootc_leafc", c(0.93,1.4, 1.33,0.93,1.4, 1.33,0.93,1.4, 1.33)),
  # list(input_hdr$stratum_def[1], "epc.alloc_crootc_stemc", c(0.27,0.5, 0.3,0.27,0.5, 0.3,0.27,0.5, 0.3)),
  # list(input_hdr$stratum_def[1], "epc.allocation_flag", c("combined","combined","combined","waring","waring","waring","dickenson","dickenson","dickenson"  ))
  # 
  # list(input_hdr$stratum_def[1], "epc.alloc_stemc_leafc", c(0.2)),
  # list(input_hdr$stratum_def[1], "epc.alloc_frootc_leafc", c(1.4)),
  # list(input_hdr$stratum_def[1], "epc.alloc_crootc_stemc", c(0.5)),
  # list(input_hdr$stratum_def[1], "epc.allocation_flag", c("waring")),
  
  # list(input_hdr$stratum_def[1], "epc.branch_turnover", 0.05), # 0.01 p301
  # list(input_hdr$stratum_def[1], "epc.max_daily_mortality", 0.005),
  # list(input_hdr$stratum_def[1], "epc.min_daily_mortality", 0.005),
  
  
  # SHRUB
  # list(input_hdr$stratum_def[2], "epc.livewood_turnover", 0.27),
  # list(input_hdr$stratum_def[2], "epc.leaf_turnover", 0.33),
  # list(input_hdr$stratum_def[2], "epc.froot_turnover", 0.45), # c(0.27, 0.6)
  # list(input_hdr$stratum_def[2], "epc.alloc_livewoodc_woodc", 0.95),
  # list(input_hdr$stratum_def[2], "epc.branch_turnover", 0.05), # 0.2 p301
  # list(input_hdr$stratum_def[2], "epc.max_daily_mortality", 0.005),
  # list(input_hdr$stratum_def[2], "epc.min_daily_mortality", 0.005),
  
  
  #   
  #   list(input_hdr$landuse_def[1], "sh_l", 1),
  #   list(input_hdr$landuse_def[1], "sh_g", 1),
  #   
  #   # ----- soils -----
  #   list(input_hdr$soil_def[1], "active_zone_z", 10.0),       # 10
  #   list(input_hdr$soil_def[1], "DOM_decay_rate", 0.05),     # 0
  #   # list(input_hdr$soil_def[1], "NO3_adsorption_rate", 0.0), # 0.0000005
  #   list(input_hdr$soil_def[1], "N_decay", 0.12),            # 1.2
  #   list(input_hdr$soil_def[1], "soil_depth", 200.0),         # 10
  # 
  #   list(input_hdr$soil_def[1], "gsurf_intercept", 0.01),  # 0.001
  #   
  #   # ============================== VARIED PARAMETERS ==============================
  # list(input_hdr$soil_def[1], "sat_to_gw_coeff", 0.01), # c(0.0, 0.4) 0.01325159 0.3132843
  # list(input_hdr$soil_def[1], "psi_air_entry", c(0.58)), # 0.18 c(0.1, 0.4) 0.218 0.18217991 0.45
  # list(input_hdr$soil_def[1], "pore_size_index", c(0.2)), # c(0.2, 0.25) 0.2525  0.2040000
  # list(input_hdr$soil_def[1], "porosity_0", 0.5), # c(0.35, 0.6)0.49110000  0.4350000
  # list(input_hdr$soil_def[1], "m", c(0.2, 0.5)), # c(0.5, 0.6) c(0.04, 0.6) .56 0.35
  # list(input_hdr$soil_def[1], "m_v", c( 0.17)), # c(0.04, 0.6)
  # list(input_hdr$soil_def[1], "Ksat_0", c(150,300)), # c(200, 325) c(50, 500) 250 228.94942248
  # list(input_hdr$soil_def[1], "Ksat_0_v", c(250)), # c(50, 500)
  #   list(input_hdr$soil_def[1], "snow_melt_Tcoef", c(0.001)),    # c(0.001, 0.02) 0.005  0.02, default is 0.05, 0.001 is better
  #   list(input_hdr$soil_def[1], "albedo", 0.03), # c(0.01, 0.28) 0.01 0.28
  #   list(input_hdr$soil_def[1], "maximum_snow_energy_deficit", -70), # -10.0000000 -70 maybe 71
  # # list(input_hdr$soil_def[1], "fs_spill", 1), # 1, 0-1
  # # list(input_hdr$soil_def[1], "fs_percolation", c(0, 0.2, 0.4, 0.6, 0.8, 1)), # 1, 0-1
  # # list(input_hdr$soil_def[1], "fs_threshold", c(0.05,0.2, 0.3)), # 0.2 m
  #   
  #   # ----- ZONE -----
  #   list(input_hdr$zone_def, "lapse_rate_precip_default", c(0.011)), #
  #   list(input_hdr$zone_def, "min_rain_temp", -4.5), #c(-5, -2) c(-5, 0) -4.93615248 -4.5
  #   list(input_hdr$zone_def, "max_snow_temp", c(5)), #c(-1, 5) c(-2, 5) 4.96597559
  #   list(input_hdr$zone_def, "lapse_rate_tmax", c(0.007)), #c(-0.0052, -0.002, 0, 0.001, 0.003, 0.0063)
  #   list(input_hdr$zone_def, "lapse_rate_tmin", c(0.007)) #c(-0.0064,-0.002, 0, 0.0012, 0.003, 0.006)
  # )
  # )
}


pars_list = list(
  # EVERGREEN
  list(input_hdr$stratum_def[1], "epc.froot_turnover", 0.26),
  # list(input_hdr$stratum_def[1], "epc.alloc_stemc_leafc", c(0.8,0.2, 1.62, 0.8,0.2, 1.62,0.8,0.2, 1.62)),
  # list(input_hdr$stratum_def[1], "epc.alloc_frootc_leafc", c(0.93,1.4, 1.33,0.93,1.4, 1.33,0.93,1.4, 1.33)),
  # list(input_hdr$stratum_def[1], "epc.alloc_crootc_stemc", c(0.27,0.5, 0.3,0.27,0.5, 0.3,0.27,0.5, 0.3)),
  # list(input_hdr$stratum_def[1], "epc.allocation_flag", c("combined","combined","combined","waring","waring","waring","dickenson","dickenson","dickenson"  ))
  
  list(input_hdr$stratum_def[1], "epc.alloc_stemc_leafc", c(0.2)),
  list(input_hdr$stratum_def[1], "epc.alloc_frootc_leafc", c(1.4)),
  list(input_hdr$stratum_def[1], "epc.alloc_crootc_stemc", c(0.5)),
  list(input_hdr$stratum_def[1], "epc.allocation_flag", c("waring")),

    # ----- soils -----
  list(input_hdr$soil_def[1], "m", c(0.15, 0.5) ), #  0.22
  list(input_hdr$soil_def[1], "m_v", c(0.22)), # 0.17 c(0.04, 0.6)
  list(input_hdr$soil_def[1], "Ksat_0",  c(150,350) ), # 264
  list(input_hdr$soil_def[1], "Ksat_0_v", 264), #250 c(50, 500)
  list(input_hdr$soil_def[1], "sat_to_gw_coeff", c(0.0, 0.1)), #  0.01
  list(input_hdr$soil_def[1], "psi_air_entry", c(0.25, 0.75)), # 0.58
  list(input_hdr$soil_def[1], "pore_size_index", c(0.15, 0.25) ), # 0.2525  0.204
  list(input_hdr$soil_def[1], "porosity_0",  c(0.35, 0.6)), # 0.49110
  # list(input_hdr$soil_def[1], "snow_melt_Tcoef", c(0.001)),    # c(0.001, 0.02) 0.005  0.02, default is 0.05, 0.001 is better
  
  # list(input_hdr$zone_def, "lapse_rate_precip_default", c(0.011)), #
  list(input_hdr$zone_def, "min_rain_temp", c(-5, -2)), #c(-5, -2) c(-5, 0) -4.93615248 -4.5
  list(input_hdr$zone_def, "max_snow_temp", c(-1, 5)), #c(-1, 5) c(-2, 5) 4.96597559
  list(input_hdr$zone_def, "lapse_rate_tmax", c(0.003, 0.007)), #c(-0.0052, -0.002, 0, 0.001, 0.003, 0.0063)
  list(input_hdr$zone_def, "lapse_rate_tmin", c(0.003, 0.007)) #c(-0.0064,-0.002, 0, 0.0012, 0.003, 0.006)
  
)


#n = 1
# n = 88
# n = 1000
n = 250

input_def_pars = lapply(pars_list, runif_sample, n)
# pars_list = fill_rep_pars(pars_list )

# def_par_allcomb = function(defpars) {
#   npars = lapply(X = defpars, FUN = function(X) {out = length(X[[3]]); return(out)})
#   if (!any(npars > 1)) {
#     cat("No pars to combine")
#     return(defpars)
#   }
#   pars_mult = lapply(defpars[npars>1], "[[", 3)
#   #pars_comb = expand.grid(as.list(unname(data.frame(pars_mult))))
#   pars_comb = expand.grid(pars_mult)
#   defpars[npars>1] = mapply(function(X, Y) {X[[3]] = Y; return(X)}, defpars[npars>1], pars_comb, SIMPLIFY = F)
#   defpars[npars==1] = lapply(pars_list[npars==1], function(X, Y) {X[[3]] = rep.int(X[[3]], Y); return(X)}, length(pars_comb[[1]]))
#   cat("Output def pars length: ", length(pars_comb[[1]]))
#   return(defpars)
# }
# input_def_pars = def_par_allcomb(pars_list)
# 
# set mz and ksatz to be same as base values for efficiency mostly
vars = lapply(input_def_pars, "[[",2)
if ("m" %in% vars) { input_def_pars[[which(vars=="m_v")]][[3]] = input_def_pars[[which(vars=="m")]][[3]] }
if ("Ksat_0" %in% vars) { input_def_pars[[which(vars=="Ksat_0_v")]][[3]] = input_def_pars[[which(vars=="Ksat_0")]][[3]] }
# 
# # copy pars for all soils
input_def_pars = dup_soil_pars(input_def_pars, input_hdr)

# input_def_pars = IOin_def_pars_simple(pars_list)

# input_def_pars = NULL

#save(input_def_pars, file = "r_obj/cal_defpars.rdata")

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1990 10 1 1", end = dates[2], output_state = F)
#input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = F)

# -------------------- Output filter --------------------
source("../R/output_aliases.R")
outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  # output_filename = "cal_vegsens_msr_basin",
  output_filename = "cal_msr_basin",
  spatial_level = "basin",
  spatial_ID = 1,
  # variables = output_carb
  variables = c(output_cal)
)

output_filter = IOin_output_filters(outfilter, file_name = "output/filters/cal_filter.yml")

# save rhessys objs
# outname = paste0("msr_cal_inputs.rdata")
# save(input_rhessys, input_hdr, input_def_pars, input_tec_data, output_filter, file = file.path("r_obj", outname))

# -------------------- RUN --------------------
write_param_table(input_def_pars)
rhout = run_rhessys_multi(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr,
  def_pars = input_def_pars,
  tec_data = input_tec_data,
  output_filter = output_filter,
  return_cmd = T,
  n_cores = 11
)
out_dir = collect_output(basename = name)

plotpdf_allvars(out_dir = out_dir,out_name = "cal_basin")

# -------------------- HPC --------------------
rhout_str = rhout_write_for_hpc(rhout, name)
filelist = list_rh_input_files(input_rhessys, input_hdr, input_def_pars, input_tec_data, output_filter)

# build copy/rsync command
# src = "~/Projects/CARB/Ward/"
dest = "/data/gpfs/assoc/firelab/william.burke/Ward/"
usr = "wburke@pronghorn.rc.unr.edu:"

copy_cmd1 = paste0("rsync -azPR ", paste0(c(filelist$all_files,rhout_str),collapse = " ") , " ", usr,dest)
dput(copy_cmd1)

writeLines(text = c(filelist$all_files,rhout_str),con = "scripts/filelist.txt")
copy_cmd2 = paste0("rsync -avvvPR --files-from=scripts/filelist.txt ./ ", usr,dest)
dput(copy_cmd2)
# all(file.exists(c(filelist$all_files,rhout_str)))

# writeLines(text = c(filelist$all_files,rhout_str),con = "scripts/filelist.txt")
transfer_dirs = unique(dirname(filelist$all_files))
transfer_dirs = transfer_dirs[!transfer_dirs %in% c("worldfiles", "flowtables", "clim", "defs", "tecfiles")]

copy_cmd3 = paste0("rsync -avvvnPR ",paste0(transfer_dirs,collapse = " ")," ", usr,dest)
dput(copy_cmd3)



# folders for many files - defs/vegX worldfiles/hdrX output/filters