library(RHESSysIOinR)
library(rhutils)
source("R/0_global_vars.R")
source("../R/output_aliases.R")
#  MSR RHESSYS CALIBRATION INPUTS

# scenario_df[1,]

# -------------------- Project/Run Name --------------------
name = paste0(site,"_cal")
# -------------------- Input RHESSys --------------------
clim = "clim/netcdf"
dates = c("1989 1 1 1", "2020 9 30 24")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.5",
  tec_file = paste0("tecfiles/",name,".tec"),
  # world_file = scenario_df[1,]$worldfiles,
  world_file = "worldfiles/BigCreek_msr90m_nc_v4.world",
  world_hdr_prefix = paste0("hdr_",name),
  # flowtable = scenario_df[1,]$flowtables,
  flowtable = "flowtables/BigCreek_msr90m_nc.flow",
  start = dates[1], end = dates[2],
  output_folder = "output/", output_prefix = name,
  commandline_options = "-g -vmort_off -ncgridinterp 10 -netcdfgrid"
  # -climrepeat -msr
)

# -------------------- Input Headers --------------------
input_hdr = IOin_hdr(
  basin = "defs/basin.def",
  hillslope = "defs/hill.def",
  zone = "defs/zone_BigCreek.def",
  soil = c("defs/soil_sandy-loam_BigCreek.def", "defs/soil_loam_BigCreek.def"),
  landuse = "defs/lu.def",
  stratum = c("defs/stratum_evergreen_BigCreek.def", # "defs/stratum_evergreen_Ward.def" "defs/veg_evergreen.def" "defs/stratum_evergreen.def"
              "defs/stratum_shrub.def", # "defs/stratum_shrub_Ward.def",
              "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  basestations = paste0(clim, ".base")
)

# --------------------  Def File Parameters --------------------
pars_list = list(
  # EVERGREEN
  list(input_hdr$stratum_def[1], "epc.allocation_flag","dickenson"), # this mattered for stemc dickenson was better
  list(input_hdr$stratum_def[1], "epc.dickenson_pa",c(0.25)), # 0.25 default
  list(input_hdr$stratum_def[1], "epc.livewood_turnover",c(0.7)), # was 0.01 
  list(input_hdr$stratum_def[1], "epc.branch_turnover", c(0.01 )), #0.03,  0.01 p301
  list(input_hdr$stratum_def[1], "epc.froot_turnover", c(0.27)),
  list(input_hdr$stratum_def[1], "epc.alloc_stemc_leafc", c(1.2)), # c(0.85, 1.6) probably doesnt matter for stemc itself
  list(input_hdr$stratum_def[1], "epc.alloc_livewoodc_woodc", c(1)),
  list(input_hdr$stratum_def[1], "epc.root_distrib_parm", c(20)),
  list(input_hdr$stratum_def[1], "epc.min_percent_leafg", c( 0.6)),
  list(input_hdr$stratum_def[1], "epc.height_to_stem_coef", c(1)),
  list(input_hdr$stratum_def[1], "epc.storage_transfer_prop", c(1)), # 0.7 1  
  list(input_hdr$stratum_def[1], "epc.gl_c", 0.0002 ), #c(0.00006, 0.0006) # seems to not matter
  list(input_hdr$stratum_def[1], "epc.proj_sla", 9), # this matters 9 was better for stemc c(2, 9)
  list(input_hdr$stratum_def[1], "epc.min_daily_mortality", 0.005), # c(0.002, 0.005)
  list(input_hdr$stratum_def[1], "epc.froot_cn", 58), # c(58, 140)
  
  # SHRUB
  # list(input_hdr$stratum_def[2], "epc.livewood_turnover", 0.2),

  # SOIL PARS
  list(input_hdr$soil_def[1], "m", c(0.1, 0.5) ), #  0.22 c(0.08) c(0.3)
  list(input_hdr$soil_def[1], "m_v", c(0.1, 0.5)), # 0.17 c(0.04, 0.6)
  list(input_hdr$soil_def[1], "Ksat_0",  c(100, 500)), # 264 c(202)
  list(input_hdr$soil_def[1], "Ksat_0_v", c(100, 500)), #250 c(50, 500)
  list(input_hdr$soil_def[1], "sat_to_gw_coeff", c(0, 0.2)), #  0.01 c(0.02)
  list(input_hdr$soil_def[1], "psi_air_entry", c(0.3, 0.8)), # 0.58 c(0.56)
  list(input_hdr$soil_def[1], "pore_size_index", c(0.1, 0.4) ), # 0.2525  0.204 c(0.21)
  list(input_hdr$soil_def[1], "porosity_0",  c(0.2, 0.8)), # 0.49110 0.53
  
  list(input_hdr$soil_def[1], "gsurf_intercept",  0.001),
  # list(input_hdr$soil_def[1], "gsurf_slope",  0.01),
  list(input_hdr$soil_def[1], "snow_melt_Tcoef", c(0.001, 0.05)),    # c(0.001, 0.02) 0.005  0.02, default is 0.05, 0.001 is better
  list(input_hdr$soil_def[1], "albedo", c(0.01, 0.4)), # c(0.01, 0.28) 0.01 0.28
  list(input_hdr$soil_def[1], "maximum_snow_energy_deficit", -70), # -10.0000000 -70 maybe 71
  list(input_hdr$soil_def[1], "fs_spill", 1), # 1, 0-1
  list(input_hdr$soil_def[1], "fs_percolation", 1), # 1, 0-1
  list(input_hdr$soil_def[1], "fs_threshold", 0.2), # 0.2 m
  
  list(input_hdr$hillslope_def[1], "gw_loss_coeff", c(0, 1)), # 0.76
  
  list(input_hdr$zone_def, "lapse_rate_precip_default", c(0.001, 0.05)), # 0.011
  list(input_hdr$zone_def, "min_rain_temp", c(-5, 0)), #c(-5, -2) c(-5, 0) -4.93615248 -4.5
  list(input_hdr$zone_def, "max_snow_temp", c(-1, 5)) #c(-1, 5) c(-2, 5) 4.96597559 2.1
  
)

n = 10
# n = 88
# n = 1000
# n = 200

input_def_pars = lapply(pars_list, runif_sample, n)
# pars_list = fill_rep_pars(pars_list )

# input_def_pars = def_par_allcomb(pars_list)
 
# set mz and ksatz to be same as base values for efficiency mostly
# vars = lapply(input_def_pars, "[[",2)
# if ("m" %in% vars) { input_def_pars[[which(vars=="m_v")]][[3]] = input_def_pars[[which(vars=="m")]][[3]] }
# if ("Ksat_0" %in% vars) { input_def_pars[[which(vars=="Ksat_0_v")]][[3]] = input_def_pars[[which(vars=="Ksat_0")]][[3]] }

# # copy pars for all soils
input_def_pars = dup_soil_pars(input_def_pars, input_hdr)

# input_def_pars = IOin_def_pars_simple(pars_list)

# input_def_pars = NULL

#save(input_def_pars, file = "r_obj/cal_defpars.rdata")

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1990 10 1 1", end = dates[2], output_state = F)
#input_tec_data = IOin_tec_std(start = "1979 10 1 1", end = dates[2], output_state = F)

# -------------------- Output filter --------------------
outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = paste0(name,"_basin"),
  spatial_level = "basin",
  spatial_ID = 1,
  variables = c(output_cal) # output_carb
)

output_filter = IOin_output_filters(outfilter, file_name = paste0("output/filters/filter_",name,".yml"))

# -------------------- HANDLE RUN IDS --------------------
runid = GetRunID(increment = T)

save(input_rhessys, input_hdr, input_def_pars, input_tec_data, output_filter, file = paste0("robj/", name,"_RunID",runid, ".RData"))
# change output files to have an ID?
output_filter = AddRunIDtoOutputFilters(output_filter,runid)

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
# out_dir = collect_output(basename = name)
# plotpdf_allvars(out_dir = out_dir,out_name = "cal_basin")

hpc = T
if (hpc) {
  # TURN ON VPN FIRST
  runscript = rhessysIO2pronghorn(rhout, name, rh_bin_replace = "/RHESSys/rhessys7.5", 
                                  dest = "/data/gpfs/assoc/firelab/william.burke/BigCreek/", usr = "wburke@pronghorn.rc.unr.edu:",
                                  input_rhessys, input_hdr, input_def_pars, input_tec_data, output_filter,
                                  transfer_method = "scp")
  
  # THIS ASSUMES WORLDFILE, FLOWTABLE, AND CLIMATE FILES ARE ALREADY MOVED
  hpccmd <- paste0(
    "wsl ssh -i ~/rsa_key wburke@pronghorn.rc.unr.edu ",
    "\"cd /data/gpfs/assoc/firelab/william.burke/BigCreek/ && ./scripts/run_dyn.sh ", runscript, "\""
  )
  system(hpccmd)
  # "Submitted batch job XXYYZZ" IS INDICATION OF SUCCESS
  
}

gethpc = T
if (gethpc) {
  library(rhutils)
  # runid = GetRunID()
  
  gethpcout <- paste0(
    "wsl ssh -i ~/rsa_key wburke@pronghorn.rc.unr.edu ",
    "\"cd /data/gpfs/assoc/firelab/william.burke/BigCreek/ && ./scripts/collect_output.sh\"")
  cmdout = system(gethpcout,intern = T)
  # "Submitted batch job XXYYZZ" IS INDICATION OF SUCCESS
  hpcoutdir = gsub("Output CSVs moved to ","", cmdout[startsWith(cmdout, "Output CSVs moved to ")])
  rsync_cmd = paste0("wsl rsync -avz -e \"ssh -i ~/rsa_key\" wburke@pronghorn.rc.unr.edu:", hpcoutdir," ./output/")
  system(rsync_cmd)
  out_dir = paste0("output/",basename(hpcoutdir))
}

# hpcout = F
# if (hpcout) {
#   
# 
#   
#   # -------------------- HPC --------------------
#   #fix output strings and write to file
#   rhout_str = rhout_write_for_hpc(rhout, name, rh_bin_replace = "/RHESSys/rhessys7.5")
#   # make list of changed files
#   filelist = list_rh_input_files(input_rhessys, input_hdr, input_def_pars, input_tec_data, output_filter)
#   
#   # build copy/rsync command
#   # src = "~/Projects/CARB/Ward/"
#   dest = "/data/gpfs/assoc/firelab/william.burke/BigCreek/"
#   usr = "wburke@pronghorn.rc.unr.edu:"
#   
#   rsync = F
#   if (rsync) {
#     copy_cmd1 = paste0("rsync -azPR ", paste0(c(filelist$all_files,rhout_str),collapse = " ") , " ", usr,dest)
#     dput(copy_cmd1)
#     
#     writeLines(text = c(filelist$all_files,rhout_str),con = "scripts/filelist.txt")
#     copy_cmd2 = paste0("rsync -avvvPR --files-from=scripts/filelist.txt ./ ", usr,dest)
#     dput(copy_cmd2)
#   }
# 
#   # SCP version
#   # files scp file.txt remote_username@10.10.0.2:/remote/directory preserves source file name
#   # whole directory scp -r /local/directory remote_username@10.10.0.2:/remote/directory
#   # -p prserves access times
#   
#   # scptxt = c("#!/bin/bash",
#   #            # paste0("scp -i ~/rsa_key -p -r ",unique(dirname(filelist$tec_files))," ", usr,dest,unique(dirname(filelist$tec_files))),
#   #            paste0("scp -i ~/rsa_key -p -r ",unique(dirname(filelist$tec_files))," ", usr,dest),
#   #            paste0("scp -i ~/rsa_key -p -r ",unique(dirname(filelist$outputfilters))," ", usr,dest,"output"),
#   #            paste0("scp -i ~/rsa_key -p -r ",unique(dirname(filelist$hdr_files))," ", usr,dest,"worldfiles"),
#   #            paste0("scp -i ~/rsa_key -p -r ","defs"," ", usr,dest),
#   #            paste0("scp -i ~/rsa_key -p ",rhout_str," ", usr,dest,"scripts"))
#   
#   # alt
#   defstr = c()
#   for (i in seq_along(unique(dirname(filelist$def_files)))) {
#     defdir = unique(dirname(filelist$def_files))[i]
#     defdirfiles = filelist$def_files[dirname(filelist$def_files) == defdir]
#     defstr = c(defstr,paste0("scp -i ~/rsa_key -p ",paste(defdirfiles,collapse =" ")," ", usr,dest,defdir))
#   }
#   
#   scptxt = 
#     c("#!/bin/bash",
#       paste0("scp -i ~/rsa_key -p ",paste(filelist$tec_files,collapse =" ")," ", usr,dest,"tecfiles/"),
#       paste0("scp -i ~/rsa_key -p ",paste(filelist$outputfilters,collapse =" ")," ", usr,dest,"output/filters"),
#       paste0("scp -i ~/rsa_key -p ",paste(filelist$hdr_files,collapse =" ")," ", file.path(usr,dest,unique(dirname(filelist$hdr_files)))),
#       defstr,
#       paste0("scp -i ~/rsa_key -p ",rhout_str," ", usr,dest,"scripts")
#     )
#   
#   scpcon = file("scripts/scpfilelist.sh","wb")
#   writeLines(text = scptxt, con = scpcon, sep = "\n")
#   close(scpcon)
#   system("wsl ./scripts/scpfilelist.sh")
#   
# }
