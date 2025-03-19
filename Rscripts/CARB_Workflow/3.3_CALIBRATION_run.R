# -------------------- Execute RHESSys Runs - calibration --------------------
library(RHESSysIOinR)
library(rhutils)
source("R/3.1_CALIBRATION_inputs.R")

# -------------------- PARALELL RUNS --------------------
write_param_table(input_def_pars)
rhout = run_rhessys_multi(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr,
  def_pars = input_def_pars,
  tec_data = input_tec_data,
  output_filter = output_filter,
  return_cmd = F,
  n_cores = 11
)
out_dir = collect_output()


# rhout2 = gsub('\"',"", as.character(gsub("bash -c ", "",unlist(rhout))))

# writeLines(unlist(rhout2),"cal_inputs.txt")
# save(rhout, file = "r_obj/std_rhout.rdata")


# -------------------- MANUAL PARALLEL RUNS --------------------
manualpara = F
if (manualpara) {
  #load("r_obj/std_rhout.rdata")
  rhout= readLines("scripts/BigCreek_cal_runcmds_2024-05-30--11-47-35_TEST.txt")
  rhout =
  
  library(parallel)
  n_cores = parallel::detectCores() - 1
  start = Sys.time()
  cl = parallel::makeCluster(n_cores)
  parallel::clusterExport(cl = cl, varlist = c("rhout"), envir = environment())
  parallel::parLapply(cl = cl, X = 201:400, fun = function(X, rhout) { system(rhout[[X]])}, rhout = rhout)
  # stop the cluster
  parallel::stopCluster(cl)
  end = Sys.time()
  end - start
}


# -------------------- SINGLE RUN --------------------
singlerun = F
if (singlerun) {
  run_rhessys_single(
    input_rhessys = input_rhessys,
    hdr_files = input_hdr,
    def_pars = input_def_pars,
    tec_data = input_tec_data,
    output_filter = output_filter,
    return_cmd = F,
    write_run_metadata = T
  )
  # write_param_table(input_def_pars)
  out_dir = collect_output()
}




