# -------------------- Execute RHESSys Runs - calibration --------------------
library(RHESSysIOinR)
library(rhutils)
load("r_obj/std_cal_inputs.rdata")

# -------------------- PARALELL RUNS --------------------
start = Sys.time()
rhout = run_rhessys_multi(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr,
  def_pars = input_def_pars,
  tec_data = input_tec_data,
  output_filter = output_filter,
  return_cmd = F
)
end = Sys.time()
end - start

write_param_table(input_def_pars)
out_dir = collect_output()
beepr::beep(2)


rhout2 = gsub('\"',"", as.character(gsub("bash -c ", "",unlist(rhout))))

writeLines(unlist(rhout2),"cal_inputs.txt")
save(rhout, file = "r_obj/std_rhout.rdata")


# -------------------- MANUAL PARALLEL RUNS --------------------
manualpara = F
if (manualpara) {
  load("r_obj/std_rhout.rdata")
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
beepr::beep()

write_param_table(input_def_pars)
out_dir = collect_output()
beepr::beep(2)





vars_defs = data.frame(variable = sapply(input_def_pars, "[[", 2), def_file = sapply(input_def_pars, "[[", 1))
param_table = cbind(vars_defs, t(sapply(input_def_pars, "[[", 3)))
names(param_table)[3:length(param_table[1,])] = paste0("run_",c(1:(length(param_table[1,])-2)))

#plotpdf_allvars_basin(out_dir = out_dir, out_name = "p301calibration")

# ----- RUN TIMES -----
# 14 runs 13 mins