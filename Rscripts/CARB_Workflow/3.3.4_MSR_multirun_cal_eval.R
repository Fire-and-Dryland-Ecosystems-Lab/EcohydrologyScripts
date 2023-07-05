# MSR_5b_multirun_cal_eval

# include multiple directories in calibration eval and sensitivity

library(RHESSysIOinR)
library(rhutils)
library(hydroGOF)
library(tidyverse)

# ================================================================================
# dirs to include in eval
dirs = paste0("output/", c(
  "rh_out_2023-03-01--14-47-10",
  "rh_out_2023-02-28--13-46-15",
  "rh_out_2023-02-27--20-35-55",
  "rh_out_2023-02-27--12-39-33",
  "rh_out_2023-02-26--15-43-41",
  "rh_out_2023-02-25--12-53-40",
  "rh_out_2023-02-24--11-43-05"))


dirlist = list()

Qobs = fread("clim/WARD_C_AT_HWY_89_NR_TAHOE_PINES_CA_1972_2022") # CHECK IF IN MM

# parse param files for each sim
itr = 0
for (i in seq_along(dirs)) {
  
  parfiles = list.files(paste0(dirs[i],"/params"), pattern = "^_.*", full.names = T)
  filenames = list.files(paste0(dirs[i],"/params"), pattern = "^_.*") # this has repeats
  par_tables_in = lapply(parfiles, read.table, header = FALSE, stringsAsFactors = FALSE, col.names=c("value","variable"))
  filebase = gsub("_\\d+","", 
                  gsub("^_", "", gsub("hillslope.params|basin.params|zone.params|stratum.params|soil.params|landuse.params","", filenames)))
  fbase_unique = unique(filebase)
  par_tables = mapply(FUN = function(X,Y,Z,n) {X$file = Y; X$path = Z; return(X)}, X = par_tables_in, 
                      Y = filebase, Z = parfiles, SIMPLIFY = F)

  filelist = lapply(fbase_unique, function(X) {parfiles[grepl(X,parfiles)] }) 
  len = sapply(filelist, length)
  # repeat def pars that weren't varied, all should be equal length now
  filelist_filled = lapply(filelist, function(X, Y) {if(length(X) == 1 & max(Y) > 1) { return(rep_len(X, max(Y))) } else {return(X)} }, len)
  filedf = as.data.frame(filelist_filled, col.names = fbase_unique)
  filedf$IDs  = as.numeric(gsub("\\D","", gsub(".*/.*/.*/", "",filedf[,which(len>1)[1]])))
  
  par_tables_list = lapply(seq_along(filelist_filled[[1]]), FUN = function(X) {
    params = par_tables[parfiles %in% filedf[X,]]
    return(rbindlist(params))
  } )

  eval_mn = cal_eval(out_dir = dirs[i], Qobs = Qobs, monthly = T)
  par_tables_list = mapply(function(X,Y,Z,nse) {X$itrID = Y; X$ID = Z; X$nse = nse; return(X)}, X = par_tables_list, 
                           Y = itr + seq_along(filelist_filled[[1]]),
                           Z = filedf$IDs, nse = unlist(eval_mn[order(filedf$IDs, decreasing = F),"NSE"], use.names = F),
                           SIMPLIFY = F)
  
  dirlist[[i]] = rbindlist(par_tables_list)
  
  itr = itr + length(filelist_filled[[1]])
}


all_pars = rbindlist(dirlist)

# ======================================== sensativity ========================================
tmp = pivot_wider(all_pars, id_cols = "itrID", names_from = c("variable","file"), values_from = "value")
#unique(all_pars$nse)
tmp = as.data.frame(tmp)
row.names(tmp) = paste0("run_",tmp[,1])
tmp = tmp[,-1]
lenunq = apply(tmp, 2, function(X) length(unique(X)))
tmp_diff = tmp[, lenunq > 1]

srcout = sensitivity::src(sapply(tmp_diff, as.numeric), unique(all_pars$nse))

# ======================================== max nse ========================================

all_pars[all_pars$nse == max(all_pars$nse),]

readpars = read.csv(list.files(paste0("output/rh_out_2023-02-27--12-39-33/params/"),pattern = "all_def_changes",full.names = T ))

table2list = function(X, Y) {
  value = unname(X[names(X)==Y])
  if (!is.na(suppressWarnings(as.numeric(value)))) {
    value = as.numeric(value)
  }
  out = list(unname(X["def_file"]), unname(X["variable"]), value)
  return(out) 
}
inputdefsread = apply(readpars,MARGIN = 1, FUN = table2list, "run_9")

