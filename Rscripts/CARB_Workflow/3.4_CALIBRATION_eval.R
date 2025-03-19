# calibration_eval
library(RHESSysIOinR)
library(rhutils)
library(hydroGOF)
library(tidyverse)
# install.packages("tidyverse")

# devtools::install_github(repo = "wburke24/rhutils", force = T)
# devtools::install_github(repo = "RHESSys/RHESSysIOinR")

# outfiles = list.files(out_dir,".csv",full.names = T)
# file.rename(outfiles, gsub("RunID11_","",outfiles))

options(scipen = 9999)

out_dir = "output/rh_out_2024-08-07--12-02-26/"
# ============================== Inputs ==============================
obs_source = "clim/BIG_C_AB_WHITES_GULCH_NR_GROVELAND_CA_1969_2024"

# cal_eval(Qsim = out_dir, Qobs = obs_source)

# ============================== Cal Metrics ==============================
Qobs = fread(obs_source) # CHECK IF IN MM
# Qobs = Qobs[Qobs$Date < "2013-1-1"]

sim_DT = get_basin_daily(out_dir)
fixQ = T
if (fixQ){
  sim_DT$streamflow = sim_DT$streamflow  + sim_DT$base_flow
}
eval = cal_eval(Qsim = sim_DT, Qobs = Qobs)
eval_mn = cal_eval(Qsim = sim_DT, Qobs = Qobs, monthly = T)

eval[order(eval$NSE, decreasing = T),]
eval_mn[order(eval_mn$NSE, decreasing = T),]

# runid = find_runID(out_dir)
runid = unique(str_extract(list.files(out_dir,".csv"),"(?<=RunID)\\d+(?=_\\d+)"))
load(list.files("robj/",runid,full.names = T)) # loads all the ioinr inputs including input_def

defchg_nse = def_changes_by_evaldf(defpar_df = defpars_list2df(input_def_pars),eval = eval_mn,stat = "NSE")
# defchg_nse = def_changes_by_eval(out_dir, eval_mn)
defchg_nse[1:nrow(defchg_nse),1:5]

defchg_nselog = def_changes_by_evaldf(defpar_df = defpars_list2df(input_def_pars),eval = eval_mn,stat = "NSElog")
defchg_nselog[1:nrow(defchg_nselog),1:5]

par_ranges_by_pctle = function(defchg_nse, pct = .95) {
  x = as.numeric(defchg_nse[1,3:ncol(defchg_nse)])
  pctl = quantile(x = x, probs = pct)
  tmp = defchg_nse[,3:ncol(defchg_nse)]
  def90 = tmp[,tmp[1,] >= pctl]
  ranges = data.frame(min = as.numeric(apply(def90,1,min)), max = as.numeric(apply(def90,1,max)))
  # signif(ranges, digits = 5)
  ranges = cbind(variable = defchg_nse[,1],ranges)
  cat("90th Percentile stat - Ranges of Pars\n------------------------------------\n")
  print(ranges)
  return(ranges)
}
ranges = par_ranges_by_pctle(defchg_nse)
ranges2 = par_ranges_by_pctle(defchg_nselog)

# ============================== Normal Plots ==============================
# if number of runs isnt huge
plotpdf_allvars(out_dir, "BigCreek_calplots", pdfwidth = 10)

# ============================== Sensitivity ==============================
sensout = defpars_sens(defpars = input_def_pars, eval_mn, stat = "NSE")
plot(sensout)

sensout2 = defpars_sens(defpars = input_def_pars, eval, stat = "RMSE")
plot(sensout2)

sensout3 = defpars_sens(defpars = input_def_pars, eval, stat = "NSElog")
plot(sensout3)

# ============================== Get Parameter Sets From Calibration Results ==============================
# just best set, based on monthly NSE
bestrun = eval_mn[order(eval_mn$NSE, decreasing = T),]$run[1]
bestrunID = as.numeric(str_extract(bestrun,"(?<=_)\\d+$"))
# defpars_list2df(input_def_pars)
input_def_pars_selected = lapply(input_def_pars,FUN = function(X,Y){X[[3]] = X[[3]][Y]; return(X)},bestrunID)
savebest = F
if (savebest) {
  save(input_def_pars_selected, file = "robj/input_def_pars_selected.RData")
}

# ================================================================================
sim_DT = get_basin_daily(out_dir)
sim_DT_mn = basin_daily2mn(sim_DT)
sim_DT_mn$month_wy = ifelse(sim_DT_mn$month >= 10, sim_DT_mn$month - 9, sim_DT_mn$month + 3)

# ================================================================================

# compare to obs
min_comb = max(min(sim_DT$date), min(Qobs$Date))
max_comb = min(max(sim_DT$date), max(Qobs$Date))
Qobs_sub = Qobs[Qobs$Date >= min_comb & Qobs$Date <= max_comb,]
sim_DT_sub = sim_DT[sim_DT$date >= min_comb & sim_DT$date <= max_comb,]
sim_DT_sub$streamflow = sim_DT_sub$streamflow * 1000
run_IDs = unique(sim_DT_sub$run)
run_IDs = run_IDs[order(as.numeric(gsub("[^0-9]+_", "",run_IDs)))] #

# sim_DT_sub_mn = aggregate(sim_DT_sub$streamflow, by = list(sim_DT_sub$run, sim_DT_sub$year, sim_DT_sub$month), FUN = mean)
sim_DT_sub_mn <- sim_DT_sub %>% group_by(run, year, month) %>% summarise_if(is.numeric, mean, na.rm=T)
sim_DT_sub_mn$year_month = as.yearmon(paste0(sim_DT_sub_mn$year,"-",sim_DT_sub_mn$month))
# names(sim_DT_sub_mn) = c("run", "year","month","streamflow")

Qobs_sub$year = as.numeric(format(Qobs_sub$Date, "%Y"))
Qobs_sub$month = as.numeric(format(Qobs_sub$Date, "%m"))
Qobs_sub_mn = aggregate(Qobs_sub$Flow_mmd, by = list(Qobs_sub$year, Qobs_sub$month), FUN = mean)
names(Qobs_sub_mn) = c("year","month","Flow_mmd")
Qobs_sub_mn$year_month = as.yearmon(paste0(Qobs_sub_mn$year,"-",Qobs_sub_mn$month))
datestr = paste0(gsub( ":", "-", sub( " ", "--", Sys.time())))

cal_plot_tsmn_obscompare = sim_DT_sub_mn %>%
  ggplot() + 
  geom_line(aes(x = year_month, y = streamflow, color = run), linewidth = 1.5) + 
  geom_line(data = Qobs_sub_mn, aes(x = year_month, y=Flow_mmd), linewidth = 1) +
  labs(title = "Streamflow timeseries", x = "Month")


ggsave(filename = paste0("plots/cal_tsmn_obscompare", datestr,".jpg"), cal_plot_tsmn_obscompare, height = 8, width = 12)

cal_plot_tsmn_obscompare_subset = sim_DT_sub_mn %>% filter(run %in% c("BigCreek_cal_8", "BigCreek_cal_47")) %>% 
  ggplot() + 
  geom_line(aes(x = year_month, y = streamflow, color = run), linewidth = 1.5) + 
  geom_line(aes(x = year_month, y = snowpack.water_equivalent_depth, color = run, linetype = "dash"), linewidth = 1.5) + 
  geom_line(data = Qobs_sub_mn, aes(x = year_month, y=Flow_mmd), linewidth = 1) +
  labs(title = "Streamflow timeseries", x = "Month")
ggsave(filename = paste0("plots/cal_tsmn_obscompare_subset", datestr,".jpg"), cal_plot_tsmn_obscompare_subset, height = 8, width = 12)

sim_DT_sub_mn %>% filter(run %in% c("BigCreek_cal_137", "BigCreek_cal_199")) %>% 
  ggplot() + 
  geom_line(aes(x = year_month, y = snowpack.water_equivalent_depth, color = run), linewidth = 1.5) 

# ================================================================================
cal_single_compare = sim_DT_sub_mn %>%
  ggplot() + 
  geom_line(aes(x = year_month, y = streamflow, color = "Sim",), linewidth = 1.5) + 
  geom_line(data = Qobs_sub_mn, aes(x = year_month, y=Flow_mmd, color = "Obs"), linewidth = 1) +
  scale_color_manual(values = c("Obs" = "black", "Sim"="red")) + 
  labs(title = "Streamflow timeseries", x = "Month", color = "Run")


datestr = paste0(gsub( ":", "-", sub( " ", "--", Sys.time())))
ggsave(filename = paste0("plots/cal_single_compare_", datestr,".jpg"), cal_single_compare, height = 8, width = 12)
# ================================================================================
