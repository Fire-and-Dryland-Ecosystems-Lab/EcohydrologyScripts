# calibration_eval
library(RHESSysIOinR)
library(rhutils)
library(hydroGOF)
library(tidyverse)
options(scipen = 9999)

out_dir = "output/Ward_msr_cal2024-04-08--12-50-05/"
# ============================== Inputs ==============================
obs_source = "clim/WARD_C_AT_HWY_89_NR_TAHOE_PINES_CA_1972_2022"

# ============================== Cal Metrics ==============================
Qobs = fread(obs_source) # CHECK IF IN MM
sim_DT = get_basin_daily(out_dir)
fixQ = T
if (fixQ){
  sim_DT$streamflow = sim_DT$streamflow  + sim_DT$base_flow
}
eval = cal_eval(sim_DT = sim_DT, Qobs = Qobs)
eval_mn = cal_eval(sim_DT = sim_DT, Qobs = Qobs, monthly = T)

eval[order(eval$NSE, decreasing = T),]
eval_mn[order(eval_mn$NSE, decreasing = T),]

defchg_nse = def_changes_by_eval(out_dir, eval_mn)
defchg_nse

par_ranges_by_pctle = function(defchg_nse, pct = .90) {
  x = as.numeric(defchg_nse[1,3:ncol(defchg_nse)])
  pctl = quantile(x = x, probs = pct)
  tmp = defchg_nse[,3:ncol(defchg_nse)]
  def90 = tmp[,tmp[1,] >= pctl]
  ranges = data.frame(min = as.numeric(apply(def90,1,min)), max = as.numeric(apply(def90,1,max)))
  # signif(ranges, digits = 5)
  ranges = cbind(variable = defchg_nse[,1],ranges)
  cat("90th Percentile NSE - Ranges of Pars\n------------------------------------\n")
  print(ranges)
  return(ranges)
}
ranges = par_ranges_by_pctle(defchg_nse)

IOin_def_pars_sobol()
sensitivity::sobol()

# ============================== Sensitivity ==============================
pars_sens(out_dir, eval_mn)

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

sim_DT_sub_mn = aggregate(sim_DT_sub$streamflow, by = list(sim_DT_sub$run, sim_DT_sub$year, sim_DT_sub$month), FUN = mean)
names(sim_DT_sub_mn) = c("run", "year","month","streamflow")
Qobs_sub$year = as.numeric(format(Qobs_sub$Date, "%Y"))
Qobs_sub$month = as.numeric(format(Qobs_sub$Date, "%m"))
Qobs_sub_mn = aggregate(Qobs_sub$Flow_mmd, by = list(Qobs_sub$year, Qobs_sub$month), FUN = mean)
names(Qobs_sub_mn) = c("year","month","Flow_mmd")
Qobs_sub_mn$year_month = as.yearmon(paste0(Qobs_sub_mn$year,"-",Qobs_sub_mn$month))
sim_DT_sub_mn$year_month = as.yearmon(paste0(sim_DT_sub_mn$year,"-",sim_DT_sub_mn$month))

cal_plot_tsmn_obscompare = sim_DT_sub_mn %>%
  ggplot() + 
  geom_line(aes(x = year_month, y = streamflow, color = run), linewidth = 1.5) + 
  geom_line(data = Qobs_sub_mn, aes(x = year_month, y=Flow_mmd), linewidth = 1) +
  labs(title = "Streamflow timeseries", x = "Month")

datestr = paste0(gsub( ":", "-", sub( " ", "--", Sys.time())))
ggsave(filename = paste0("plots/cal_tsmn_obscompare", datestr,".jpg"), cal_plot_tsmn_obscompare, height = 8, width = 12)


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
