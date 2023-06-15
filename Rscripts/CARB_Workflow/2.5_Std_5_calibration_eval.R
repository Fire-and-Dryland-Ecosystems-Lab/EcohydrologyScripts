# calibration_eval

library(RHESSysIOinR)
library(rhutils)
library(hydroGOF)
source("R/misc_utils.R")

load("r_obj/std_cal_inputs.rdata")

#out_dir = "output/rh_out_2022-11-19--11-54-24/"

Qobs = fread("clim/WARD_C_AT_HWY_89_NR_TAHOE_PINES_CA_1972_2022")
# CHECK IF IN MM

eval = cal_eval(out_dir = out_dir, Qobs = Qobs)
eval_mn = cal_eval(out_dir = out_dir, Qobs = Qobs, monthly = T)

eval[order(eval$NSE, decreasing = T),]


inputpars = get_param_table(input_def_pars)

max(eval$NSE)
max(eval$R2)



# --------
pars = t(inputpars[,-c(1,2)])
colnames(pars) = paste0(inputpars$variable,"__",inputpars$def_file)
lenunq = apply(pars, 2, function(X) length(unique(X)))
inputpars_diff = pars[, lenunq > 1]

srcout = sensitivity::src(inputpars_diff, eval$NSE)


# --------
run_nsemax = eval[which.max(eval$NSE), "run"] # run 60

pars_nsemax = inputpars_diff[gsub("\\D","",rownames(inputpars_diff)) == gsub("\\D","",run_nsemax) , ]

readpars = read.csv("output/rh_out_2022-12-08--09-39-06/params/all_def_changes_2022-12-08_09.39.06.params",)

table2list = function(X, Y) {
  value = unname(X[names(X)==Y])
  if (!is.na(suppressWarnings(as.numeric(value)))) {
    value = as.numeric(value)
  }
  out = list(unname(X["def_file"]), unname(X["variable"]), value)
  return(out) 
}
inputdefsread = apply(readpars,MARGIN = 1, FUN = table2list, "run_60")

dput(inputdefsread)

# -------- PARAMETER SETS
# NSE 0.25:  sat_to_gw_coeff__defs/sandyloam_test.def   psi_air_entry__defs/sandyloam_test.def               m__defs/sandyloam_test.def 
#             0.01325159                               0.18217991                               0.68355381 
#   m_v__defs/sandyloam_test.def          Ksat_0__defs/sandyloam_test.def        Ksat_0_v__defs/sandyloam_test.def 
#   0.68355381                             228.94942248                             228.94942248 
#   min_snow_temp__defs/zone.def             max_snow_temp__defs/zone.def
#  -4.93615248                               4.9659755
# --------
sim_DT = get_basin_daily(out_dir)
sim_DT_mn = basin_daily2mn(sim_DT)
sim_DT_mn$month_wy = ifelse(sim_DT_mn$month >= 10, sim_DT_mn$month - 9, sim_DT_mn$month + 3)


#%>% filter(run == "cal_1")

cal_plot = sim_DT_mn %>% group_by(run, month_wy) %>% summarise(Streamflow = mean(streamflow)) %>% 
  ggplot() + 
  geom_line(aes(x = month_wy, y = Streamflow, color = run), linewidth = 1.5) + xlab("Month")



plotpdf_allvars_basin(out_dir = "./output/", out_name = "cal")


# take top 10 runs based on nse
nse_runs = eval[1:10,'run']

vars_defs = data.frame(variable = sapply(input_def_pars, "[[", 2), def_file = sapply(input_def_pars, "[[", 1))
param_table = cbind(vars_defs, t(sapply(input_def_pars, "[[", 3)))
names(param_table)[3:length(param_table[1,])] = paste0("run_",c(1:(length(param_table[1,])-2)))

nse_runs_pars = cbind(vars_defs, param_table[,paste0("run_", nse_runs) ])
nse_runs_pars$min = apply(nse_runs_pars[,3:12],1,min)
nse_runs_pars$max = apply(nse_runs_pars[,3:12],1,max)


keepruns = nse_runs[1]
nse_def_pars = input_def_pars
nse_def_pars = lapply(nse_def_pars, function(X,Y) {X[[3]] = X[[3]][Y]; return(X) }, keepruns)
save(nse_def_pars, file = "r_obj/nse_def_pars1.rdata")


ggplot() + geom_line(data = DT, aes(x = date, y = streamflow*1000 ), color = "blue") + geom_line(data = Qobs, aes(x = date, y = mm))

sum(Qobs$mm)
sum(DT[DT$run == i, streamflow])

length(Qobs$mm)
length(DT[DT$run == i, streamflow])

DT$et = DT$evaporation + DT$transpiration_unsat_zone + DT$transpiration_sat_zone

DT %>% group_by(wy) %>% summarise(et)


# pdfname = file.path("plots", paste0(gsub(".pdf","", out_name), gsub( ":", ".", sub( " ", "_", Sys.time())), ".pdf"  ))
# pdf(pdfname)
# for (i in seq_along(vars)) {
#   tmpplot = DT %>% ggplot() + aes(x = year_month, color = as.factor(run), linetype =as.factor(run)) + aes_string(y = vars[i]) + geom_line() + ggtitle(vars[i])
#   plot(tmpplot)
# }
# dev.off()