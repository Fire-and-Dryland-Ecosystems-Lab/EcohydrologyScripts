
library(RHESSysIOinR)
library(rhutils)

# WARD CREEK STANDARD RHESSYS SOIL SPINUP

# -------------------- Project/Run Name --------------------
name = "Ward_std_veggrow"
# -------------------- Input RHESSys --------------------
#clim = "clim/ward_500y"
clim = "clim/netcdf"
# dates = read_clim(clim,dates_out = T)
# clim_repeat(clim, "clim/ward_1000y", 1000, "years")

dates = c("1980 9 1 1", "2479 9 30 24")
# dates = c("1980 9 1 1", "2020 9 30 24")

input_rhessys = IOin_rhessys_input(
  version = "../bin/rhessys7.4netcdf",
  tec_file = "tecfiles/ward.tec",
  world_file = "worldfiles/Ward90m_nc_1000yr_soilspin_reset.world",
  world_hdr_prefix = name,
  flowtable = "flowtables/Ward90m.flow",
  start = dates[1],
  end = dates[2],
  output_folder = "output/",
  output_prefix = name,
  commandline_options = "-g -netcdfgrid -ncgridinterp 11 -vmort_off -climrepeat"
)
 
# -------------------- Input Headers --------------------
input_hdr = IOin_hdr(
  basin = "defs/basin.def",
  hillslope = "defs/hill.def",
  zone = "defs/zone.def",
  soil = c("defs/soil_sandy-loam.def", "defs/soil_loam.def"),
  landuse = "defs/lu.def",
  stratum = c("defs/stratum_evergreen.def", "defs/stratum_shrub.def", "defs/stratum_grass.def", "defs/stratum_nonveg.def"),
  basestations = paste0(clim, ".base")
)

# --------------------  Def File Parameters --------------------
pars_list = list(
  # EVERGREEN
  list("defs/stratum_evergreen.def", "epc.froot_turnover", 0.3),
  list("defs/stratum_evergreen.def", "epc.livewood_turnover", 0.015), # was 0.01
  list("defs/stratum_evergreen.def", "epc.alloc_livewoodc_woodc", 0.55),
  list("defs/stratum_evergreen.def", "epc.alloc_stemc_leafc", 0.8),
  list("defs/stratum_evergreen.def", "epc.branch_turnover", 0.05),
  list("defs/stratum_evergreen.def", "epc.max_daily_mortality", 0.005),
  list("defs/stratum_evergreen.def", "epc.min_daily_mortality", 0.005),
  # SHRUB
  list("defs/stratum_shrub.def", "epc.livewood_turnover", 0.27),
  list("defs/stratum_shrub.def", "epc.leaf_turnover", 0.33),
  list("defs/stratum_shrub.def", "epc.froot_turnover", 0.45), # c(0.27, 0.6)
  list("defs/stratum_shrub.def", "epc.alloc_livewoodc_woodc", 0.95),
  list("defs/stratum_shrub.def", "epc.branch_turnover", 0.05),
  list("defs/stratum_shrub.def", "epc.max_daily_mortality", 0.005),
  list("defs/stratum_shrub.def", "epc.min_daily_mortality", 0.005)
)

pars_list = fill_rep_pars(pars_list)

input_def_pars = IOin_def_pars_simple(pars_list)
# input_def_pars = NULL

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start = "1980 10 1 1", end = dates[2], output_state = T)

# -------------------- Output filter --------------------
source("../R/output_aliases.R")

outfilter = build_output_filter(
  timestep = "daily",
  output_format = "csv",
  output_path = "output",
  output_filename = "growveg_basin",
  spatial_level = "basin",
  spatial_ID = 1,
  variables = c(output_vars_minimal,"stratum.cs.totalc", "stratum.cs.leafc", "stratum.cs.cpool", "stratum.cs.live_stemc", "stratum.cs.dead_stemc", "stratum.cs.live_crootc", "stratum.cs.dead_crootc",
  "stratum.cs.frootc", "stratum.ns.npool")
)
output_filter = IOin_output_filters(outfilter, file_name = "./output/filters/filter.yml")

# -------------------- Write Params --------------------
# if there is > 1 def file input
if (max(sapply(input_def_pars, function(X) length(X[[3]]))) > 1) {
  write_param_table(input_def_pars)
}

# -------------------- Run --------------------
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

# run_rhessys_multi(
#   input_rhessys = input_rhessys,
#   hdr_files = input_hdr,
#   def_pars = input_def_pars,
#   tec_data = input_tec_data,
#   output_filter = output_filter
# )

out_dir = collect_output()

plotpdf_allvars(out_dir, "growveg")


# ============================== PLOTS ==============================
plots = F
if (plots) {
  # ".//output/rh_out_2022-08-28--14-14-40"
  DT = get_basin_daily(out_dir = out_dir)
  DTm = basin_daily2mn(DT)
  
  names(DTm)
  
  theme_set(ggdark::dark_mode(cowplot::theme_cowplot()))
  
  
  ggplot(DTm) + 
    geom_line(aes(x = year_month, y = cs.dead_stemc, color = run))
  
  ggplot(DTm) + 
    geom_line(aes(x = year_month, y = cs.dead_crootc, color = run))
  
  
  p1 = ggplot(DTm) + 
    geom_line(aes(x = year_month, y = soil_cs.totalc)) +
    ylim(0,24) + 
    labs(title = "Soil Carbon", x = "Date", y = "Soil Carbon (Kg/m^2)") 
  p2 = ggplot(DTm) + 
    geom_line(aes(x = year_month, y = soil_ns.totaln)) +
    ylim(0,5) + 
    labs(title = "Soil Nitrogen", x = "Date", y = "Soil Nitrogen (Kg/m^2)")
  p3 = ggplot(DTm) + 
    geom_line(aes(x = year_month, y = lai)) +
    ylim(0,3) + 
    labs(title = "Leaf Area Index", x = "Date", y = "Leaf Area Index")
  p4 = ggplot(DTm) + 
    geom_line(aes(x = year_month, y = cs.leafc + cs.live_stemc)) +
    ylim(0,4) + 
    labs(title = "Live Leaf and Stem Carbon", x = "Date", y = "Live Leaf and Stem Carbon (Kg/m^2)")
  
  pg = plot_grid(p1, p2, p3, p4, ncol = 1, align="v")
  ggsave("plots/init_cnlai.jpg", pg, width = 10, height = 12)
}

# ============================== DATA EXPLORATION ==============================
dataexplore = F
if (dataexplore) {
  library(tidyverse)
  
  dtavg = DTm[, lapply(.SD, mean), by=c("run"), .SDcols = names(DTm)[!names(DTm) %in% c("run", "year", "month", "year_month","ym_ind")] ]
  
  # acending, so sorting by lowest dead stem/croot
  deadstem = dtavg[order(dtavg$cs.dead_stemc)]$run
  deadcroot = dtavg[order(dtavg$cs.dead_crootc)]$run
  deadcroot[1:5] %in% deadstem[1:5]
  
  newruns = as.numeric(gsub("growveg_","",deadcroot[1:5]))
  # [1] 14 28 10 40 32
  
  defpars_subset = function(defpars, run_ids) {
    newdefpars = lapply(defpars, function(X, Y) {X[[3]] = X[[3]][Y];return(X)}, run_ids)
    return(newdefpars)
  }
  
  newdefpars = defpars_subset(defpars, newruns)
  
  write_param_table(newdefpars)
  save(defpars, file = "r_obj/defpars_vegsens.rdata")
  
  
  get_param_table(input_def_pars)
}

