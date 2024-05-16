# MSR hillslope plots -- SOIL DEBUG

library(RHESSysIOinR)
library(rhutils)
library(tidyverse)
library(fasstr)
library(cowplot)
library(zoo)
library(terra)
source("R/0_global_vars.R")

treat_date = "1999-8-1"
treat_yrmn = as.yearmon(treat_date)

# --------------------------------------------------------------
# -------------------- Spatial Data --------------------
# --------------------------------------------------------------
hillmap = trim(rast("preprocessing/spatial90m/subbasins.tif"))
hillct = length(unique(values(hillmap)))
my_palette <- viridis(hillct)
plot(hillmap,type = "classes", col = sample(my_palette))



summary(as.factor(values(hillmap)))
streams = trim(rast("preprocessing/spatial90m/streams.tif"))
plot(streams)

# --------------------------------------------------------------
# -------------------- Data Ingest --------------------
# --------------------------------------------------------------

load("r_obj/patch_monthly/allruns_hill_monthly.rda") # hill_mn_DT

hill_mn_DT$Treatment = "NA"
hill_mn_DT[yr_mn < treat_yrmn, "Treatment"] = "Pre Treatment"
hill_mn_DT[yr_mn >= treat_yrmn, "Treatment"] = "Post Treatment"
hill_mn_DT$Treatment = factor(hill_mn_DT$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))
hill_mn_DT$sm = (hill_mn_DT$rz_storage + hill_mn_DT$unsat_storage)/hill_mn_DT$sat_deficit

hill_year = hill_mn_DT %>% group_by(wy, hillID, Scenario) %>% summarise_if(is.numeric, mean, na.rm=T)
hill_year = hill_year %>% mutate(across(c(streamflow, et, trans, PET,psn ,theta, rz_storage, unsat_storage, sat_deficit,pcp, rain_thru, snow_thru,sm), ~ . * 365))

# --------------------------------------------------------------
# -------------------- Plots --------------------
# --------------------------------------------------------------

# average across hills
tmp = hill_mn_DT %>% group_by(yr_mn) %>% summarize(unsat_storage = mean(unsat_storage)) %>%
  ggplot() +
  aes(x = yr_mn, y = unsat_storage) +
  geom_line(linewidth = 1)

tmp2 = hill_mn_DT %>% group_by(yr_mn) %>% summarize(sat_deficit = mean(sat_deficit)) %>%
  ggplot() +
  aes(x = yr_mn, y = sat_deficit) +
  geom_line(linewidth = 1)

tmp3 = hill_mn_DT %>% group_by(yr_mn) %>% summarize(rz_storage = mean(rz_storage)) %>%
  ggplot() +
  aes(x = yr_mn, y = rz_storage) +
  geom_line(linewidth = 1)

soilvars = plot_grid(tmp, tmp2, tmp3)


# include all hills, boxplots
tmp = hill_mn_DT %>%
  ggplot() +
  aes(x = hillID, y = unsat_storage, group = hillID) +
  geom_boxplot() +
  stat_summary(geom = 'text', label = sort(unique(hill_mn_DT$hillID)), fun = max, vjust = -1) +
  ggtitle("Unsat Storage")


tmp2 = hill_mn_DT %>% 
  ggplot() +
  aes(x = hillID, y = sat_deficit, group = hillID) +
  geom_boxplot()+
  stat_summary(geom = 'text', label = sort(unique(hill_mn_DT$hillID)), fun = max, vjust = -1) +
  ggtitle("Sat Deficit Storage")

tmp3 = hill_mn_DT %>% 
  ggplot() +
  aes(x = hillID, y = rz_storage, group = hillID) +
  geom_boxplot() +
  stat_summary(geom = 'text', label = sort(unique(hill_mn_DT$hillID)), fun = max, vjust = -1) +
  ggtitle("Root Zone Storage")

tmpdf = freq(hillmap)

tmp4 = tmpdf %>% ggplot() +
  aes(x = value, y = count) +
  geom_col() +
  geom_text(aes(label = value), vjust = -0.5, size = 4) + ggtitle("Subbasin area by ID")

soilvarsbox = plot_grid(tmp, tmp2, tmp3, tmp4)

ggsave("plots/HillSoilMoistComparison_boxplots.jpg", soilvarsbox, height = 8, width = 12)

