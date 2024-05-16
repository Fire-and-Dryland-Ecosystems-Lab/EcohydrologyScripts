# MSR hillslope plots

library(RHESSysIOinR)
library(rhutils)
library(tidyverse)
library(fasstr)
library(cowplot)
library(zoo)
source("R/0_global_vars.R")

treat_date = "1999-8-1"
treat_yrmn = as.yearmon(treat_date)

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

data.merge = merge(hill_year[hill_year$Scenario != "Baseline", ], hill_year[hill_year$Scenario == "Baseline",], by=c("wy", "hillID", "basinID"), all.x=T)

data.merge = data.merge %>% mutate(streamflow_diff = streamflow.x - streamflow.y) %>%
  mutate(et_diff = et.x - et.y) %>%
  mutate(rz_s_diff = theta.x - theta.y)

aridity = hill_year[hill_year$Scenario == "Baseline",] %>% group_by(hillID, Scenario) %>%
  summarise_if(is.numeric, mean, na.rm = T) %>%
  dplyr::select(PET, et, rain_thru, theta, hillID, Scenario) %>%
  as.data.frame() %>% mutate(aridity = PET/rain_thru) %>% mutate(rz_s = theta/365) %>%
  dplyr::select(-c(Scenario, theta))

summary(aridity)

# --------------------------------------------------------------
# -------------------- Plots --------------------
# --------------------------------------------------------------

data.merge.first = data.merge %>% dplyr::filter(wy >=1999 & wy <=2002) %>%
  group_by(Scenario.x, hillID) %>% summarise_if(is.numeric, mean, na.rm=T) %>%
  as.data.frame() %>% 
  mutate(streamflow_diff = streamflow.x - streamflow.y) %>%
  mutate(et_diff = et.x - et.y) %>%
  mutate(rz_s_diff = theta.x - theta.y)

data.merge.first2 = merge(data.merge.first, aridity, by="hillID", all.x=T)

arid1 = ggplot(data.merge.first2) +
  aes(x = aridity, y = streamflow_diff) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Scenario.x), nrow = 1) + 
  geom_hline(yintercept = 0) + ggtitle("First 3 years (1999 - 2002)") + ylim(-1,7)

rzs1 = ggplot(data.merge.first2) +
  aes(x = rz_s, y = streamflow_diff) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Scenario.x), nrow = 1) +
  geom_hline(yintercept = 0) + ggtitle("First 3 years (1999 - 2002)")+ ylim(0,7)

# ------------------------------

data.merge.first = data.merge %>% dplyr::filter(wy >=2003 & wy <=2009) %>%
  group_by(Scenario.x, hillID) %>% summarise_if(is.numeric, mean, na.rm=T) %>%
  as.data.frame() %>% 
  mutate(streamflow_diff = streamflow.x - streamflow.y) %>%
  mutate(et_diff = et.x - et.y) %>%
  mutate(rz_s_diff = theta.x - theta.y)

data.merge.first2 = merge(data.merge.first, aridity, by="hillID", all.x=T)

arid2 = ggplot(data.merge.first2) +
  aes(x = aridity, y = streamflow_diff) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Scenario.x), nrow = 1) + 
  geom_hline(yintercept = 0) + ggtitle("Next 5 years (2003 - 2009)")+ ylim(-1,7)

rzs2 = ggplot(data.merge.first2) +
  aes(x = rz_s, y = streamflow_diff) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Scenario.x), nrow = 1) + 
  geom_hline(yintercept = 0) + ggtitle("Next 5 years (2003 - 2009)")+ ylim(0,7)

# ------------------------------
year.start = 2011
year.end = 2020

data.merge.first = data.merge %>% dplyr::filter(wy >=year.start & wy <=year.end) %>%
  group_by(Scenario.x, hillID) %>% summarise_if(is.numeric, mean, na.rm=T) %>%
  as.data.frame() %>% 
  mutate(streamflow_diff = streamflow.x - streamflow.y) %>%
  mutate(et_diff = et.x - et.y) %>%
  mutate(rz_s_diff = theta.x - theta.y)

data.merge.first2 = merge(data.merge.first, aridity, by="hillID", all.x=T)

arid3 = ggplot(data.merge.first2) +
  aes(x = aridity, y = streamflow_diff) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Scenario.x), nrow = 1) +
  geom_hline(yintercept = 0) + ggtitle("Next 10 years (2011 - 2020)")+ ylim(-1,7)

rzs3 = ggplot(data.merge.first2) +
  aes(x = rz_s, y = streamflow_diff) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Scenario.x), nrow = 1) + 
  geom_hline(yintercept = 0) + ggtitle("Next 10 years (2011 - 2020)")+ ylim(0,7)

# --------------------
# year.start = 2021
# year.end = 2030
# 
# data.merge.first = data.merge %>% dplyr::filter(wy >=year.start & wy <=year.end) %>%
#   group_by(Scenario.x, hillID) %>% summarise_if(is.numeric, mean, na.rm=T) %>%
#   as.data.frame() %>% 
#   mutate(streamflow_diff = streamflow.x - streamflow.y) %>%
#   mutate(et_diff = et.x - et.y) %>%
#   mutate(rz_s_diff = theta.x - theta.y)
# 
# data.merge.first2 = merge(data.merge.first, aridity, by="hillID", all.x=T)
# 
# ggplot(data.merge.first2) +
#   aes(x = aridity, y = streamflow_diff) +
#   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
#   theme_minimal() +
#   facet_wrap(vars(Scenario.x), nrow = 1) + 
#   geom_hline(yintercept = 0) + ggtitle("Next 10 years (2021 - 2030)")
# 
# ggplot(data.merge.first2) +
#   aes(x = rz_s, y = streamflow_diff) +
#   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
#   theme_minimal() +
#   facet_wrap(vars(Scenario.x), nrow = 1) + 
  # geom_hline(yintercept = 0) + ggtitle("Next 10 years (2021 - 2030)")


aridplots = plot_grid(arid1, arid2, arid3, ncol=1)

rzsplots = plot_grid(rzs1, rzs2, rzs3, ncol=1)

ggsave("figures/aridplots.jpg",aridplots, height = 12, width = 8, scale = 0.75)
ggsave("figures/rzs.jpg",rzsplots, height = 12, width = 8, scale = 0.75)

tmpsm = hill_mn_DT %>% group_by(year, run) %>% summarise("Mean Soil Moisture" = mean(sm))


ggplot(tmpsm) +
  aes(x = year, y = `Mean Soil Moisture`, color = run) +
  geom_line(linewidth = 1)

hill_mn_DT %>%
  filter(run %in% "Baseline") %>%
  ggplot() +
  aes(x = hillID, y = sm, group = hillID) +
  geom_boxplot() +
  facet_wrap(vars(Treatment))

hill_mn_DT %>%
  filter(hillID == 6 & run == "Baseline") %>%
  ggplot() +
  aes(x = yr_mn, y = sm) +
  geom_line(linewidth = 1)


hill_mn_DT %>%
  filter(hillID == 52 & run == "Baseline") %>%
  ggplot() +
  aes(x = yr_mn, y = sm) +
  geom_line(linewidth = 1)


summary(hillvegid_mn_DT[hillvegid_mn_DT$hillID == 6])
summary(hillvegid_mn_DT[hillvegid_mn_DT$hillID == 52])


hillvegid_mn_DT %>%
  filter(hillID == 52 & run == "Baseline") %>%
  ggplot() +
  aes(y = veg_parm_ID) +
  geom_histogram()


ggplot(hillvegid_mn_DT) +
  aes(x = "", y = cs.live_stemc, fill = run) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(veg_parm_ID))



library(ggplot2)

ggplot(hill_mn_DT) +
 aes(x = soilmoist, y = streamflow) +
 geom_point(shape = "circle", size = 1.5, colour = "#112446") +
 theme_minimal() +
 facet_wrap(vars(run))


# Soil moisture pre treatment
tmpsm = hill_mn_DT %>% group_by(hillID) %>% filter(yr_mn < "Mar 2005") %>% summarise("Mean Soil Moisture" = mean(soilmoist))
hill_mn_DT = merge(hill_mn_DT, tmpsm, by = "hillID")
# ET/PET pre treatment
tmpet = hill_mn_DT %>% group_by(hillID) %>% filter(yr_mn < "Mar 2005") %>% summarise("AET/PET" = mean(et/PET))
hill_mn_DT = merge(hill_mn_DT, tmpet, by = "hillID")

hill_mn_DT %>%
ggplot() +
  aes(x = `AET/PET`, y = streamflow) +
  geom_point(shape = "circle", size = 1.5, colour = Treatment) +
  theme_minimal() +
  facet_wrap(vars(run))


hill_5yr_DT = hill_mn_DT %>% filter(year<2005) %>% group_by(hillID, run) %>% summarise("streamflow" = mean(streamflow))
hill_5yr_DT$Treatment = "Pre Treatment"
hill_5yr_DT2 = hill_mn_DT %>% filter(year>=2005 & year < 2010) %>% group_by(hillID, run) %>% summarise("streamflow" = mean(streamflow))
hill_5yr_DT2$Treatment = "Post Treatment"
hill_5yr_DTcomb = rbind(hill_5yr_DT, hill_5yr_DT2)

hill_5yr_DTcomb = merge(hill_5yr_DTcomb, tmpet, by = "hillID")
hill_5yr_DTcomb = merge(hill_5yr_DTcomb, tmpsm, by = "hillID")


library(ggplot2)

ggplot(hill_5yr_DTcomb) +
 aes(x = `AET/PET`, y = streamflow, colour = Treatment) +
 geom_point(shape = "circle", 
 size = 2) +
 scale_color_hue(direction = 1) +
 theme_minimal() +
 facet_wrap(vars(run))

ggplot(hill_5yr_DTcomb) +
  aes(x = `Mean Soil Moisture`, y = streamflow, colour = Treatment) +
  geom_point(shape = "circle", 
             size = 2) +
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(run))



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


tmp = hill_mn_DT %>%
  ggplot() +
  aes(x = yr_mn, y = unsat_storage, group = hillID, color = hillID) +
  geom_line(linewidth = 1)

tmp2 = hill_mn_DT %>% 
  ggplot() +
  aes(x = yr_mn, y = sat_deficit, group = hillID, color = hillID) +
  geom_line(linewidth = 1)

tmp3 = hill_mn_DT %>% 
  ggplot() +
  aes(x = yr_mn, y = rz_storage, group = hillID, color = hillID) +
  geom_line(linewidth = 1)

soilvars2 = plot_grid(tmp, tmp2, tmp3)


tmp = hill_mn_DT %>%
  ggplot() +
  aes(x = yr_mn, y = unsat_storage, group = yr_mn) +
  geom_boxplot()

tmp2 = hill_mn_DT %>% 
  ggplot() +
  aes(x = yr_mn, y = sat_deficit, group = yr_mn) +
  geom_boxplot()

tmp3 = hill_mn_DT %>% 
  ggplot() +
  aes(x = yr_mn, y = rz_storage, group = yr_mn) +
  geom_boxplot()

soilvarsbox = plot_grid(tmp, tmp2, tmp3)

