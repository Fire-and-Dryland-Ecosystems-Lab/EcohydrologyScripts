# MSR_6_plot

library(RHESSysIOinR)
library(rhutils)
library(tidyverse)
library(fasstr)
library(cowplot)
library(zoo)

options(scipen = 999)

out_dir = "output/rh_out_2023-08-28--17-25-31"
# treat_date = "2005-03-1"
treat_date = "1999-8-1"
treat_yrmn = as.yearmon(treat_date)

# --------------------------------------------------------------
# -------------------- Add vars and aggregate  --------------------
# --------------------------------------------------------------
treat_year = as.POSIXlt(treat_date)$year + 1900

DT = get_basin_daily(out_dir)

DT$run[DT$run == paste0(site,"_msrbaseline")] = "Baseline"
DT$run[DT$run == paste0(site,"_msrthin2")] = "Prescribed Fire"
DT$run[DT$run == paste0(site,"_msrthin3")] = "Heavy Thinning"
DT$run[DT$run == paste0(site,"_msrthin4")] = "Moderate Thinning"
DT$run[DT$run == paste0(site,"_msrthin5")] = "Mastication"
DT$run[DT$run == paste0(site,"_msrthin6")] = "Clearcut"

DT_mn = basin_daily2mn(DT)

vars = names(DT)[!names(DT) %in% c("day","month","year","basinID", "date","run", "wy","yd")]
DT_yr = DT %>% group_by(run, wy) %>% summarise(across(all_of(vars), mean))

DT$Treatment = "NA"
DT[date < treat_date, "Treatment"] = "Pre Treatment"
DT[date >= treat_date, "Treatment"] = "Post Treatment"
DT$Treatment = factor(DT$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

# add water year month
DT_mn$month_wy = ifelse(DT_mn$month >= 10, DT_mn$month - 9, DT_mn$month + 3)
wy_mn = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")
DT_mn$month_wy_str = factor(wy_mn[DT_mn$month_wy], levels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"))

DT_mn$Treatment = "NA"
DT_mn[year_month < as.yearmon(treat_date), "Treatment"] = "Pre Treatment"
DT_mn[year_month >= as.yearmon(treat_date), "Treatment"] = "Post Treatment"
DT_mn$Treatment = factor(DT_mn$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

# DT_yr$Treatment = "NA"
# DT_yr[wy <= treat_year, "Treatment"] = "Pre Treatment"
# DT_yr[wy > treat_year, "Treatment"] = "Post Treatment"
# DT_yr$Treatment = factor(DT_yr$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

# streamflow timing
DT_Q_time = calc_annual_flow_timing(data = DT, dates = "date", values = "streamflow", groups = "run")
DT_Q_time$Treatment = "NA"
DT_Q_time[DT_Q_time$Year < treat_year, "Treatment"] = "Pre Treatment"
DT_Q_time[DT_Q_time$Year >= treat_year, "Treatment"] = "Post Treatment"
DT_Q_time$Treatment = factor(DT_Q_time$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

DT_yr$Evapotranspiration = DT_yr$evaporation + DT_yr$evaporation_surf + DT_yr$transpiration_sat_zone + DT_yr$transpiration_unsat_zone

DT_yr = DT_yr %>% rename("NPP" = "cs.net_psn","LAI"= "lai","Plant Carbon" = "totalc", "Streamflow"="streamflow","Stem Carbon" = "cs.live_stemc","Root Carbon" = "cs.live_crootc")

DT_yr$Evapotranspiration = DT_yr$Evapotranspiration * 365
DT_yr$NPP = DT_yr$NPP * 365
DT_yr$`Soil Moisture` = ((DT_yr$rz_storage + DT_yr$unsat_storage )/DT_yr$sat_deficit)
DT_yr$Streamflow = DT_yr$Streamflow * 365

# --------------------------------------------------------------
# -------------------- Carbon --------------------
# --------------------------------------------------------------

vars1 = c("Plant Carbon", "Stem Carbon", "Root Carbon", "LAI")

plotdata = DT_yr %>% select(eval(vars1), run, wy) %>% pivot_longer(cols = eval(vars1), names_to = "variable")
tmpplot = ggplot(data=plotdata, aes(x=wy, y=value, col=run)) + # here water year
  geom_line(aes(linetype= run), lwd=1.3) +
  facet_wrap(~variable, scales="free") +
  scale_x_continuous(breaks=scales::pretty_breaks(n = 8)) + xlab("Water Year") + ggtitle(site)

ggsave("figures/quad_plant.jpg", tmpplot, height = 8, width = 12, scale = 0.75)


# --------------------------------------------------------------
# -------------------- Streamflow --------------------
# --------------------------------------------------------------

vars1 = c("Streamflow", "Soil Moisture", "Evapotranspiration", "NPP")

plotdata = DT_yr %>% select(eval(vars1), run, wy) %>% pivot_longer(cols = eval(vars1), names_to = "variable")
tmpplot = ggplot(data=plotdata, aes(x=wy, y=value, col=run)) + # here water year
  geom_line(aes(linetype= run), lwd=1.3) +
  facet_wrap(~variable, scales="free") +
  scale_x_continuous(breaks=scales::pretty_breaks(n = 8)) + xlab("Water Year") + ggtitle(site)

ggsave("figures/quad_streamflow.jpg", tmpplot, height = 8, width = 12, scale = 0.75)

# 
# 
# Q_prepost = DT_mn %>% 
#   group_by(run, Treatment, month_wy) %>% 
#   summarise(streamflow = mean(streamflow)) %>%
#   ggplot() + aes(x = month_wy, y = streamflow*1000, color = as.factor(run), linetype = Treatment) + 
#   geom_line(linewidth=1) + ggtitle("Average Monthly Streamflow") + ylab("Streamflow (mm)") + 
#   labs(color = "Scenario") +
#   scale_x_continuous(labels = wy_mn, breaks = c(1:12)) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 315))
# 
# Q_prepost2 = DT_mn %>% 
#   group_by(run, Treatment, month_wy) %>%
#   filter(year_month < treat_yrmn) %>%
#   summarise(streamflow = mean(streamflow)) %>%
#   ggplot() + aes(x = month_wy, y = streamflow*1000, color = as.factor(run), linetype = Treatment) + 
#   geom_line(linewidth=1) + ggtitle("Average Monthly Streamflow") + ylab("Streamflow (mm)") + 
#   labs(color = "Scenario") +
#   scale_x_continuous(labels = wy_mn, breaks = c(1:12)) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 315))
# 
# Q_prepost2DOY = DT %>% 
#   group_by(run, Treatment, yd) %>%
#   filter(date < "2010-03-1") %>%
#   summarise(streamflow = mean(streamflow)) %>%
#   ggplot() + aes(x = yd, y = streamflow*1000, color = as.factor(run), linetype = Treatment) + 
#   geom_line(linewidth=1) + ggtitle("Average Monthly Streamflow") + ylab("Streamflow (mm)") + 
#   labs(color = "Scenario") +
#   scale_x_continuous(labels = wy_mn, breaks = c(1:12)) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 315))
# 
# 
# Q_timeseries = DT %>% 
#   filter(date >= "2004-01-1" & date < "2015-01-1") %>%
#   ggplot() + aes(x = date, y = streamflow*1000, color = as.factor(run)) + 
#   geom_line(linewidth=1) + ggtitle("Streamflow Following Treatment") + ylab("Streamflow (mm)") + 
#   geom_vline(xintercept = as.numeric(as.Date("1989-10-01")), alpha = 0.5, linewidth=1) +
#   labs(color = "Scenario") + xlab("Date")
# 
# qtiming_line = DT_Q_time %>%
#   filter(Treatment %in% "Post Treatment") %>%
#   ggplot() +
#   aes(x = Year, y = DoY_50pct_TotalQ, colour = run) +
#   geom_line(linewidth=1) +
#   scale_color_hue(direction = 1) +
#   ggtitle("Streamflow Timing") +
#   ylab("Day of Year (50th Percentile)") +
#   labs(color = "Scenario")
# 
# 
# 
# 
# Q_cumsum = DT %>% 
#   ggplot() + aes(x = date, y = `Cumulative Streamflow`*1000, color = as.factor(run)) + 
#   geom_line(linewidth=1) + ggtitle("Streamflow") + ylab("Streamflow (mm)") + 
#   labs(color = "Scenario") + xlab("Date")
# 
# qtiming_box = DT_Q_time %>%
#   filter(Treatment %in% "Post Treatment") %>%
#   ggplot() +
#   aes(x = "", y = DoY_50pct_TotalQ, fill = run) +
#   geom_boxplot() +
#   labs(color = "Scenario") +
#   ggtitle("Streamflow Timing") +
#   scale_fill_hue(direction = 1) +
#   ylab("Day of Year (50th Percentile)") +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# qtiming_density = DT_Q_time %>%
#   filter(Treatment %in% "Post Treatment") %>%
#   ggplot() +
#   aes(x = DoY_50pct_TotalQ, fill = run) +
#   geom_density(adjust = 1L) +
#   scale_fill_hue(direction = 1) +
#   ggtitle("Streamflow Timing") +
#   ylab("Day of Year (50th Percentile)") +
#   labs(fill = "Scenario")
# 
# 
# leg = g_legend(Q_prepost)
# Q_quad = plot_grid(Q_prepost + theme(legend.position="none"),
#                    Q_timeseries+ theme(legend.position="none"),
#                    qtiming_line + theme(legend.position="none"),
#                    leg, align = "v", nrow = 2)
# 
# ggsave(filename = "plots/Q_quad.png", plot = Q_quad, width = 12, height = 9, dpi = 500, scale = 0.75)
# 
# # --------------------------------------------------------------
# # -------------------- Carbon --------------------
# # --------------------------------------------------------------
# 
# totalc = DT_mn %>% filter(year_month >= "Feb 2005") %>% 
#   ggplot() +
#   aes(x = ym_ind, y = totalc, colour = run) +
#   geom_line(size = 1) +
#   ggtitle("Total Carbon") +
#   scale_color_hue(direction = 1)
# leafc  = DT_mn %>% filter(year_month >= "Feb 2005") %>% 
#   ggplot() +
#   aes(x = ym_ind, y = cs.leafc, colour = run) +
#   geom_line(size = 1) +
#   ggtitle("Leaf Carbon") +
#   scale_color_hue(direction = 1)
# stemc = DT_mn %>% filter(year_month >= "Feb 2005") %>% 
#   ggplot() +
#   aes(x = ym_ind, y = cs.dead_stemc+cs.live_stemc, colour = run) +
#   geom_line(size = 1) +
#   ggtitle("Stem Carbon") +
#   scale_color_hue(direction = 1)
# rootc = DT_mn %>% filter(year_month >= "Feb 2005") %>% 
#   ggplot() +
#   aes(x = ym_ind, y = cs.live_crootc+cs.dead_crootc+cs.frootc, colour = run) +
#   geom_line(size = 1) +
#   ggtitle("Root Carbon") +
#   scale_color_hue(direction = 1)
# 
# leg = g_legend(totalc)
# 
# carbon_quad = plot_grid(plot_grid(totalc + theme(legend.position="none"),
#                    leafc+ theme(legend.position="none"),
#                    stemc + theme(legend.position="none"),
#                    rootc + theme(legend.position="none"), align = "v", nrow = 2),leg, ncol = 2,rel_widths = c(4,1))
# 
# ggsave(filename = "plots/carbon_quad.png", plot = carbon_quad, width = 12, height = 9, dpi = 500, scale = 0.75)
# 
# psn = ggplot(DT_mn) +
#   aes(x = ym_ind, y = cs.net_psn, colour = run) +
#   geom_line(size = 1) +
#   scale_color_hue(direction = 1)


