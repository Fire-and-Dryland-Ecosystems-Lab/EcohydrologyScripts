# MSR_6_plot

library(RHESSysIOinR)
library(rhutils)
library(tidyverse)
library(fasstr)
library(cowplot)

out_dir = "output/rh_out_2023-06-10--04-55-35/"
treat_date = "2005-03-1"

DT = get_basin_daily(out_dir)
DT_mn = basin_daily2mn(DT)

# problem_vars = c("cs.net_psn", "Kstar_potential_both", "cdf.psn_to_cpool")
treat_year = as.POSIXlt(treat_date)$year + 1900

DT$Treatment = "NA"
DT[date < treat_date, "Treatment"] = "Pre Treatment"
DT[date >= treat_date, "Treatment"] = "Post Treatment"
DT$Treatment = factor(DT$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

DT = DT %>% group_by(run) %>% mutate(`Cumulative Streamflow` = cumsum(streamflow))

DT$run[DT$run == " Ward_msrbaseline"] = "Baseline"
DT$run[DT$run == "Ward_msrthin2"] = "Prescribed Fire"
DT$run[DT$run == "Ward_msrthin3"] = "Heavy Thinning"
DT$run[DT$run == "Ward_msrthin4"] = "Moderate Thinning"
DT$run[DT$run == "Ward_msrthin5"] = "Mastication"
DT$run[DT$run == "Ward_msrthin6"] = "Clearcut"

# thin_3yr = DT %>% group_by(run) %>% filter(date > "1989-10-1" & date <= "1992-10-1") %>% summarise(sum(streamflow))
# thin_3yr$`sum(streamflow)`[2]/thin_3yr$`sum(streamflow)`[1]
# thin_after3yr = DT %>% group_by(run) %>% filter(date > "1992-10-1") %>% summarise(sum(streamflow))
# thin_after3yr$`sum(streamflow)`[2]/thin_after3yr$`sum(streamflow)`[1]

# streamflow timing
DT_Q_time = calc_annual_flow_timing(data = DT, dates = "date", values = "streamflow", groups = "run")
DT_Q_time$Treatment = "NA"
DT_Q_time[DT_Q_time$Year < treat_year, "Treatment"] = "Pre Treatment"
DT_Q_time[DT_Q_time$Year >= treat_year, "Treatment"] = "Post Treatment"
DT_Q_time$Treatment = factor(DT_Q_time$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

# add water year month
DT_mn$month_wy = ifelse(DT_mn$month >= 10, DT_mn$month - 9, DT_mn$month + 3)
wy_mn = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")
DT_mn$month_wy_str = factor(wy_mn[DT_mn$month_wy], levels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"))

DT_mn$Treatment = "NA"
DT_mn[year_month < "Oct 2005", "Treatment"] = "Pre Treatment"
DT_mn[year_month >= "Oct 2005", "Treatment"] = "Post Treatment"
DT_mn$Treatment = factor(DT_mn$Treatment, levels = c("Pre Treatment", "Post Treatment", "NA"))

DT_mn$run[DT_mn$run == " Ward_msrbaseline"] = "Baseline"
DT_mn$run[DT_mn$run == "Ward_msrthin2"] = "Prescribed Fire"
DT_mn$run[DT_mn$run == "Ward_msrthin3"] = "Heavy Thinning"
DT_mn$run[DT_mn$run == "Ward_msrthin4"] = "Moderate Thinning"
DT_mn$run[DT_mn$run == "Ward_msrthin5"] = "Mastication"
DT_mn$run[DT_mn$run == "Ward_msrthin6"] = "Clearcut"

theme_set(theme_cowplot())

Q_prepost = DT_mn %>% 
  group_by(run, Treatment, month_wy) %>% 
  summarise(streamflow = mean(streamflow)) %>%
  ggplot() + aes(x = month_wy, y = streamflow*1000, color = as.factor(run), linetype = Treatment) + 
  geom_line(linewidth=1) + ggtitle("Average Monthly Streamflow") + ylab("Streamflow (mm)") + 
  labs(color = "Scenario") +
  scale_x_continuous(labels = wy_mn, breaks = c(1:12)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 315))


Q_prepost2 = DT_mn %>% 
  group_by(run, Treatment, month_wy) %>%
  filter(year_month < "Mar 2010") %>%
  summarise(streamflow = mean(streamflow)) %>%
  ggplot() + aes(x = month_wy, y = streamflow*1000, color = as.factor(run), linetype = Treatment) + 
  geom_line(linewidth=1) + ggtitle("Average Monthly Streamflow") + ylab("Streamflow (mm)") + 
  labs(color = "Scenario") +
  scale_x_continuous(labels = wy_mn, breaks = c(1:12)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 315))


Q_prepost2DOY = DT %>% 
  group_by(run, Treatment, yd) %>%
  filter(date < "2010-03-1") %>%
  summarise(streamflow = mean(streamflow)) %>%
  ggplot() + aes(x = yd, y = streamflow*1000, color = as.factor(run), linetype = Treatment) + 
  geom_line(linewidth=1) + ggtitle("Average Monthly Streamflow") + ylab("Streamflow (mm)") + 
  labs(color = "Scenario") +
  scale_x_continuous(labels = wy_mn, breaks = c(1:12)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 315))


DT_mn %>%
 filter(year >= 2000L & year <= 2010L) %>%
 ggplot() +
 aes(x = "", y = lai, fill = run) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_minimal() +
 facet_wrap(vars(Treatment))




Q_timeseries = DT %>% 
  filter(date >= "2004-01-1" & date < "2015-01-1") %>%
  ggplot() + aes(x = date, y = streamflow*1000, color = as.factor(run)) + 
  geom_line(linewidth=1) + ggtitle("Streamflow Following Treatment") + ylab("Streamflow (mm)") + 
  geom_vline(xintercept = as.numeric(as.Date("1989-10-01")), alpha = 0.5, size = 1) +
  labs(color = "Scenario") + xlab("Date")

qtiming_line = DT_Q_time %>%
  filter(Treatment %in% "Post Treatment") %>%
  ggplot() +
  aes(x = Year, y = DoY_50pct_TotalQ, colour = run) +
  geom_line(linewidth=1) +
  scale_color_hue(direction = 1) +
  ggtitle("Streamflow Timing") +
  ylab("Day of Year (50th Percentile)") +
  labs(color = "Scenario")


leg = g_legend(Q_prepost)

#rel_heights = c(3,3,1,1),

Q_quad = plot_grid(Q_prepost + theme(legend.position="none"),
                   Q_timeseries+ theme(legend.position="none"),
                   qtiming_line + theme(legend.position="none"),
                   leg, align = "v", nrow = 2)

# ggsave(filename = "plots/Q_quad.png", plot = Q_quad, width = 12, height = 9, dpi = 500, scale = 0.75)


library(dplyr)
library(ggplot2)


library(dplyr)
library(ggplot2)

DT %>%
 filter(year >= 2005L & year <= 2007L) %>%
 ggplot() +
 aes(x = date, y = `Cumulative Streamflow`, colour = run) +
 geom_line() +
 scale_color_hue(direction = 1) +
 theme_minimal()




Q_cumsum = DT %>% 
  ggplot() + aes(x = date, y = `Cumulative Streamflow`*1000, fill = as.factor(run)) + 
  geom_line(linewidth=1) + ggtitle("Streamflow") + ylab("Streamflow (mm)") + 
  labs(color = "Scenario") + xlab("Date")

qtiming_box = DT_Q_time %>%
  filter(Treatment %in% "Post Treatment") %>%
  ggplot() +
  aes(x = "", y = DoY_50pct_TotalQ, fill = run) +
  geom_boxplot() +
  labs(color = "Scenario") +
  ggtitle("Streamflow Timing") +
  scale_fill_hue(direction = 1) +
  ylab("Day of Year (50th Percentile)") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

qtiming_density = DT_Q_time %>%
  filter(Treatment %in% "Post Treatment") %>%
  ggplot() +
  aes(x = DoY_50pct_TotalQ, fill = run) +
  geom_density(adjust = 1L) +
  scale_fill_hue(direction = 1) +
  ggtitle("Streamflow Timing") +
  ylab("Day of Year (50th Percentile)") +
  labs(fill = "Scenario")


ggplot(DT_mn) +
  aes(x = ym_ind, y = totalc, colour = run) +
  geom_line() +
  scale_color_hue(direction = 1) 

ggplot(DT_mn) +
  aes(x = ym_ind, y = streamflow, colour = run) +
  geom_line() +
  scale_color_hue(direction = 1) 

ggplot(DT_mn) +
  aes(x = ym_ind, y = evaporation, colour = run) +
  geom_line() +
  scale_color_hue(direction = 1) 

ggplot(DT_mn) +
  aes(x = ym_ind, y = cs.net_psn, colour = run) +
  geom_line() +
  scale_color_hue(direction = 1)


ggplot(DT_mn) +
  aes(x = "", y = streamflow, fill = run) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal()



ggplot(DT_mn) +
 aes(x = "", y = totalc, fill = run, colour = Treatment) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 scale_color_hue(direction = 1) +
 theme_minimal()


ggplot(DT_mn) +
  aes(x = "", y = cs.net_psn, fill = run) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal()

