# spatial_get_basin_points

library(tidyverse)
library(terra)
library(sf)
library(cowplot)

basin = trim(rast("preprocessing/spatial90m/basin.tif"))
basin = project(basin, "EPSG:3310")
plot(basin)
basin_agg = aggregate(basin, fact = 8)
plot(basin_agg)

basin_points = as.points(basin_agg)
plot(basin_points)
basin_points

points_df = as.data.frame(geom(basin_points)[,c("geom","x","y")])
points_df

writetable = F
if (writetable) {
  write.csv(points_df, file = "../data/Ward_grid_epsg3310.csv")
}

DEM = trim(rast("preprocessing/spatial90m/dem.tif"))
DEM = project(DEM, "EPSG:3310")
plot(DEM)
DEM_agg = aggregate(DEM, fact = 8)
plot(DEM_agg)

DEM_points = as.points(DEM_agg)
plot(DEM_points)
DEM_points
DEM_df = as.data.frame(geom(DEM_points)[,c("geom","x","y")])
DEM_df$elevation = unlist(unname(values(DEM_points)))

# USE THOSE POINTS
input_ml_gridmet = "../data/20230427_will_results.csv"
gridmet_ml_df = read.csv(input_ml_gridmet)
gridmet_ml_df$X = NULL
gridmet_ml_df = gridmet_ml_df[gridmet_ml_df$quantile == 50,]

dfmerge = merge(gridmet_ml_df, DEM_df)
dfmerge$yearmon = zoo::as.yearmon(paste0(dfmerge$year,"-",dfmerge$month))

# dumb tidy doesnt play well with lm
df_pr = dfmerge[dfmerge$variable == "pr",]
df_tmin = dfmerge[dfmerge$variable == "tmin",]
df_tmax = dfmerge[dfmerge$variable == "tmax",]

df_ym = data.frame(yearmon= unique(dfmerge$yearmon))
df_ym$pr_int = NA
df_ym$pr_slope = NA
df_ym$tmin_int = NA
df_ym$tmin_slope = NA
df_ym$tmax_int = NA
df_ym$tmax_slope = NA


for (i in seq_along(df_ym$yearmon) ) {
  coefs = df_pr %>% filter(yearmon == df_ym$yearmon[i]) %>% lm(prediction~elevation, data = .) %>% .$coefficients %>% unname()
  df_ym$pr_int[i] = coefs[1]
  df_ym$pr_slope[i] = coefs[2]
  coefs = df_tmin %>% filter(yearmon == df_ym$yearmon[i]) %>% lm(prediction~elevation, data = .) %>% .$coefficients %>% unname()
  df_ym$tmin_int[i] = coefs[1]
  df_ym$tmin_slope[i] = coefs[2]
  coefs = df_tmax %>% filter(yearmon == df_ym$yearmon[i]) %>% lm(prediction~elevation, data = .) %>% .$coefficients %>% unname()
  df_ym$tmax_int[i] = coefs[1]
  df_ym$tmax_slope[i] = coefs[2]
}

df_ym$month = as.integer(format(as.POSIXlt(df_ym$yearmon), "%m"))
df_ym$year = as.integer(format(as.POSIXlt(df_ym$yearmon), "%Y"))

avg_pr_slope <- df_ym %>% summarise(avg_pr_slope = mean(pr_slope))
avg_tmin_slope <- df_ym %>% summarise(avg_tmin_slope = mean(tmin_slope))
avg_tmax_slope <- df_ym %>% summarise(avg_tmax_slope = mean(tmax_slope))

fig_pr = df_ym %>%
ggplot() +
  #geom_point(aes(x = month, y = pr_slope,color = factor(year)) ) +
  geom_line(aes(x = month, y = pr_slope,color = factor(year)), linewidth = 1 ) +
  geom_text(data = avg_pr_slope, aes(x = max(df_ym$month), y = max(df_ym$pr_slope), 
                                     label = paste("Avg pr_slope:", round(avg_pr_slope, 5))), hjust = 1, vjust = 1) +
  labs(title = "Precipitation Lapse Rates", x = "Month", y="Reg Slope" ) +
  theme(legend.title=element_blank())

fig_tmin = df_ym %>%
  ggplot() +
  #geom_point(aes(x = month, y = pr_slope,color = factor(year)) ) +
  geom_line(aes(x = month, y = tmin_slope,color = factor(year)), linewidth = 1 ) +
  geom_text(data = avg_tmin_slope, aes(x = max(df_ym$month), y = max(df_ym$tmin_slope), 
                                     label = paste("Avg tmin_slope:", round(avg_tmin_slope, 5))), hjust = 1, vjust = 1) +
  labs(title = "Min Temperature Lapse Rates", x = "Month", y="Reg Slope" ) +
  theme(legend.title=element_blank())

fig_tmax = df_ym %>%
  ggplot() +
  #geom_point(aes(x = month, y = pr_slope,color = factor(year)) ) +
  geom_line(aes(x = month, y = tmax_slope,color = factor(year)), linewidth = 1 ) +
  geom_text(data = avg_tmax_slope, aes(x = max(df_ym$month), y = max(df_ym$tmax_slope), 
                                       label = paste("Avg tmax_slope:", round(avg_tmax_slope, 5))), hjust = 1, vjust = 1) +
  labs(title = "Min Temperature Lapse Rates", x = "Month", y="Reg Slope" ) +
  theme(legend.title=element_blank())

fig_lapse = plot_grid(fig_pr,fig_tmin,fig_tmax, ncol = 1)

ggsave("plots/lapse_rates_gridmetml.jpg", fig_lapse, width = 10, height = 12)
