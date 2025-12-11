library(WtRegDO)
library(tidyverse)
library(patchwork)

# SAPDC July -------------------------------------------------------------

data(SAPDC)
tz <- 'America/Jamaica'
lat <- 31.39
long <- -81.28
winls <- list(3.075, 0.1, 0.1)

tomod <- SAPDC[grepl('\\-07\\-', SAPDC$DateTimeStamp),]

resnosine <- wtreg(tomod, tz = tz, lat = lat, long = long, wins = winls)
reswtsine <- wtreg(tomod, tz = tz, lat = lat, long = long, wins = winls, sine = TRUE)

toplo <- list(
    `No Sine` = resnosine,
    `With Sine` = reswtsine
  ) |> 
  enframe() |> 
  unnest('value')

p1 <- ggplot(toplo, aes(x = DateTimeStamp, y = DO_obs)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_line(aes(y = DO_prd, color = name), size = 1) +
  labs(
    title = 'Weighted Regression with and without Sine Wave',
    subtitle = 'Model predicted',
    y = 'Dissolved Oxygen (mg/L)',
    x = NULL,
    color = 'Model'
  ) +
  theme_minimal()

p2 <- ggplot(toplo, aes(x = DateTimeStamp, y = DO_obs)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_line(aes(y = DO_nrm, color = name), size = 1) +
  labs(
    subtitle = 'Model normalized',
    y = 'Dissolved Oxygen (mg/L)',
    x = NULL,
    color = 'Model'
  ) +
  theme_minimal()

p <- p1 + p2 + plot_layout(ncol = 1, guides = 'collect', axes = 'collect_y') &
  theme(legend.position = 'bottom')

png('~/Desktop/SAPDC_wtreg_sine_vs_nosine.png', width = 10, height = 8, units = 'in', res = 300)
print(p)
dev.off()

# ML 2016 July -----------------------------------------------------------

# import from file on Desktop
ML2016 <- read.csv('~/Desktop/ML_2016.csv') |> 
  mutate(
    DateTimeStamp = lubridate::mdy_hm(DateTimeStamp, tz = 'America/Jamaica')
  ) |> 
  filter(!is.na(Tide))

tz <- 'America/Jamaica'
lat <- 31.39
long <- -81.28
winls <- list(3.075, 0.1, 0.1)

tomod <- ML2016[grepl('\\-07\\-', ML2016$DateTimeStamp),]

resnosine <- wtreg(tomod, tz = tz, lat = lat, long = long, wins = winls)
reswtsine <- wtreg(tomod, tz = tz, lat = lat, long = long, wins = winls, sine = TRUE)

toplo <- list(
    `No Sine` = resnosine,
    `With Sine` = reswtsine
  ) |> 
  enframe() |> 
  unnest('value')

p1 <- ggplot(toplo, aes(x = DateTimeStamp, y = DO_obs)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_line(aes(y = DO_prd, color = name), size = 1) +
  labs(
    title = 'Weighted Regression with and without Sine Wave',
    subtitle = 'Model predicted',
    y = 'Dissolved Oxygen (mg/L)',
    x = NULL,
    color = 'Model'
  ) +
  theme_minimal()

p2 <- ggplot(toplo, aes(x = DateTimeStamp, y = DO_obs)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_line(aes(y = DO_nrm, color = name), size = 1) +
  labs(
    subtitle = 'Model normalized',
    y = 'Dissolved Oxygen (mg/L)',
    x = NULL,
    color = 'Model'
  ) +
  theme_minimal()

p <- p1 + p2 + plot_layout(ncol = 1, guides = 'collect', axes = 'collect_y') &
  theme(legend.position = 'bottom')

png('~/Desktop/ML2016_wtreg_sine_vs_nosine.png', width = 10, height = 8, units = 'in', res = 300)
print(p)
dev.off()