library(WtRegDO)
library(doSNOW)
library(foreach)
library(tidyverse)
library(patchwork)
library(plotly)

outpth <- '~/Desktop/wtreg_out'

# run grid models for wtreg -------------------------------------------------------------------

dir.create(outpth, showWarnings = FALSE)

# import from file on Desktop
dat <- read.csv('~/Desktop/ML_2016.csv') |> 
  mutate(
    DateTimeStamp = lubridate::mdy_hm(DateTimeStamp, tz = 'America/Jamaica')
  ) |> 
  filter(!is.na(Tide))

# setup window width grid
n <- 5
grd <- expand.grid(days = seq(0.1, 12, length.out = n),
                   hrs  = seq(0.1, 12, length.out = n),
                   tide = seq(0.1, 1, length.out = n)
)

# metadata for the location
tz <- 'America/Jamaica'
lat <- 31.39
long <- -81.28

# setup parallel backend with doSNOW
ncores <- parallel::detectCores()
cl <- makeCluster(ncores - 10)
registerDoSNOW(cl)

# setup progress bar
pb <- txtProgressBar(max = nrow(grd), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# parallel loop using foreach with progress bar
res <- foreach(i = 1:nrow(grd),
               .options.snow = opts,
               .packages = c('WtRegDO', 'dplyr', 'tibble'),
               .export = c('dat', 'tz', 'lat', 'long')
               ) %dopar% {

  # get windows
  wins <- list(grd$days[i], grd$hrs[i], grd$tide[i])

  # weighted regression - turn off internal parallelization
  wtreg_res <- wtreg(dat, parallel = FALSE, wins = wins,
                    tz = tz, lat = lat, long = long)

  # save wtreg output to csv
  outfl <- file.path(outpth, 
    paste0(
      'wtreg_days', round(grd$days[i], 3),
      '_hrs', round(grd$hrs[i], 3),
      '_tide', round(grd$tide[i], 3),
      '.csv'
    )
  )
  write.csv(wtreg_res, file = outfl, row.names = FALSE)
                 
  # estimate ecosystem metabolism using observed DO time series
  metab_obs <- ecometab(wtreg_res, DO_var = 'DO_obs', tz = tz,
                       lat = lat, long = long)

  # estimate ecosystem metabolism using detided DO time series
  metab_dtd <- ecometab(wtreg_res, DO_var = 'DO_nrm', tz = tz,
                       lat = lat, long = long)

  # return results for this iteration
  list(
      obs = metab_obs,
      dtd = metab_dtd
    ) |>
    enframe() |>
    bind_cols(grd[i,])

}

# clean up
close(pb)
stopCluster(cl)

save(res, file = '~/Desktop/res.RData')

# eval observed, predicted, and detided DO ----------------------------------------------------

# get filenames
fls <- list.files(outpth, full.names = TRUE)

# pick which set of results to view
days <- 0.1
hrs <- 0.1
tide <- 0.1

# import the desired file, convert DateTimeStamp to POSIXct
toview <- fls[grepl(paste0('days', days, '_hrs', hrs, '_tide', tide, '\\.csv$'), fls)] |> 
  read.csv() |> 
  mutate(
    DateTimeStamp = case_when(
      !grepl('\\:', DateTimeStamp) ~ paste0(DateTimeStamp, ' 00:00:00'),
      TRUE ~ DateTimeStamp
    ),
    DateTimeStamp = lubridate::ymd_hms(DateTimeStamp, tz = 'America/Jamaica')
  )

# use plotly to show observed do (DO_obs) as points, predicted DO (DO_prd) as line, detided DO (DO_nrm) as line
plot_ly(toview, x = ~ DateTimeStamp) |> 
  add_markers(y = ~ DO_obs, name = 'Observed DO', marker = list(color = 'blue', size = 2)) |> 
  add_lines(y = ~ DO_prd, name = 'Predicted DO', line = list(color = 'red', width = 1)) |> 
  add_lines(y = ~ DO_nrm, name = 'Detided DO', line = list(color = 'green', width = 1)) |> 
  layout(yaxis = list(title = 'DO (mg/L)'),
         xaxis = list(title = 'DateTimeStamp'))

# eval monthy metab for observed and detided -----------------------------

load(file = '~/Desktop/res.RData')

# pick which set of results to view
days <- 0.1
hrs <- 0.1
tide <- 0.1

toview <- do.call('rbind', res) |>
  filter(days == !!days & hrs == !!hrs & tide == !!tide) |> 
  mutate(
    value = map(value, WtRegDO:::aggregate.metab, by = 'months')
  ) |> 
  unnest(cols = value)

ggplot(toview, aes(x = Date, y = val, color = name)) +
  geom_point(position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 5)) + 
  facet_wrap(~ Estimate, ncol = 1, scales = 'free_y') +
  labs(title = paste('Monthly metabolism for observed vs detided DO, days =', days, ', hrs =', hrs, ', tide =', tide),
       x = NULL, y = 'Metabolism (mmol O2 m-2 d-1)', color = 'DO type') +
  theme_minimal()

# eval metab results using criteria -----------------------------------------------------------

load(file = '~/Desktop/res.RData')

toeval <- do.call('rbind', res) |>
  group_by(days, hrs, tide) |>
  nest() |>
  mutate(
    objall = purrr::map(data, function(x){
      objfun(x$value[[1]], x$value[[2]])
    }),
    objnomean = purrr::map(data, function(x){
      objfun(x$value[[1]], x$value[[2]], vls = c('sdPg', 'anomPg', 'sdRt', 'anomRt'))
    }),
    objanom = purrr::map(data, function(x){
      objfun(x$value[[1]], x$value[[2]], vls = c('anomPg', 'anomRt'))
    }),
    objmean = purrr::map(data, function(x){
      objfun(x$value[[1]], x$value[[2]], vls = c('meanPg', 'meanRt'))
    }),
    objsd = purrr::map(data, function(x){
      objfun(x$value[[1]], x$value[[2]], vls = c('sdPg', 'sdRt'))
    })
  ) |>
  select(-data)

# plot objective fun results
objplo_fun <- function(toeval, obj, ttl){

  toplo <- toeval |>
    unnest(!!obj) |> 
    rename(obj = !!obj) |> 
    select(days, hrs, tide, obj) |> 
    ungroup()

  out <- ggplot(toplo, aes(x = days, y = hrs)) +
    geom_tile(aes(fill = obj), color = 'grey') +
    geom_tile(data = toplo |> filter(obj == min(obj)), fill = NA, color = 'red', width = 3, height = 3, linewidth = 1) +
    scale_fill_viridis_c() +
    scale_x_continuous(breaks = round(unique(toeval$days), 2), expand = c(0, 0)) +
    scale_y_continuous(breaks = round(unique(toeval$hrs), 2), expand = c(0, 0)) +
    facet_wrap(~ tide) +
    labs(title = paste('Objective function for metabolism from observed vs detided DO, ', ttl),
        subtitle = 'Grid search over window widths (days, hours, tide) for weighted regression',
        x = 'Days window', y = 'Hours window',
        fill = 'Objective\nfunction') +
    theme_minimal()
  
  return(out)
  
}

# get window width combo for minimum objective function, picks only one window comb if multiple min
minobj_fun <- function(res, toeval, obj){

  minobj <- toeval |>
    unnest(!!obj) |> 
    rename(obj = !!obj) |> 
    select(days, hrs, tide, obj) |> 
    ungroup() |> 
    pull(obj) |> 
    which.min()

  out <- res[[minobj[1]]]

  return(out)

}

# compare metab eval criteria for lowest obj, picks only one window comb if multiple min
checkeval_fun <- function(res, toeval, obj){

  tochk <- minobj_fun(res, toeval, obj)

  out <- list(
      obs = meteval(tochk[tochk$name == 'obs', 'value'][[1]][[1]])$cmp,
      dtd = meteval(tochk[tochk$name == 'dtd', 'value'][[1]][[1]])$cmp
    ) |> 
    enframe() |> 
    unnest('value')
  
  return(out)

}

# plot metab obs v dtd for lowest obj, picks only one window comb if multiple min
metplo_fun <- function(res, toeval, obj, ttl){

  tochk <- minobj_fun(res, toeval, obj)

  # observed DO plot
  p1 <- plot(tochk[tochk$name == 'obs', 'value'][[1]][[1]], by = 'days') + 
    labs(title = 'Observed DO metabolism')
  p2 <- plot(tochk[tochk$name == 'dtd', 'value'][[1]][[1]], by = 'days') + 
    coord_cartesian(ylim = range(p1$data$val, na.rm = T)) +
    labs(title = paste('Detided DO metabolism, ', ttl))

  out <- p1 + p2 + plot_layout(ncol = 1)

  return(out)

}

# all criteria
crit <- 'objall'
ttl <- 'all criteria'
objplo_fun(toeval, crit, ttl)
minobj_fun(res, toeval, crit)
checkeval_fun(res, toeval, crit)
metplo_fun(res, toeval, crit, ttl)

# anom and sd criteria
crit <- 'objnomean'
ttl <- 'no mean criteria'
objplo_fun(toeval, crit, ttl)
minobj_fun(res, toeval, crit)
checkeval_fun(res, toeval, crit)
metplo_fun(res, toeval, crit, ttl)

# anom criteria only
crit <- 'objanom'
ttl <- 'anomalous criteria only'
objplo_fun(toeval, crit, ttl)
minobj_fun(res, toeval, crit)
checkeval_fun(res, toeval, crit)
metplo_fun(res, toeval, crit, ttl)

# mean criteria only
crit <- 'objmean'
ttl <- 'mean criteria only'
objplo_fun(toeval, crit, ttl)
minobj_fun(res, toeval, crit)
checkeval_fun(res, toeval, crit)
metplo_fun(res, toeval, crit, ttl)

# sd criteria only
crit <- 'objsd'
ttl <- 'sd criteria only'
objplo_fun(toeval, crit, ttl)
minobj_fun(res, toeval, crit)
checkeval_fun(res, toeval, crit)
                 metplo_fun(res, toeval, crit, ttl)