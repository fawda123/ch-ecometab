# this script will import multiple year files for Sapelo, combine, and run WtRegDO on year, window with combos

source("R/f_calcWan09.R")
source("R/ecometab_HW.R")

# ---- Libraries ----
# remotes::install_github("fawda123/WtRegDO")
# Load plyr before tidyverse to avoid conflicts, and let furrr handle parallelism
library(plyr)        # ddply, join (used in monthly evals)
library(tidyverse)   # dplyr, readr, tidyr, lubridate
library(WtRegDO)
library(future)
library(furrr)
library(lubridate)   # explicitly included for clarity

# ---- Read input & prep ----

# file paths
fls <- list.files("data", pattern = "^Sapelo", full.names = TRUE)

# import as list, drop NA values, check for DateTimeStamp column
dat_list <- lapply(fls, function(f) {
  cat("\n📄 Loading data from:", f, "\n")
  fl <- readr::read_csv(f, show_col_types = FALSE) %>%
    tidyr::drop_na()

  if (!("DateTimeStamp" %in% names(fl))) {
    stop("❌ No 'DateTimeStamp' column found in ", basename(f))
  }

  cat("\n✅ Files loaded successfully:", basename(f), "\n")
  cat("Columns detected:", paste(names(fl), collapse = ", "), "\n")
  
  return(fl)
  
})

# combine, add Year column
dat <- bind_rows(dat_list) %>% 
  mutate(
    Year = year(DateTimeStamp)
  )

# ---- Detect and correct timezone interpretation ----
dt_col <- dat$DateTimeStamp
cat("\n⏱  Checking DateTimeStamp format...\n")

if (inherits(dt_col, "POSIXt")) {
  tz_current <- attr(dt_col, "tzone")
  if (is.null(tz_current) || tz_current == "" || tz_current == "UTC") {
    cat("🕓 Detected POSIXct in UTC — reinterpreting as local EST (America/Jamaica)...\n")
    # The times are *actually* local (EST) but labeled UTC, so relabel without shifting
    dat$DateTimeStamp <- lubridate::force_tz(dt_col, "America/Jamaica")
  } else {
    cat("🕓 Detected POSIXct with tz =", tz_current, " — keeping as is.\n")
  }
} else if (is.character(dt_col)) {
  if (grepl("T\\d{2}:\\d{2}:\\d{2}Z$", dt_col[1])) {
    cat("🌎 Detected ISO8601 UTC format — parsing and converting to local EST...\n")
    dat$DateTimeStamp <- lubridate::ymd_hms(dt_col, tz = "UTC")
    dat$DateTimeStamp <- lubridate::with_tz(dat$DateTimeStamp, "America/Jamaica")
  } else if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}", dt_col[1])) {
    cat("📅 Detected m/d/yyyy format — parsing as local EST (America/Jamaica)...\n")
    dat$DateTimeStamp <- lubridate::mdy_hm(dt_col, tz = "America/Jamaica")
  } else {
    stop("❌ Could not determine timestamp format. Check CSV manually.")
  }
} else {
  stop("❌ Unrecognized DateTimeStamp type: ", class(dt_col)[1])
}

# ---- Sanity checks ----
if (any(is.na(dat$DateTimeStamp))) {
  warning("⚠️ Some timestamps could not be parsed — check CSV format.")
}

dup_count <- sum(duplicated(dat$DateTimeStamp))
if (dup_count > 0) {
  warning("⚠️ Found ", dup_count, " duplicated timestamps.")
}

cat("✅ Final DateTimeStamp class:", class(dat$DateTimeStamp)[1], "\n")
cat("🗓  First timestamp:", format(min(dat$DateTimeStamp, na.rm = TRUE), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("🗓  Last timestamp: ", format(max(dat$DateTimeStamp, na.rm = TRUE), "%Y-%m-%d %H:%M:%S %Z"), "\n")

# ---- Directories (ABSOLUTE PATHS to avoid "no output" when using parallel workers) ----
out_dir <- normalizePath("results",
                         winslash = "/", mustWork = TRUE)
setwd(out_dir)

# ---- Assign cleaned data ----
dat_nonan <- dat

# ---- Settings ----
tz   <- 'America/Jamaica'
lat  <- 31.42
long <- -81.30

# ---- Outputs ----
log_path     <- file.path(out_dir, "sapelo_metabolism_model_log.txt")
summary_path <- file.path(out_dir, "sapelo_summary_metabolism_results.csv")

readr::write_lines("", log_path)  # initialize log

summary_template <- tibble::tibble(
  obs_Pg_mean = numeric(), obs_Pg_sd = numeric(), obs_Pg_pct_neg = numeric(),
  obs_Rt_mean = numeric(), obs_Rt_sd = numeric(), obs_Rt_pct_pos = numeric(),
  obs_NEM_mean = numeric(), obs_NEM_sd = numeric(),
  dtd_Pg_mean = numeric(), dtd_Pg_sd = numeric(), dtd_Pg_pct_neg = numeric(),
  dtd_Rt_mean = numeric(), dtd_Rt_sd = numeric(), dtd_Rt_pct_pos = numeric(),
  dtd_NEM_mean = numeric(), dtd_NEM_sd = numeric(),
  DOcor_mean = numeric(), Pgcor_mean = numeric(), Rtcor_mean = numeric(),
  meanPg_mon = numeric(), sdPg_mon = numeric(), anomPg_mon = numeric(),
  meanRt_mon = numeric(), sdRt_mon = numeric(), anomRt_mon = numeric(),
  yr = numeric(), w1 = numeric(), w2 = numeric(), w3 = numeric(), tag = character()
)
readr::write_csv(summary_template, summary_path)

# ---- Parallel plan ----
plan(multisession, workers = max(1, parallel::detectCores() - 1))

# helper for logging
log_msg <- function(msg) {
  timestamp <- as.character(Sys.time())
  readr::write_lines(paste(timestamp, "-", msg), log_path, append = TRUE)
}

# ---- Compute OBS metabolism ONCE ----
# metabolism for the entire time series
metab_obs <- ecometab_HW(dat_nonan, DO_var = 'DO_obs', tz = tz, lat = lat, long = long,
                         depth_vec = 4.4, replacemet = TRUE,
                         gasex = 'Wanninkhof', gasave = 'instant')

# Save for reference
readr::write_csv(metab_obs, file.path(out_dir, "sapelo_metab_obs_reference.csv"))

# add Year for subsetting
metab_obs$Year <- year(metab_obs$Date)

# ---- Per-combination worker ----
process_combo <- function(yr, w1, w2, w3) {
  combo_tag <- paste0("yr_", yr, "_w1_", w1, "_w2_", w2, "_w3_", w3)
  log_msg(paste("Starting", combo_tag))
  
  # get the year-specific observed metabolism summary
  summary_obs_yr <- metab_obs[metab_obs$Year == yr, ]
  summary_obs_yr <- dplyr::summarise(
    summary_obs_yr,
    obs_Pg_mean = mean(Pg, na.rm = TRUE),
    obs_Pg_sd   = sd(Pg, na.rm = TRUE),
    obs_Pg_pct_neg = mean(Pg < 0, na.rm = TRUE) * 100,
    obs_Rt_mean   = mean(Rt, na.rm = TRUE),
    obs_Rt_sd     = sd(Rt, na.rm = TRUE),
    obs_Rt_pct_pos = mean(Rt > 0, na.rm = TRUE) * 100,
    obs_NEM_mean = mean(NEM, na.rm = TRUE),
    obs_NEM_sd   = sd(NEM, na.rm = TRUE)
  )

  safe_row <- function() {
    out <- dplyr::bind_cols(summary_obs_yr, tibble::tibble(
      dtd_Pg_mean = NA_real_, dtd_Pg_sd = NA_real_, dtd_Pg_pct_neg = NA_real_,
      dtd_Rt_mean = NA_real_, dtd_Rt_sd = NA_real_, dtd_Rt_pct_pos = NA_real_,
      dtd_NEM_mean = NA_real_, dtd_NEM_sd = NA_real_,
      DOcor_mean = NA_real_, Pgcor_mean = NA_real_, Rtcor_mean = NA_real_,
      meanPg_mon = NA_real_, sdPg_mon = NA_real_, anomPg_mon = NA_real_,
      meanRt_mon = NA_real_, sdRt_mon = NA_real_, anomRt_mon = NA_real_,
      yr = yr, w1 = w1, w2 = w2, w3 = w3, tag = combo_tag
    ))
    out %>% dplyr::select(tidyselect::all_of(names(summary_template)))
  }
  
  tryCatch({
    # Weighted regression and detided metabolism for THIS combo

    # subset input data by year
    dat_nonan_yr <- dat_nonan[dat_nonan$Year == yr, ]
    
    # run weighted regression
    wtreg_res <- wtreg(dat_nonan_yr, wins = list(w1, w2, w3),
                      progress = FALSE, tz = tz, lat = lat, long = long,
                      parallel = FALSE)

    # calculate metabolism from weighted regression output
    metab_dtd <- ecometab_HW(wtreg_res, DO_var = 'DO_nrm', tz = tz, lat = lat, long = long,
                            depth_vec = 4.4, replacemet = TRUE,
                            gasex = 'Wanninkhof', gasave = 'instant')
        
    # Save full outputs (absolute paths => safe in parallel workers)
    readr::write_csv(wtreg_res, file.path(out_dir, paste0("wtreg_res_sapelo_", combo_tag, ".csv")))
    readr::write_csv(metab_dtd, file.path(out_dir, paste0("metab_dtd_sapelo_", combo_tag, ".csv")))
    
    # ---- Monthly & correlation analysis (integrated from meteval.r) ----
    toeval    <- stats::na.omit(metab_dtd)
    rawdat    <- attr(metab_dtd, 'rawdat')
    depth_val <- attr(metab_dtd, 'depth_val')
    DO_var    <- 'DO_nrm'
    
    if (is.null(rawdat) || is.null(depth_val) || is.null(DO_var)) {
      log_msg(paste("WARN: missing rawdat/depth_val/DO_var for", combo_tag))
      DOcor_mean <- Pgcor_mean <- Rtcor_mean <- NA_real_
      meanPg_mon <- sdPg_mon <- anomPg_mon <- NA_real_
      meanRt_mon <- sdRt_mon <- anomRt_mon <- NA_real_
    } else {
      rawdat$month <- strftime(rawdat$Date, '%m')
      
      DOcor <- plyr::ddply(
        rawdat, .variables = c('month'),
        .fun = function(x) suppressWarnings(cor.test(x[, DO_var], x[, depth_val])$estimate)
      )
      names(DOcor)[2] <- 'DOcor'
      
      tide_rngs <- plyr::ddply(rawdat, .variables = c('metab_date'), .fun = function(x) {
        sunrise <- mean(diff(x[x$solar_period %in% 'sunrise', 'Tide']), na.rm = TRUE)
        sunset  <- mean(diff(x[x$solar_period %in% 'sunset',  'Tide']), na.rm = TRUE)
        if (is.infinite(sunrise)) sunrise <- NA_real_
        if (is.infinite(sunset))  sunset  <- NA_real_
        daytot  <- mean(diff(x$Tide, na.rm = TRUE))
        c(daytot, sunrise, sunset)
      })
      names(tide_rngs) <- c('Date','daytot','sunrise','sunset')
      
      toeval <- merge(toeval, tide_rngs, by = 'Date', all.x = TRUE)
      toeval$month <- strftime(toeval$Date, '%m')
      
      metcor <- plyr::ddply(
        toeval, .variables = c('month'),
        .fun = function(x) with(x, c(
          Pgcor = suppressWarnings(try(cor.test(Pg, sunrise)$estimate, silent = TRUE)),
          Rtcor = suppressWarnings(try(cor.test(Rt, sunset)$estimate, silent = TRUE))
        ))
      )
      
      PgRtmon <- plyr::ddply(
        toeval, .variables = c('month'),
        .fun = function(x) with(x, c(
          meanPg = mean(Pg, na.rm = TRUE),
          sdPg   = sd(Pg,   na.rm = TRUE),
          anomPg = 100 * sum(Pg <= 0) / nrow(x),
          meanRt = mean(Rt, na.rm = TRUE),
          sdRt   = sd(Rt,   na.rm = TRUE),
          anomRt = 100 * sum(Rt >= 0) / nrow(x)
        ))
      )
      
      mos <- plyr::join(DOcor, metcor, by = 'month')
      mos <- plyr::join(mos, PgRtmon, by = 'month')
      
      # Save per-month table for this combo
     readr::write_csv(mos, file.path(out_dir, paste0("sapelo_", combo_tag, ".csv")))
      
      # Monthly means for the summary row
      DOcor_mean <- mean(mos$DOcor, na.rm = TRUE)
      Pgcor_mean <- mean(mos$Pgcor, na.rm = TRUE)
      Rtcor_mean <- mean(mos$Rtcor, na.rm = TRUE)
      
      meanPg_mon <- mean(mos$meanPg, na.rm = TRUE)
      sdPg_mon   <- mean(mos$sdPg,   na.rm = TRUE)
      anomPg_mon <- mean(mos$anomPg, na.rm = TRUE)
      
      meanRt_mon <- mean(mos$meanRt, na.rm = TRUE)
      sdRt_mon   <- mean(mos$sdRt,   na.rm = TRUE)
      anomRt_mon <- mean(mos$anomRt, na.rm = TRUE)
    }
    
    # ---- DTD summaries for this combo ----
    summary_dtd <- toeval %>% dplyr::summarise(
      dtd_Pg_mean = mean(Pg, na.rm = TRUE),
      dtd_Pg_sd   = sd(Pg,   na.rm = TRUE),
      dtd_Pg_pct_neg = mean(Pg < 0, na.rm = TRUE) * 100,
      dtd_Rt_mean = mean(Rt, na.rm = TRUE),
      dtd_Rt_sd   = sd(Rt,   na.rm = TRUE),
      dtd_Rt_pct_pos = mean(Rt > 0, na.rm = TRUE) * 100,
      dtd_NEM_mean = mean(NEM, na.rm = TRUE),
      dtd_NEM_sd   = sd(NEM,   na.rm = TRUE)
    )
    
    # ---- Assemble single summary row ----
    summary_row <- dplyr::bind_cols(summary_obs_yr, summary_dtd) %>%
      dplyr::mutate(
        DOcor_mean = DOcor_mean, Pgcor_mean = Pgcor_mean, Rtcor_mean = Rtcor_mean,
        meanPg_mon = meanPg_mon, sdPg_mon = sdPg_mon, anomPg_mon = anomPg_mon,
        meanRt_mon = meanRt_mon, sdRt_mon = sdRt_mon, anomRt_mon = anomRt_mon,
        yr = yr, w1 = w1, w2 = w2, w3 = w3, tag = combo_tag
      ) %>%
      dplyr::select(tidyselect::all_of(names(summary_template)))
    
    log_msg(paste("Finished", combo_tag))
    return(summary_row)
    
  }, error = function(e) {
    log_msg(paste("ERROR in", combo_tag, ":", conditionMessage(e)))
    return(safe_row())
  })
}

# ---- Build grid and run in parallel ----

yr_values <- unique(dat_nonan$Year)
w1_values <- seq(0.1,12, length.out=5)
w2_values <- seq(0.1,12, length.out=5)
w3_values <- seq(0.1, 1.0, length.out = 5)

param_grid <- tidyr::expand_grid(yr = yr_values, w1 = w1_values, w2 = w2_values, w3 = w3_values)
all_rows <- furrr::future_pmap(param_grid, ~process_combo(..1, ..2, ..3, ..4),
                               .options = furrr::furrr_options(seed = TRUE))
summary_all <- dplyr::bind_rows(all_rows)

# Write once to avoid parallel-write issues
readr::write_csv(summary_all, summary_path, append = TRUE)

# ---- Graceful shutdown ----
plan(sequential)


