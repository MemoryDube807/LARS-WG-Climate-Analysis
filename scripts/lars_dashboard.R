###############################################################################
# LARS-WG 8 projections for Mount Makulu: change from the observed baseline
#
# Reads the observed 1981-2010 station record and the 20 LARS-WG runs
# (5 GCMs x 2 SSPs x 2 periods), works out how each run differs from the
# observed climate, summarises the five models for each scenario, and writes
# dashboard_data.json plus three report figures.
#
# Run from the repository root:  Rscript scripts/lars_dashboard.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(jsonlite)
  library(ggplot2)
})

# ============================================================================
# 1. Settings
# ============================================================================
# Paths follow the LARS-WG install folder, with this repository in Output/
SCEN_DIR   <- "."                          # the 20 *WG.dat / *WG.st runs
OBS_FILE   <- "../Data/mmku.dat"           # observed daily record, 1981-2010
TST_FILE   <- "../Sitebase/mmku.tst"       # LARS-WG calibration test report
PVAL_FILE  <- "../Skill.PValues.txt"       # p-values of the calibration tests
OUT_JSON   <- "dashboard_data.json"
FIG_DIR    <- "figures"

MISSING    <- -99
HOT_DAY    <- 35    # daily maximum at or above this (degC) counts as a hot day
WARM_NIGHT <- 20    # daily minimum at or above this (degC) counts as a warm night
HEAVY_RAIN <- 30    # daily rain at or above this (mm) counts as a heavy rain day
SEASONS    <- c("NDJ", "FMA", "MJJ", "ASO")

SITE <- list(name = "Mount Makulu", code = "mmku",
             lat = -15.55, lon = 28.25, alt = 1200)

# ============================================================================
# 2. Readers
# ============================================================================
month_of_jday <- rep(1:12, times = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

read_obs <- function(path) {
  d <- read.table(path, col.names = c("year", "month", "day", "t_min", "t_max", "rain"))
  d[d == MISSING] <- NA
  d
}

# LARS-WG writes 365-day years numbered 1 to 100 with columns YEAR JDAY MIN MAX RAIN
read_wg <- function(path) {
  d <- read.table(path, col.names = c("year", "jday", "t_min", "t_max", "rain"))
  d$month <- month_of_jday[d$jday]
  d
}

parse_run_name <- function(path) {
  m <- regmatches(basename(path),
                  regexec("^mmku_(.+)\\[mmku,(ssp[0-9]+),([0-9]{4}-[0-9]{4})\\]WG\\.dat$",
                          basename(path)))[[1]]
  list(gcm = m[2], ssp = m[3], period = m[4])
}

assign_season <- function(month) {
  case_when(
    month %in% c(11, 12, 1) ~ "NDJ",
    month %in% c(2, 3, 4)   ~ "FMA",
    month %in% c(5, 6, 7)   ~ "MJJ",
    month %in% c(8, 9, 10)  ~ "ASO"
  )
}

# ============================================================================
# 3. Summaries for one daily series (observed or generated)
# ============================================================================
# Monthly rain totals only count when the month has no missing days. Temperatures use
# every available day. Annual and seasonal values are built from the monthly
# climatology, so a missing dry-season month in one year doesn't bias them.
summarise_series <- function(d) {
  monthly_totals <- d %>%
    group_by(year, month) %>%
    summarise(rain = sum(rain), .groups = "drop") %>%
    filter(!is.na(rain))

  clim <- d %>%
    group_by(month) %>%
    summarise(t_max = mean(t_max, na.rm = TRUE),
              t_min = mean(t_min, na.rm = TRUE), .groups = "drop") %>%
    left_join(monthly_totals %>% group_by(month) %>%
                summarise(rain = mean(rain, na.rm = TRUE), .groups = "drop"),
              by = "month") %>%
    arrange(month)

  seasonal <- clim %>%
    mutate(season = assign_season(month)) %>%
    group_by(season) %>%
    summarise(rain = sum(rain), t_max = mean(t_max), t_min = mean(t_min), .groups = "drop")
  annual <- tibble(season = "Year", rain = sum(clim$rain),
                   t_max = mean(clim$t_max), t_min = mean(clim$t_min))
  seasonal <- bind_rows(seasonal, annual)

  # Yearly counts, using only years with no missing temperature (or rain) data
  yearly <- d %>%
    group_by(year) %>%
    summarise(
      ok_t  = all(!is.na(t_max)) && all(!is.na(t_min)),
      ok_r  = all(!is.na(rain)),
      hot   = sum(t_max >= HOT_DAY, na.rm = TRUE),
      warm  = sum(t_min >= WARM_NIGHT, na.rm = TRUE),
      heavy = sum(rain >= HEAVY_RAIN, na.rm = TRUE),
      .groups = "drop")

  # Rainy-season totals, October to April, so each season stays in one piece
  rainy <- d %>%
    filter(month %in% c(10:12, 1:4)) %>%
    mutate(season_year = ifelse(month >= 10, year, year - 1)) %>%
    group_by(season_year) %>%
    summarise(rain = sum(rain), n = n(), .groups = "drop") %>%
    filter(n >= 212, !is.na(rain))   # 212 days from 1 Oct to 30 Apr (213 in leap years)

  list(
    monthly  = clim,
    seasonal = seasonal,
    hot_days     = mean(yearly$hot[yearly$ok_t]),
    warm_nights  = mean(yearly$warm[yearly$ok_t]),
    heavy_days   = mean(yearly$heavy[yearly$ok_r]),
    rainy_totals = rainy$rain
  )
}

# ============================================================================
# 4. Observed baseline and the 20 runs
# ============================================================================
obs <- summarise_series(read_obs(OBS_FILE))
cat(sprintf("Observed baseline: %.0f mm a year, Tmax %.2f, Tmin %.2f, %d rainy seasons\n",
            obs$seasonal$rain[obs$seasonal$season == "Year"],
            mean(obs$monthly$t_max), mean(obs$monthly$t_min), length(obs$rainy_totals)))

files <- list.files(SCEN_DIR, pattern = "^mmku_.*WG\\.dat$", full.names = TRUE)
stopifnot(length(files) == 20)

runs <- lapply(files, function(f) {
  meta <- parse_run_name(f)
  s <- summarise_series(read_wg(f))
  base <- obs$seasonal
  change <- s$seasonal %>%
    left_join(base, by = "season", suffix = c("", "_obs")) %>%
    transmute(season,
              rain_pct = 100 * (rain - rain_obs) / rain_obs,
              rain = rain - rain_obs,
              t_max = t_max - t_max_obs,
              t_min = t_min - t_min_obs)
  c(meta, list(
    monthly = s$monthly, seasonal = s$seasonal, change = change,
    hot_days = s$hot_days, warm_nights = s$warm_nights, heavy_days = s$heavy_days,
    rainy_q = unname(quantile(s$rainy_totals, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  ))
})

# ============================================================================
# 5. Ensemble summary: median and range across the five models
# ============================================================================
change_tbl <- bind_rows(lapply(runs, function(r)
  r$change %>% mutate(gcm = r$gcm, ssp = r$ssp, period = r$period)))

ensemble <- change_tbl %>%
  group_by(ssp, period, season) %>%
  summarise(across(c(rain, rain_pct, t_max, t_min),
                   list(med = median, min = min, max = max)),
            n_drier = sum(rain < 0), .groups = "drop")

cat("\nAnnual change from 1981-2010, median of five models (range):\n")
ensemble %>% filter(season == "Year") %>%
  mutate(txt = sprintf("%s %s  rain %+.0f mm (%+.0f to %+.0f), %d of 5 drier   Tmax %+.2f (%+.2f to %+.2f)   Tmin %+.2f (%+.2f to %+.2f)",
                       ssp, period, rain_med, rain_min, rain_max, n_drier,
                       t_max_med, t_max_min, t_max_max, t_min_med, t_min_min, t_min_max)) %>%
  pull(txt) %>% cat(sep = "\n")

ind_tbl <- bind_rows(lapply(runs, function(r)
  tibble(gcm = r$gcm, ssp = r$ssp, period = r$period, hot_days = r$hot_days,
         warm_nights = r$warm_nights, heavy_days = r$heavy_days,
         rainy_median = r$rainy_q[3], rainy_iqr = r$rainy_q[4] - r$rainy_q[2])))
cat("\nYearly counts, median of five models:\n")
print(ind_tbl %>% group_by(ssp, period) %>%
        summarise(across(c(hot_days, warm_nights, heavy_days, rainy_median, rainy_iqr), median), .groups = "drop"),
      width = 120)
cat(sprintf("Observed: %.1f hot days, %.1f warm nights, %.1f heavy rain days a year\n",
            obs$hot_days, obs$warm_nights, obs$heavy_days))

# ============================================================================
# 6. Calibration and validation against the 1981-2010 record
# ============================================================================
# LARS-WG was calibrated on the observed 1981-2010 series (site analysis), then
# validated by generating a synthetic series of the same length and testing it
# against the observations (QTest). The results sit in the .tst report:
# Kolmogorov-Smirnov tests on daily distributions and spell lengths, t-tests on
# monthly means and F-tests on the year-to-year variance of monthly means.
tst <- gsub("\r", "", readLines(TST_FILE, warn = FALSE))

num_row <- function(l) suppressWarnings(as.numeric(strsplit(trimws(l), "\t")[[1]]))

# Monthly block: obs mean, obs sd, gen mean, gen sd, t, p(t), F, p(F)
read_monthly <- function(header) {
  i <- grep(header, tst, fixed = TRUE)[1]
  r <- lapply(tst[(i + 1):(i + 8)], num_row)
  list(obs_mean = r[[1]], obs_sd = r[[2]], gen_mean = r[[3]], gen_sd = r[[4]],
       p_t = r[[6]], p_f = r[[8]])
}
# Daily extremes block: obs median, 95th percentile, maximum, then the same for gen
read_extremes <- function(header) {
  i <- grep(header, tst, fixed = TRUE)[1]
  r <- lapply(tst[(i + 1):(i + 6)], num_row)
  list(obs_median = r[[1]], obs_p95 = r[[2]], obs_max = r[[3]],
       gen_median = r[[4]], gen_p95 = r[[5]], gen_max = r[[6]])
}
# KS block: one line per season or month, p-value last. Months with no rain have no test.
read_ks <- function(header) {
  i <- grep(header, tst, fixed = TRUE)[1]
  j <- i + 1
  p <- c()
  while (j <= length(tst) && !startsWith(tst[j], "[")) {
    v <- num_row(tst[j]); v <- v[!is.na(v)]
    if (length(v) >= 3) p <- c(p, tail(v, 1))
    j <- j + 1
  }
  p
}

val_monthly <- list(
  rain  = read_monthly("[RAIN monthly mean & sd"),
  t_min = read_monthly("[MIN monthly mean & sd"),
  t_max = read_monthly("[MAX monthly mean & sd"))
val_extremes <- list(
  rain  = read_extremes("[RAIN daily maxima"),
  t_max = read_extremes("[MAX daily maxima"))
ks_spells <- read_ks("[KS-test for seasonal wet/dry SERIES")
ks_rain   <- read_ks("[KS-test for daily RAIN")
ks_min    <- read_ks("[KS-test for daily MIN")
ks_max    <- read_ks("[KS-test for daily MAX")

# Rain variance is only testable in months with rain in every year (not Jul and Aug)
rain_f <- val_monthly$rain$p_f[val_monthly$rain$obs_sd > 0]
test_tbl <- tibble(
  check = c("Wet and dry spell lengths", "Daily rainfall amounts", "Monthly mean rainfall",
            "Year-to-year variation in monthly rainfall",
            "Daily minimum temperature", "Monthly mean minimum temperature",
            "Year-to-year variation in monthly minimum",
            "Daily maximum temperature", "Monthly mean maximum temperature",
            "Year-to-year variation in monthly maximum"),
  test = c("KS", "KS", "t", "F", "KS", "t", "F", "KS", "t", "F"),
  p = list(ks_spells, ks_rain, val_monthly$rain$p_t, rain_f,
           ks_min, val_monthly$t_min$p_t, val_monthly$t_min$p_f,
           ks_max, val_monthly$t_max$p_t, val_monthly$t_max$p_f)) %>%
  mutate(n = lengths(p), passed = sapply(p, function(x) sum(x >= 0.01))) %>%
  select(-p)

# LARS-WG's own headline count: the KS and t-tests above, as listed in Skill.PValues.txt
pv <- read.table(PVAL_FILE, header = TRUE, sep = "\t", skip = 1, check.names = FALSE)
pcols <- 8:(8 + 4 + 4 + 12 * 6 - 1)   # spells, rain KS and t, Tmin KS and t, Tmax KS and t
pvals <- as.numeric(pv[1, pcols])

validation <- list(
  monthly = val_monthly,
  extremes = val_extremes,
  tests = test_tbl,
  tests_total = length(pvals), tests_failed = sum(pvals < 0.01),
  # Spread of October-April totals in the observations against the median spread
  # of the 20 runs, which share the calibrated parameters. Each generated month
  # is drawn largely on its own, so whole wet or dry seasons are rarer than in reality.
  rainy_iqr_obs = unname(diff(quantile(obs$rainy_totals, c(0.25, 0.75)))),
  rainy_iqr_runs = median(sapply(runs, function(r) r$rainy_q[4] - r$rainy_q[2]))
)
cat(sprintf("\nValidation: %d of %d LARS-WG tests failed at the 1%% level\n",
            validation$tests_failed, validation$tests_total))
print(test_tbl, width = 120)
cat(sprintf("Spread of Oct-Apr totals (interquartile range): observed %.0f mm, runs %.0f mm\n",
            validation$rainy_iqr_obs, validation$rainy_iqr_runs))

# ============================================================================
# 7. Export for the dashboard
# ============================================================================
out <- list(
  site = SITE,
  thresholds = list(hot_day = HOT_DAY, warm_night = WARM_NIGHT,
                    heavy_rain = HEAVY_RAIN),
  observed = list(
    monthly = obs$monthly %>% select(rain, t_max, t_min) %>% as.list(),
    seasonal = obs$seasonal,
    hot_days = obs$hot_days, warm_nights = obs$warm_nights, heavy_days = obs$heavy_days,
    rainy_q = unname(quantile(obs$rainy_totals, c(0.1, 0.25, 0.5, 0.75, 0.9))),
    n_rainy_seasons = length(obs$rainy_totals)),
  runs = lapply(runs, function(r) list(
    gcm = r$gcm, ssp = r$ssp, period = r$period,
    monthly = r$monthly %>% select(rain, t_max, t_min) %>% as.list(),
    seasonal = r$seasonal, change = r$change,
    hot_days = r$hot_days, warm_nights = r$warm_nights, heavy_days = r$heavy_days,
    rainy_q = r$rainy_q)),
  ensemble = ensemble,
  validation = validation
)
write_json(out, OUT_JSON, auto_unbox = TRUE, digits = 4, dataframe = "rows")
cat("Saved", OUT_JSON, "\n")

# ============================================================================
# 8. Report figures (small multiples)
# ============================================================================
dir.create(FIG_DIR, showWarnings = FALSE)
season_names <- c(NDJ = "Nov to Jan", FMA = "Feb to Apr", MJJ = "May to Jul",
                  ASO = "Aug to Oct", Year = "Whole year")
run_levels <- c("SSP2-4.5, 2021 to 2040", "SSP2-4.5, 2031 to 2050",
                "SSP5-8.5, 2021 to 2040", "SSP5-8.5, 2031 to 2050")
run_colours <- c("#E08E2B", "#E08E2B", "#B2182B", "#B2182B")
run_shapes  <- c(21, 19, 21, 19)    # open for 2021-2040, filled for 2031-2050
names(run_colours) <- names(run_shapes) <- run_levels
model_levels <- c("Median of five", rev(sort(unique(change_tbl$gcm))))   # median plots at the bottom

label_run <- function(ssp, period) {
  paste0(ifelse(ssp == "ssp245", "SSP2-4.5", "SSP5-8.5"), ", ",
         sub("-", " to ", period))
}

# Projected change for every model plus the five-model median
dot_df <- bind_rows(
  change_tbl,
  change_tbl %>% group_by(ssp, period, season) %>%
    summarise(across(c(rain, rain_pct, t_max, t_min), median), .groups = "drop") %>%
    mutate(gcm = "Median of five")) %>%
  mutate(run = factor(label_run(ssp, period), levels = run_levels),
         gcm = factor(gcm, levels = model_levels),
         season = factor(season_names[season], levels = season_names))

dot_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold"),
        plot.title.position = "plot")
dot_layers <- list(
  geom_vline(xintercept = 0, colour = "grey45"),
  geom_point(aes(colour = run, shape = run), size = 2.6, stroke = 1.1),
  scale_colour_manual(values = run_colours, name = NULL),
  scale_shape_manual(values = run_shapes, name = NULL),
  guides(colour = guide_legend(nrow = 2)))

# 8a. Temperature: one panel per season, daily maximum above daily minimum
p_temp <- dot_df %>%
  select(gcm, run, season, t_max, t_min) %>%
  pivot_longer(c(t_max, t_min), names_to = "var", values_to = "change") %>%
  mutate(var = factor(ifelse(var == "t_max", "Daily maximum", "Daily minimum"),
                      levels = c("Daily maximum", "Daily minimum"))) %>%
  ggplot(aes(x = change, y = gcm)) + dot_layers +
  facet_grid(var ~ season) +
  labs(x = "Projected change from the observed 1981 to 2010 mean (\u00b0C)", y = NULL,
       title = "Projected temperature change at Mount Makulu, by season and model",
       subtitle = "Each row is a climate model run through LARS-WG; the bottom row is the median of the five") +
  dot_theme
ggsave(file.path(FIG_DIR, "temp_change_by_model.png"), p_temp, width = 12, height = 6.5, dpi = 200)

# 8b. Rainfall: one panel per season (May to July holds about 3 mm and is left out)
p_rain <- dot_df %>%
  filter(season != "May to Jul") %>%
  ggplot(aes(x = rain, y = gcm)) + dot_layers +
  facet_wrap(~ season, nrow = 1) +
  labs(x = "Projected change from the observed 1981 to 2010 total (mm)", y = NULL,
       title = "Projected rainfall change at Mount Makulu, by season and model",
       subtitle = "Each row is a climate model run through LARS-WG; the bottom row is the median of the five") +
  dot_theme
ggsave(file.path(FIG_DIR, "rain_change_by_model.png"), p_rain, width = 12, height = 4.8, dpi = 200)

# 8c. Monthly rainfall, observed against projected 2031-2050, one panel per model
obs_month <- tibble(month = 1:12, rain = obs$monthly$rain)
run_month <- bind_rows(lapply(runs, function(r)
  tibble(gcm = r$gcm, ssp = r$ssp, period = r$period, month = 1:12, rain = r$monthly$rain)))
run_month <- bind_rows(
  run_month,
  run_month %>% group_by(ssp, period, month) %>%
    summarise(rain = median(rain), .groups = "drop") %>% mutate(gcm = "Median of five")) %>%
  filter(period == "2031-2050") %>%
  mutate(run = ifelse(ssp == "ssp245", "SSP2-4.5, 2031 to 2050", "SSP5-8.5, 2031 to 2050"),
         gcm = factor(gcm, levels = c(sort(unique(change_tbl$gcm)), "Median of five")))
obs_panel <- tidyr::crossing(gcm = factor(levels(run_month$gcm), levels = levels(run_month$gcm)),
                             obs_month)
p_month <- ggplot() +
  geom_col(data = obs_panel, aes(x = month, y = rain, fill = "Observed 1981 to 2010"), width = 0.75) +
  geom_line(data = run_month, aes(x = month, y = rain, colour = run), linewidth = 0.9) +
  geom_point(data = run_month, aes(x = month, y = rain, colour = run), size = 1.3) +
  facet_wrap(~ gcm, ncol = 3) +
  scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1)) +
  scale_fill_manual(values = c("Observed 1981 to 2010" = "#C9D1D8"), name = NULL) +
  scale_colour_manual(values = c("SSP2-4.5, 2031 to 2050" = "#E08E2B",
                                 "SSP5-8.5, 2031 to 2050" = "#B2182B"), name = NULL) +
  labs(x = NULL, y = "Mean monthly rainfall (mm)",
       title = "Monthly rainfall at Mount Makulu, observed 1981 to 2010 and projected for 2031 to 2050") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.title.position = "plot")
ggsave(file.path(FIG_DIR, "rain_monthly_by_model.png"), p_month, width = 11, height = 6.5, dpi = 200)

# 8d. Validation: observed and generated 1981-2010 monthly means, with +/- 1 standard
#     deviation of the year-to-year values as a band
val_df <- bind_rows(lapply(names(val_monthly), function(v) {
  m <- val_monthly[[v]]
  bind_rows(
    tibble(var = v, month = 1:12, source = "Observed 1981 to 2010", mean = m$obs_mean, sd = m$obs_sd),
    tibble(var = v, month = 1:12, source = "Generated by LARS-WG", mean = m$gen_mean, sd = m$gen_sd))
})) %>%
  mutate(var = factor(c(rain = "Rainfall (mm per month)", t_max = "Daily maximum (\u00b0C)",
                        t_min = "Daily minimum (\u00b0C)")[var],
                      levels = c("Rainfall (mm per month)", "Daily maximum (\u00b0C)", "Daily minimum (\u00b0C)")))
p_val <- ggplot(val_df, aes(x = month, y = mean, colour = source, fill = source)) +
  geom_ribbon(aes(ymin = ifelse(var == "Rainfall (mm per month)", pmax(mean - sd, 0), mean - sd), ymax = mean + sd),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.4) +
  facet_wrap(~ var, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1)) +
  scale_colour_manual(values = c("Observed 1981 to 2010" = "#4A5C66", "Generated by LARS-WG" = "#2A6F97"), name = NULL) +
  scale_fill_manual(values = c("Observed 1981 to 2010" = "#4A5C66", "Generated by LARS-WG" = "#2A6F97"), name = NULL) +
  labs(x = NULL, y = NULL,
       title = "Validation: the calibrated generator against the observed 1981 to 2010 record",
       subtitle = "Lines are monthly means; bands show one standard deviation of the year-to-year values") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.title.position = "plot")
ggsave(file.path(FIG_DIR, "validation_monthly.png"), p_val, width = 12, height = 4.6, dpi = 200)

cat("Saved figures to", FIG_DIR, "\n")
