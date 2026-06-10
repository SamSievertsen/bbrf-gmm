# =============================================================================
# 01_build_dataset.R
# Build the analytic person-day dataset for the GMM from raw BBRF sources.
#
# This is a clean, config-driven re-implementation of the data wrangle
# (originally by Dani Y. Del Rubin & Scott A. Jones). It reproduces the substantive
# logic, i.e., joining the daily EMA survey, the daily go/no-go summaries, and
# actigraphy data,but with modular helpers, dynamic output naming, and a single
# config switch for the actigraphy algorithm (Sadeh <-> Cole-Kripke).
#
# STATUS: v1. The numerical/EM core (gmm_em.R) is unit-tested; THIS script has
# not yet been run against the real files and will need a local pass to confirm
# exact column names, encodings, and the relevant (e.g., no-actigraphy) edge cases. 
# The sleep/nap timing correction (Section 4) is the highest-risk logic and is the 
# first thing to validate against known cases.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr);
  library(tidyr);
  library(readr);
  library(stringr);
  library(lubridate);
  library(hms);
  library(zoo);
  library(purrr);
  library(glue);
  library(here)
})

source(here::here("src", "utils.R"))

cfg <- load_config()
set.seed(cfg$gmm$seed)
options(scipen = 999)

#0.1 Resolve source paths from config
RAW <- cfg$paths$data_raw
actig_dir <- file.path(RAW, cfg$sources$actigraphy_root, cfg$actigraphy$algorithm)


## 1. Loaders (one small, testable function per source) ##

#1.1 Load the long-format EMA survey (one row per completed survey/person-day)
load_ema_long <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    
    #1.1.1 Standardize the known REDCap-side typo so downstream names are clean
    dplyr::rename_with(~ stringr::str_replace(., "^ACITVESUICIDALTHOUGHTS$",
                                              "ACTIVESUICIDALTHOUGHTS"))
}

#1.2 Load and stack the per-subject daily go/no-go summaries.
#    Subject ID is inferred from the file name; rows carry the join keys
#    sessionid (-> instance_id) and subjectrspid (-> rsp_id).
load_ema_gonogo <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) stop("No EMA go/no-go files matched '", pattern, "' in ", dir)
  purrr::map_dfr(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE) %>%
      dplyr::mutate(source_file = basename(f))
  })
}

#1.3 Load and stack per-subject actigraphy. Skips the metadata header rows,
#    derives the subject ID from the file name, and concatenates multi-part
#    exports (e.g., "SleepBD27 part 1" + "part 2") into a single subject.
load_actigraphy <- function(dir, skip) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No actigraphy CSVs found in ", dir)
  purrr::map_dfr(files, function(f) {
    
    #1.3.1 Subject = filename minus extension, multi-part suffixes, and algo tag
    subj <- basename(f) %>%
      stringr::str_remove("\\.csv$") %>%
      stringr::str_remove("(?i)\\s*part\\s*\\d+") %>%
      stringr::str_remove("(?i)\\s*Sadeh") %>%
      stringr::str_trim()
    readr::read_csv(f, skip = skip, show_col_types = FALSE) %>%
      dplyr::mutate(subject = subj)
  })
}


## 2. Actigraphy cleaning + EMA-day alignment ##

#2.1 Standardize actigraphy column names and build in-bed midpoint timestamps.
#    Each actigraphy row is a sleep period; we align it to the EMA "day" of the
#    SUBSEQUENT wake period using the midpoint of the in-bed interval.
clean_actigraphy <- function(actig) {
  actig %>%
    dplyr::rename(
      sleep_onsetdate = `Onset Date`,
      sleep_onsettime = `Onset Time`,
      sleep_latency = `Latency`,
      sleep_totmoves = `Total Counts`,
      sleep_eff = `Efficiency`,
      sleep_totalmin = `Total Sleep Time (TST)`,
      awakenings = `Number of Awakenings`,
      wakelengths = `Average Awakening Length`,
      sleep_fragidx = `Sleep Fragmentation Index`,
      bedtime_totalmin = `Total Minutes in Bed`) %>%
    dplyr::group_by(subject) %>%
    dplyr::arrange(`In Bed Date`, `In Bed Time`, .by_group = TRUE) %>%
    dplyr::mutate(
      in_bed_dt = lubridate::mdy_hms(paste(`In Bed Date`, `In Bed Time`), tz = "UTC"),
      out_bed_dt = lubridate::mdy_hms(paste(`Out Bed Date`, `Out Bed Time`), tz = "UTC"),
      in_bed_mid = in_bed_dt + (out_bed_dt - in_bed_dt) / 2,
      next_in_bed_mid = dplyr::lead(in_bed_mid)) %>%
    dplyr::ungroup()
}

#2.2 Assign each actigraphy row the EMA `day` whose survey timestamp falls within
#    that sleep period's [in_bed_mid, next_in_bed_mid) wake interval. Vectorized
#    per subject; non-overlapping intervals guarantee at most one match.
align_actigraphy_to_day <- function(actig, ema) {
  actig$day <- NA_integer_
  subjects <- intersect(unique(actig$subject), unique(ema$subject))
  for (subj in subjects) {
    a_idx <- which(actig$subject == subj)
    e_sub <- ema[ema$subject == subj, c("day", "actual_start_local")]
    e_sub <- e_sub[!is.na(e_sub$actual_start_local), ]
    for (i in a_idx) {
      iv <- lubridate::interval(actig$in_bed_mid[i], actig$next_in_bed_mid[i])
      hit <- which(e_sub$actual_start_local %within% iv)
      if (length(hit) >= 1L) actig$day[i] <- e_sub$day[hit[1L]]
    }
  }
  actig
}


## 3. EMA survey cleaning ##

#3.0 Parse a free-text clock value to hms. The EMA records times on an explicit
#    12-hour clock WITH AM/PM (e.g., "12:00 AM", "11:00 PM"), so we parse that
#    format directly -- NO heuristic AM/PM disambiguation is needed or wanted.
#    Sentinels (e.g., "True") -> NA; a 24-hour fallback catches any stragglers
#    so unexpected entries are never silently dropped. Vectorized over a column.
parse_clock <- function(t) {
  
  #3.0.1 Non-character (already hms/numeric/difftime): coerce straight through
  if (!is.character(t)) {
    return(tryCatch(hms::as_hms(t),
                    error = function(e) hms::as_hms(rep(NA_real_, length(t)))))
  }
  
  #3.0.2 Blank out known sentinels / empties before parsing
  t[t %in% c("", "True", "False", "NA", "NaN")] <- NA_character_
  
  #3.0.3 Helper: parse with a given format and return seconds-since-midnight
  to_secs <- function(x, fmt) {
    p <- strptime(x, format = fmt, tz = "UTC")
    p$hour * 3600 + p$min * 60 + p$sec
  }
  
  #3.0.4 Primary parse: 12-hour clock with AM/PM (case-insensitive; trims spaces)
  secs <- to_secs(toupper(trimws(t)), "%I:%M %p")
  
  #3.0.5 Defensive 24-hour fallbacks for any non-sentinel value that didn't match
  missed <- is.na(secs) & !is.na(t)
  
  if (any(missed)) secs[missed] <- to_secs(t[missed], "%H:%M:%S")
  
  missed <- is.na(secs) & !is.na(t)
  
  if (any(missed)) secs[missed] <- to_secs(t[missed], "%H:%M")
  
  hms::as_hms(secs)
}

#3.1 Hours-from-anchor for a continuous "how late to bed" measure. Anchoring at
#    15:00 avoids the midnight wrap that makes very-early and very-late onset
#    look identical on a raw clock. Vectorized + NA-safe via parse_clock()
hours_from_anchor <- function(t, anchor = "15:00:00") {
  secs <- as.numeric(parse_clock(t))
  (secs - as.numeric(hms::as_hms(anchor))) %% (24 * 3600) / 3600
}

#3.2 Clean the EMA survey: recode binary SI items to 0/1, build SI composites,
#    correct sleep/wake times, and derive continuous timing features.
clean_ema <- function(ema) {
  ema %>%
    dplyr::mutate(
      
      #3.3.1 Parse self-report sleep onset time (explicit AM/PM clock)
      sleeptime_parsed = parse_clock(SLEEPTIME),
      hrs_to_sleep_ema = hours_from_anchor(sleeptime_parsed),
      
      #3.2.2 SI binary items: instrument codes 1/2 -> 0/1
      WISHTOBEDEAD = WISHTOBEDEAD - 1,
      ACTIVESUICIDALTHOUGHTS = ACTIVESUICIDALTHOUGHTS - 1,
      SUICIDEATTEMPT = SUICIDEATTEMPT - 1,
      SELFHARM = SELFHARM - 1,
      
      #3.2.3 SI composites
      si_severity = WISHTOBEDEAD + ACTIVESUICIDALTHOUGHTS + SUICIDEATTEMPT,
      si_any = pmax(WISHTOBEDEAD, ACTIVESUICIDALTHOUGHTS, SUICIDEATTEMPT),
      
      #3.2.4 Sleep self-report composite (items SLEEP2-4 are 1-anchored)
      sleep_disturbance = SLEEP2 + SLEEP3 + SLEEP4 - 3) %>%
        dplyr::rename(
          depression = Depression,
          anxiety = Anxiety,
          sleep_quality = Sleep1,
          si_bother = SICont,
          controllability = Control_Slider)
}


## 4. Assemble the analytic person-day table ##

#4.1 Orchestrate: load -> crosswalk -> clean -> dedupe -> join -> engineer -> save.
build_dataset <- function(cfg) {
  
  #4.1.1 Load sources
  ema <- load_ema_long(file.path(RAW, cfg$sources$ema_long_csv))
  gonogo <- load_ema_gonogo(file.path(RAW, cfg$sources$ema_gonogo_dir),
                            cfg$sources$ema_gonogo_glob)
  actig <- load_actigraphy(actig_dir, cfg$actigraphy$header_skip)
  
  #4.1.2 Participant ID crosswalk. The go/no-go table is the only source carrying
  #       BOTH the platform id (subjectrspid == EMA `id`) and the canonical study
  #       label (subject, e.g. "SleepBD09"); we use it to map id -> subject. Must 
  #       be 1:1 at the participant level
  id_xwalk <- gonogo %>%
    dplyr::distinct(subjectrspid, subject) %>%
    dplyr::filter(!is.na(subjectrspid), !is.na(subject))
  
  if (any(duplicated(id_xwalk$subjectrspid)))
    stop("ID crosswalk is not 1:1: a subjectrspid maps to multiple subject labels.")
  
  #4.1.3 Clean EMA -> PERSON-DAY BACKBONE, then attach canonical `subject` up front
  #       (EMA `id` == go/no-go `subjectrspid`) so every downstream key uses it
  ema_clean <- clean_ema(ema) %>%
    dplyr::left_join(id_xwalk, by = c("id" = "subjectrspid"))
  
  n_unmapped <- sum(is.na(ema_clean$subject))
  
  if (n_unmapped > 0L)
    message("Note: ", n_unmapped, " EMA row(s) had no id->subject match ",
            "(EMA ids absent from the go/no-go files); they carry NA subject.")
  
  #4.1.3.1 Enforce one survey row per (subject, day); keep the most-complete row
  dup_ema <- ema_clean %>% dplyr::filter(!is.na(subject)) %>%
    dplyr::count(subject, day) %>% dplyr::filter(n > 1L)
  
  if (nrow(dup_ema) > 0L) {
    message("Collapsing ", sum(dup_ema$n - 1L),
            " duplicate EMA person-day row(s) to the most-complete record.")
    
    ema_clean <- ema_clean %>%
      dplyr::mutate(.n_obs = rowSums(!is.na(dplyr::across(dplyr::everything())))) %>%
      dplyr::group_by(subject, day) %>%
      dplyr::slice_max(.n_obs, n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.n_obs)
  }
  
  #4.1.4 Attach daily go/no-go METRICS on the per-instance key (instance + id).
  #      Select keys + metrics only so we never re-merge `subject`; many-to-one
  #      guards against fan-out
  gonogo_metrics <- gonogo %>%
    dplyr::select(sessionid, subjectrspid, meanrt, rt_cv, commission, errorrate) %>%
    dplyr::distinct(sessionid, subjectrspid, .keep_all = TRUE)
  
  ema_task <- ema_clean %>%
    dplyr::left_join(gonogo_metrics,
                     by = c("instance_id" = "sessionid", "id" = "subjectrspid"),
                     relationship = "many-to-one") %>%
    dplyr::mutate(
      meanrt = dplyr::if_else(meanrt > 1.5, NA_real_, meanrt),
      rt_cv  = dplyr::if_else(rt_cv  > 2.0, NA_real_, rt_cv))
  
  #4.1.5 Actigraphy: align to EMA day, keep ONE primary (longest-TST) sleep period
  #       per (subject, day), drop unmatched periods, LEFT-join (many-to-one)
  actig_clean <- clean_actigraphy(actig)
  
  actig_only <- setdiff(unique(actig_clean$subject), unique(ema_clean$subject))
  
  if (length(actig_only) > 0L)
    message("Note: ", length(actig_only), " actigraphy label(s) have no EMA match ",
            "(format mismatch would drop their sleep): ",
            paste(utils::head(actig_only, 5), collapse = ", "))
  
  actig_dayed <- align_actigraphy_to_day(actig_clean, ema_task) %>%
    dplyr::filter(!is.na(day)) %>%
    dplyr::group_by(subject, day) %>%
    dplyr::slice_max(sleep_totalmin, n = 1L, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  daily <- ema_task %>%
    dplyr::left_join(actig_dayed, by = c("subject", "day"),
                     relationship = "many-to-one")
  
  #4.1.6 Derived sleep-timing + empirical-logit cognition feature
  daily <- daily %>%
    dplyr::mutate(
      hrs_to_sleep = hours_from_anchor(sleep_onsettime),
      commission_elogit = empirical_logit(commission, n_trials = cfg$cognition$n_nogo_trials))
  
  #4.1.7 Within-person variability scaffolding at rolling window from config
  w  <- cfg$dynamics$ews_window
  cv <- function(v) stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE)
  daily <- daily %>%
    dplyr::group_by(subject) %>%
    dplyr::arrange(day, .by_group = TRUE) %>%
    dplyr::mutate(
      sleepeff_movingcv = zoo::rollapply(sleep_eff,    w, cv, fill = NA, align = "right"),
      sleeptime_movingcv = zoo::rollapply(hrs_to_sleep, w, cv, fill = NA, align = "right")) %>%
    dplyr::ungroup()
  
  #4.1.8 Sample filters: protocol exclusions + drop onboarding day 1
  daily <- daily %>% dplyr::filter(!subject %in% cfg$sample$exclude_ids)
  if (isTRUE(cfg$sample$exclude_day_1)) daily <- daily %>% dplyr::filter(day != 1L)
  
  #4.1.9 Backbone integrity: exactly one row per (subject, day)
  n_rows <- nrow(daily); n_keys <- dplyr::n_distinct(daily[c("subject", "day")])
  
  message(glue::glue(
    "Backbone integrity: {n_rows} rows / {n_keys} unique (subject, day) -> ",
    "{ifelse(n_rows == n_keys, 'OK (one row per person-day)', 'STILL DUPLICATED')}"))
  
  daily
}


## 5. Run + persist ##

#5.1 Build and write with a dynamic, informative filename
if (sys.nframe() == 0L) {
  daily <- build_dataset(cfg)

  out_csv <- file.path(cfg$paths$data_processed,
                       out_name("ema_daily", "csv", algo = cfg$actigraphy$algorithm))
  
  dir.create(cfg$paths$data_processed, recursive = TRUE, showWarnings = FALSE)
  
  write_if_any(daily, out_csv)

  #5.1.1 Quick size + missingness sanity check
  message(glue::glue(
    "Built analytic table: {nrow(daily)} person-days across ",
    "{dplyr::n_distinct(daily$subject)} subjects."))
}
