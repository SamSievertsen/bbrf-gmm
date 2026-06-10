# =============================================================================
# utils.R
# Shared, dependency-light utilities used across the wrangle, analysis, and
# tests: config loading, feature transforms, defensive IO, and a decluttered
# ggplot theme aligned with Storytelling-with-Data principles.
# =============================================================================

suppressPackageStartupMessages({
  library(here);
  library(yaml)
})


## 1. Configuration ##

#1.1 Load config.yml and resolve all repo-relative paths against here::here().
#    A config.local.yml, if present, is shallow-merged on top to allow
#    machine-specific overrides without touching the tracked config.
load_config <- function(path = here::here("config.yml"), profile = "default") {
  cfg <- yaml::read_yaml(path)[[profile]]

  #1.1.1 Optional local overrides (untracked)
  local_path <- here::here("config.yml")
  if (file.exists(local_path)) {
    local_cfg <- yaml::read_yaml(local_path)[[profile]]
    cfg <- utils::modifyList(cfg, local_cfg)
  }

  #1.1.2 Resolve top-level paths to absolute, repo-anchored paths
  cfg$paths <- lapply(cfg$paths, function(p) here::here(p))
  cfg
}


## 2. Feature transforms ##

#2.1 Empirical-logit transform for a bounded proportion in [0, 1].
#    Adds a 0.5/(n+1) style cushion so exact 0 and 1 map to finite values.
#    `n_trials` is the number of go/no-go trials the proportion is based on; if
#    unknown, a constant pseudo-count (default 0.5) is used.
empirical_logit <- function(p, n_trials = NULL, eps = 0.5) {
  if (is.null(n_trials)) {
    num <- p + eps
    den <- (1 - p) + eps
  } else {
    num <- p * n_trials + eps
    den <- (1 - p) * n_trials + eps
  }
  log(num / den)
}

#2.2 Column-wise z-score that preserves NA and stores center/scale as attributes
#    (so new data can be projected onto the same standardization later).
zscore_matrix <- function(X) {
  X <- as.matrix(X)
  ctr <- colMeans(X, na.rm = TRUE)
  scl <- apply(X, 2L, stats::sd, na.rm = TRUE)
  Z <- sweep(sweep(X, 2L, ctr, "-"), 2L, scl, "/")
  attr(Z, "center") <- ctr
  attr(Z, "scale")  <- scl
  Z
}

#2.3 Person-mean center each feature within subject (within-person "state" view).
#    Returns the centered matrix; used for the complementary within-person GMM
#    that strips between-person trait differences.
person_center <- function(X, person_id) {
  X <- as.matrix(X)
  out <- X
  for (j in seq_len(ncol(X))) {
    pm <- tapply(X[, j], person_id, mean, na.rm = TRUE)
    out[, j] <- X[, j] - pm[as.character(person_id)]
  }
  out
}


## 3. Defensive IO ##

#3.1 Write a data frame only if it is non-empty; otherwise log and skip.
write_if_any <- function(df, path) {
  if (is.data.frame(df) && nrow(df) > 0L && ncol(df) > 0L) {
    readr::write_csv(df, path)
    message("Wrote ", basename(path), " (", nrow(df), " rows).")
  } else {
    message("Skipping write for ", basename(path), " (empty).")
  }
}

#3.2 Build a dynamic, informative output filename: <stem>_<algo>_<date>.<ext>.
out_name <- function(stem, ext = "csv", algo = NULL, date = Sys.Date()) {
  parts <- c(stem, algo, format(as.Date(date), "%Y%m%d"))
  parts <- parts[!is.null(parts) & nzchar(parts)]
  paste0(paste(parts, collapse = "_"), ".", ext)
}


## 4. Visualization theme ##

#4.1 Decluttered ggplot theme: no chartjunk, muted axes, readable type,
#    left-aligned title. Mirrors the seaborn/SWD aesthetic from prior work,
#    adapted to ggplot2 since this project is in R. For interactive versions,
#    wrap any ggplot in plotly::ggplotly().
theme_swd <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = ggplot2::element_text(hjust = 0, colour = "grey30"),
      axis.title = ggplot2::element_text(colour = "grey20"),
      axis.text = ggplot2::element_text(colour = "grey35"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey88", linewidth = 0.4),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(colour = "grey20"),
      plot.title.position = "plot")
}

#4.2 Restrained categorical palette for latent states (color used purposefully)
state_palette <- function(n) {
  base <- c("#2C7FB8", "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#666666")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

#4.3 theme_hw(): primary project theme with bold title, viridis-friendly, bottom legend
theme_hw <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, colour = "grey40"),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"))
}

## 5. Presentation & model-interpretation helpers ##

#5.1 nice_tbl(): clean kable + kableExtra styling
nice_tbl <- function(x, caption = NULL, digits = 3, stripe = FALSE) {
  opts <- c("hover", "condensed", "responsive")
  if (stripe) opts <- c("striped", opts)
  knitr::kable(x, caption = caption, digits = digits) |>
    kableExtra::kable_styling(bootstrap_options = opts, full_width = FALSE)
}

#5.2 gmm_state_profile(): un-standardize component means back to NATIVE units so
#    each state reads interpretably (slider points, %, minutes, logits). Requires
#    X built by zscore_matrix() (carries center/scale attributes
gmm_state_profile <- function(fit, X) {
  ctr <- attr(X, "center"); scl <- attr(X, "scale")
  if (is.null(ctr) || is.null(scl))
    stop("X lacks center/scale attributes; build it with zscore_matrix().")
  native <- sweep(sweep(fit$mu, 2L, scl, "*"), 2L, ctr, "+")
  colnames(native) <- colnames(X)
  out <- tibble::as_tibble(round(native, 2))
  out$state <- factor(seq_len(nrow(out)))
  out$weight <- round(fit$weights, 3)
  out$n_days <- as.integer(round(fit$weights * nrow(X)))
  dplyr::relocate(out, state, weight, n_days)
}

#5.3 gmm_feature_discrimination(): how much each feature separates the states.
#    between = weighted variance of component means; within = weighted mean of
#    component variances (both standardized, so comparable). Ratio >> 1 drives
#    separation; near 0 = ~constant across states (a demotion candidate)
gmm_feature_discrimination <- function(fit, X) {
  w <- fit$weights; K <- fit$K
  grand   <- colSums(w * fit$mu)
  between <- colSums(w * sweep(fit$mu, 2L, grand, "-")^2)
  within  <- Reduce(`+`, lapply(seq_len(K), function(k) w[k] * diag(fit$Sigma[[k]])))
  tibble::tibble(feature = colnames(X),
                 between = round(between, 3),
                 within  = round(within, 3),
                 discrimination = round(between / within, 3)) |>
    dplyr::arrange(dplyr::desc(discrimination))
}
