# BBRF GMM: Latent Daily States & Suicide Risk

A from-scratch Gaussian Mixture Model (GMM), fit by Expectation-Maximization
(EM), applied to intensive longitudinal data from the BBRF study of adolescents
with bipolar disorder. The model discovers latent daily *psychophysiological
states* from continuous mood, actigraphy-sleep, and go/no-go cognition features;
those states are then characterized by their suicidal-ideation profiles and
analyzed dynamically across each person's 30 days

This is the final project for an ML course, so the learning algorithm is
**implemented entirely from scratch** (`R/gmm_em.R`): the Gaussian log-density,
E-step responsibilities, M-step updates, observed-data log-likelihood, BIC model
selection, and multi-restart initialization are all hand-rolled. No mixture-model
package (e.g., mclust, mixtools, flexmix, etc.) is used

## Scientific framing

Daily snapshots (one per person-day) are pooled into a single feature matrix and
clustered into K multivariate-Gaussian states. Each component's covariance
encodes how the features co-vary *within* a state (the "how does variability
relate" question); the per-person sequence of state assignments over days
recovers the dynamical/trajectory view (state dwell, switching, and Scheffer-style
early-warning signals before transitions into high-SI states)

## Data safety (read before touching data)

This study contains identifiable, IRB-protected daily suicide/self-harm reports
and wearable data. **No participant-level data may ever be committed.** The entire
`data/`, `results/`, and `figures/` trees are git-ignored. Place downloaded study
files under `data/raw/` locally; the pipeline writes the analytic table to
`data/processed/` (also ignored). Commit only code and the config

## Repository structure

```
bbrf-gmm/
├── config.yml              # single source of truth (paths, features, GMM settings)
├── src/
│   ├── gmm_em.R            # FROM-SCRATCH GMM/EM (unit-tested) + BIC / AIC / LIC + helpers
│   ├── utils.R             # config loader, transforms, IO, SWD ggplot theme
├── analysis/
│   └── 01_build_dataset.R  # clean, config-driven wrangle (run locally, v1)
│   └── 02_gmm_states.Rmd   # 3-stage analysis driver (EDA -> GMM -> dynamics)
├── tests/
│   └── test_gmm_em.R       # simulation-based recovery test for the EM
├── data/  results/  figures/   # all git-ignored (PHI)
```

## Quick start

```bash
# 1. Install R package dependencies (one time)
Rscript -e 'install.packages(c("here","yaml","dplyr","tidyr","readr","stringr",
  "lubridate","hms","zoo","purrr","glue","ggplot2","rmarkdown", "kableExtra"))'

# 2. Verify the from-scratch EM (no data needed; pure base R)
Rscript tests/test_gmm_em.R

# 3. Place raw study files under data/raw/, then build the analytic table
Rscript analysis/01_build_dataset.R

# 4. Render the analysis
Rscript -e 'rmarkdown::render("analysis/02_gmm_states.Rmd")'
```

## Status

- `src/gmm_em.R` + `tests/test_gmm_em.R`: **verified** (recovers known means,
  weights, and the true K on simulated data)
- `analysis/01_build_dataset.R`: **v1, not yet run on the real files.** Confirm exact
  column names/encodings and the SleepBD27 (two-part) / SleepBD08 (no-actigraphy)
  edge cases on first local run; validate the sleep/nap timing correction first
