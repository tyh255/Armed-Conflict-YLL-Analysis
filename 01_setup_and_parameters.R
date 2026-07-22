# ---- 0. Packages -----------------------------------------------------------
required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "purrr", "tibble",
  "zoo",                       # na.approx
  "MASS",                      # mvrnorm for ITS bootstrap
  "nlme",                      # gls AR(1) REML for ITS (Prais-Winsten-like)
  "tidysynth",                 # synthetic control
  "ggplot2", "scales",         # plotting
  "patchwork",                  # multi-panel figures
  "HonestDiD", "did", "fixest"
)
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
invisible(lapply(required_pkgs, install_if_missing))

# ---- 1. Paths --------------------------------------------------------------
# DATA_DIR root can be overridden without editing source by setting the
# environment variable CVD_DATA_DIR. The default matches the prior pipeline.
DATA_DIR     <- Sys.getenv("CVD_DATA_DIR",
                           unset = "~/Documents/Conflict Vaccination DALYs")
COVERAGE_DIR <- file.path(DATA_DIR, "Vaccination Coverage Rates")
YLL_DATA_DIR <- file.path(DATA_DIR, "YLL Data")
WHO_DIR      <- Sys.getenv("CVD_WHO_DIR",
                           unset = "~/Documents/Conflict Structural Equation/Data")
OUT_DIR      <- file.path(DATA_DIR, "Output_New")
FIG_DIR      <- file.path(OUT_DIR, "figures")
TAB_DIR      <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Country / conflict periods (locked to draft Table 1) --------------
conflict_info <- tibble::tribble(
  ~country,        ~ISO3, ~conflict_start, ~conflict_end,
  "Burkina Faso",  "BFA", 2022L,           2024L,
  "Ethiopia",      "ETH", 2012L,           2024L,
  "Iraq",          "IRQ", 2013L,           2022L,
  "Myanmar",       "MMR", 2021L,           2024L,
  "Nigeria",       "NGA", 2010L,           2024L,
  "Pakistan",      "PAK", 2004L,           2015L,
  "Somalia",       "SOM", 2006L,           2024L,
  "Sri Lanka",     "LKA", 2006L,           2009L,
  "Syria",         "SYR", 2011L,           2024L,
  "Yemen",         "YEM", 2014L,           2024L
) %>%
  mutate(
    preconflict_year   = conflict_start - 1L,
    # 3-year pre-conflict window for the baseline (robust to single-year noise)
    preconflict_window = purrr::map(preconflict_year, ~ (.x - 2L):.x)
  )

# ---- 3. Vaccine -> disease links -------------------------------------------
vaccine_disease_links <- tibble::tribble(
  ~vaccine, ~disease,
  "BCG",    "Tuberculosis",
  "MCV1",   "Measles",
  "DTP3",   "Diphtheria",
  "DTP3",   "Pertussis",
  "DTP3",   "Tetanus"
)

# ---- 4. Monte Carlo & discounting settings ---------------------------------
N_SIM         <- 10000        # Monte Carlo draws per country-vaccine-disease
# 3% continuous discounting, no age weighting: the GBD2010 / WHO-CHOICE
# convention. (GBD2019 uses 0% and no age weighting.) Undiscounted YLLs are
# reported alongside discounted throughout, so the choice is transparent.
DISCOUNT_RATE <- 0.03
CATCHUP_RATES <- c(0, 0.25, 0.50, 0.75, 1.0)   # catch-up scenarios
INC_LOGNORM_SIGMA <- 0.20     # fallback log-space SD when a GBD row lacks a UI
# Seed used for the SHARED structural-parameter draws (CFR/VE/multiplier).
# Kept separate from the year-specific noise seed in 04 so the two streams are
# independent and individually reproducible.
SHARED_DRAW_SEED <- 2024L

# ---- 4a. Methodological switches (red-team remediation) --------------------
# MEASLES_MULT_MODE: whether the conflict measles-CFR multiplier (~x2.43 mean)
#   is applied. Folding it into the point estimate inflated the measles arm and
#   accounted for ~40% of the global YLL total, so it is now OFF by default and
#   intended to be reported only as an explicit, labelled scenario. Set 'on'
#   (or env CVD_MEASLES_MULT=on) to reinstate it as a sensitivity, NOT headline.
MEASLES_MULT_MODE <- Sys.getenv("CVD_MEASLES_MULT", unset = "off")
stopifnot(MEASLES_MULT_MODE %in% c("off", "on"))

# TB_TX_WEIGHTING: apply treatment-access weighting to the TB CFR. The Jenkins
#   2017 value (0.436) is the UNTREATED paediatric CFR; applying it to every
#   incident case assumes zero treatment access. With weighting ON (default) the
#   effective CFR is p_treat*CFR_treated + (1-p_treat)*CFR_untreated, with both
#   the treated CFR and the treatment-coverage fraction drawn as shared
#   structural parameters. Set 'off' (CVD_TB_TX=off) to recover the untreated-CFR
#   upper bound. NOTE: TB_treat_coverage (params table) is a literature-anchored
#   PLACEHOLDER and should be calibrated to country-level paediatric TB treatment
#   coverage (WHO Global TB Report) before publication.
TB_TX_WEIGHTING <- identical(Sys.getenv("CVD_TB_TX", unset = "on"), "on")

# ---- 4b. Age-band & tetanus-arm configuration (recommendations #7, #8) -----
AGE_BAND_MODE <- Sys.getenv("CVD_AGE_BAND_MODE", unset = "infant")
stopifnot(AGE_BAND_MODE %in% c("infant", "under5_tb_diph"))

.under5_diseases <- c("Tuberculosis", "Diphtheria")
.use_under5 <- function(disease)
  AGE_BAND_MODE == "under5_tb_diph" && disease %in% .under5_diseases

# IHME incidence age_name to query per disease (consumed by 02).
disease_age_band <- function(disease)
  if (.use_under5(disease)) "Under 5" else "<1 year"   # verify exact age_name string

# Population denominator column per disease in covariates_panel (consumed by 04).
disease_denominator_col <- function(disease)
  if (.use_under5(disease)) "Pop_Under5" else "Birth_Cohort"

TETANUS_ARM <- Sys.getenv("CVD_TETANUS_ARM", unset = "post_neonatal")
stopifnot(TETANUS_ARM %in% c("neonatal", "post_neonatal", "drop"))
if (TETANUS_ARM == "neonatal") {
  warning(
    "TETANUS_ARM='neonatal': applying the Seale 2013 neonatal CFR and the ",
    "maternal-TT VE. This is ONLY valid if step 02 has been switched to load ",
    "GBD cause 'Neonatal tetanus'. With the default non-neonatal incidence ",
    "still in place, this applies a neonatal CFR to non-neonatal incidence - ",
    "a regime mismatch. Switch step 02 first.",
    call. = FALSE
  )
}
if (TETANUS_ARM == "drop") {
  vaccine_disease_links <- vaccine_disease_links %>% dplyr::filter(disease != "Tetanus")
  message("TETANUS_ARM='drop': DTP3 -> Tetanus arm excluded from the analysis.")
}

# ---- 5. Parameter table with distributions ---------------------------------
params <- tibble::tribble(
  ~parameter,                ~population,                        ~mean,    ~lo,     ~hi,    ~dist,          ~source,
  # ----- Case fatality ratios -----
  "CFR_TB",                  "Ages 0-4",                         0.436,    0.368,   0.506,  "beta",         "Jenkins 2017",
  "CFR_TB",                  "Ages 5-14",                        0.149,    0.115,   0.191,  "beta",         "Jenkins 2017",
  # Treated paediatric TB CFR (used by the TB_TX_WEIGHTING blend; ~an order of
  # magnitude below the untreated 0.436). Verified against Jenkins 2017 treated
  # arm (Table S1: 2.0%, 95% CI 0.5-7.4); note the wide, right-skewed interval.
  "CFR_TB",                  "Ages 0-4, treated",                0.020,    0.005,   0.074,  "beta",         "Jenkins 2017 (treated paediatric TB)",
  "CFR_Measles",             "Ages <1, community",               0.0303,   0.0289,  0.0316, "beta",         "Sbarra 2023",
  "CFR_Measles",             "Ages 1-4, community",              0.0163,   0.0158,  0.0168, "beta",         "Sbarra 2023",
  "CFR_Diphtheria",          "Ages 0-4",                         0.435,    0.403,   0.467,  "beta",         "Truelove 2020",
  "CFR_Diphtheria",          "Ages 5-19",                        0.232,    0.230,   0.263,  "beta",         "Truelove 2020",
  "CFR_Diphtheria",          "Ages 20+",                         0.290,    0.288,   0.292,  "beta",         "Truelove 2020",
  "CFR_Pertussis",           "Ages 0-1",                         0.037,    0.020,   0.060,  "beta",         "Crowcroft 2003",
  "CFR_Pertussis",           "Ages 1-4",                         0.010,    0.000,   0.020,  "beta",         "Crowcroft 2003",
  "CFR_Pertussis",           "Ages 5+",                          0.000,    0.000,   0.001,  "beta",         "Crowcroft 2003",
  "CFR_Tetanus",             "Ages 0-1",                         0.33,     0.24,    0.43,   "beta",         "Patel & Mehta 1999",
  "CFR_Tetanus",             "Neonatal",                         0.64,     0.55,    0.72,   "beta",         "Seale 2013 (with care)",
  # ----- Vaccine effectiveness -----
  "VE_BCG",                  "Prevent TB death",                 0.71,     0.47,    0.84,   "beta",         "Colditz 1994",
  "VE_MCV1",                 "First dose >=9 months",            0.83,     0.76,    0.88,   "beta",         "Nic Lochlainn 2019",
  "VE_MCV1_M12",             "First dose >=12 months (Syria)",   0.92,     0.84,    0.95,   "beta",         "Uzicanin/WHO 2017",
  "VE_DTP3_Diphtheria",      "3 doses, clinical disease",        0.87,     0.68,    0.97,   "beta",         "Truelove 2020",
  "VE_DTP3_Pertussis_wP",    "3 doses, acellular",               0.74,     0.55,    0.85,   "beta",         "Simondon 1997 (3 doses, acellular)",
  "VE_DTP3_Tetanus",         "3 doses, clinical disease",        0.97,     0.90,    0.99,   "beta",         "WHO/CDC",
  "VE_TT_maternal",          "2+ doses (pregnant women)",        0.94,     0.80,    0.98,   "beta",         "Blencowe 2010",
  # ----- Measles conflict CFR multiplier (percentile-matched lognormal) -----
  # lo/hi are the 2.5/97.5 percentiles; point estimate shown is the MEDIAN
  # (sqrt(lo*hi) = 2.24). Implied mean ~2.43. Draws floored at 1.0.
  "Measles_conflict_mult",   "Conflict-context multiplier",      2.24,     1.0,     5.0,    "lognormal_ci", "Lam 2010; Salama 2001; NEJM 2025",
  # ----- Catch-up campaign reach (shared structural draw; recommendation #3) -----
  "Catchup_reach",           "Campaign reach (frac zero-dose)",  0.66,     0.35,    0.90,   "beta",         "Portnoy 2018",
  # ----- TB treatment-access weighting (TB_TX_WEIGHTING) --------------------
  "TB_treat_coverage",       "Frac paediatric TB treated (conflict)", 0.45, 0.25,   0.65,   "beta",         "WHO Global TB Report (paediatric Tx coverage)"
)

is_supplementary_row <- function(parameter, population) {
  (parameter == "CFR_Tetanus" & population == "Neonatal") |
    (parameter == "VE_TT_maternal")
}

# ---- 6. Distribution-fitting helpers ---------------------------------------
fit_beta <- function(mean_val, lo, hi) {
  if (mean_val <= 0 || mean_val >= 1)
    stop("Beta mean must be in (0,1); got ", mean_val)
  sd_val  <- (hi - lo) / 3.92
  var_val <- sd_val^2
  # Cap variance below mean*(1-mean) to keep alpha, beta positive
  max_var <- mean_val * (1 - mean_val) * 0.999
  var_val <- min(var_val, max_var)
  if (var_val <= 0) var_val <- 1e-6
  common  <- (mean_val * (1 - mean_val) / var_val) - 1
  alpha <- mean_val * common
  beta  <- (1 - mean_val) * common
  list(alpha = alpha, beta = beta, sd = sd_val)
}

fit_lognormal <- function(mean_val, lo, hi) {
  if (mean_val <= 0) stop("Lognormal mean must be > 0; got ", mean_val)
  lo_safe <- max(lo, mean_val * 0.01)   # guard log(0)
  sigma   <- (log(hi) - log(lo_safe)) / 3.92
  if (sigma <= 0) sigma <- 0.01
  mu      <- log(mean_val) - sigma^2 / 2
  list(mu = mu, sigma = sigma)
}

fit_lognormal_ci <- function(lo, hi) {
  if (lo <= 0 || hi <= lo) stop("lognormal_ci needs 0 < lo < hi; got ",
                                lo, ", ", hi)
  sigma <- (log(hi) - log(lo)) / 3.92
  mu    <- (log(lo) + log(hi)) / 2          # = log(median)
  list(mu = mu, sigma = sigma)
}

# Generic draw function
draw_param <- function(n, mean_val, lo, hi, dist) {
  if (dist == "beta") {
    f <- fit_beta(mean_val, lo, hi)
    rbeta(n, f$alpha, f$beta)
  } else if (dist == "lognormal") {
    f <- fit_lognormal(mean_val, lo, hi)
    rlnorm(n, f$mu, f$sigma)
  } else if (dist == "lognormal_ci") {
    f <- fit_lognormal_ci(lo, hi)
    x <- rlnorm(n, f$mu, f$sigma)
    pmax(x, lo)                              # floor at lo (2.5th pct); >= 1 here
  } else stop("Unknown distribution: ", dist)
}

# Lookup parameter by name (and optional population qualifier)
get_param <- function(param_name, pop = NULL) {
  p <- params %>% dplyr::filter(parameter == param_name)
  if (!is.null(pop)) p <- p %>% dplyr::filter(population == pop)
  if (nrow(p) == 0) stop("Parameter not found: ", param_name, " / ", pop)
  if (nrow(p) > 1)  stop("Parameter ambiguous: ", param_name, " / ", pop)
  p
}

draw_named <- function(n, param_name, pop = NULL) {
  p <- get_param(param_name, pop)
  draw_param(n, p$mean, p$lo, p$hi, p$dist)
}

# ---- 7. Per-country MCV1 schedule (M9 vs M12) ------------------------------
mcv1_schedule <- tibble::tribble(
  ~country,        ~mcv1_param,
  "Burkina Faso",  "VE_MCV1",
  "Ethiopia",      "VE_MCV1",
  "Iraq",          "VE_MCV1",
  "Myanmar",       "VE_MCV1",
  "Nigeria",       "VE_MCV1",
  "Pakistan",      "VE_MCV1",
  "Somalia",       "VE_MCV1",
  "Sri Lanka",     "VE_MCV1",
  "Syria",         "VE_MCV1_M12",
  "Yemen",         "VE_MCV1"
)

# ---- 8. Disease -> CFR / VE row selectors ---------------------------------
cfr_population_for_disease <- function(disease) {
  switch(disease,
         "Tuberculosis" = "Ages 0-4",
         "Measles"      = "Ages <1, community",
         "Diphtheria"   = "Ages 0-4",
         "Pertussis"    = "Ages 0-1",
         "Tetanus"      = if (TETANUS_ARM == "neonatal")
           "Neonatal"
         else
           "Ages 0-1")
}

cfr_param_for_disease <- function(disease) {
  switch(disease,
         "Tuberculosis" = "CFR_TB",
         "Measles"      = "CFR_Measles",
         "Diphtheria"   = "CFR_Diphtheria",
         "Pertussis"    = "CFR_Pertussis",
         "Tetanus"      = "CFR_Tetanus")
}

ve_param_for_link <- function(vaccine, disease, country) {
  if (vaccine == "BCG") {
    "VE_BCG"
  } else if (vaccine == "MCV1") {
    mcv1_schedule$mcv1_param[mcv1_schedule$country == country]
  } else if (vaccine == "DTP3" && disease == "Diphtheria") {
    "VE_DTP3_Diphtheria"
  } else if (vaccine == "DTP3" && disease == "Pertussis") {
    "VE_DTP3_Pertussis_wP"
  } else if (vaccine == "DTP3" && disease == "Tetanus") {
    if (TETANUS_ARM == "neonatal") "VE_TT_maternal" else "VE_DTP3_Tetanus"
  } else {
    stop("Unknown vaccine-disease combination: ", vaccine, "/", disease)
  }
}

# ---- 8b. SHARED (common-random-number) structural parameter draws ---------

shared_key <- function(param, pop = NULL) {
  if (is.null(pop) || (length(pop) == 1 && is.na(pop)))
    paste0(param, "||*")
  else
    paste0(param, "||", pop)
}

# Distinct (parameter, population) pairs actually used by the pipeline.
build_shared_param_registry <- function(conflict_df = conflict_info,
                                        vd = vaccine_disease_links) {
  diseases <- unique(vd$disease)
  cfr_reg <- tibble::tibble(
    param = vapply(diseases, cfr_param_for_disease, character(1)),
    pop   = vapply(diseases, cfr_population_for_disease, character(1))
  ) %>% dplyr::distinct(param, pop)
  
  ve_names <- character(0)
  for (k in seq_len(nrow(vd))) {
    for (cn in conflict_df$country) {
      ve_names <- c(ve_names, ve_param_for_link(vd$vaccine[k], vd$disease[k], cn))
    }
  }
  ve_reg <- tibble::tibble(param = unique(ve_names), pop = NA_character_)
  
  mult_reg <- tibble::tibble(param = "Measles_conflict_mult", pop = NA_character_)
  
  # Campaign reach (recommendation #3): shared across all countries' campaigns.
  reach_reg <- tibble::tibble(param = "Catchup_reach", pop = NA_character_)

  tb_tx_reg <- if (isTRUE(get0("TB_TX_WEIGHTING", ifnotfound = FALSE)))
    tibble::tibble(param = c("CFR_TB", "TB_treat_coverage"),
                   pop   = c("Ages 0-4, treated", NA_character_))
  else tibble::tibble(param = character(0), pop = character(0))
  
  dplyr::bind_rows(cfr_reg, ve_reg, mult_reg, reach_reg, tb_tx_reg) %>%
    dplyr::distinct(param, pop)
}

# Draw every shared parameter once, under an isolated, reproducible RNG stream.
draw_shared_params <- function(n_sim = N_SIM, seed = SHARED_DRAW_SEED,
                               conflict_df = conflict_info,
                               vd = vaccine_disease_links) {
  reg <- build_shared_param_registry(conflict_df, vd)
  
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  
  out  <- vector("list", nrow(reg))
  keys <- character(nrow(reg))
  for (i in seq_len(nrow(reg))) {
    pop_i   <- if (is.na(reg$pop[i])) NULL else reg$pop[i]
    p       <- get_param(reg$param[i], pop_i)
    out[[i]] <- draw_param(n_sim, p$mean, p$lo, p$hi, p$dist)
    keys[i]  <- shared_key(reg$param[i], reg$pop[i])
  }
  names(out) <- keys
  attr(out, "n_sim") <- n_sim
  attr(out, "seed")  <- seed
  message("Shared parameter draws built: ", nrow(reg),
          " parameters x ", n_sim, " draws (seed ", seed, ").")
  out
}

get_shared_draws <- function(shared, param, pop = NULL) {
  key <- shared_key(param, pop)
  if (!key %in% names(shared)) {
    stop("Shared draws not found for '", key, "'. Build them with ",
         "draw_shared_params() and pass the result into the simulator.")
  }
  shared[[key]]
}

# ---- 9. Print summary -----------------------------------------------------
message("Setup loaded:")
message("  Countries: ", paste(conflict_info$country, collapse = ", "))
message("  Vaccines: BCG, MCV1, DTP3")
message("  Diseases: TB, Measles, Diphtheria, Pertussis, Tetanus")
message("  Monte Carlo draws: ", N_SIM)
message("  Discount rate: ", DISCOUNT_RATE, " (undiscounted reported alongside)")
message("  Tetanus arm: ", TETANUS_ARM,
        if (TETANUS_ARM == "neonatal")
          " (Seale CFR + maternal-TT VE; step-02 MUST load neonatal incidence)"
        else if (TETANUS_ARM == "post_neonatal")
          " (Patel & Mehta CFR + DTP3 VE; non-neonatal <1 incidence)"
        else "")
message("  Output directory: ", OUT_DIR)
