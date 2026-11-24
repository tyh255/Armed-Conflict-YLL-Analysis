# ---- packages (versions pinned in renv.lock; a missing package stops the run) ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr); library(tibble)
  library(zoo); library(MASS); library(nlme)
  library(tidysynth); library(augsynth); library(did); library(HonestDiD)
  library(countrycode); library(digest)
  library(ggplot2); library(scales); library(patchwork); library(ggridges); library(ragg)
})
select <- dplyr::select; filter <- dplyr::filter   # MASS masks select

# ---- paths (relative to repo root; run from root or set CVD_ROOT) ----
ROOT    <- normalizePath(Sys.getenv("CVD_ROOT", unset = "."), mustWork = TRUE)
RAW_DIR <- file.path(ROOT, "data", "raw")
GBD_DIR <- file.path(ROOT, "data", "gbd")
OUT_DIR <- file.path(ROOT, Sys.getenv("CVD_OUT", unset = "out"))
FIG_DIR <- file.path(OUT_DIR, "figs")
TAB_DIR <- file.path(OUT_DIR, "tables")
for (d in c(OUT_DIR, FIG_DIR, TAB_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- seeds (single block; values unchanged from the working pipeline) ----
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
SEED_SHARED  <- 2024L       # shared structural draws (CFR, VE, multiplier, reach, TB tx)
SEED_GAP     <- 20240601L   # every gap estimator (ITS bootstrap, SC, SC+cov, ASCM, DiD)
SEED_MC      <- 42L         # YLL Monte Carlo (incidence + gap noise); all re-runs reuse it
SEED_PHI_SIM <- 20261111L   # 13_phi_sim.R individual-level check

# ---- conflicts (onset expert-coded, Methods; end = last WUENIC year 2023 if ongoing) ----
conflict_info <- tibble::tribble(
  ~country,        ~ISO3, ~conflict_start, ~conflict_end,
  "Burkina Faso",  "BFA", 2022L,           2023L,
  "Ethiopia",      "ETH", 2012L,           2023L,
  "Iraq",          "IRQ", 2013L,           2022L,
  "Myanmar",       "MMR", 2021L,           2023L,
  "Nigeria",       "NGA", 2010L,           2023L,
  "Pakistan",      "PAK", 2004L,           2015L,
  "Somalia",       "SOM", 2006L,           2023L,
  "Sri Lanka",     "LKA", 2006L,           2009L,
  "Syria",         "SYR", 2011L,           2023L,
  "Yemen",         "YEM", 2014L,           2023L
) %>%
  mutate(preconflict_year   = conflict_start - 1L,
         preconflict_window = purrr::map(preconflict_year, ~ (.x - 2L):.x))

# vaccine -> disease (tetanus = post-neonatal infant tetanus via DTP3)
vaccine_disease_links <- tibble::tribble(
  ~vaccine, ~disease,
  "BCG",    "Tuberculosis",
  "MCV1",   "Measles",
  "DTP3",   "Diphtheria",
  "DTP3",   "Pertussis",
  "DTP3",   "Tetanus")

# ---- MC / discounting ----
N_SIM             <- 10000
DISCOUNT_RATE     <- 0.03                          # continuous, no age weighting
CATCHUP_RATES     <- c(0, 0.25, 0.50, 0.75, 1.0)   # lambda grid
INC_LOGNORM_SIGMA <- 0.20                          # fallback log-sd if a GBD row has no UI

# ---- switches Methods needs ----
# measles conflict CFR multiplier: off in headline, on = sensitivity
MEASLES_MULT_MODE <- Sys.getenv("CVD_MEASLES_MULT", unset = "off")
stopifnot(MEASLES_MULT_MODE %in% c("off", "on"))
# TB CFR = p_treat*CFR_treated + (1-p_treat)*CFR_untreated; off = untreated CFR (upper bound)
TB_TX_WEIGHTING <- identical(Sys.getenv("CVD_TB_TX", unset = "on"), "on")
# incidence basis: I_u = I / (1 - phi*c_obs*VE_I)
#   unvaccinated_timed (headline) phi by schedule | unvaccinated phi=1 (upper) | population no conversion (lower)
INC_BASIS <- Sys.getenv("CVD_INC_BASIS", unset = "unvaccinated_timed")
stopifnot(INC_BASIS %in% c("unvaccinated", "unvaccinated_timed", "population"))
UNVAX_DENOM_FLOOR <- 0.02   # numerical guard on (1 - phi*c*VE_I); caps factor at 50x

# ---- parameter table (Supplementary Table 1) ----
# beta: method of moments, SD = (hi-lo)/3.92 | lognormal_ci: lo/hi = 2.5/97.5 pct, floored at lo
params <- tibble::tribble(
  ~parameter,                ~population,                        ~mean,    ~lo,     ~hi,    ~dist,          ~source,
  "CFR_TB",                  "Ages 0-4",                         0.436,    0.368,   0.506,  "beta",         "Jenkins 2017",
  "CFR_TB",                  "Ages 5-14",                        0.149,    0.115,   0.191,  "beta",         "Jenkins 2017",
  "CFR_TB",                  "Ages 0-4, treated",                0.020,    0.005,   0.074,  "beta",         "Jenkins 2017 (treated paediatric TB)",
  "CFR_Measles",             "Ages <1, community",               0.0303,   0.0289,  0.0316, "beta",         "Sbarra 2023",
  "CFR_Measles",             "Ages 1-4, community",              0.0163,   0.0158,  0.0168, "beta",         "Sbarra 2023",
  "CFR_Diphtheria",          "Ages 0-4",                         0.435,    0.403,   0.467,  "beta",         "Truelove 2020",
  "CFR_Diphtheria",          "Ages 5-19",                        0.232,    0.230,   0.263,  "beta",         "Truelove 2020",
  "CFR_Diphtheria",          "Ages 20+",                         0.290,    0.288,   0.292,  "beta",         "Truelove 2020",
  "CFR_Pertussis",           "Ages 0-1",                         0.037,    0.020,   0.060,  "beta",         "Crowcroft 2003",
  "CFR_Pertussis",           "Ages 1-4",                         0.010,    0.000,   0.020,  "beta",         "Crowcroft 2003",
  "CFR_Pertussis",           "Ages 5+",                          0.000,    0.000,   0.001,  "beta",         "Crowcroft 2003",   # reference only, never drawn
  "CFR_Tetanus",             "Ages 0-1",                         0.33,     0.24,    0.43,   "beta",         "Patel & Mehta 1999",
  "VE_BCG",                  "Prevent TB death",                 0.71,     0.47,    0.84,   "beta",         "Colditz 1994",
  "VE_BCG_disease",          "Prevent TB disease (trials)",      0.51,     0.30,    0.66,   "beta",         "Colditz 1994 (trials)",  # divisor only, quantile-mapped
  "VE_MCV1",                 "First dose >=9 months",            0.83,     0.76,    0.88,   "beta",         "Nic Lochlainn 2019",
  "VE_MCV1_M12",             "First dose >=12 months (Syria)",   0.92,     0.84,    0.95,   "beta",         "Uzicanin/WHO 2017",
  "VE_DTP3_Diphtheria",      "3 doses, clinical disease",        0.87,     0.68,    0.97,   "beta",         "Truelove 2020",
  "VE_DTP3_Pertussis_wP",    "3 doses, acellular",               0.74,     0.55,    0.85,   "beta",         "Simondon 1997 (3 doses, acellular)",  # name kept; value is acellular
  "VE_DTP3_Tetanus",         "3 doses, clinical disease",        0.97,     0.90,    0.99,   "beta",         "WHO/CDC",
  "Measles_conflict_mult",   "Conflict-context multiplier",      2.24,     1.0,     5.0,    "lognormal_ci", "Lam 2010; Salama 2001; NEJM 2025",
  "Catchup_reach",           "Campaign reach (frac zero-dose)",  0.66,     0.35,    0.90,   "beta",         "Portnoy 2018",
  "TB_treat_coverage",       "Frac paediatric TB treated (conflict)", 0.45, 0.25,   0.65,   "beta",         "WHO Global TB Report (paediatric Tx coverage)"
)

# ---- distribution helpers ----
fit_beta <- function(mean_val, lo, hi) {
  if (mean_val <= 0 || mean_val >= 1) stop("Beta mean must be in (0,1); got ", mean_val)
  sd_val  <- (hi - lo) / 3.92
  var_val <- min(sd_val^2, mean_val * (1 - mean_val) * 0.999)
  if (var_val <= 0) var_val <- 1e-6
  common <- (mean_val * (1 - mean_val) / var_val) - 1
  list(alpha = mean_val * common, beta = (1 - mean_val) * common, sd = sd_val)
}
fit_lognormal_ci <- function(lo, hi) {
  if (lo <= 0 || hi <= lo) stop("lognormal_ci needs 0 < lo < hi")
  list(mu = (log(lo) + log(hi)) / 2, sigma = (log(hi) - log(lo)) / 3.92)
}
draw_param <- function(n, mean_val, lo, hi, dist) {
  if (dist == "beta") {
    f <- fit_beta(mean_val, lo, hi); rbeta(n, f$alpha, f$beta)
  } else if (dist == "lognormal_ci") {
    f <- fit_lognormal_ci(lo, hi); pmax(rlnorm(n, f$mu, f$sigma), lo)
  } else stop("Unknown distribution: ", dist)
}
get_param <- function(param_name, pop = NULL) {
  p <- params %>% dplyr::filter(parameter == param_name)
  if (!is.null(pop)) p <- p %>% dplyr::filter(population == pop)
  if (nrow(p) != 1) stop("Parameter not unique/found: ", param_name, " / ", pop)
  p
}

# ---- MCV1 schedule (Syria M12) ----
mcv1_schedule <- tibble::tibble(country = conflict_info$country,
                                mcv1_param = ifelse(conflict_info$country == "Syria",
                                                    "VE_MCV1_M12", "VE_MCV1"))

# ---- disease -> CFR / VE selectors (infant cohort: youngest band) ----
cfr_population_for_disease <- function(disease)
  switch(disease, "Tuberculosis" = "Ages 0-4", "Measles" = "Ages <1, community",
         "Diphtheria" = "Ages 0-4", "Pertussis" = "Ages 0-1", "Tetanus" = "Ages 0-1")
cfr_param_for_disease <- function(disease)
  switch(disease, "Tuberculosis" = "CFR_TB", "Measles" = "CFR_Measles",
         "Diphtheria" = "CFR_Diphtheria", "Pertussis" = "CFR_Pertussis", "Tetanus" = "CFR_Tetanus")
ve_param_for_link <- function(vaccine, disease, country) {
  if (vaccine == "BCG") return("VE_BCG")
  if (vaccine == "MCV1") return(mcv1_schedule$mcv1_param[mcv1_schedule$country == country])
  if (vaccine == "DTP3") return(switch(disease, "Diphtheria" = "VE_DTP3_Diphtheria",
                                       "Pertussis" = "VE_DTP3_Pertussis_wP",
                                       "Tetanus" = "VE_DTP3_Tetanus"))
  stop("Unknown vaccine-disease: ", vaccine, "/", disease)
}

# ---- incidence-basis helpers ----
# phi = share of infancy after dose protects (uniform hazard): 1 - (age + 2 wk)/52 wk
PHI_SCHEDULE <- c(BCG = 1, DTP3 = 1 - 16 / 52, MCV1 = 1 - 41 / 52, MCV1_M12 = 0)
unvax_phi <- function(vaccine, country) {
  if (INC_BASIS != "unvaccinated_timed") return(1)
  key <- vaccine
  if (vaccine == "MCV1" &&
      mcv1_schedule$mcv1_param[mcv1_schedule$country == country] == "VE_MCV1_M12") key <- "MCV1_M12"
  unname(PHI_SCHEDULE[key])
}
# VE against disease: identity except BCG (same-quantile map 0.71 death -> 0.51 disease; no RNG)
ve_incidence_draws <- function(ve_draws, vaccine) {
  if (vaccine != "BCG") return(ve_draws)
  pd <- get_param("VE_BCG");         fd <- fit_beta(pd$mean, pd$lo, pd$hi)
  pv <- get_param("VE_BCG_disease"); fi <- fit_beta(pv$mean, pv$lo, pv$hi)
  stats::qbeta(stats::pbeta(ve_draws, fd$alpha, fd$beta), fi$alpha, fi$beta)
}

# ---- shared (common-random-number) structural draws ----
shared_key <- function(param, pop = NULL)
  if (is.null(pop) || (length(pop) == 1 && is.na(pop))) paste0(param, "||*") else paste0(param, "||", pop)

build_shared_param_registry <- function(conflict_df = conflict_info, vd = vaccine_disease_links) {
  diseases <- unique(vd$disease)
  cfr_reg <- tibble::tibble(param = vapply(diseases, cfr_param_for_disease, character(1)),
                            pop   = vapply(diseases, cfr_population_for_disease, character(1))) %>%
    dplyr::distinct(param, pop)
  ve_names <- character(0)
  for (k in seq_len(nrow(vd))) for (cn in conflict_df$country)
    ve_names <- c(ve_names, ve_param_for_link(vd$vaccine[k], vd$disease[k], cn))
  tb_tx <- if (TB_TX_WEIGHTING)
    tibble::tibble(param = c("CFR_TB", "TB_treat_coverage"), pop = c("Ages 0-4, treated", NA_character_))
  else tibble::tibble(param = character(0), pop = character(0))
  dplyr::bind_rows(cfr_reg,
                   tibble::tibble(param = unique(ve_names), pop = NA_character_),
                   tibble::tibble(param = "Measles_conflict_mult", pop = NA_character_),
                   tibble::tibble(param = "Catchup_reach", pop = NA_character_),
                   tb_tx) %>%
    dplyr::distinct(param, pop)
}

# isolated RNG stream (restores global state on exit)
draw_shared_params <- function(n_sim = N_SIM, seed = SEED_SHARED,
                               conflict_df = conflict_info, vd = vaccine_disease_links) {
  reg <- build_shared_param_registry(conflict_df, vd)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  out <- vector("list", nrow(reg)); keys <- character(nrow(reg))
  for (i in seq_len(nrow(reg))) {
    p <- get_param(reg$param[i], if (is.na(reg$pop[i])) NULL else reg$pop[i])
    out[[i]] <- draw_param(n_sim, p$mean, p$lo, p$hi, p$dist)
    keys[i]  <- shared_key(reg$param[i], reg$pop[i])
  }
  names(out) <- keys
  out
}
get_shared_draws <- function(shared, param, pop = NULL) {
  key <- shared_key(param, pop)
  if (!key %in% names(shared)) stop("Shared draws not found for '", key, "'.")
  shared[[key]]
}

message("Setup: ", nrow(conflict_info), " countries | N_SIM ", N_SIM, " | INC_BASIS ", INC_BASIS,
        " | out ", OUT_DIR)
