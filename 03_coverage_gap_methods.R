# ============================================================================
# 03_COVERAGE_GAP_METHODS.R   (revised)
# ----------------------------------------------------------------------------
# Three methods for estimating conflict-attributable coverage gaps:
#
#   (1) gap_baseline   - observed vs. 3-year pre-conflict mean
#   (2) gap_synthetic  - tidysynth SC, full-window-complete donor pool,
#                        RMSPE-filtered placebo SD, optional lag predictors
#   (3) gap_its        - ITS with logit link, Jensen-corrected back-transform,
#                        parametric bootstrap that includes RESIDUAL variance,
#                        capped counterfactual with cap-binding diagnostic
#
# Each returns a long-format tibble with at least:
#   country, ISO3, vaccine, year, coverage_obs, coverage_cf, gap, gap_sd, method
# (SC/ITS also carry method-specific diagnostic columns.)
#
# Units: coverage and gap in percentage points (0-100), gap_sd in pp.
# estimate_all_gaps() runs all three, stacks, and adds proportion versions.
#
# ---------------------------------------------------------------------------
# SIGN CONVENTION (consistent across all methods):
#     gap = counterfactual - observed   ==>  positive gap = coverage shortfall
#
# KEY ASSUMPTIONS (read before trusting any number):
#   * Coverage point estimates are treated as known (WUENIC's own uncertainty
#     is NOT propagated).
#   * The three methods estimate DIFFERENT counterfactuals (flat / donor-
#     weighted / trend-extrapolated). They are a triangulation, not three
#     independent estimates of one quantity; their gap_sd are not on the same
#     inferential footing.
#   * No conflict spillover onto donor / comparator countries (SUTVA).
#   * Baseline anchors on preconflict_year; SC and ITS anchor on
#     conflict_start. estimate_all_gaps() warns if these disagree.
#   * Method-specific assumptions are documented at each function.
#
# NOT YET IMPLEMENTED (candidate next steps, left as comments where relevant):
#   * Augmented SC (Ben-Michael et al. 2021) ridge bias-correction for
#     imperfect pre-fit; penalized SC (Abadie & L'Hour 2021).
#   * Conformal SC inference (Chernozhukov, Wuthrich & Zhu 2021) as an
#     alternative to placebo permutation.
#   * Multiple imputation to propagate WUENIC coverage uncertainty.
#
# IMPLEMENTED (recommendation #6): ITS autocorrelation handling now selects OLS
#   (short pre-segments) vs GLS AR(1) by REML / Prais-Winsten-like (longer
#   segments) by series length; Newey-West HAC is demoted to an explicit option.
#   Coefficient draws use a multivariate-t (df = residual df) to widen short-
#   series intervals. (Bottomley 2023; Turner 2021.)
#
# LITERATURE BASIS for the choices below:
#   SC  - Abadie, Diamond & Hainmueller (2010); Abadie (2021, JEL); epidemiology
#         tutorial Bouttell/Bonell-style guidance: placebo inference restricted
#         to comparable-fit donors, RMSPE-ratio permutation p-value.
#   ITS - Bottomley et al. (2023, Stat Med); Turner et al. (2021, BMC MRM):
#         account for autocorrelation, but prefer OLS/REML in short segments and
#         GLS AR(1)/Prais-Winsten in longer ones over Newey-West HAC; do not
#         trust Durbin-Watson in short series.
# ============================================================================
# augsynth/did are used (guarded by requireNamespace) in 03b. Load them SOFTLY
# here so a missing optional package does not crash sourcing of the pipeline;
# 03 itself does not require them, and 03b self-skips the methods that do.
for (.p in c("augsynth", "did"))
  if (requireNamespace(.p, quietly = TRUE))
    suppressPackageStartupMessages(library(.p, character.only = TRUE))

# ---- Global reproducibility default ---------------------------------------
DEFAULT_SEED <- 20240601L

# ---- ITS autocorrelation policy (recommendation #6) -----------------------
# Series-length threshold: pre-segments with >= ITS_GLS_MIN_N points use GLS
# AR(1) by REML (Prais-Winsten-like); shorter segments use OLS with no
# autocorrelation adjustment. This follows the simulation evidence that, for
# short series, OLS / REML are better-calibrated than Newey-West HAC, which is
# unstable and tends to under-cover (Bottomley et al. 2023; Turner et al. 2021).
ITS_GLS_MIN_N <- 12L

# ITS extrapolation horizon cap (red-team remediation). The previous unbounded
# default (max_horizon = Inf) let the within-unit logit trend be projected the
# full length of long conflicts (e.g. Nigeria ~15 yr, Ethiopia ~13 yr), which
# BOTH over-attributed (Nigeria, large fabricated gaps) and saturated to spurious
# zeros (Ethiopia, Iraq: cf overtakes observed -> gap <= 0). Post-onset years at
# or beyond this horizon are now DROPPED from the ITS gap, so only the credible
# near-onset window enters the YLL chain. Long-conflict period totals are
# therefore truncated to this window (surfaced via n_valid_years); set the env
# var CVD_ITS_MAX_HORIZON=Inf to recover the prior behaviour. Aligned with
# reliable_horizon (8): years already flagged horizon_unreliable are excluded.
ITS_MAX_HORIZON <- {
  .h <- Sys.getenv("CVD_ITS_MAX_HORIZON", unset = "8")
  if (identical(tolower(.h), "inf")) Inf else as.integer(.h)
}

# Multivariate-t coefficient draws (df = residual df): heavier tails than the
# normal mvrnorm, to reflect estimated-variance uncertainty in short ITS
# segments where normal predictive intervals under-cover (Turner et al. 2021).
mvt_draws <- function(n, mu, Sigma, df) {
  z <- MASS::mvrnorm(n, mu = rep(0, length(mu)), Sigma = Sigma)
  if (is.matrix(z) && is.finite(df) && df > 0)
    z <- z * sqrt(df / stats::rchisq(n, df = df))
  sweep(z, 2, mu, "+")
}

# ---- Helpers --------------------------------------------------------------
preconflict_mean_coverage <- function(coverage_df, ctry, vacc, preconf_year) {
  yrs <- (preconf_year - 2L):preconf_year
  vals <- coverage_df %>%
    dplyr::filter(country == ctry, vaccine == vacc, year %in% yrs) %>%
    dplyr::pull(coverage)
  if (length(vals) == 0) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

preconflict_sd_coverage <- function(coverage_df, ctry, vacc, preconf_year) {
  yrs <- (preconf_year - 2L):preconf_year
  vals <- coverage_df %>%
    dplyr::filter(country == ctry, vaccine == vacc, year %in% yrs) %>%
    dplyr::pull(coverage)
  if (length(vals) < 2) return(2.0)
  max(stats::sd(vals, na.rm = TRUE), 2.0)
}

# Standard error of the pre-conflict mean (uncertainty in the counterfactual
# itself, as opposed to natural year-to-year dispersion). Exposed so callers
# can choose which uncertainty notion to use downstream.
preconflict_se_mean <- function(coverage_df, ctry, vacc, preconf_year) {
  yrs <- (preconf_year - 2L):preconf_year
  vals <- coverage_df %>%
    dplyr::filter(country == ctry, vaccine == vacc, year %in% yrs) %>%
    dplyr::pull(coverage)
  vals <- vals[!is.na(vals)]
  if (length(vals) < 2) return(NA_real_)
  stats::sd(vals) / sqrt(length(vals))
}

# ---- METHOD 1: Simple baseline -------------------------------------------
# ASSUMPTIONS: counterfactual coverage stays FLAT at the 3-year pre-conflict
# mean (no secular trend); gap_sd reflects natural pre-period dispersion of
# coverage (floored at 2pp), held constant across conflict years. We also
# carry se_cf (SE of the counterfactual mean) for callers who want estimate
# uncertainty rather than dispersion.
gap_baseline <- function(coverage_df, conflict_df,
                         vaccines = c("BCG", "MCV1", "DTP3")) {
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (vacc in vaccines) {
      pre_mu <- preconflict_mean_coverage(coverage_df, ci$country, vacc, ci$preconflict_year)
      pre_sd <- preconflict_sd_coverage(coverage_df, ci$country, vacc, ci$preconflict_year)
      pre_se <- preconflict_se_mean(coverage_df, ci$country, vacc, ci$preconflict_year)
      if (is.na(pre_mu)) next
      ctry_yrs <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc,
                      year >= ci$conflict_start, year <= ci$conflict_end) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs = coverage)
      if (nrow(ctry_yrs) == 0) next
      ctry_yrs <- ctry_yrs %>%
        dplyr::mutate(
          coverage_cf = pre_mu,
          gap         = pre_mu - coverage_obs,
          gap_sd      = pre_sd,
          se_cf       = pre_se,
          method      = "Baseline"
        )
      out[[length(out) + 1]] <- ctry_yrs
    }
  }
  dplyr::bind_rows(out)
}

# ============================================================================
# DONOR-POOL SCOPE  (F1: SUTVA neighbour exclusion + region/income matching)
# ----------------------------------------------------------------------------
# Both peer reviews flag that the donor pool is any non-conflict country
# worldwide (no neighbour screen, no comparability screen), which (a) admits
# spillover-contaminated neighbours (SUTVA) and (b) seats richer, flat-trend,
# high-coverage donors off the treated unit's support. Scope is now selectable;
# the DEFAULT ("global") reproduces the previous pool exactly, so this is
# non-breaking. Flip the whole pipeline with CVD_DONOR_SCOPE without editing any
# call site, or pass scope= per call. Used by SC, SC+cov, AugSC and DiD controls.
# ============================================================================

# Land/near-maritime neighbours of the 10 treated units (ISO3). Excluding these
# removes the donors most exposed to cross-border spillover (refugee flows,
# regional campaign disruption, shared outbreaks). Source: UN geoscheme / CIA WFB.
CONFLICT_NEIGHBORS <- list(
  BFA = c("MLI","NER","BEN","TGO","GHA","CIV"),
  ETH = c("ERI","DJI","SOM","KEN","SSD","SDN"),
  IRQ = c("SYR","TUR","IRN","KWT","SAU","JOR"),
  MMR = c("BGD","IND","CHN","LAO","THA"),
  NGA = c("BEN","NER","TCD","CMR"),
  PAK = c("AFG","IRN","IND","CHN"),
  SOM = c("ETH","DJI","KEN","YEM"),
  LKA = c("IND","MDV"),
  SYR = c("TUR","IRQ","JOR","LBN","ISR"),
  YEM = c("SAU","OMN","DJI","SOM"))

# Pipeline-wide default scope (env-overridable):
#   "global" (default) | "no_neighbors" | "region" | "income" | "region_income"
DONOR_SCOPE <- Sys.getenv("CVD_DONOR_SCOPE", unset = "global")

# Build (ISO3 -> region, income_tier) once. Region via {countrycode} (World Bank
# 7-region; a defensible, available proxy for "same WHO region" -- swap a WHO map
# in here if you have one). Income tier derived from econ_panel GDP per capita
# using World Bank FY24 thresholds, so it needs no extra file. Assign the result
# to `.donor_meta` in 00 after 02b/03 are sourced.
build_donor_meta <- function(coverage_df,
                             econ = if (exists("econ_panel")) econ_panel else NULL,
                             gdp_max_year = NULL) {
  iso <- sort(unique(coverage_df$ISO3))
  region <- stats::setNames(rep(NA_character_, length(iso)), iso)
  if (requireNamespace("countrycode", quietly = TRUE))
    region[] <- countrycode::countrycode(iso, "iso3c", "region", warn = FALSE)
  tier <- stats::setNames(rep(NA_character_, length(iso)), iso)
  if (!is.null(econ) && "GDP_pc" %in% names(econ)) {
    g <- econ %>% dplyr::filter(!is.na(GDP_pc))
    if (!is.null(gdp_max_year)) g <- g %>% dplyr::filter(year <= gdp_max_year)
    g <- g %>% dplyr::group_by(ISO3) %>%
      dplyr::slice_max(year, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
    band <- cut(g$GDP_pc, c(-Inf, 1145, 4515, 14005, Inf),
                labels = c("LIC","LMIC","UMIC","HIC"))
    tier[g$ISO3] <- as.character(band)
  }
  tibble::tibble(ISO3 = iso, region = unname(region[iso]),
                 income_tier = unname(tier[iso]))
}

# Apply a scope to candidate donor ISO3s relative to ONE treated unit. Degrades
# gracefully: a dimension the metadata can't resolve is skipped, not emptied.
.apply_donor_scope <- function(cand_iso3, treated_iso3,
                               scope = DONOR_SCOPE,
                               meta = if (exists(".donor_meta")) .donor_meta else NULL,
                               neighbors = CONFLICT_NEIGHBORS) {
  if (identical(scope, "global") || is.null(treated_iso3)) return(cand_iso3)
  keep <- setdiff(cand_iso3, neighbors[[treated_iso3]])   # SUTVA drop for any non-global scope
  if (scope %in% c("region","region_income") && !is.null(meta)) {
    tr <- meta$region[match(treated_iso3, meta$ISO3)]
    if (!is.na(tr)) keep <- intersect(keep,
                                      meta$ISO3[!is.na(meta$region) & meta$region == tr])
  }
  if (scope %in% c("income","region_income") && !is.null(meta)) {
    ti <- meta$income_tier[match(treated_iso3, meta$ISO3)]
    if (!is.na(ti)) keep <- intersect(keep,
                                      meta$ISO3[!is.na(meta$income_tier) & meta$income_tier == ti])
  }
  keep
}

# Pooled-design variant (DiD): no single treated unit, so take the UNION of the
# per-cohort scoped pools (drops every cohort's neighbours; keeps any donor in
# at least one treated region/tier).
did_scope_union <- function(donor_iso3, treated_iso3_vec, scope = DONOR_SCOPE,
                            meta = if (exists(".donor_meta")) .donor_meta else NULL,
                            neighbors = CONFLICT_NEIGHBORS) {
  if (identical(scope, "global")) return(donor_iso3)
  Reduce(union, lapply(treated_iso3_vec,
                       function(t) .apply_donor_scope(donor_iso3, t, scope, meta, neighbors)))
}

# ---- SC/ASCM long-horizon flag (F3) ---------------------------------------
# Donor-comparison gaps grow monotonically with event time (ASCM ~9->33pp by
# year 17, all on Somalia's 19-yr tail), unlike the flat design-based methods.
# Mirror ITS: carry horizon_unreliable, and optionally DROP beyond a cap so the
# YLL chain isn't inflated by extrapolated long-horizon donor gaps. Default keeps
# all years (flag only) so behaviour is unchanged unless CVD_SC_MAX_HORIZON is set.
SC_RELIABLE_HORIZON <- as.integer(Sys.getenv("CVD_SC_RELIABLE_HORIZON", unset = "8"))
SC_MAX_HORIZON <- { .h <- Sys.getenv("CVD_SC_MAX_HORIZON", unset = "Inf")
if (identical(tolower(.h), "inf")) Inf else as.integer(.h) }

flag_sc_horizon <- function(gap_tbl, conflict_start,
                            reliable = SC_RELIABLE_HORIZON, max_h = SC_MAX_HORIZON) {
  if (nrow(gap_tbl) == 0) return(gap_tbl)
  gap_tbl <- gap_tbl %>%
    dplyr::mutate(event_time = year - conflict_start,
                  horizon_unreliable = event_time >= reliable)
  if (is.finite(max_h)) gap_tbl <- gap_tbl %>% dplyr::filter(event_time < max_h)
  dplyr::select(gap_tbl, -event_time)
}

# ---- METHOD 2: Synthetic Control via tidysynth ---------------------------
# Donor pool: non-conflict countries with full pre-period coverage data and
# meaningful variation (not stuck at the ceiling for the entire window), then
# the optional comparability/SUTVA scope above.
build_donor_pool <- function(coverage_df, vacc, exclude_iso3, pre_window,
                             treated_iso3 = NULL, scope = DONOR_SCOPE,
                             meta = if (exists(".donor_meta")) .donor_meta else NULL,
                             neighbors = CONFLICT_NEIGHBORS) {
  cov_ceiling <- if (vacc == "BCG") 99.5 else 98
  cand <- coverage_df %>%
    dplyr::filter(vaccine == vacc, year %in% pre_window,
                  !ISO3 %in% exclude_iso3) %>%
    dplyr::group_by(country, ISO3) %>%
    dplyr::summarise(
      n_obs    = dplyr::n(),
      mean_cov = mean(coverage, na.rm = TRUE),
      sd_cov   = stats::sd(coverage, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::filter(n_obs == length(pre_window),
                  mean_cov >= 30, mean_cov <= cov_ceiling,
                  !is.na(sd_cov))
  ok_iso <- .apply_donor_scope(cand$ISO3, treated_iso3, scope, meta, neighbors)
  cand %>% dplyr::filter(ISO3 %in% ok_iso) %>% dplyr::pull(country)
}

# ASSUMPTIONS / CHANGES vs. original:
#   * Donors must be complete over the FULL pre+post window (was: pre only).
#     A donor missing a post year would otherwise propagate NA into synth_y.
#   * Placebos are filtered by pre-period RMSPE before the gap SD is computed
#     (Abadie-style): placebos that fit their own pre-period poorly produce
#     spurious post gaps and inflate gap_sd. rmspe_ratio_cap controls the cut.
#   * Optional lagged-outcome predictors (use_lag_predictors) strengthen the
#     match beyond a single pre-window mean.
#   * Output carries pre_rmspe (treated pre-period fit) and n_good_placebos
#     as diagnostics.
gap_synthetic <- function(coverage_df, conflict_df,
                          vaccines = c("BCG", "MCV1", "DTP3"),
                          pre_len = 7L, min_donors = 10L,
                          rmspe_ratio_cap = 5,
                          use_lag_predictors = TRUE,
                          seed = DEFAULT_SEED) {
  set.seed(seed)
  conflict_iso3 <- conflict_df$ISO3
  out <- list()
  
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    # Shorter pre-window recovers donor pools for late-onset / data-sparse
    # conflicts, at the cost of matching quality -- a genuine trade-off.
    pre_window  <- (ci$conflict_start - pre_len):(ci$conflict_start - 1L)
    post_window <- ci$conflict_start:ci$conflict_end
    
    for (vacc in vaccines) {
      message("  SC: ", ci$country, " | ", vacc)
      donors <- build_donor_pool(coverage_df, vacc, conflict_iso3, pre_window,
                                 treated_iso3 = ci$ISO3, scope = DONOR_SCOPE)
      if (length(donors) < min_donors) {
        message("    skipping (too few donors: ", length(donors), ")")
        next
      }
      
      # Balanced panel of treated + candidate donors over the full window
      panel <- coverage_df %>%
        dplyr::filter(vaccine == vacc,
                      country %in% c(ci$country, donors),
                      year %in% c(pre_window, post_window)) %>%
        dplyr::select(country, year, coverage) %>%
        tidyr::complete(country, year) %>%
        dplyr::arrange(country, year)
      
      # Treated must be complete in pre-period (needed to fit).
      treated_pre_ok <- panel %>%
        dplyr::filter(country == ci$country, year %in% pre_window) %>%
        dplyr::summarise(ok = all(!is.na(coverage))) %>% dplyr::pull(ok)
      if (!isTRUE(treated_pre_ok)) {
        message("    skipping (treated has NA in pre-window)")
        next
      }
      
      # Donors must be complete over the FULL window (else NA leaks into synth).
      donor_full_ok <- panel %>%
        dplyr::filter(country != ci$country) %>%
        dplyr::group_by(country) %>%
        dplyr::summarise(ok = all(!is.na(coverage)), .groups = "drop") %>%
        dplyr::filter(ok) %>% dplyr::pull(country)
      if (length(donor_full_ok) < min_donors) {
        message("    skipping (insufficient full-window donors: ",
                length(donor_full_ok), ")")
        next
      }
      keep  <- c(ci$country, donor_full_ok)
      panel <- panel %>% dplyr::filter(country %in% keep)
      
      # Lag-predictor years (start / middle / end of pre-window), de-duplicated.
      y1 <- pre_window[1]
      y2 <- pre_window[ceiling(length(pre_window) / 2)]
      y3 <- pre_window[length(pre_window)]
      lag_years <- unique(c(y1, y2, y3))
      
      # Fit synthetic control
      sc_obj <- tryCatch({
        sc <- panel %>%
          tidysynth::synthetic_control(
            outcome = coverage, unit = country, time = year,
            i_unit = ci$country, i_time = ci$conflict_start,
            generate_placebos = TRUE
          ) %>%
          tidysynth::generate_predictor(
            time_window = pre_window,
            cov_mean = mean(coverage, na.rm = TRUE)
          )
        if (use_lag_predictors) {
          for (ly in lag_years) {
            sc <- tidysynth::generate_predictor(
              sc, time_window = ly,
              !!paste0("lag_", ly) := mean(coverage, na.rm = TRUE)
            )
          }
        }
        sc %>%
          tidysynth::generate_weights(
            optimization_window = pre_window,
            margin_ipop = 0.02, sigf_ipop = 7, bound_ipop = 6
          ) %>%
          tidysynth::generate_control()
      }, error = function(e) {
        message("    tidysynth error: ", conditionMessage(e)); NULL
      })
      if (is.null(sc_obj)) next
      
      # Single grab of treated + placebos; split by .placebo.
      full_sc <- tryCatch(
        tidysynth::grab_synthetic_control(sc_obj, placebo = TRUE),
        error = function(e) { message("    grab error: ", conditionMessage(e)); NULL }
      )
      if (is.null(full_sc) || nrow(full_sc) == 0) next
      if (!all(c(".id", ".placebo", "time_unit", "real_y", "synth_y") %in% names(full_sc))) {
        message("    unexpected grab schema; skipping"); next
      }
      full_sc <- dplyr::rename(full_sc, year = time_unit)
      
      # Treated gap (.placebo == 0)
      gap_tbl <- full_sc %>%
        dplyr::filter(.placebo == 0, year %in% post_window) %>%
        dplyr::transmute(
          country = ci$country, ISO3 = ci$ISO3, vaccine = vacc, year = year,
          coverage_obs = real_y, coverage_cf = synth_y,
          gap = synth_y - real_y
        )
      if (nrow(gap_tbl) == 0) next
      
      # Pre-period RMSPE per unit -> filter poorly-fitting placebos.
      rmspe_tbl <- full_sc %>%
        dplyr::filter(year %in% pre_window) %>%
        dplyr::group_by(.id, .placebo) %>%
        dplyr::summarise(
          pre_rmspe = sqrt(mean((real_y - synth_y)^2, na.rm = TRUE)),
          .groups = "drop"
        )
      treated_rmspe <- rmspe_tbl %>%
        dplyr::filter(.placebo == 0) %>% dplyr::pull(pre_rmspe)
      treated_rmspe <- if (length(treated_rmspe) == 1 && is.finite(treated_rmspe)) treated_rmspe else NA_real_
      
      good_ids <- rmspe_tbl %>%
        dplyr::filter(.placebo == 1, is.finite(pre_rmspe),
                      is.na(treated_rmspe) | pre_rmspe <= rmspe_ratio_cap * treated_rmspe) %>%
        dplyr::pull(.id)
      
      placebo_tbl <- full_sc %>%
        dplyr::filter(.placebo == 1, .id %in% good_ids, year %in% post_window) %>%
        dplyr::mutate(placebo_gap = synth_y - real_y) %>%
        dplyr::group_by(year) %>%
        dplyr::summarise(gap_sd = stats::sd(placebo_gap, na.rm = TRUE),
                         .groups = "drop")
      
      # Field-standard SC inference (Abadie et al. 2010; Abadie 2021): the
      # post/pre RMSPE RATIO, ranked against placebos. This corrects for
      # pre-period fit quality (a big post gap means little if the pre fit was
      # already poor) in a way a raw placebo SD does not. The permutation
      # p-value is the treated unit's rank among ALL units' ratios. One-sided
      # (shortfall) is the honest test here given a directional hypothesis.
      post_rmspe_tbl <- full_sc %>%
        dplyr::filter(year %in% post_window) %>%
        dplyr::group_by(.id, .placebo) %>%
        dplyr::summarise(post_rmspe = sqrt(mean((real_y - synth_y)^2, na.rm = TRUE)),
                         .groups = "drop")
      ratio_tbl <- rmspe_tbl %>%
        dplyr::left_join(post_rmspe_tbl, by = c(".id", ".placebo")) %>%
        dplyr::mutate(rmspe_ratio = post_rmspe / pre_rmspe) %>%
        dplyr::filter(is.finite(rmspe_ratio))
      treated_ratio <- ratio_tbl %>%
        dplyr::filter(.placebo == 0) %>% dplyr::pull(rmspe_ratio)
      treated_ratio <- if (length(treated_ratio) == 1) treated_ratio else NA_real_
      # rank-based permutation p-value over treated + all placebos
      perm_p <- if (is.na(treated_ratio)) NA_real_ else
        mean(ratio_tbl$rmspe_ratio >= treated_ratio, na.rm = TRUE)
      
      gap_tbl <- gap_tbl %>%
        dplyr::left_join(placebo_tbl, by = "year") %>%
        dplyr::mutate(
          gap_sd          = ifelse(is.na(gap_sd) | gap_sd == 0, 5.0, gap_sd),
          pre_rmspe       = treated_rmspe,
          rmspe_ratio     = treated_ratio,
          perm_pvalue     = perm_p,
          n_good_placebos = length(good_ids),
          n_donors_total  = length(donor_full_ok),
          method          = "Synthetic Control"
        ) %>%
        dplyr::select(country, ISO3, vaccine, year,
                      coverage_obs, coverage_cf, gap, gap_sd,
                      pre_rmspe, rmspe_ratio, perm_pvalue,
                      n_good_placebos, n_donors_total, method)
      gap_tbl <- flag_sc_horizon(gap_tbl, ci$conflict_start)   # F3: long-horizon flag
      out[[length(out) + 1]] <- gap_tbl
    }
  }
  dplyr::bind_rows(out)
}

# ---- METHOD 3: ITS with logit link ---------------------------------------
# Pre-period: logit(coverage/100) = b0 + b1*time
# Counterfactual for conflict year t: plogis(b0 + b1*t)
# Jensen correction: average plogis(draws), not plogis(mean).
#
# LITERATURE ALIGNMENT (ITS best practice):
#   * Annual coverage series are autocorrelated; OLS that ignores this yields
#     standard errors that are too small and CIs that under-cover (Bottomley
#     et al. 2023, Stat Med; Turner et al. 2021, BMC Med Res Methodol). For the
#     SHORT pre-segments here, those same papers find OLS / REML better-
#     calibrated than Newey-West HAC, which is unstable in short series. We
#     therefore select by length (vcov_type = "auto"): OLS for n_pre <
#     ITS_GLS_MIN_N, else GLS AR(1) by REML (Prais-Winsten-like). "NW" and a
#     forced "OLS"/"GLS" remain as explicit options. Coefficient draws use a
#     multivariate-t (df = residual df) for honest short-sample tails.
#   * Durbin-Watson is unreliable in short series (Turner et al. 2021); we flag
#     short segments via short_series / n_pre rather than relying on it.
#
# ASSUMPTIONS / CHANGES vs. original:
#   * The pre-conflict logit-linear trend would continue absent conflict.
#   * Forecast bootstrap adds RESIDUAL variance (include_resid_var), so cf_sd
#     is a predictive interval, not just a coefficient interval.
#   * Coefficient covariance reflects serial correlation via length-based method
#     selection (OLS for short segments, GLS AR(1)-REML for longer); draws use a
#     multivariate-t, which widens cf_sd appropriately in short series.
#   * NO one-sided level cap. The logit link already bounds the counterfactual
#     to (0,100); a prior cap at pre_max+2 was asymmetric (could only lower the
#     cf) and pinned the counterfactual below observed for improving countries
#     over long horizons. Instead, horizon is handled explicitly: years beyond
#     `reliable_horizon` are flagged (horizon_unreliable) and years beyond
#     `max_horizon` (default Inf) are dropped. cf_saturated reports how often the
#     extrapolated cf reaches the ceiling -- the honest signal that a no-control
#     forecast is overshooting and the comparison-group methods should be trusted.
gap_its <- function(coverage_df, conflict_df,
                    vaccines = c("BCG", "MCV1", "DTP3"),
                    n_pre_min = 6, n_boot = 1000, pre_lookback = 15L,
                    include_resid_var = TRUE,
                    reliable_horizon = 8L, max_horizon = ITS_MAX_HORIZON,
                    vcov_type = c("auto", "OLS", "GLS", "NW"),
                    seed = DEFAULT_SEED) {
  vcov_type <- match.arg(vcov_type)
  set.seed(seed)
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (vacc in vaccines) {
      message("  ITS: ", ci$country, " | ", vacc)
      ctry <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc) %>%
        dplyr::arrange(year)
      pre  <- ctry %>% dplyr::filter(year < ci$conflict_start,
                                     year >= ci$conflict_start - pre_lookback)
      post <- ctry %>% dplyr::filter(year >= ci$conflict_start,
                                     year <= ci$conflict_end)
      if (nrow(pre) < n_pre_min) {
        message("    skipping (only ", nrow(pre), " pre years)")
        next
      }
      
      df_fit <- pre %>%
        dplyr::mutate(
          p       = pmin(pmax(coverage, 0.5), 99.5) / 100,
          logit_p = stats::qlogis(p),
          time    = year - ci$conflict_start
        )
      fit <- tryCatch(stats::lm(logit_p ~ time, data = df_fit),
                      error = function(e) NULL)
      if (is.null(fit) || is.na(stats::coef(fit)[2])) {
        message("    fit failed"); next
      }
      
      coefs <- stats::coef(fit)
      n_pre <- nrow(pre)
      # Autocorrelation handling by series length (recommendation #6):
      #   short pre-segments  -> OLS, no autocorrelation adjustment;
      #   longer pre-segments  -> GLS AR(1) by REML (Prais-Winsten-like).
      # Newey-West is demoted to an explicit option ("NW"); it is unstable and
      # under-covers in short series (Bottomley 2023; Turner 2021).
      method_used <- vcov_type
      if (vcov_type == "auto")
        method_used <- if (n_pre < ITS_GLS_MIN_N) "OLS" else "GLS"
      
      V           <- stats::vcov(fit)
      sigma_resid <- summary(fit)$sigma            # residual SD in logit space
      resid_df    <- fit$df.residual
      
      if (method_used == "GLS" && requireNamespace("nlme", quietly = TRUE) &&
          n_pre >= 5) {
        gls_fit <- tryCatch(
          nlme::gls(logit_p ~ time, data = df_fit,
                    correlation = nlme::corAR1(form = ~ time),
                    method = "REML"),
          error = function(e) NULL)
        if (!is.null(gls_fit)) {
          coefs       <- stats::coef(gls_fit)
          V           <- as.matrix(stats::vcov(gls_fit))
          sigma_resid <- gls_fit$sigma
          resid_df    <- n_pre - 2L
        } else {
          method_used <- "OLS(GLS-failed)"
        }
      } else if (method_used == "NW" &&
                 requireNamespace("sandwich", quietly = TRUE)) {
        # Retained only as an explicit comparison option, not the default.
        V <- tryCatch(sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE),
                      error = function(e) stats::vcov(fit))
      }
      
      # Multivariate-t coefficient draws widen intervals to reflect estimated-
      # variance uncertainty in short segments; ITS CIs still tend to under-
      # cover (Turner 2021), which is reported as a caveat in the diagnostics.
      beta_draws  <- mvt_draws(n_boot, mu = coefs, Sigma = V,
                               df = max(resid_df, 1))
      
      # Extrapolation horizon. ITS forecasts a within-unit pre-trend forward, so
      # its credibility decays with distance from the pre-period. We do NOT cap
      # the counterfactual level (the previous one-sided cap at pre_max+2 could
      # only push the counterfactual DOWN: for an improving country the logit
      # trend clears pre_max almost at once, so the cap pinned the cf near the
      # pre-conflict level while observed coverage kept rising past it, producing
      # a spurious negative gap -- e.g. Ethiopia, 2012-2024). The logit link
      # already bounds the cf to (0, 100); instead we (a) flag years beyond
      # `reliable_horizon` as horizon_unreliable, and (b) optionally drop years
      # beyond `max_horizon` (default Inf = keep all, flagged) so long-conflict
      # ITS estimates can be down-weighted or excluded rather than fabricated.
      post_years_all <- post$year
      keep_h <- (post_years_all - ci$conflict_start) < max_horizon
      post   <- post[keep_h, , drop = FALSE]
      post_years <- post$year
      if (length(post_years) == 0) {
        message("    skipping (no post years within max_horizon)"); next
      }
      
      cf_mean   <- numeric(length(post_years))
      cf_sd     <- numeric(length(post_years))
      cf_satur  <- numeric(length(post_years))
      for (k in seq_along(post_years)) {
        t_k <- post_years[k] - ci$conflict_start
        lin <- beta_draws[, 1] + beta_draws[, 2] * t_k
        # Predictive draw: coefficient + residual uncertainty, then Jensen
        # back-transform (average plogis(draws), not plogis(mean)).
        eps <- if (include_resid_var) stats::rnorm(n_boot, 0, sigma_resid) else 0
        p_draws    <- stats::plogis(lin + eps)
        cf_satur[k] <- mean(p_draws > 0.98)   # how often the cf saturates the ceiling
        cf_mean[k] <- mean(p_draws) * 100
        cf_sd[k]   <- stats::sd(p_draws) * 100
      }
      
      gap_tbl <- tibble::tibble(
        country      = ci$country,
        ISO3         = ci$ISO3,
        vaccine      = vacc,
        year         = post_years,
        coverage_obs = post$coverage,
        coverage_cf  = cf_mean,
        gap          = cf_mean - post$coverage,
        gap_sd       = cf_sd,
        cf_saturated = cf_satur,                       # cf pinned near 100% (overshoot signal)
        slope_logit_per_yr = unname(coefs[2]),         # pre-trend slope (logit/yr)
        horizon_unreliable = (post_years - ci$conflict_start) >= reliable_horizon,
        resid_df     = resid_df,
        n_pre        = n_pre,
        short_series = n_pre < ITS_GLS_MIN_N,  # below the GLS threshold: fragile
        vcov_type    = method_used,            # resolved method (OLS / GLS / NW)
        method       = "ITS"
      )
      out[[length(out) + 1]] <- gap_tbl
    }
  }
  dplyr::bind_rows(out)
}

# ---- MASTER: estimate all three methods ----------------------------------
estimate_all_gaps <- function(coverage_df, conflict_df,
                              vaccines = c("BCG", "MCV1", "DTP3"),
                              seed = DEFAULT_SEED) {
  # Harmonization check: Baseline uses preconflict_year; SC/ITS use
  # conflict_start. Warn if they disagree so gaps aren't compared across
  # mismatched pre/post boundaries.
  if ("preconflict_year" %in% names(conflict_df)) {
    mismatch <- conflict_df %>%
      dplyr::filter(preconflict_year != conflict_start - 1L)
    if (nrow(mismatch) > 0) {
      warning("preconflict_year != conflict_start - 1 for ",
              nrow(mismatch), " conflict row(s); Baseline and SC/ITS ",
              "are anchored differently for these.")
    }
  }
  
  message("Estimating coverage gaps with three methods...")
  message("== Method 1: Baseline ==")
  b <- gap_baseline(coverage_df, conflict_df, vaccines)
  message("== Method 2: Synthetic Control ==")
  s <- gap_synthetic(coverage_df, conflict_df, vaccines, seed = seed)
  message("== Method 3: ITS ==")
  i <- gap_its(coverage_df, conflict_df, vaccines, seed = seed)
  
  dplyr::bind_rows(b, s, i) %>%
    dplyr::mutate(
      gap_sd            = pmax(gap_sd, 1.0),     # floor SD at 1pp
      gap_prop          = gap / 100,
      gap_sd_prop       = gap_sd / 100,
      coverage_obs_prop = coverage_obs / 100,
      coverage_cf_prop  = coverage_cf / 100
    )
}

# ============================================================================
# DIAGNOSTIC / SENSITIVITY TESTS
# ----------------------------------------------------------------------------
# These are validation tools, not part of the production estimate. Run them
# after estimate_all_gaps() to sanity-check the counterfactuals.
# ============================================================================

# ---- D1: Placebo-in-time --------------------------------------------------
# Pretend each conflict happened entirely BEFORE its real start, in a window
# where no conflict occurred. A well-behaved method should return gaps ~ 0
# there. Large systematic placebo gaps => the counterfactual model is picking
# up trend/noise, not conflict. (SC omitted here: expensive; rerun gap_synthetic
# on the shifted conflict_df if needed.)
diag_placebo_in_time <- function(coverage_df, conflict_df,
                                 method = c("baseline", "its"),
                                 vaccines = c("BCG", "MCV1", "DTP3"),
                                 seed = DEFAULT_SEED) {
  method <- match.arg(method)
  dur <- pmax(conflict_df$conflict_end - conflict_df$conflict_start, 1L)
  fake <- conflict_df %>%
    dplyr::mutate(
      .dur             = dur,
      conflict_end     = conflict_start - 1L,
      conflict_start   = conflict_start - .dur,
      preconflict_year = conflict_start - 1L
    ) %>%
    dplyr::select(-.dur)
  
  res <- if (method == "baseline") {
    gap_baseline(coverage_df, fake, vaccines)
  } else {
    gap_its(coverage_df, fake, vaccines, seed = seed)
  }
  if (nrow(res) == 0) return(tibble::tibble())
  
  res %>%
    dplyr::summarise(
      method        = method,
      n             = dplyr::n(),
      mean_placebo_gap   = mean(gap, na.rm = TRUE),
      median_placebo_gap = stats::median(gap, na.rm = TRUE),
      sd_placebo_gap     = stats::sd(gap, na.rm = TRUE),
      # share of placebo gaps that are "large" relative to a 5pp tolerance
      frac_abs_gt_5pp    = mean(abs(gap) > 5, na.rm = TRUE)
    )
}

# ---- D2: Cross-method agreement -------------------------------------------
# Pairwise correlation, bias (method A - method B), and mean abs difference on
# matched country-vaccine-year cells. High disagreement flags model-driven
# (not data-driven) gaps.
diag_cross_method <- function(all_gaps) {
  w <- all_gaps %>%
    dplyr::mutate(m = dplyr::case_when(
      method %in% c("Baseline", "Baseline (3yr)") ~ "baseline",
      method == "Synthetic Control" ~ "sc",
      method == "ITS"               ~ "its",
      TRUE                          ~ method
    )) %>%
    dplyr::select(country, vaccine, year, m, gap) %>%
    tidyr::pivot_wider(names_from = m, values_from = gap)
  
  pairs <- list(c("baseline", "sc"), c("baseline", "its"), c("sc", "its"))
  res <- lapply(pairs, function(p) {
    if (!all(p %in% names(w))) return(NULL)
    d <- w[stats::complete.cases(w[, p]), , drop = FALSE]
    if (nrow(d) < 3) return(NULL)
    tibble::tibble(
      pair          = paste(p[1], "vs", p[2]),
      n             = nrow(d),
      pearson_r     = stats::cor(d[[p[1]]], d[[p[2]]]),
      bias_A_minus_B = mean(d[[p[1]]] - d[[p[2]]]),
      mean_abs_diff = mean(abs(d[[p[1]]] - d[[p[2]]]))
    )
  })
  dplyr::bind_rows(res)
}

# ---- D3: ITS residual / autocorrelation report ----------------------------
# Surfaces pre-period fit health: residual df (extrapolation risk), residual
# sigma (forecast width), lag-1 autocorrelation and Durbin-Watson (independence
# assumption). NOTE (Turner et al. 2021, BMC MRM): the DW test is unreliable in
# short series, so we report it but key the recommendation off series length and
# the AR(1) magnitude rather than DW alone.
diag_its_residuals <- function(coverage_df, conflict_df,
                               vaccines = c("BCG", "MCV1", "DTP3"),
                               pre_lookback = 15L) {
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (vacc in vaccines) {
      pre <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc,
                      year < ci$conflict_start,
                      year >= ci$conflict_start - pre_lookback) %>%
        dplyr::arrange(year)
      if (nrow(pre) < 3) next
      df_fit <- pre %>%
        dplyr::mutate(p = pmin(pmax(coverage, 0.5), 99.5) / 100,
                      logit_p = stats::qlogis(p),
                      time = year - ci$conflict_start)
      fit <- tryCatch(stats::lm(logit_p ~ time, data = df_fit),
                      error = function(e) NULL)
      if (is.null(fit)) next
      r   <- stats::residuals(fit)
      ac1 <- if (length(r) > 2) stats::cor(r[-1], r[-length(r)]) else NA_real_
      dw  <- if (length(r) > 2) sum(diff(r)^2) / sum(r^2) else NA_real_
      n_pre <- nrow(pre)
      # Recommendation per ITS literature (Bottomley 2023; Turner 2021): for
      # short pre-segments use OLS (no autocorrelation adjustment) or REML; for
      # longer segments use GLS AR(1) / Prais-Winsten. gap_its now selects OLS
      # vs GLS AR(1)-REML by series length (threshold ITS_GLS_MIN_N); Newey-West
      # is an option, not the default. All ITS intervals tend to under-cover.
      rec <- if (n_pre < ITS_GLS_MIN_N)
        "short series: OLS/REML preferred (no HAC); CIs under-cover"
      else if (!is.na(ac1) && abs(ac1) > 0.3)
        "autocorrelation present: GLS AR(1)/Prais-Winsten (gap_its default at this length)"
      else "GLS AR(1)/OLS adequate"
      out[[length(out) + 1]] <- tibble::tibble(
        country            = ci$country,
        vaccine            = vacc,
        n_pre              = n_pre,
        resid_df           = fit$df.residual,
        slope_logit_per_yr = unname(stats::coef(fit)[2]),
        resid_sigma_logit  = summary(fit)$sigma,
        resid_acf1         = ac1,
        durbin_watson      = dw,
        flag_short_series  = n_pre < ITS_GLS_MIN_N,
        flag_autocorr      = !is.na(ac1) && abs(ac1) > 0.3,
        recommendation     = rec
      )
    }
  }
  dplyr::bind_rows(out)
}

# ---- D4: Gap distribution summary -----------------------------------------
# Per-method shape check: a large share of negative gaps (coverage above the
# counterfactual) or implausibly large gaps warrants scrutiny.
diag_gap_summary <- function(all_gaps) {
  all_gaps %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n            = dplyr::n(),
      mean_gap     = mean(gap, na.rm = TRUE),
      median_gap   = stats::median(gap, na.rm = TRUE),
      sd_gap       = stats::sd(gap, na.rm = TRUE),
      frac_negative = mean(gap < 0, na.rm = TRUE),
      frac_gt_30pp = mean(gap > 30, na.rm = TRUE),
      mean_gap_sd  = mean(gap_sd, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- D5: ITS pre-window sensitivity ---------------------------------------
# Re-estimate ITS under several lookback lengths and compare per-cell gaps.
# Large swings mean the extrapolation is window-sensitive (fragile).
diag_its_window_sensitivity <- function(coverage_df, conflict_df,
                                        vaccines = c("BCG", "MCV1", "DTP3"),
                                        lookbacks = c(10L, 15L, 20L),
                                        seed = DEFAULT_SEED) {
  runs <- lapply(lookbacks, function(L) {
    g <- gap_its(coverage_df, conflict_df, vaccines,
                 pre_lookback = L, seed = seed)
    if (nrow(g) == 0) return(NULL)
    g %>% dplyr::transmute(country, vaccine, year, lookback = L, gap)
  })
  runs <- dplyr::bind_rows(runs)
  if (nrow(runs) == 0) return(tibble::tibble())
  runs %>%
    tidyr::pivot_wider(names_from = lookback, values_from = gap,
                       names_prefix = "lb_") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      gap_range = {
        v <- dplyr::c_across(dplyr::starts_with("lb_"))
        v <- v[!is.na(v)]
        if (length(v) < 2) NA_real_ else max(v) - min(v)
      }
    ) %>%
    dplyr::ungroup()
}

# ---- Diagnostics runner ---------------------------------------------------
run_all_diagnostics <- function(coverage_df, conflict_df, all_gaps,
                                seed = DEFAULT_SEED) {
  list(
    placebo_baseline   = diag_placebo_in_time(coverage_df, conflict_df, "baseline", seed = seed),
    placebo_its        = diag_placebo_in_time(coverage_df, conflict_df, "its", seed = seed),
    cross_method       = diag_cross_method(all_gaps),
    its_residuals      = diag_its_residuals(coverage_df, conflict_df),
    gap_summary        = diag_gap_summary(all_gaps),
    its_window_sens    = diag_its_window_sensitivity(coverage_df, conflict_df, seed = seed)
  )
}

# ---- Example usage (not executed on source) -------------------------------
if (FALSE) {
  library(dplyr); library(tidyr); library(tibble)
  # requires: tidysynth, MASS, nlme; suggests: sandwich (NW option only)
  all_gaps <- estimate_all_gaps(coverage_df, conflict_df)
  diags    <- run_all_diagnostics(coverage_df, conflict_df, all_gaps)
  
  print(diags$gap_summary)        # shape / negative-gap share by method
  print(diags$cross_method)       # do the three methods agree?
  print(diags$placebo_baseline)   # should be ~0 if counterfactual is sound
  print(diags$placebo_its)        # should be ~0
  print(diags$its_residuals)      # df / autocorrelation flags + recommendation
  print(diags$its_window_sens)    # is ITS lookback-sensitive?
  
  # SC effects with field-standard inference: filter to credible permutation
  # p-values rather than reading every gap as causal.
  all_gaps %>%
    dplyr::filter(method == "Synthetic Control") %>%
    dplyr::distinct(country, vaccine, rmspe_ratio, perm_pvalue, n_good_placebos) %>%
    dplyr::arrange(perm_pvalue) %>%
    print(n = 50)
  
  # Compare ITS uncertainty across autocorrelation handling. Default ("auto")
  # picks OLS (short) or GLS AR(1)-REML (longer) by series length.
  its_auto <- gap_its(coverage_df, conflict_df)                    # length-based
  its_gls  <- gap_its(coverage_df, conflict_df, vcov_type = "GLS") # force GLS AR(1)
  its_ols  <- gap_its(coverage_df, conflict_df, vcov_type = "OLS") # force OLS
}