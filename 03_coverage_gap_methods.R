# 03_gaps.R | in: coverage_long, conflict_info, econ_panel | out: estimate_all_gaps_plus() -> gap_df (8 estimators) | Methods: Estimation of conflict-associated coverage loss
# gap = counterfactual - observed (pp); positive = shortfall. Each estimator reseeds with SEED_GAP,
# so estimators are order-independent (and can run in parallel with identical results).

METHOD_ITS <- "ITS"; METHOD_SC <- "Synthetic Control"; METHOD_SC_COV <- "SC + covariates"
METHOD_AUGSC <- "Augmented SC"; METHOD_DID <- "DiD (CS)"
BASELINE_WINDOWS <- c(1L, 3L, 5L)            # 3-yr = headline baseline
baseline_label   <- function(w) paste0("Baseline (", w, "yr)")
METHOD_BASELINE_LABELS <- baseline_label(BASELINE_WINDOWS)
ALL_METHODS <- c(METHOD_BASELINE_LABELS, METHOD_ITS, METHOD_SC, METHOD_SC_COV, METHOD_AUGSC, METHOD_DID)

ITS_GLS_MIN_N      <- 12L   # >= 12 pre-years: GLS AR(1) by REML; else OLS
ITS_MAX_HORIZON    <- 8L    # ITS years >= 8 after onset dropped
SC_RELIABLE_HORIZON <- 8L   # donor-based years >= 8 after onset flagged (not dropped)

# multivariate-t coefficient draws (df = residual df)
mvt_draws <- function(n, mu, Sigma, df) {
  z <- MASS::mvrnorm(n, mu = rep(0, length(mu)), Sigma = Sigma)
  if (is.matrix(z) && is.finite(df) && df > 0) z <- z * sqrt(df / stats::rchisq(n, df = df))
  sweep(z, 2, mu, "+")
}

# ---- donor scope: "global" (headline) | "region_income" (robustness: no neighbours, same WB region + income tier) ----
DONOR_SCOPE <- Sys.getenv("CVD_DONOR_SCOPE", unset = "global")
stopifnot(DONOR_SCOPE %in% c("global", "region_income"))
CONFLICT_NEIGHBORS <- list(
  BFA = c("MLI","NER","BEN","TGO","GHA","CIV"), ETH = c("ERI","DJI","SOM","KEN","SSD","SDN"),
  IRQ = c("SYR","TUR","IRN","KWT","SAU","JOR"), MMR = c("BGD","IND","CHN","LAO","THA"),
  NGA = c("BEN","NER","TCD","CMR"),             PAK = c("AFG","IRN","IND","CHN"),
  SOM = c("ETH","DJI","KEN","YEM"),             LKA = c("IND","MDV"),
  SYR = c("TUR","IRQ","JOR","LBN","ISR"),       YEM = c("SAU","OMN","DJI","SOM"))

# ISO3 -> WB region (countrycode) + income tier (latest GDP pc, WB FY24 thresholds)
build_donor_meta <- function(coverage_df, econ = econ_panel) {
  iso <- sort(unique(coverage_df$ISO3))
  region <- stats::setNames(countrycode::countrycode(iso, "iso3c", "region", warn = FALSE), iso)
  tier <- stats::setNames(rep(NA_character_, length(iso)), iso)
  g <- econ %>% dplyr::filter(!is.na(GDP_pc)) %>% dplyr::group_by(ISO3) %>%
    dplyr::slice_max(year, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
  tier[g$ISO3] <- as.character(cut(g$GDP_pc, c(-Inf, 1145, 4515, 14005, Inf),
                                   labels = c("LIC","LMIC","UMIC","HIC")))
  tibble::tibble(ISO3 = iso, region = unname(region[iso]), income_tier = unname(tier[iso]))
}

# restrict candidate donors for one treated unit (an unresolved dimension is skipped)
.apply_donor_scope <- function(cand_iso3, treated_iso3, scope = DONOR_SCOPE,
                               meta = .donor_meta, neighbors = CONFLICT_NEIGHBORS) {
  if (identical(scope, "global") || is.null(treated_iso3)) return(cand_iso3)
  keep <- setdiff(cand_iso3, neighbors[[treated_iso3]])
  tr <- meta$region[match(treated_iso3, meta$ISO3)]
  if (!is.na(tr)) keep <- intersect(keep, meta$ISO3[!is.na(meta$region) & meta$region == tr])
  ti <- meta$income_tier[match(treated_iso3, meta$ISO3)]
  if (!is.na(ti)) keep <- intersect(keep, meta$ISO3[!is.na(meta$income_tier) & meta$income_tier == ti])
  keep
}
# pooled DiD: union of per-cohort scoped pools
did_scope_union <- function(donor_iso3, treated_iso3_vec, scope = DONOR_SCOPE) {
  if (identical(scope, "global")) return(donor_iso3)
  Reduce(union, lapply(treated_iso3_vec, function(t) .apply_donor_scope(donor_iso3, t, scope)))
}

flag_sc_horizon <- function(gap_tbl, conflict_start, reliable = SC_RELIABLE_HORIZON) {
  if (nrow(gap_tbl) == 0) return(gap_tbl)
  gap_tbl %>% dplyr::mutate(horizon_unreliable = (year - conflict_start) >= reliable)
}

# ---- Baseline (w-yr pre-conflict mean, held flat) ----
.gap_baseline_w <- function(coverage_df, conflict_df, vaccines, w) {
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (vacc in vaccines) {
      yrs  <- (ci$preconflict_year - (w - 1L)):ci$preconflict_year
      vals <- coverage_df %>% dplyr::filter(country == ci$country, vaccine == vacc, year %in% yrs) %>%
        dplyr::pull(coverage)
      if (length(vals) == 0) next
      pre_mu <- mean(vals, na.rm = TRUE)
      pre_sd <- if (length(vals) < 2) 2.0 else max(stats::sd(vals, na.rm = TRUE), 2.0)
      ctry_yrs <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc,
                      year >= ci$conflict_start, year <= ci$conflict_end) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs = coverage)
      if (nrow(ctry_yrs) == 0) next
      out[[length(out) + 1]] <- ctry_yrs %>%
        dplyr::mutate(coverage_cf = pre_mu, gap = pre_mu - coverage_obs, gap_sd = pre_sd)
    }
  }
  dplyr::bind_rows(out)
}

# ---- Synthetic control donor pool: non-conflict, complete pre-window, mean 30%..98% (99.5% BCG) ----
build_donor_pool <- function(coverage_df, vacc, exclude_iso3, pre_window, treated_iso3 = NULL) {
  cov_ceiling <- if (vacc == "BCG") 99.5 else 98
  cand <- coverage_df %>%
    dplyr::filter(vaccine == vacc, year %in% pre_window, !ISO3 %in% exclude_iso3) %>%
    dplyr::group_by(country, ISO3) %>%
    dplyr::summarise(n_obs = dplyr::n(), mean_cov = mean(coverage, na.rm = TRUE),
                     sd_cov = stats::sd(coverage, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(n_obs == length(pre_window), mean_cov >= 30, mean_cov <= cov_ceiling, !is.na(sd_cov))
  ok_iso <- .apply_donor_scope(cand$ISO3, treated_iso3)
  cand %>% dplyr::filter(ISO3 %in% ok_iso) %>% dplyr::pull(country)
}

# shared post-fit block for SC / SC+cov: treated gap, placebo-SD (RMSPE-filtered), RMSPE-ratio p
.sc_gap_table <- function(full_sc, ci, vacc, pre_window, post_window, rmspe_ratio_cap,
                          n_donors, method_label) {
  full_sc <- dplyr::rename(full_sc, year = time_unit)
  gap_tbl <- full_sc %>%
    dplyr::filter(.placebo == 0, year %in% post_window) %>%
    dplyr::transmute(country = ci$country, ISO3 = ci$ISO3, vaccine = vacc, year = year,
                     coverage_obs = real_y, coverage_cf = synth_y, gap = synth_y - real_y)
  if (nrow(gap_tbl) == 0) return(NULL)
  rmspe_tbl <- full_sc %>% dplyr::filter(year %in% pre_window) %>%
    dplyr::group_by(.id, .placebo) %>%
    dplyr::summarise(pre_rmspe = sqrt(mean((real_y - synth_y)^2, na.rm = TRUE)), .groups = "drop")
  treated_rmspe <- rmspe_tbl %>% dplyr::filter(.placebo == 0) %>% dplyr::pull(pre_rmspe)
  treated_rmspe <- if (length(treated_rmspe) == 1 && is.finite(treated_rmspe)) treated_rmspe else NA_real_
  good_ids <- rmspe_tbl %>%
    dplyr::filter(.placebo == 1, is.finite(pre_rmspe),
                  is.na(treated_rmspe) | pre_rmspe <= rmspe_ratio_cap * treated_rmspe) %>%
    dplyr::pull(.id)
  placebo_tbl <- full_sc %>%
    dplyr::filter(.placebo == 1, .id %in% good_ids, year %in% post_window) %>%
    dplyr::mutate(placebo_gap = synth_y - real_y) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(gap_sd = stats::sd(placebo_gap, na.rm = TRUE), .groups = "drop")
  post_rmspe_tbl <- full_sc %>% dplyr::filter(year %in% post_window) %>%
    dplyr::group_by(.id, .placebo) %>%
    dplyr::summarise(post_rmspe = sqrt(mean((real_y - synth_y)^2, na.rm = TRUE)), .groups = "drop")
  ratio_tbl <- rmspe_tbl %>% dplyr::left_join(post_rmspe_tbl, by = c(".id", ".placebo")) %>%
    dplyr::mutate(rmspe_ratio = post_rmspe / pre_rmspe) %>% dplyr::filter(is.finite(rmspe_ratio))
  treated_ratio <- ratio_tbl %>% dplyr::filter(.placebo == 0) %>% dplyr::pull(rmspe_ratio)
  treated_ratio <- if (length(treated_ratio) == 1) treated_ratio else NA_real_
  perm_p <- if (is.na(treated_ratio)) NA_real_ else mean(ratio_tbl$rmspe_ratio >= treated_ratio, na.rm = TRUE)
  gap_tbl %>%
    dplyr::left_join(placebo_tbl, by = "year") %>%
    dplyr::mutate(gap_sd = ifelse(is.na(gap_sd) | gap_sd == 0, 5.0, gap_sd),
                  pre_rmspe = treated_rmspe, rmspe_ratio = treated_ratio, perm_pvalue = perm_p,
                  n_good_placebos = length(good_ids), n_donors_total = n_donors,
                  method = method_label) %>%
    dplyr::select(country, ISO3, vaccine, year, coverage_obs, coverage_cf, gap, gap_sd,
                  pre_rmspe, rmspe_ratio, perm_pvalue, n_good_placebos, n_donors_total, method) %>%
    flag_sc_horizon(ci$conflict_start)
}

# ---- Synthetic control (tidysynth): 7-yr pre-window; predictors = pre mean + 3 lagged values ----
# A cell whose QP fit errors is skipped (logged), as in the working pipeline.
gap_synthetic <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                          pre_len = 7L, min_donors = 10L, rmspe_ratio_cap = 5, seed = SEED_GAP) {
  set.seed(seed)
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    pre_window  <- (ci$conflict_start - pre_len):(ci$conflict_start - 1L)
    post_window <- ci$conflict_start:ci$conflict_end
    for (vacc in vaccines) {
      message("  SC: ", ci$country, " | ", vacc)
      donors <- build_donor_pool(coverage_df, vacc, conflict_df$ISO3, pre_window, treated_iso3 = ci$ISO3)
      if (length(donors) < min_donors) { message("    skipping (too few donors: ", length(donors), ")"); next }
      panel <- coverage_df %>%
        dplyr::filter(vaccine == vacc, country %in% c(ci$country, donors),
                      year %in% c(pre_window, post_window)) %>%
        dplyr::select(country, year, coverage) %>%
        tidyr::complete(country, year) %>% dplyr::arrange(country, year)
      treated_pre_ok <- panel %>% dplyr::filter(country == ci$country, year %in% pre_window) %>%
        dplyr::summarise(ok = all(!is.na(coverage))) %>% dplyr::pull(ok)
      if (!isTRUE(treated_pre_ok)) { message("    skipping (treated has NA in pre-window)"); next }
      donor_full_ok <- panel %>% dplyr::filter(country != ci$country) %>%
        dplyr::group_by(country) %>% dplyr::summarise(ok = all(!is.na(coverage)), .groups = "drop") %>%
        dplyr::filter(ok) %>% dplyr::pull(country)
      if (length(donor_full_ok) < min_donors) {
        message("    skipping (insufficient full-window donors: ", length(donor_full_ok), ")"); next }
      panel <- panel %>% dplyr::filter(country %in% c(ci$country, donor_full_ok))
      lag_years <- unique(c(pre_window[1], pre_window[ceiling(length(pre_window) / 2)],
                            pre_window[length(pre_window)]))
      sc_obj <- tryCatch({
        sc <- panel %>%
          tidysynth::synthetic_control(outcome = coverage, unit = country, time = year,
                                       i_unit = ci$country, i_time = ci$conflict_start,
                                       generate_placebos = TRUE) %>%
          tidysynth::generate_predictor(time_window = pre_window, cov_mean = mean(coverage, na.rm = TRUE))
        for (ly in lag_years)
          sc <- tidysynth::generate_predictor(sc, time_window = ly,
                                              !!paste0("lag_", ly) := mean(coverage, na.rm = TRUE))
        sc %>%
          tidysynth::generate_weights(optimization_window = pre_window,
                                      margin_ipop = 0.02, sigf_ipop = 7, bound_ipop = 6) %>%
          tidysynth::generate_control()
      }, error = function(e) { message("    tidysynth error: ", conditionMessage(e)); NULL })
      if (is.null(sc_obj)) next
      full_sc <- tidysynth::grab_synthetic_control(sc_obj, placebo = TRUE)
      if (nrow(full_sc) == 0) next
      g <- .sc_gap_table(full_sc, ci, vacc, pre_window, post_window, rmspe_ratio_cap,
                         length(donor_full_ok), METHOD_SC)
      if (!is.null(g)) out[[length(out) + 1]] <- g
    }
  }
  dplyr::bind_rows(out)
}

# ---- SC + covariates: as SC plus pre-window means of log GDP pc, schooling, HDI ----
gap_synthetic_cov <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                              econ = econ_panel, predictors = ECON_PREDICTORS,
                              pre_len = 7L, min_donors = 10L, rmspe_ratio_cap = 5, seed = SEED_GAP) {
  set.seed(seed)
  econ_use <- econ %>% dplyr::select(dplyr::all_of(c("country", "year", intersect(predictors, names(econ))))) %>%
    dplyr::distinct()
  preds <- intersect(predictors, names(econ_use))
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    pre_window  <- (ci$conflict_start - pre_len):(ci$conflict_start - 1L)
    post_window <- ci$conflict_start:ci$conflict_end
    for (vacc in vaccines) {
      message("  SC+cov: ", ci$country, " | ", vacc)
      donors <- build_donor_pool(coverage_df, vacc, conflict_df$ISO3, pre_window, treated_iso3 = ci$ISO3)
      if (length(donors) < min_donors) { message("    skipping (too few donors: ", length(donors), ")"); next }
      panel <- coverage_df %>%
        dplyr::filter(vaccine == vacc, country %in% c(ci$country, donors),
                      year %in% c(pre_window, post_window)) %>%
        dplyr::select(country, year, coverage) %>%
        dplyr::left_join(econ_use, by = c("country", "year")) %>%
        tidyr::complete(country, year) %>% dplyr::arrange(country, year)
      treated_pre_ok <- panel %>% dplyr::filter(country == ci$country, year %in% pre_window) %>%
        dplyr::summarise(ok = all(!is.na(coverage))) %>% dplyr::pull(ok)
      if (!isTRUE(treated_pre_ok)) { message("    skipping (treated NA in pre-window)"); next }
      # donors: outcome complete over full window + >=1 non-missing value per covariate in pre-window
      out_ok_ctry <- panel %>% dplyr::filter(country != ci$country) %>% dplyr::group_by(country) %>%
        dplyr::summarise(out_ok = all(!is.na(coverage)), .groups = "drop")
      cov_ok_ctry <- panel %>% dplyr::filter(country != ci$country, year %in% pre_window) %>%
        dplyr::group_by(country) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(preds), ~ any(!is.na(.x)), .names = "ok_{.col}"),
                         .groups = "drop")
      cov_ok_ctry$cov_ok <- rowSums(as.matrix(cov_ok_ctry[paste0("ok_", preds)])) == length(preds)
      donor_ok <- out_ok_ctry %>%
        dplyr::left_join(cov_ok_ctry[c("country", "cov_ok")], by = "country") %>%
        dplyr::mutate(cov_ok = !is.na(cov_ok) & cov_ok) %>%
        dplyr::filter(out_ok, cov_ok) %>% dplyr::pull(country)
      if (length(donor_ok) < min_donors) {
        message("    skipping (insufficient covariate-complete donors: ", length(donor_ok), ")"); next }
      panel <- panel %>% dplyr::filter(country %in% c(ci$country, donor_ok))
      lag_years <- unique(c(pre_window[1], pre_window[ceiling(length(pre_window) / 2)],
                            pre_window[length(pre_window)]))
      sc_obj <- tryCatch({
        sc <- panel %>%
          tidysynth::synthetic_control(outcome = coverage, unit = country, time = year,
                                       i_unit = ci$country, i_time = ci$conflict_start,
                                       generate_placebos = TRUE) %>%
          tidysynth::generate_predictor(time_window = pre_window, cov_mean = mean(coverage, na.rm = TRUE))
        for (ly in lag_years)
          sc <- tidysynth::generate_predictor(sc, time_window = ly,
                                              !!paste0("lag_", ly) := mean(coverage, na.rm = TRUE))
        if ("log_GDP_pc" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window, p_log_gdp = mean(log_GDP_pc, na.rm = TRUE))
        if ("sch_eys" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window, p_eys = mean(sch_eys, na.rm = TRUE))
        if ("sch_mys" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window, p_mys = mean(sch_mys, na.rm = TRUE))
        if ("hdi" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window, p_hdi = mean(hdi, na.rm = TRUE))
        sc %>%
          tidysynth::generate_weights(optimization_window = pre_window,
                                      margin_ipop = 0.02, sigf_ipop = 7, bound_ipop = 6) %>%
          tidysynth::generate_control()
      }, error = function(e) { message("    tidysynth error: ", conditionMessage(e)); NULL })
      if (is.null(sc_obj)) next
      full_sc <- tidysynth::grab_synthetic_control(sc_obj, placebo = TRUE)
      if (nrow(full_sc) == 0) next
      g <- .sc_gap_table(full_sc, ci, vacc, pre_window, post_window, rmspe_ratio_cap,
                         length(donor_ok), METHOD_SC_COV)
      if (!is.null(g)) out[[length(out) + 1]] <- g
    }
  }
  dplyr::bind_rows(out)
}

# ---- Augmented SC (ridge ASCM, augsynth): adaptive window = years with treated obs and >= min_donors donors ----
gap_augsynth <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                         pre_len = 10L, min_pre = 5L, min_donors = 8L, progfunc = "ridge", seed = SEED_GAP) {
  set.seed(seed)
  conflict_iso3 <- conflict_df$ISO3
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    pre_window  <- (ci$conflict_start - pre_len):(ci$conflict_start - 1L)
    post_window <- ci$conflict_start:ci$conflict_end
    for (vacc in vaccines) {
      message("  augSC: ", ci$country, " | ", vacc)
      vacc_obs <- coverage_df %>% dplyr::filter(vaccine == vacc, !is.na(coverage))
      treated_years <- sort(vacc_obs$year[vacc_obs$country == ci$country])
      donors_by_year <- vacc_obs %>% dplyr::filter(!ISO3 %in% conflict_iso3) %>% dplyr::count(year, name = "n_donor")
      good_years <- donors_by_year$year[donors_by_year$n_donor >= min_donors]
      pre_obs  <- sort(intersect(intersect(pre_window,  treated_years), good_years))
      post_obs <- sort(intersect(intersect(post_window, treated_years), good_years))
      if (length(pre_obs) > pre_len) pre_obs <- utils::tail(pre_obs, pre_len)
      if (length(pre_obs) < min_pre || length(post_obs) < 1) {
        message("    skipping (realised window: ", length(pre_obs), " pre / ", length(post_obs), " post yrs)"); next }
      analysis_years <- c(pre_obs, post_obs)
      donors <- vacc_obs %>% dplyr::filter(!ISO3 %in% conflict_iso3, year %in% analysis_years) %>%
        dplyr::count(country, name = "n") %>% dplyr::filter(n == length(analysis_years)) %>% dplyr::pull(country)
      if (!identical(DONOR_SCOPE, "global")) {
        .di <- vacc_obs %>% dplyr::filter(country %in% donors) %>% dplyr::distinct(country, ISO3)
        donors <- .di$country[.di$ISO3 %in% .apply_donor_scope(.di$ISO3, ci$ISO3)]
      }
      if (length(donors) < min_donors) { message("    skipping (balanced donors: ", length(donors), ")"); next }
      panel <- vacc_obs %>%
        dplyr::filter(country %in% c(ci$country, donors), year %in% analysis_years) %>%
        dplyr::select(country, year, coverage) %>%
        dplyr::mutate(trt = as.integer(country == ci$country & year >= ci$conflict_start))
      asyn <- tryCatch(augsynth::augsynth(coverage ~ trt, unit = country, time = year, data = panel,
                                          progfunc = progfunc, scm = TRUE),
                       error = function(e) { message("    augsynth error: ", conditionMessage(e)); NULL })
      if (is.null(asyn)) next
      # jackknife+ CI; fall back to the default summary if it errors (as in the working pipeline)
      sm <- tryCatch(summary(asyn, inf_type = "jackknife+"), error = function(e) NULL)
      if (is.null(sm) || is.null(sm$att)) {
        sm <- tryCatch(summary(asyn), error = function(e) NULL)
        if (is.null(sm) || is.null(sm$att)) next
      }
      att <- sm$att
      tcol  <- intersect(c("Time", "time"), names(att))[1]
      ecol  <- intersect(c("Estimate", "estimate"), names(att))[1]
      locol <- intersect(c("lower_bound", "lower", "Lower"), names(att))[1]
      hicol <- intersect(c("upper_bound", "upper", "Upper"), names(att))[1]
      if (any(is.na(c(tcol, ecol)))) next
      eff <- tibble::tibble(year = att[[tcol]], est = att[[ecol]],
                            lo = if (!is.na(locol)) att[[locol]] else NA_real_,
                            hi = if (!is.na(hicol)) att[[hicol]] else NA_real_) %>%
        dplyr::filter(year %in% post_window)
      if (nrow(eff) == 0) next
      obs <- coverage_df %>% dplyr::filter(country == ci$country, vaccine == vacc, year %in% post_window) %>%
        dplyr::select(year, coverage_obs = coverage)
      out[[length(out) + 1]] <- eff %>%
        dplyr::left_join(obs, by = "year") %>%
        dplyr::mutate(country = ci$country, ISO3 = ci$ISO3, vaccine = vacc, gap = -est,
                      coverage_cf = coverage_obs + (-est),
                      gap_sd = ifelse(is.finite(hi) & is.finite(lo), (hi - lo) / 3.92, NA_real_),
                      method = METHOD_AUGSC) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs, coverage_cf, gap, gap_sd, method) %>%
        flag_sc_horizon(ci$conflict_start)
    }
  }
  dplyr::bind_rows(out)
}

# ---- Callaway-Sant'Anna staggered DiD (did): never-treated controls, outcome regression,
#      universal base period, multiplier bootstrap; post-end treated years censored ----
gap_did_cs <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                       donor_pre_len = 8L, min_year_frac = 0.8, seed = SEED_GAP) {
  set.seed(seed)
  conflict_iso3 <- conflict_df$ISO3
  win_years <- (min(conflict_df$conflict_start) - donor_pre_len):max(conflict_df$conflict_end)
  n_win <- length(win_years)
  out_gaps <- list(); es_list <- list(); overall_list <- list(); attgt_list <- list()
  for (vacc in vaccines) {
    message("  DiD(CS): ", vacc)
    treated <- coverage_df %>%
      dplyr::filter(vaccine == vacc, ISO3 %in% conflict_iso3, year %in% win_years) %>%
      dplyr::select(country, ISO3, year, coverage) %>%
      dplyr::left_join(conflict_df %>% dplyr::select(country, conflict_start, conflict_end), by = "country") %>%
      dplyr::filter(year <= conflict_end) %>%
      dplyr::rename(first_treat = conflict_start) %>% dplyr::select(-conflict_end)
    donor_ok <- coverage_df %>%
      dplyr::filter(vaccine == vacc, !ISO3 %in% conflict_iso3, year %in% win_years, !is.na(coverage)) %>%
      dplyr::group_by(country, ISO3) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n >= ceiling(min_year_frac * n_win)) %>% dplyr::pull(country)
    if (!identical(DONOR_SCOPE, "global")) {
      .di <- coverage_df %>% dplyr::filter(country %in% donor_ok) %>% dplyr::distinct(country, ISO3)
      donor_ok <- .di$country[.di$ISO3 %in% did_scope_union(.di$ISO3, conflict_df$ISO3)]
    }
    controls <- coverage_df %>%
      dplyr::filter(vaccine == vacc, country %in% donor_ok, year %in% win_years) %>%
      dplyr::select(country, ISO3, year, coverage) %>% dplyr::mutate(first_treat = 0L)
    panel <- dplyr::bind_rows(treated, controls) %>%
      dplyr::mutate(id = as.integer(factor(country))) %>% dplyr::filter(!is.na(coverage))
    if (dplyr::n_distinct(panel$id[panel$first_treat == 0]) < 5) {
      message("    skipping ", vacc, " (too few never-treated controls)"); next }
    att <- did::att_gt(yname = "coverage", tname = "year", idname = "id", gname = "first_treat",
                       xformla = ~1, data = panel, control_group = "nevertreated", est_method = "reg",
                       allow_unbalanced_panel = TRUE, base_period = "universal",
                       bstrap = TRUE, biters = 1000, print_details = FALSE)
    attgt_list[[vacc]] <- att
    .grp <- att$group; .yr <- att$t; .att <- att$att; .se <- att$se   # locals: a column named att would mask the object
    gt <- tibble::tibble(group = .grp, year = .yr, att = .att, se = .se)
    for (i in seq_len(nrow(conflict_df))) {
      ci <- conflict_df[i, ]
      yrs <- ci$conflict_start:ci$conflict_end
      sub <- gt %>% dplyr::filter(group == ci$conflict_start, year %in% yrs)
      if (nrow(sub) == 0) next
      obs <- coverage_df %>% dplyr::filter(country == ci$country, vaccine == vacc, year %in% yrs) %>%
        dplyr::select(year, coverage_obs = coverage)
      out_gaps[[length(out_gaps) + 1]] <- sub %>%
        dplyr::left_join(obs, by = "year") %>%
        dplyr::mutate(country = ci$country, ISO3 = ci$ISO3, vaccine = vacc, gap = -att,
                      coverage_cf = coverage_obs + (-att),
                      gap_sd = ifelse(is.finite(se) & se > 0, se, NA_real_), method = METHOD_DID) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs, coverage_cf, gap, gap_sd, method)
    }
    dyn <- did::aggte(att, type = "dynamic", na.rm = TRUE)
    es_list[[length(es_list) + 1]] <- tibble::tibble(vaccine = vacc, event_time = dyn$egt,
                                                     att = dyn$att.egt, se = dyn$se.egt)
    ov <- did::aggte(att, type = "group", na.rm = TRUE)
    overall_list[[length(overall_list) + 1]] <- tibble::tibble(vaccine = vacc, overall_att = ov$overall.att,
                                                               overall_se = ov$overall.se)
  }
  res <- dplyr::bind_rows(out_gaps)
  attr(res, "event_study") <- dplyr::bind_rows(es_list)
  attr(res, "overall")     <- dplyr::bind_rows(overall_list)
  attr(res, "cs_attgt")    <- attgt_list
  res
}

# ---- ITS: logit-linear pre-trend (<= 15 pre-years), OLS (<12) or GLS AR(1)-REML, mvt coefficient
#      draws + residual noise, per-draw back-transform; years >= ITS_MAX_HORIZON dropped ----
gap_its <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                    n_pre_min = 6, n_boot = 1000, pre_lookback = 15L,
                    reliable_horizon = 8L, max_horizon = ITS_MAX_HORIZON, seed = SEED_GAP) {
  set.seed(seed)
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (vacc in vaccines) {
      message("  ITS: ", ci$country, " | ", vacc)
      ctry <- coverage_df %>% dplyr::filter(country == ci$country, vaccine == vacc) %>% dplyr::arrange(year)
      pre  <- ctry %>% dplyr::filter(year < ci$conflict_start, year >= ci$conflict_start - pre_lookback)
      post <- ctry %>% dplyr::filter(year >= ci$conflict_start, year <= ci$conflict_end)
      if (nrow(pre) < n_pre_min) { message("    skipping (only ", nrow(pre), " pre years)"); next }
      df_fit <- pre %>% dplyr::mutate(p = pmin(pmax(coverage, 0.5), 99.5) / 100,
                                      logit_p = stats::qlogis(p), time = year - ci$conflict_start)
      fit <- stats::lm(logit_p ~ time, data = df_fit)
      if (is.na(stats::coef(fit)[2])) { message("    fit failed"); next }
      coefs <- stats::coef(fit); n_pre <- nrow(pre)
      method_used <- if (n_pre < ITS_GLS_MIN_N) "OLS" else "GLS"
      V <- stats::vcov(fit); sigma_resid <- summary(fit)$sigma; resid_df <- fit$df.residual
      if (method_used == "GLS") {
        # GLS can fail to converge on a short/flat series -> OLS (as in the working pipeline)
        gls_fit <- tryCatch(nlme::gls(logit_p ~ time, data = df_fit,
                                      correlation = nlme::corAR1(form = ~ time), method = "REML"),
                            error = function(e) NULL)
        if (!is.null(gls_fit)) {
          coefs <- stats::coef(gls_fit); V <- as.matrix(stats::vcov(gls_fit))
          sigma_resid <- gls_fit$sigma; resid_df <- n_pre - 2L
        } else method_used <- "OLS(GLS-failed)"
      }
      beta_draws <- mvt_draws(n_boot, mu = coefs, Sigma = V, df = max(resid_df, 1))
      post <- post[(post$year - ci$conflict_start) < max_horizon, , drop = FALSE]
      post_years <- post$year
      if (length(post_years) == 0) { message("    skipping (no post years within max_horizon)"); next }
      cf_mean <- cf_sd <- cf_satur <- numeric(length(post_years))
      for (k in seq_along(post_years)) {
        lin <- beta_draws[, 1] + beta_draws[, 2] * (post_years[k] - ci$conflict_start)
        p_draws <- stats::plogis(lin + stats::rnorm(n_boot, 0, sigma_resid))
        cf_satur[k] <- mean(p_draws > 0.98)
        cf_mean[k]  <- mean(p_draws) * 100
        cf_sd[k]    <- stats::sd(p_draws) * 100
      }
      out[[length(out) + 1]] <- tibble::tibble(
        country = ci$country, ISO3 = ci$ISO3, vaccine = vacc, year = post_years,
        coverage_obs = post$coverage, coverage_cf = cf_mean, gap = cf_mean - post$coverage,
        gap_sd = cf_sd, cf_saturated = cf_satur, slope_logit_per_yr = unname(coefs[2]),
        horizon_unreliable = (post_years - ci$conflict_start) >= reliable_horizon,
        resid_df = resid_df, n_pre = n_pre, short_series = n_pre < ITS_GLS_MIN_N,
        vcov_type = method_used, method = METHOD_ITS)
    }
  }
  dplyr::bind_rows(out)
}

# ---- driver: all estimators (optionally forked in parallel; results identical) ----
# cores for forked estimation: set in 00_run.R (CORES / CVD_CORES); 1 = sequential (default, all OSes)
N_CORES <- if (.Platform$OS.type == "windows" || is.na(parallel::detectCores())) 1L else
  as.integer(Sys.getenv("CVD_CORES", unset = "1"))

# one thunk per estimator, in the fixed stacking order (baselines 1/3/5, ITS, SC, SC+cov, ASCM, DiD)
.gap_jobs <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                      methods = ALL_METHODS, seed = SEED_GAP) {
  force(coverage_df); force(conflict_df); force(vaccines); force(methods); force(seed)  # thunks must not see later loop values
  jobs <- list()
  for (w in BASELINE_WINDOWS) {
    lab <- baseline_label(w)
    if (lab %in% methods) jobs[[lab]] <- local({ w <- w; lab <- lab; function() {
      gw <- .gap_baseline_w(coverage_df, conflict_df, vaccines, w)
      if (nrow(gw) > 0) gw %>% dplyr::mutate(method = lab) else NULL } })
  }
  if (METHOD_ITS    %in% methods) jobs[[METHOD_ITS]]    <- function() gap_its(coverage_df, conflict_df, vaccines, seed = seed)
  if (METHOD_SC     %in% methods) jobs[[METHOD_SC]]     <- function() gap_synthetic(coverage_df, conflict_df, vaccines, seed = seed)
  if (METHOD_SC_COV %in% methods) jobs[[METHOD_SC_COV]] <- function() gap_synthetic_cov(coverage_df, conflict_df, vaccines, seed = seed)
  if (METHOD_AUGSC  %in% methods) jobs[[METHOD_AUGSC]]  <- function() gap_augsynth(coverage_df, conflict_df, vaccines, seed = seed)
  if (METHOD_DID    %in% methods) jobs[[METHOD_DID]]    <- function() gap_did_cs(coverage_df, conflict_df, vaccines, seed = seed)
  jobs
}

# run thunks (forked if n_cores > 1); stop on any failure
.run_jobs <- function(jobs, n_cores = N_CORES) {
  res <- if (n_cores > 1 && length(jobs) > 1)
    parallel::mclapply(jobs, function(f) f(), mc.cores = min(n_cores, length(jobs)), mc.preschedule = FALSE)
  else lapply(jobs, function(f) f())
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("estimation failed for: ", paste(names(res)[bad], collapse = ", "), "\n",
                     paste(unlist(res[bad]), collapse = "\n"))
  res
}

# stack pieces -> gap_df (sd floored at 1 pp; proportions); DiD event study kept as attributes
.combine_gaps <- function(pieces) {
  did_res <- pieces[[METHOD_DID]]
  stacked <- dplyr::bind_rows(pieces)
  if (nrow(stacked) == 0) stop("no methods produced gaps.")
  res <- stacked %>%
    dplyr::mutate(gap_sd = pmax(gap_sd, 1.0), gap_prop = gap / 100, gap_sd_prop = gap_sd / 100,
                  coverage_obs_prop = coverage_obs / 100, coverage_cf_prop = coverage_cf / 100)
  attr(res, "event_study") <- if (!is.null(did_res)) attr(did_res, "event_study")
  attr(res, "overall_did") <- if (!is.null(did_res)) attr(did_res, "overall")
  attr(res, "cs_attgt")    <- if (!is.null(did_res)) attr(did_res, "cs_attgt")
  res
}

estimate_all_gaps_plus <- function(coverage_df, conflict_df, vaccines = c("BCG", "MCV1", "DTP3"),
                                   methods = ALL_METHODS, seed = SEED_GAP, n_cores = N_CORES)
  .combine_gaps(.run_jobs(.gap_jobs(coverage_df, conflict_df, vaccines, methods, seed), n_cores))

# ---- complete-case total YLL per estimator (cells present under every method), draw level ----
aggregate_global_by_method_draws <- function(yll_draws, methods = NULL, catchup_focus = 0,
                                             framing_focus = "campaign_topup") {
  d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  if (!is.null(methods)) d <- d %>% dplyr::filter(method %in% methods)
  d <- d[vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1)), , drop = FALSE]
  if (nrow(d) == 0) return(tibble::tibble())
  n_methods <- dplyr::n_distinct(d$method)
  complete <- d %>% dplyr::distinct(disease, country, method) %>%
    dplyr::count(disease, country, name = "nm") %>% dplyr::filter(nm == n_methods) %>%
    dplyr::select(disease, country)
  d %>% dplyr::semi_join(complete, by = c("disease", "country")) %>%
    dplyr::group_by(method) %>%
    dplyr::group_modify(~ {
      s_un <- colSums(do.call(rbind, .x$d_undisc)); s_di <- colSums(do.call(rbind, .x$d_disc))
      tibble::tibble(n_cells = nrow(.x),
                     yll_undisc_mean = mean(s_un),
                     yll_undisc_lo = stats::quantile(s_un, 0.025, names = FALSE),
                     yll_undisc_hi = stats::quantile(s_un, 0.975, names = FALSE),
                     yll_disc_mean = mean(s_di),
                     yll_disc_lo = stats::quantile(s_di, 0.025, names = FALSE),
                     yll_disc_hi = stats::quantile(s_di, 0.975, names = FALSE))
    }) %>% dplyr::ungroup()
}
