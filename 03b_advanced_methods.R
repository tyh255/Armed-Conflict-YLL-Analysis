stopifnot(exists("build_donor_pool"), exists("DEFAULT_SEED"))

# Method labels (single source of truth; 05 palette keys off these).
METHOD_BASELINE <- "Baseline"
METHOD_ITS      <- "ITS"
METHOD_SC       <- "Synthetic Control"
METHOD_SC_COV   <- "SC + covariates"
METHOD_AUGSC    <- "Augmented SC"
METHOD_DID      <- "DiD (CS)"

BASELINE_WINDOWS         <- as.integer(strsplit(
  Sys.getenv("CVD_BASELINE_WINDOWS", unset = "1,3,5"), ",")[[1]])
BASELINE_HEADLINE_WINDOW <- 3L
baseline_label         <- function(w) paste0("Baseline (", w, "yr)")
METHOD_BASELINE_LABELS <- baseline_label(BASELINE_WINDOWS)

ADVANCED_METHODS <- c(METHOD_SC_COV, METHOD_AUGSC, METHOD_DID)
ALL_METHODS      <- c(METHOD_BASELINE_LABELS, METHOD_ITS, METHOD_SC,
                      METHOD_SC_COV, METHOD_AUGSC, METHOD_DID)

# ---------------------------------------------------------------------------
# METHOD 4: Synthetic Control with economic/development covariates
# ---------------------------------------------------------------------------
gap_synthetic_cov <- function(coverage_df, conflict_df,
                              vaccines = c("BCG", "MCV1", "DTP3"),
                              econ = if (exists("econ_panel")) econ_panel else NULL,
                              predictors = if (exists("ECON_PREDICTORS")) ECON_PREDICTORS else character(0),
                              pre_len = 7L, min_donors = 10L,
                              rmspe_ratio_cap = 5,
                              use_lag_predictors = TRUE,
                              seed = DEFAULT_SEED) {
  if (is.null(econ) || length(predictors) == 0) {
    message("  SC+cov: econ_panel / ECON_PREDICTORS unavailable -> skipping ",
            "(source 02b_econ_covariates.R with the GDP/HDI files present).")
    return(tibble::tibble())
  }
  if (!requireNamespace("tidysynth", quietly = TRUE)) {
    message("  SC+cov: package 'tidysynth' not installed -> skipping."); return(tibble::tibble())
  }
  set.seed(seed)
  conflict_iso3 <- conflict_df$ISO3
  econ_keep <- c("country", "year", intersect(predictors, names(econ)))
  econ_use  <- econ %>% dplyr::select(dplyr::all_of(econ_keep)) %>% dplyr::distinct()
  preds     <- intersect(predictors, names(econ_use))
  out <- list()
  
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    pre_window  <- (ci$conflict_start - pre_len):(ci$conflict_start - 1L)
    post_window <- ci$conflict_start:ci$conflict_end
    
    for (vacc in vaccines) {
      message("  SC+cov: ", ci$country, " | ", vacc)
      donors <- build_donor_pool(coverage_df, vacc, conflict_iso3, pre_window,
                                 treated_iso3 = ci$ISO3, scope = DONOR_SCOPE)
      if (length(donors) < min_donors) {
        message("    skipping (too few donors: ", length(donors), ")"); next
      }
      
      panel <- coverage_df %>%
        dplyr::filter(vaccine == vacc, country %in% c(ci$country, donors),
                      year %in% c(pre_window, post_window)) %>%
        dplyr::select(country, year, coverage) %>%
        dplyr::left_join(econ_use, by = c("country", "year")) %>%
        tidyr::complete(country, year) %>%
        dplyr::arrange(country, year)
      
      treated_pre_ok <- panel %>%
        dplyr::filter(country == ci$country, year %in% pre_window) %>%
        dplyr::summarise(ok = all(!is.na(coverage))) %>% dplyr::pull(ok)
      if (!isTRUE(treated_pre_ok)) { message("    skipping (treated NA in pre-window)"); next }
      
    
      out_ok_ctry <- panel %>%
        dplyr::filter(country != ci$country) %>%
        dplyr::group_by(country) %>%
        dplyr::summarise(out_ok = all(!is.na(coverage)), .groups = "drop")
      cov_ok_ctry <- panel %>%
        dplyr::filter(country != ci$country, year %in% pre_window) %>%
        dplyr::group_by(country) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(preds),
                                       ~ any(!is.na(.x)), .names = "ok_{.col}"),
                         .groups = "drop")
      ok_cols <- paste0("ok_", preds)
      cov_ok_ctry$cov_ok <-
        rowSums(as.matrix(cov_ok_ctry[ok_cols])) == length(preds)
      donor_ok <- out_ok_ctry %>%
        dplyr::left_join(cov_ok_ctry[c("country", "cov_ok")], by = "country") %>%
        dplyr::mutate(cov_ok = !is.na(cov_ok) & cov_ok) %>%
        dplyr::filter(out_ok, cov_ok) %>% dplyr::pull(country)
      if (length(donor_ok) < min_donors) {
        message("    skipping (insufficient covariate-complete donors: ",
                length(donor_ok), ")"); next
      }
      panel <- panel %>% dplyr::filter(country %in% c(ci$country, donor_ok))
      
      y1 <- pre_window[1]; y2 <- pre_window[ceiling(length(pre_window) / 2)]
      y3 <- pre_window[length(pre_window)]
      lag_years <- unique(c(y1, y2, y3))
      
      sc_obj <- tryCatch({
        sc <- panel %>%
          tidysynth::synthetic_control(
            outcome = coverage, unit = country, time = year,
            i_unit = ci$country, i_time = ci$conflict_start,
            generate_placebos = TRUE) %>%
          tidysynth::generate_predictor(time_window = pre_window,
                                        cov_mean = mean(coverage, na.rm = TRUE))
        if (use_lag_predictors)
          for (ly in lag_years)
            sc <- tidysynth::generate_predictor(
              sc, time_window = ly,
              !!paste0("lag_", ly) := mean(coverage, na.rm = TRUE))
        # Covariate predictors (explicit per-covariate calls; literal column
        # names so tidysynth's NSE evaluates them in the panel context).
        if ("log_GDP_pc" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window,
                                              p_log_gdp = mean(log_GDP_pc, na.rm = TRUE))
        if ("sch_eys" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window,
                                              p_eys = mean(sch_eys, na.rm = TRUE))
        if ("sch_mys" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window,
                                              p_mys = mean(sch_mys, na.rm = TRUE))
        if ("hdi" %in% preds)
          sc <- tidysynth::generate_predictor(sc, time_window = pre_window,
                                              p_hdi = mean(hdi, na.rm = TRUE))
        sc %>%
          tidysynth::generate_weights(optimization_window = pre_window,
                                      margin_ipop = 0.02, sigf_ipop = 7,
                                      bound_ipop = 6) %>%
          tidysynth::generate_control()
      }, error = function(e) { message("    tidysynth error: ", conditionMessage(e)); NULL })
      if (is.null(sc_obj)) next
      
      full_sc <- tryCatch(tidysynth::grab_synthetic_control(sc_obj, placebo = TRUE),
                          error = function(e) NULL)
      if (is.null(full_sc) || nrow(full_sc) == 0) next
      if (!all(c(".id", ".placebo", "time_unit", "real_y", "synth_y") %in% names(full_sc))) next
      full_sc <- dplyr::rename(full_sc, year = time_unit)
      
      gap_tbl <- full_sc %>%
        dplyr::filter(.placebo == 0, year %in% post_window) %>%
        dplyr::transmute(country = ci$country, ISO3 = ci$ISO3, vaccine = vacc, year = year,
                         coverage_obs = real_y, coverage_cf = synth_y, gap = synth_y - real_y)
      if (nrow(gap_tbl) == 0) next
      
      rmspe_tbl <- full_sc %>%
        dplyr::filter(year %in% pre_window) %>%
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
      
      post_rmspe_tbl <- full_sc %>%
        dplyr::filter(year %in% post_window) %>%
        dplyr::group_by(.id, .placebo) %>%
        dplyr::summarise(post_rmspe = sqrt(mean((real_y - synth_y)^2, na.rm = TRUE)), .groups = "drop")
      ratio_tbl <- rmspe_tbl %>%
        dplyr::left_join(post_rmspe_tbl, by = c(".id", ".placebo")) %>%
        dplyr::mutate(rmspe_ratio = post_rmspe / pre_rmspe) %>%
        dplyr::filter(is.finite(rmspe_ratio))
      treated_ratio <- ratio_tbl %>% dplyr::filter(.placebo == 0) %>% dplyr::pull(rmspe_ratio)
      treated_ratio <- if (length(treated_ratio) == 1) treated_ratio else NA_real_
      perm_p <- if (is.na(treated_ratio)) NA_real_ else
        mean(ratio_tbl$rmspe_ratio >= treated_ratio, na.rm = TRUE)
      
      gap_tbl <- gap_tbl %>%
        dplyr::left_join(placebo_tbl, by = "year") %>%
        dplyr::mutate(gap_sd = ifelse(is.na(gap_sd) | gap_sd == 0, 5.0, gap_sd),
                      pre_rmspe = treated_rmspe, rmspe_ratio = treated_ratio,
                      perm_pvalue = perm_p, n_good_placebos = length(good_ids),
                      n_donors_total = length(donor_ok), method = METHOD_SC_COV) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs, coverage_cf, gap,
                      gap_sd, pre_rmspe, rmspe_ratio, perm_pvalue,
                      n_good_placebos, n_donors_total, method)
      gap_tbl <- flag_sc_horizon(gap_tbl, ci$conflict_start)   # F3: long-horizon flag
      out[[length(out) + 1]] <- gap_tbl
    }
  }
  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# METHOD 5: Augmented Synthetic Control (ridge), {augsynth}
# ---------------------------------------------------------------------------
gap_augsynth <- function(coverage_df, conflict_df,
                         vaccines = c("BCG", "MCV1", "DTP3"),
                         econ = if (exists("econ_panel")) econ_panel else NULL,
                         predictors = if (exists("ECON_PREDICTORS")) ECON_PREDICTORS else character(0),
                         pre_len = 10L, min_pre = 5L, min_donors = 8L,
                         progfunc = "ridge", seed = DEFAULT_SEED) {
  if (!requireNamespace("augsynth", quietly = TRUE)) {
    message("  augSC: package 'augsynth' not installed (GitHub: ebenmichael/augsynth) -> skipping.")
    return(tibble::tibble())
  }
  set.seed(seed)
  conflict_iso3 <- conflict_df$ISO3
  preds <- if (!is.null(econ)) intersect(predictors, names(econ)) else character(0)
  out <- list()
  
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    pre_window  <- (ci$conflict_start - pre_len):(ci$conflict_start - 1L)
    post_window <- ci$conflict_start:ci$conflict_end
    
    for (vacc in vaccines) {
      message("  augSC: ", ci$country, " | ", vacc)
      vacc_obs <- coverage_df %>%
        dplyr::filter(vaccine == vacc, !is.na(coverage))
      
    
      treated_years <- sort(vacc_obs$year[vacc_obs$country == ci$country])
      donors_by_year <- vacc_obs %>%
        dplyr::filter(!ISO3 %in% conflict_iso3) %>%
        dplyr::count(year, name = "n_donor")
      good_years <- donors_by_year$year[donors_by_year$n_donor >= min_donors]
      
      pre_obs  <- sort(intersect(intersect(pre_window,  treated_years), good_years))
      post_obs <- sort(intersect(intersect(post_window, treated_years), good_years))
      if (length(pre_obs) > pre_len) pre_obs <- utils::tail(pre_obs, pre_len)  # most recent
      if (length(pre_obs) < min_pre || length(post_obs) < 1) {
        message("    skipping (realised window: ", length(pre_obs), " pre / ",
                length(post_obs), " post yrs; need >=", min_pre, "/1)"); next }
      analysis_years <- c(pre_obs, post_obs)
      
      
      donors <- vacc_obs %>%
        dplyr::filter(!ISO3 %in% conflict_iso3, year %in% analysis_years) %>%
        dplyr::count(country, name = "n") %>%
        dplyr::filter(n == length(analysis_years)) %>% dplyr::pull(country)
      if (!identical(DONOR_SCOPE, "global")) {            # F1: SUTVA/comparability scope
        .di <- vacc_obs %>% dplyr::filter(country %in% donors) %>%
          dplyr::distinct(country, ISO3)
        donors <- .di$country[.di$ISO3 %in%
                                .apply_donor_scope(.di$ISO3, ci$ISO3, scope = DONOR_SCOPE)]
      }
      if (length(donors) < min_donors) {
        message("    skipping (balanced donors: ", length(donors), ")"); next }
      
      panel <- vacc_obs %>%
        dplyr::filter(country %in% c(ci$country, donors), year %in% analysis_years) %>%
        dplyr::select(country, year, coverage) %>%
        dplyr::mutate(trt = as.integer(country == ci$country & year >= ci$conflict_start))
      fml <- coverage ~ trt
      asyn <- tryCatch(
        augsynth::augsynth(fml, unit = country, time = year, data = panel,
                           progfunc = progfunc, scm = TRUE),
        error = function(e) { message("    augsynth error: ", conditionMessage(e)); NULL })
      if (is.null(asyn)) next
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
      
      obs <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc, year %in% post_window) %>%
        dplyr::select(year, coverage_obs = coverage)
      gap_tbl <- eff %>%
        dplyr::left_join(obs, by = "year") %>%
        dplyr::mutate(
          country = ci$country, ISO3 = ci$ISO3, vaccine = vacc,
          gap = -est,                          # effect on coverage is negative
          coverage_cf = coverage_obs + (-est), # cf = obs + gap
          gap_sd = ifelse(is.finite(hi) & is.finite(lo), (hi - lo) / 3.92, NA_real_),
          method = METHOD_AUGSC) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs, coverage_cf,
                      gap, gap_sd, method)
      gap_tbl <- flag_sc_horizon(gap_tbl, ci$conflict_start)   # F3: long-horizon flag
      out[[length(out) + 1]] <- gap_tbl
    }
  }
  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# METHOD 6: Callaway & Sant'Anna (2021) staggered DiD, {did}
# ---------------------------------------------------------------------------
gap_did_cs <- function(coverage_df, conflict_df,
                       vaccines = c("BCG", "MCV1", "DTP3"),
                       econ = if (exists("econ_panel")) econ_panel else NULL,
                       predictors = if (exists("ECON_PREDICTORS")) ECON_PREDICTORS else character(0),
                       donor_pre_len = 8L, min_year_frac = 0.8,
                       control_group = "nevertreated",
                       use_covariates = FALSE, seed = DEFAULT_SEED) {
  if (!requireNamespace("did", quietly = TRUE)) {
    message("  DiD(CS): package 'did' not installed -> skipping."); return(tibble::tibble())
  }
  set.seed(seed)
  conflict_iso3 <- conflict_df$ISO3
  preds <- if (use_covariates && !is.null(econ)) intersect(predictors, names(econ)) else character(0)
  
  win_start <- min(conflict_df$conflict_start) - donor_pre_len
  win_end   <- max(conflict_df$conflict_end)
  win_years <- win_start:win_end
  n_win     <- length(win_years)
  
  out_gaps <- list(); es_list <- list(); overall_list <- list()
  # Retain the full att_gt objects per vaccine so the HonestDiD module (07) can
  # re-aggregate the dynamic event study with min_e/max_e and read the influence
  # function. Stored as an attribute; adds a few MB to the gap_df RDS.
  attgt_list <- list()
  
  for (vacc in vaccines) {
    message("  DiD(CS): ", vacc)
    # Treated conflict countries with their onset.
    treated <- coverage_df %>%
      dplyr::filter(vaccine == vacc, ISO3 %in% conflict_iso3, year %in% win_years) %>%
      dplyr::select(country, ISO3, year, coverage) %>%
      dplyr::left_join(conflict_df %>%
                         dplyr::select(country, conflict_start, conflict_end),
                       by = "country") %>%
      
      dplyr::filter(year <= conflict_end) %>%
      dplyr::rename(first_treat = conflict_start) %>%
      dplyr::select(-conflict_end)
    # Never-treated donor controls: non-conflict, present in >= min_year_frac of
    # the window with non-missing coverage.
    donor_ok <- coverage_df %>%
      dplyr::filter(vaccine == vacc, !ISO3 %in% conflict_iso3, year %in% win_years,
                    !is.na(coverage)) %>%
      dplyr::group_by(country, ISO3) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n >= ceiling(min_year_frac * n_win)) %>%
      dplyr::pull(country)
    if (!identical(DONOR_SCOPE, "global")) {             # F1: pooled SUTVA/comparability scope
      .di <- coverage_df %>% dplyr::filter(country %in% donor_ok) %>%
        dplyr::distinct(country, ISO3)
      donor_ok <- .di$country[.di$ISO3 %in%
                                did_scope_union(.di$ISO3, conflict_df$ISO3, scope = DONOR_SCOPE)]
    }
    controls <- coverage_df %>%
      dplyr::filter(vaccine == vacc, country %in% donor_ok, year %in% win_years) %>%
      dplyr::select(country, ISO3, year, coverage) %>%
      dplyr::mutate(first_treat = 0L)
    
    panel <- dplyr::bind_rows(treated, controls)
    if (length(preds) > 0)
      panel <- panel %>% dplyr::left_join(
        econ %>% dplyr::select(dplyr::all_of(c("country", "year", preds))),
        by = c("country", "year"))
    # did needs an integer unit id and (for DR covariates) complete covariate
    # rows; drop covariate-missing rows only when covariates are requested.
    panel <- panel %>%
      dplyr::mutate(id = as.integer(factor(country))) %>%
      dplyr::filter(!is.na(coverage))
    use_cov_here <- length(preds) > 0
    if (use_cov_here) {
      before <- nrow(panel)
      panel  <- panel %>% dplyr::filter(dplyr::if_all(dplyr::all_of(preds), ~ !is.na(.)))
      if (nrow(panel) < 0.6 * before) {
        message("    covariate rows too sparse (kept ", nrow(panel), "/", before,
                "); falling back to no-covariate DiD for ", vacc, ".")
        use_cov_here <- FALSE
        panel <- dplyr::bind_rows(treated, controls) %>%
          dplyr::mutate(id = as.integer(factor(country))) %>% dplyr::filter(!is.na(coverage))
      }
    }
    if (dplyr::n_distinct(panel$id[panel$first_treat == 0]) < 5) {
      message("    skipping ", vacc, " (too few never-treated controls)"); next }
    
    xformla   <- if (use_cov_here) stats::as.formula(paste("~", paste(preds, collapse = " + "))) else ~1
    
    est_method <- "reg"
    run_attgt <- function(xf, em) tryCatch(
      did::att_gt(yname = "coverage", tname = "year", idname = "id",
                  gname = "first_treat", xformla = xf, data = panel,
                  control_group = control_group, est_method = em,
                  allow_unbalanced_panel = TRUE, base_period = "universal",
                  bstrap = TRUE, biters = 1000, print_details = FALSE),
      error = function(e) { message("    att_gt error: ", conditionMessage(e)); NULL })
    att <- run_attgt(xformla, est_method)
   
    if (is.null(att) && use_cov_here) {
      message("    retrying ", vacc, " without covariates (unconditional CS).")
      att <- run_attgt(~1, "reg")
    }
    if (is.null(att)) next
    attgt_list[[vacc]] <- att
    
    
    .grp <- att$group; .yr <- att$t; .att <- att$att; .se <- att$se
    gt <- tibble::tibble(group = .grp, year = .yr, att = .att, se = .se)
    # Map ATT(g,t) -> per conflict-country-year gap within its window.
    for (i in seq_len(nrow(conflict_df))) {
      ci <- conflict_df[i, ]
      yrs <- ci$conflict_start:ci$conflict_end
      sub <- gt %>% dplyr::filter(group == ci$conflict_start, year %in% yrs)
      if (nrow(sub) == 0) next
      obs <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc, year %in% yrs) %>%
        dplyr::select(year, coverage_obs = coverage)
      g <- sub %>%
        dplyr::left_join(obs, by = "year") %>%
        dplyr::mutate(country = ci$country, ISO3 = ci$ISO3, vaccine = vacc,
                      gap = -att, coverage_cf = coverage_obs + (-att),
                      gap_sd = ifelse(is.finite(se) & se > 0, se, NA_real_),
                      method = METHOD_DID) %>%
        dplyr::select(country, ISO3, vaccine, year, coverage_obs, coverage_cf,
                      gap, gap_sd, method)
      out_gaps[[length(out_gaps) + 1]] <- g
    }
    # Event study (dynamic) + overall ATT for the comparison figure.
    dyn <- tryCatch(did::aggte(att, type = "dynamic", na.rm = TRUE),
                    error = function(e) NULL)
    if (!is.null(dyn))
      es_list[[length(es_list) + 1]] <- tibble::tibble(
        vaccine = vacc, event_time = dyn$egt, att = dyn$att.egt, se = dyn$se.egt)
    ov <- tryCatch(did::aggte(att, type = "group", na.rm = TRUE), error = function(e) NULL)
    if (!is.null(ov))
      overall_list[[length(overall_list) + 1]] <- tibble::tibble(
        vaccine = vacc, overall_att = ov$overall.att, overall_se = ov$overall.se)
  }
  
  res <- dplyr::bind_rows(out_gaps)
  if (nrow(res) == 0) {
    message("  DiD(CS): produced 0 gap rows. Likely causes: att_gt errored for ",
            "every vaccine (see 'att_gt error:' lines above), too few never-treated ",
            "controls, or no ATT(g,t) matched a conflict cohort/window. DiD will be ",
            "absent from gap_df, YLL, and figures.")
  } else {
    message("  DiD(CS): produced ", nrow(res), " gap rows across ",
            dplyr::n_distinct(res$country), " countries.")
  }
  attr(res, "event_study") <- dplyr::bind_rows(es_list)
  attr(res, "overall")     <- dplyr::bind_rows(overall_list)
  attr(res, "cs_attgt")    <- attgt_list
  res
}

# ---------------------------------------------------------------------------
# Sun & Abraham (2021) event study companion, {fixest}
# ---------------------------------------------------------------------------
did_event_study_sa <- function(coverage_df, conflict_df, vaccine = "MCV1",
                               donor_pre_len = 8L, min_year_frac = 0.8) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    message("  SA event study: package 'fixest' not installed -> skipping.")
    return(tibble::tibble())
  }
  conflict_iso3 <- conflict_df$ISO3
  win_years <- (min(conflict_df$conflict_start) - donor_pre_len):max(conflict_df$conflict_end)
  n_win <- length(win_years)
  treated <- coverage_df %>%
    dplyr::filter(vaccine == !!vaccine, ISO3 %in% conflict_iso3, year %in% win_years) %>%
    dplyr::left_join(conflict_df %>% dplyr::select(country, conflict_start), by = "country") %>%
    dplyr::transmute(country, year, coverage, first_treat = conflict_start)
  donor_ok <- coverage_df %>%
    dplyr::filter(vaccine == !!vaccine, !ISO3 %in% conflict_iso3, year %in% win_years, !is.na(coverage)) %>%
    dplyr::count(country) %>% dplyr::filter(n >= ceiling(min_year_frac * n_win)) %>% dplyr::pull(country)
  controls <- coverage_df %>%
    dplyr::filter(vaccine == !!vaccine, country %in% donor_ok, year %in% win_years) %>%
    dplyr::transmute(country, year, coverage, first_treat = 10000L)  # fixest: never-treated as +Inf cohort
  panel <- dplyr::bind_rows(treated, controls) %>% dplyr::filter(!is.na(coverage))
  m <- tryCatch(
    fixest::feols(coverage ~ fixest::sunab(first_treat, year) | country + year, data = panel),
    error = function(e) { message("    sunab error: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(tibble::tibble())
  ct <- as.data.frame(stats::coeftable(m))
  ct$term <- rownames(ct)
  ct %>%
    dplyr::filter(grepl("year::", term)) %>%
    dplyr::transmute(event_time = as.integer(sub(".*year::(-?\\d+).*", "\\1", term)),
                     att = .data[["Estimate"]], se = .data[["Std. Error"]],
                     vaccine = vaccine)
}

# ---------------------------------------------------------------------------
# Baseline-window sensitivity
# ---------------------------------------------------------------------------
gap_baseline_windows <- function(coverage_df, conflict_df,
                                 vaccines = c("BCG", "MCV1", "DTP3"),
                                 windows = c(3L, 5L)) {
  res <- lapply(windows, function(w) {
    cf <- conflict_df %>%
      dplyr::mutate(preconflict_window = purrr::map(preconflict_year,
                                                    ~ (.x - (w - 1L)):.x))
    g <- .gap_baseline_w(coverage_df, cf, vaccines, w)
    if (nrow(g) == 0) return(NULL)
    g %>% dplyr::mutate(method = paste0("Baseline (", w, "yr)"))
  })
  dplyr::bind_rows(res)
}

# Window-parameterised baseline (the production gap_baseline() is fixed at 3yr).
.gap_baseline_w <- function(coverage_df, conflict_df, vaccines, w) {
  out <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (vacc in vaccines) {
      yrs <- (ci$preconflict_year - (w - 1L)):ci$preconflict_year
      vals <- coverage_df %>%
        dplyr::filter(country == ci$country, vaccine == vacc, year %in% yrs) %>%
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

# ---------------------------------------------------------------------------
# Extended driver: original 3 methods + the 3 advanced ones, finalised.
# ---------------------------------------------------------------------------
.finalize_gaps <- function(stacked) {
  stacked %>%
    dplyr::mutate(gap_sd = pmax(gap_sd, 1.0),
                  gap_prop = gap / 100, gap_sd_prop = gap_sd / 100,
                  coverage_obs_prop = coverage_obs / 100,
                  coverage_cf_prop  = coverage_cf / 100)
}

estimate_all_gaps_plus <- function(coverage_df, conflict_df,
                                   vaccines = c("BCG", "MCV1", "DTP3"),
                                   methods = ALL_METHODS, seed = DEFAULT_SEED) {
  pieces <- list(); es <- NULL; overall <- NULL
  
  bw_requested <- intersect(METHOD_BASELINE_LABELS, methods)
  if (METHOD_BASELINE %in% methods)
    bw_requested <- union(bw_requested, baseline_label(BASELINE_HEADLINE_WINDOW))
  if (length(bw_requested) > 0) {
    for (w in BASELINE_WINDOWS) {
      lab <- baseline_label(w)
      if (!lab %in% bw_requested) next
      message("== ", lab, " ==")
      gw <- .gap_baseline_w(coverage_df, conflict_df, vaccines, w)
      if (nrow(gw) > 0)
        pieces[[length(pieces) + 1]] <- gw %>% dplyr::mutate(method = lab)
    }
  }
  if (METHOD_ITS %in% methods) {
    message("== Method 2: ITS ==")
    pieces[[length(pieces) + 1]] <- gap_its(coverage_df, conflict_df, vaccines, seed = seed)
  }
  if (METHOD_SC %in% methods) {
    message("== Method 3: Synthetic Control ==")
    pieces[[length(pieces) + 1]] <- gap_synthetic(coverage_df, conflict_df, vaccines, seed = seed)
  }
  if (METHOD_SC_COV %in% methods) {
    message("== Method 4: SC + covariates ==")
    pieces[[length(pieces) + 1]] <- gap_synthetic_cov(coverage_df, conflict_df, vaccines, seed = seed)
  }
  if (METHOD_AUGSC %in% methods) {
    message("== Method 5: Augmented SC ==")
    pieces[[length(pieces) + 1]] <- gap_augsynth(coverage_df, conflict_df, vaccines, seed = seed)
  }
  did_res <- NULL
  cs_attgt <- NULL
  if (METHOD_DID %in% methods) {
    message("== Method 6: DiD (Callaway-Sant'Anna) ==")
    did_res <- gap_did_cs(coverage_df, conflict_df, vaccines, seed = seed)
    pieces[[length(pieces) + 1]] <- did_res
    es <- attr(did_res, "event_study"); overall <- attr(did_res, "overall")
    cs_attgt <- attr(did_res, "cs_attgt")
  }
  stacked <- dplyr::bind_rows(pieces)
  if (nrow(stacked) == 0) stop("estimate_all_gaps_plus: no methods produced gaps.")
  res <- .finalize_gaps(stacked)
  attr(res, "event_study") <- es
  attr(res, "overall_did") <- overall
  attr(res, "cs_attgt")    <- cs_attgt
  res
}

# ---------------------------------------------------------------------------
# Global (all-country, all-disease) YLL total per METHOD, draw-level, CRN-correct
# ---------------------------------------------------------------------------
aggregate_global_by_method_draws <- function(yll_draws, methods = NULL,
                                             catchup_focus = 0,
                                             framing_focus = "campaign_topup") {
  if (is.null(yll_draws) || nrow(yll_draws) == 0) return(tibble::tibble())
  d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  if (!is.null(methods)) d <- d %>% dplyr::filter(method %in% methods)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(tibble::tibble())
  
  n_methods <- dplyr::n_distinct(d$method)
  complete  <- d %>%
    dplyr::distinct(disease, country, method) %>%
    dplyr::count(disease, country, name = "nm") %>%
    dplyr::filter(nm == n_methods) %>% dplyr::select(disease, country)
  dropped <- d %>% dplyr::distinct(disease, country) %>%
    dplyr::anti_join(complete, by = c("disease", "country"))
  if (nrow(dropped) > 0)
    message("Figure-4 global-by-method: complete-case excludes ", nrow(dropped),
            " country-disease cells absent from >=1 of ", n_methods, " methods.")
  
  d %>%
    dplyr::semi_join(complete, by = c("disease", "country")) %>%
    dplyr::group_by(method) %>%
    dplyr::group_modify(~ {
      M_un <- do.call(rbind, .x$d_undisc); M_di <- do.call(rbind, .x$d_disc)
      s_un <- colSums(M_un); s_di <- colSums(M_di)
      tibble::tibble(
        n_cells = nrow(M_un),
        yll_undisc_mean = mean(s_un),
        yll_undisc_lo = stats::quantile(s_un, 0.025, names = FALSE),
        yll_undisc_hi = stats::quantile(s_un, 0.975, names = FALSE),
        yll_disc_mean = mean(s_di),
        yll_disc_lo = stats::quantile(s_di, 0.025, names = FALSE),
        yll_disc_hi = stats::quantile(s_di, 0.975, names = FALSE))
    }) %>% dplyr::ungroup()
}
