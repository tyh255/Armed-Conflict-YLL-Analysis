# ---- Incidence column lookup ----------------------------------------------
incidence_col <- function(disease) {
  switch(disease,
         "Tuberculosis" = "TB_Inc_Rate",
         "Measles"      = "Measles_Inc_Rate",
         "Diphtheria"   = "Diphtheria_Inc_Rate",
         "Pertussis"    = "Pertussis_Inc_Rate",
         "Tetanus"      = "Tetanus_Inc_Rate")
}

# Per-row GBD lognormal sigma column (paired with incidence_col)
incidence_sigma_col <- function(disease) {
  switch(disease,
         "Tuberculosis" = "TB_Inc_Sigma",
         "Measles"      = "Measles_Inc_Sigma",
         "Diphtheria"   = "Diphtheria_Inc_Sigma",
         "Pertussis"    = "Pertussis_Inc_Sigma",
         "Tetanus"      = "Tetanus_Inc_Sigma")
}

draw_gap_beta <- function(n, mean_val, sd_val) {
  m <- min(max(mean_val, 1e-9), 1 - 1e-9)
  max_sd <- sqrt(m * (1 - m)) * 0.999
  s <- min(sd_val, max_sd)
  if (!is.finite(s) || s <= 0) return(rep(m, n))
  var_v  <- s^2
  common <- m * (1 - m) / var_v - 1
  a <- m * common; b <- (1 - m) * common
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0)
    return(pmin(pmax(stats::rnorm(n, mean_val, sd_val), 0), 1))
  stats::rbeta(n, a, b)
}

# ============================================================================
# 04_CATCHUP_REVISED.R
# ----------------------------------------------------------------------------

CATCHUP_REACH <- list(mean = 0.66, lo = 0.35, hi = 0.90, dist = "beta",
                      source = "Portnoy 2018")


CATCHUP_IMMUNE_DISCOUNT_MEASLES <- 0.0


CATCHUP_FRAMINGS <- c("routine_recovery", "campaign_topup")

catchup_gap_factors <- function(framing, lambda, gap_mean_vec,
                                coverage_obs_vec, coverage_cf_vec,
                                reach_draws, immune_discount = 0) {
  n <- length(gap_mean_vec)
  if (lambda == 0 || n == 0) return(rep(list(1), n))
  
  if (framing == "campaign_topup") {
    
    close <- pmin(pmax(lambda * reach_draws * (1 - immune_discount), 0), 1)
    return(rep(list(pmax(1 - close, 0)), n))
  }
  
  if (framing == "routine_recovery") {
    
    if (anyNA(c(coverage_obs_vec[1], coverage_cf_vec[1]))) return(rep(list(1), n))
    cov_path <- numeric(n); cov_path[1] <- coverage_obs_vec[1]
    C_pre <- coverage_cf_vec[1]
    if (n >= 2) {
      for (k in 2:n) {
        cov_path[k] <- cov_path[k - 1] + lambda * (C_pre - cov_path[k - 1])
        if (!is.na(coverage_obs_vec[k])) cov_path[k] <- max(cov_path[k], coverage_obs_vec[k])
      }
    }
    recovered_gap <- pmax(coverage_cf_vec - cov_path, 0)
    # Express recovery as a multiplicative factor on the raw per-year gap mean.
    fac <- ifelse(gap_mean_vec > 1e-9, recovered_gap / gap_mean_vec, 1)
    return(as.list(pmin(pmax(fac, 0), 1)))
  }
  stop("Unknown framing: ", framing)
}

# ---- Core MC simulator: single country x vaccine x disease x method -------

simulate_yll_one <- function(country_iso, country_name, vaccine, disease,
                             method_name, gap_data, covariates, shared,
                             n_sim = N_SIM, catchup_grid = CATCHUP_RATES) {
  g <- gap_data %>%
    dplyr::filter(country == country_name, vaccine == !!vaccine,
                  method == method_name) %>%
    dplyr::arrange(year)
  if (nrow(g) == 0) return(NULL)
  
  cov_yrs <- covariates %>%
    dplyr::filter(ISO3 == country_iso, year %in% g$year) %>%
    dplyr::arrange(year) %>%
    dplyr::select(-dplyr::any_of("country"))
  if (nrow(cov_yrs) == 0) return(NULL)
  g <- g %>% dplyr::left_join(cov_yrs, by = c("ISO3", "year"))
  
  inc_col <- incidence_col(disease)
  if (!inc_col %in% names(g)) return(NULL)
  sig_col <- incidence_sigma_col(disease)
  has_sig_col <- sig_col %in% names(g)
  
  denom_col <- disease_denominator_col(disease)
  if (!denom_col %in% names(g)) {
    message("    skipping: denominator '", denom_col, "' absent for ",
            disease, " (AGE_BAND_MODE='", AGE_BAND_MODE, "')")
    return(NULL)
  }
  
  cfr_pop  <- cfr_population_for_disease(disease)
  cfr_draws_country <- get_shared_draws(shared, cfr_param_for_disease(disease), cfr_pop)
  ve_draws_country  <- get_shared_draws(shared, ve_param_for_link(vaccine, disease, country_name))
  
  cfr_effective_country <- cfr_draws_country
  
  if (disease == "Measles" &&
      identical(get0("MEASLES_MULT_MODE", ifnotfound = "off"), "on")) {
    mult_draws <- get_shared_draws(shared, "Measles_conflict_mult")
    cfr_effective_country <- pmin(cfr_draws_country * mult_draws, 1)
  }
 
  if (disease == "Tuberculosis" &&
      isTRUE(get0("TB_TX_WEIGHTING", ifnotfound = FALSE))) {
    cfr_tb_treated <- get_shared_draws(shared, "CFR_TB", "Ages 0-4, treated")
    p_tb_treat     <- get_shared_draws(shared, "TB_treat_coverage")
    cfr_effective_country <- p_tb_treat * cfr_tb_treated +
      (1 - p_tb_treat) * cfr_draws_country
  }
  
  reach_draws <- get_shared_draws(shared, "Catchup_reach")
  imm_discount <- if (disease == "Measles") CATCHUP_IMMUNE_DISCOUNT_MEASLES else 0
  
  summary_rows <- list()
  draw_rows    <- list()
  n_years <- nrow(g)
  
  for (lambda in catchup_grid) {
    for (framing in CATCHUP_FRAMINGS) {
      
      gap_factors <- catchup_gap_factors(
        framing          = framing,
        lambda           = lambda,
        gap_mean_vec     = g$gap_prop,
        coverage_obs_vec = g$coverage_obs_prop,
        coverage_cf_vec  = g$coverage_cf_prop,
        reach_draws      = reach_draws,
        immune_discount  = imm_discount
      )
      
      yll_disc_mat   <- matrix(0, nrow = n_years, ncol = n_sim)
      yll_undisc_mat <- matrix(0, nrow = n_years, ncol = n_sim)
      deaths_mat     <- matrix(0, nrow = n_years, ncol = n_sim)
      valid_year     <- logical(n_years)
      
      for (yi in seq_len(n_years)) {
        inc_val    <- g[[inc_col]][yi]
        N_val      <- g[[denom_col]][yi]
        # L = remaining life expectancy used to value each death. Resolved via
        # the life-expectancy layer (01b): GBD reference table (headline) by
        # default, or a local-LE sensitivity. For the infant cohort this is
        # ex at age 0 (~88.87) in reference mode. g$Life_Expectancy is the
        # World Bank value, used only when LIFE_TABLE_MODE == "local_wb".
        L_val      <- resolve_life_expectancy(disease, country_iso, g$year[yi],
                                              local_wb_le = g$Life_Expectancy[yi])
        gap_mean   <- g$gap_prop[yi]
        gap_sd_yi  <- g$gap_sd_prop[yi]
        if (anyNA(c(inc_val, N_val, L_val, gap_mean, gap_sd_yi)) ||
            N_val <= 0 || L_val <= 0) {
          next
        }
        valid_year[yi] <- TRUE
        
        # If the central conflict-attributable gap is <= 0 (coverage held or
        # improved), there is no attributable excess: contribute exactly 0
        # rather than letting clamped symmetric noise invent positive burden.
        if (gap_mean <= 0) {
          # matrices already 0 for this row; leave as valid-but-zero.
          next
        }
        
        inc_sigma <- if (has_sig_col) g[[sig_col]][yi] else NA_real_
        if (is.na(inc_sigma) || inc_sigma <= 0) inc_sigma <- INC_LOGNORM_SIGMA
        inc_draws <- stats::rlnorm(
          n_sim,
          meanlog = log(max(inc_val, 1e-12)) - inc_sigma^2 / 2,
          sdlog   = inc_sigma
        )
        
        
        gap_sd_yr <- max(gap_sd_yi, 1e-6)
        gap_raw   <- draw_gap_beta(n_sim, gap_mean, gap_sd_yr)
        
        # Apply catch-up multiplicatively (scalar or length-n_sim vector).
        gap_draws <- gap_raw * gap_factors[[yi]]
        gap_draws <- pmin(pmax(gap_draws, 0), 1)
        
        disc <- (1 - exp(-DISCOUNT_RATE * L_val)) / (DISCOUNT_RATE * L_val)
        
        deaths     <- N_val * inc_draws * cfr_effective_country * gap_draws * ve_draws_country
        yll_undisc <- deaths * L_val
        yll_disc   <- yll_undisc * disc
        
        deaths_mat[yi, ]     <- deaths
        yll_undisc_mat[yi, ] <- yll_undisc
        yll_disc_mat[yi, ]   <- yll_disc
      }
      
      # ---- Yearly summary rows --------------------------------------------
      for (yi in seq_len(n_years)) {
        if (!valid_year[yi]) {
          summary_rows[[length(summary_rows) + 1]] <- tibble::tibble(
            country = country_name, ISO3 = country_iso, vaccine = vaccine,
            disease = disease, method = method_name,
            year = g$year[yi], catchup = lambda, framing = framing,
            deaths_mean = NA_real_, deaths_lo = NA_real_, deaths_hi = NA_real_,
            yll_undisc_mean = NA_real_, yll_undisc_lo = NA_real_, yll_undisc_hi = NA_real_,
            yll_disc_mean   = NA_real_, yll_disc_lo   = NA_real_, yll_disc_hi   = NA_real_
          )
          next
        }
        summary_rows[[length(summary_rows) + 1]] <- tibble::tibble(
          country  = country_name, ISO3 = country_iso, vaccine = vaccine,
          disease  = disease, method = method_name,
          year = g$year[yi], catchup = lambda, framing = framing,
          deaths_mean     = mean(deaths_mat[yi, ], na.rm = TRUE),
          deaths_lo       = stats::quantile(deaths_mat[yi, ], 0.025, names = FALSE, na.rm = TRUE),
          deaths_hi       = stats::quantile(deaths_mat[yi, ], 0.975, names = FALSE, na.rm = TRUE),
          yll_undisc_mean = mean(yll_undisc_mat[yi, ], na.rm = TRUE),
          yll_undisc_lo   = stats::quantile(yll_undisc_mat[yi, ], 0.025, names = FALSE, na.rm = TRUE),
          yll_undisc_hi   = stats::quantile(yll_undisc_mat[yi, ], 0.975, names = FALSE, na.rm = TRUE),
          yll_disc_mean   = mean(yll_disc_mat[yi, ], na.rm = TRUE),
          yll_disc_lo     = stats::quantile(yll_disc_mat[yi, ], 0.025, names = FALSE, na.rm = TRUE),
          yll_disc_hi     = stats::quantile(yll_disc_mat[yi, ], 0.975, names = FALSE, na.rm = TRUE)
        )
      }
      
      # ---- Period total ----------------------------------------------------
      if (!any(valid_year)) {
        summary_rows[[length(summary_rows) + 1]] <- tibble::tibble(
          country = country_name, ISO3 = country_iso, vaccine = vaccine,
          disease = disease, method = method_name,
          year = NA_integer_, catchup = lambda, framing = framing,
          deaths_mean = NA_real_, deaths_lo = NA_real_, deaths_hi = NA_real_,
          yll_undisc_mean = NA_real_, yll_undisc_lo = NA_real_, yll_undisc_hi = NA_real_,
          yll_disc_mean   = NA_real_, yll_disc_lo   = NA_real_, yll_disc_hi   = NA_real_,
          n_valid_years = 0L, n_years = n_years
        )
        next
      }
      total_deaths     <- colSums(deaths_mat,     na.rm = TRUE)
      total_yll_undisc <- colSums(yll_undisc_mat, na.rm = TRUE)
      total_yll_disc   <- colSums(yll_disc_mat,   na.rm = TRUE)
      
      summary_rows[[length(summary_rows) + 1]] <- tibble::tibble(
        country  = country_name, ISO3 = country_iso, vaccine = vaccine,
        disease  = disease, method = method_name,
        year = NA_integer_, catchup = lambda, framing = framing,
        deaths_mean     = mean(total_deaths, na.rm = TRUE),
        deaths_lo       = stats::quantile(total_deaths, 0.025, names = FALSE, na.rm = TRUE),
        deaths_hi       = stats::quantile(total_deaths, 0.975, names = FALSE, na.rm = TRUE),
        yll_undisc_mean = mean(total_yll_undisc, na.rm = TRUE),
        yll_undisc_lo   = stats::quantile(total_yll_undisc, 0.025, names = FALSE, na.rm = TRUE),
        yll_undisc_hi   = stats::quantile(total_yll_undisc, 0.975, names = FALSE, na.rm = TRUE),
        yll_disc_mean   = mean(total_yll_disc, na.rm = TRUE),
        yll_disc_lo     = stats::quantile(total_yll_disc, 0.025, names = FALSE, na.rm = TRUE),
        yll_disc_hi     = stats::quantile(total_yll_disc, 0.975, names = FALSE, na.rm = TRUE),
        # Expose how many conflict years actually contributed, so a "total" over
        # a partially-covered period is visible rather than silently partial.
        n_valid_years = sum(valid_year), n_years = n_years
      )
      
      draw_rows[[length(draw_rows) + 1]] <- tibble::tibble(
        country = country_name, ISO3 = country_iso, vaccine = vaccine,
        disease = disease, method = method_name,
        catchup = lambda, framing = framing,
        n_valid_years = sum(valid_year), n_years = n_years,
        d_undisc = list(total_yll_undisc),
        d_disc   = list(total_yll_disc)
      )
    }
  }
  list(summary = dplyr::bind_rows(summary_rows),
       draws   = dplyr::bind_rows(draw_rows))
}

# ---- Top-level driver -----------------------------------------------------
run_full_yll <- function(gap_data, covariates, conflict_df = conflict_info,
                         vaccine_disease = vaccine_disease_links,
                         n_sim = N_SIM, seed = 42,
                         methods = NULL) {
  set.seed(seed)
  
  if (is.null(methods)) methods <- unique(gap_data$method)
  message("YLL methods: ", paste(methods, collapse = ", "))
  
  shared <- draw_shared_params(n_sim, conflict_df = conflict_df,
                               vd = vaccine_disease)
  
  all_summary <- list()
  all_draws   <- list()
  for (i in seq_len(nrow(conflict_df))) {
    ci <- conflict_df[i, ]
    for (j in seq_len(nrow(vaccine_disease))) {
      vd <- vaccine_disease[j, ]
      for (m in methods) {
        message("YLL MC: ", ci$country, " | ", vd$vaccine, " -> ",
                vd$disease, " | ", m)
        res <- simulate_yll_one(
          country_iso  = ci$ISO3,
          country_name = ci$country,
          vaccine      = vd$vaccine,
          disease      = vd$disease,
          method_name  = m,
          gap_data     = gap_data,
          covariates   = covariates,
          shared       = shared,
          n_sim        = n_sim
        )
        if (!is.null(res)) {
          all_summary[[length(all_summary) + 1]] <- res$summary
          if (!is.null(res$draws) && nrow(res$draws) > 0)
            all_draws[[length(all_draws) + 1]] <- res$draws
        }
      }
    }
  }
  list(summary = dplyr::bind_rows(all_summary),
       draws   = dplyr::bind_rows(all_draws))
}

# ---- Aggregation helpers --------------------------------------------------
# Country-period totals are already computed inside simulate_yll_one as rows
# with year == NA. Pull them out and rename for clarity.
aggregate_country_yll <- function(yll_df) {
  yll_df %>%
    dplyr::filter(is.na(year)) %>%
    dplyr::rename(
      yll_disc_total_mean   = yll_disc_mean,
      yll_disc_total_lo     = yll_disc_lo,
      yll_disc_total_hi     = yll_disc_hi,
      yll_undisc_total_mean = yll_undisc_mean,
      yll_undisc_total_lo   = yll_undisc_lo,
      yll_undisc_total_hi   = yll_undisc_hi,
      deaths_total_mean     = deaths_mean,
      deaths_total_lo       = deaths_lo,
      deaths_total_hi       = deaths_hi
    ) %>%
    dplyr::select(-year)
}

# Yearly rows only (for trend plots)
yearly_yll <- function(yll_df) {
  yll_df %>% dplyr::filter(!is.na(year))
}


aggregate_disease_yll <- function(yll_draws) {
  if (is.null(yll_draws) || nrow(yll_draws) == 0) return(tibble::tibble())
  keep <- vapply(yll_draws$d_undisc,
                 function(x) !is.null(x) && length(x) > 0, logical(1))
  yll_draws <- yll_draws[keep, , drop = FALSE]
  yll_draws %>%
    dplyr::group_by(disease, vaccine, method, catchup, framing) %>%
    dplyr::group_modify(~ {
      M_un <- do.call(rbind, .x$d_undisc)   # n_country x n_sim
      M_di <- do.call(rbind, .x$d_disc)
      s_un <- colSums(M_un)                  # length n_sim (global-by-group)
      s_di <- colSums(M_di)
      tibble::tibble(
        n_countries     = nrow(M_un),
        yll_undisc_mean = mean(s_un),
        yll_undisc_lo   = stats::quantile(s_un, 0.025, names = FALSE),
        yll_undisc_hi   = stats::quantile(s_un, 0.975, names = FALSE),
        yll_disc_mean   = mean(s_di),
        yll_disc_lo     = stats::quantile(s_di, 0.025, names = FALSE),
        yll_disc_hi     = stats::quantile(s_di, 0.975, names = FALSE)
      )
    }) %>%
    dplyr::ungroup()
}

aggregate_method_comparison_draws <- function(yll_draws, catchup_focus = 0,
                                              framing_focus = "campaign_topup") {
  if (is.null(yll_draws) || nrow(yll_draws) == 0) return(tibble::tibble())
  d <- yll_draws %>%
    dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(tibble::tibble())
  
  n_methods <- dplyr::n_distinct(d$method)
  complete  <- d %>%
    dplyr::distinct(disease, country, method) %>%
    dplyr::count(disease, country, name = "nm") %>%
    dplyr::filter(nm == n_methods) %>%
    dplyr::select(disease, country)
  
  dropped <- d %>%
    dplyr::distinct(disease, country) %>%
    dplyr::anti_join(complete, by = c("disease", "country"))
  if (nrow(dropped) > 0) {
    message("Table 2 complete-case (draw-level): excluding ", nrow(dropped),
            " country-disease cells absent from >=1 of ", n_methods, " methods.")
  }
  
  d %>%
    dplyr::semi_join(complete, by = c("disease", "country")) %>%
    dplyr::group_by(disease, method) %>%
    dplyr::group_modify(~ {
      M_un <- do.call(rbind, .x$d_undisc)
      M_di <- do.call(rbind, .x$d_disc)
      s_un <- colSums(M_un)
      s_di <- colSums(M_di)
      tibble::tibble(
        n_countries     = nrow(M_un),
        yll_undisc_mean = mean(s_un),
        yll_undisc_lo   = stats::quantile(s_un, 0.025, names = FALSE),
        yll_undisc_hi   = stats::quantile(s_un, 0.975, names = FALSE),
        yll_disc_mean   = mean(s_di),
        yll_disc_lo     = stats::quantile(s_di, 0.025, names = FALSE),
        yll_disc_hi     = stats::quantile(s_di, 0.975, names = FALSE)
      )
    }) %>%
    dplyr::ungroup()
}
