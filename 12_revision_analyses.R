# ============================================================================
# 12_REVISION_ANALYSES.R   (supplementary; PRIMARY spec unchanged)
# ----------------------------------------------------------------------------
# Additive revision analyses requested by review. NONE alters the headline
# (Tables 1-2, Figs 1-4). Each function self-skips on missing input. Five blocks:
#
#   (R5)  Temporal concentration / peak single-year burden.
#         Reviewer: cumulative + per-capita both hide acute one-cohort collapses
#         (Myanmar) vs chronic erosion (Nigeria). We add peak single-year coverage
#         gap (model-free, from gap_df) and peak single-year YLL (from the retained
#         yearly summary), plus a concentration index = peak-year / cumulative.
#
#   (R1)  COVID-19 confounding of the conflict window.
#         (a) Mean-level pandemic-cohort share: fraction of each country's
#             composite YLL accruing in the 2020-2021 (and 2020-2022) cohort-years
#             -- "how much of each country's gap survives" -- with NO re-run.
#         (b) Draw-level COVID-excluded headline: one re-run of run_full_yll on a
#             gap_df with 2020-2021 cohort-years dropped, re-assembled through the
#             SAME composite as the headline. Opt-in (re-run); CRN-matched seed.
#
#   (R2)  Country dominance / leave-one-out. Composite headline recomputed with
#         each country removed (Nigeria-excluded is the headline row), exact and
#         draw-level (the composite is a sum of per-country draw vectors).
#
#   (R10) Estimator-selection sensitivity of the composite. The 10-country
#         composite recomputed with each estimator placed FIRST in the priority
#         (best-available fills the rest), exposing how much of the total is an
#         artefact of the SC->DiD->baseline selection rule rather than the data.
#
#   (R7)  Convention corners. The four-corner table (undiscounted/discounted x
#         reference/national life table) computed EXACTLY at the draw level from
#         the existing d_undisc/d_disc draws (infant deaths use a single age, so
#         the national-LE discount factor is analytic), giving a co-equal
#         "conservative-conventions" headline beside the maximal one.
#
#   (R3)  [opt-in re-run] Pertussis incidence-uncertainty stress. Incidence is
#         ALREADY lognormally sampled from the GBD 95% UI (02_load_data.R); this
#         re-runs the headline with the pertussis log-sigma floored at the wider
#         published modelling range (deaths ~38k-670k, ~18-fold) to show the
#         decomposition is stable to the shakiest input.
#
# Source AFTER 04/05/05c and 09 (uses yll_raw, gap_df, covariates_panel,
# conflict_info, assemble_headline_totals, run_full_yll, ref_ex_at_age,
# FRAMING_HEADLINE, DISCOUNT_RATE, save_tab/OUT_DIR).
# ============================================================================

# ---- Shared small helpers --------------------------------------------------

# Reference ex0 (headline life table), matching 09's .REF_EX0 fallback.
.rev_ref_ex0 <- function() {
  if (exists("ref_ex_at_age")) ref_ex_at_age(0) else 88.8718951
}

.rev_framing <- function() if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup"
.rev_rate    <- function() if (exists("DISCOUNT_RATE")) DISCOUNT_RATE else 0.03

# Draw-vector -> c(mean, lo, hi). NA-safe for an absent vector.
.rev_q <- function(v) {
  if (is.null(v) || (length(v) == 1 && is.na(v)))
    return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  c(mean = mean(v),
    lo = stats::quantile(v, 0.025, names = FALSE),
    hi = stats::quantile(v, 0.975, names = FALSE))
}

# Sum a list of equal-length draw vectors element-wise (draw index preserved).
.rev_sum_draws <- function(vlist) {
  vlist <- vlist[vapply(vlist, function(x) !is.null(x) && length(x) > 0, logical(1))]
  if (length(vlist) == 0) return(NULL)
  colSums(do.call(rbind, vlist))
}

# Composite cell PICK (rows, not just the total): mirrors assemble_headline_totals
# / .composite_cmap EXACTLY (same anchor_priority, same slice_min), but returns
# the picked draw rows so we can leave-one-out, rescale, and restrict the yearly
# summary to the same cells. Filtered to the headline (catchup=0, headline framing).
.rev_composite_pick <- function(draws,
                                anchor_priority = c("Synthetic Control",
                                                    "DiD (CS)", "Baseline (3yr)"),
                                catchup_focus = 0,
                                framing_focus = .rev_framing()) {
  d <- draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(d)
  d %>% dplyr::filter(method %in% anchor_priority) %>%
    dplyr::mutate(prio = match(method, anchor_priority)) %>%
    dplyr::group_by(disease, country) %>%
    dplyr::slice_min(prio, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-prio)
}

.rev_save <- function(tbl, file) {
  if (is.null(tbl) || nrow(tbl) == 0) return(invisible(NULL))
  if (exists("save_tab")) save_tab(tbl, file)
  else if (exists("OUT_DIR"))
    utils::write.csv(tbl, file.path(OUT_DIR, file), row.names = FALSE)
  invisible(tbl)
}

# ============================================================================
# (R5)  TEMPORAL CONCENTRATION / PEAK SINGLE-YEAR BURDEN
# ============================================================================
# Peak single-year coverage gap is model-free (max over post-onset year x antigen
# of the percentage-point gap, taken from gap_df). Peak single-year YLL uses the
# retained yearly summary, restricted to the SAME composite cells as the headline,
# summed across diseases within country-year. Peak-year YLL is a point estimate
# (the per-year draw vectors are not retained; the MC interval lives on the
# cumulative total) -- stated as such. The concentration index peak/cumulative
# separates acute one-cohort collapses from chronic multi-year erosion.
revision_temporal_concentration <- function(yll_summary, gap_df,
                                            yll_draws,
                                            conflict_df = conflict_info,
                                            framing_focus = .rev_framing()) {
  if (is.null(yll_summary) || nrow(yll_summary) == 0) {
    message("  (R5) temporal concentration: empty summary -> skip."); return(invisible(NULL))
  }
  pick <- .rev_composite_pick(yll_draws, framing_focus = framing_focus)
  if (nrow(pick) == 0) { message("  (R5): empty composite pick -> skip."); return(invisible(NULL)) }
  cell_map <- pick %>% dplyr::distinct(country, disease, method)
  
  # --- peak & cumulative YLL from the yearly summary (composite cells) --------
  yr <- yll_summary %>%
    dplyr::filter(!is.na(year), catchup == 0, framing == framing_focus) %>%
    dplyr::semi_join(cell_map, by = c("country", "disease", "method")) %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(yll_year = sum(yll_undisc_mean, na.rm = TRUE), .groups = "drop")
  
  peak_yll <- yr %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      cumulative_yll   = sum(yll_year, na.rm = TRUE),
      peak_year_yll    = max(yll_year, na.rm = TRUE),
      peak_yll_year    = year[which.max(yll_year)],
      n_years_active   = sum(yll_year > 0, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(
      concentration_index = ifelse(cumulative_yll > 0,
                                   peak_year_yll / cumulative_yll, NA_real_))
  
  # --- peak single-year coverage gap (percentage points), model-free ---------
  # Restricted to ONE interpretable estimator (the flat 3-year baseline = the raw
  # drop from the pre-conflict mean, the construct in Fig. 1b) rather than the max
  # across estimators, which would silently pick whichever method gives the
  # largest gap in any year and overstate the acute peak.
  base_pref <- c("Baseline (3yr)", "Baseline (5yr)", "Baseline (1yr)")
  base_method <- intersect(base_pref, unique(gap_df$method))
  gap_base <- if (length(base_method) > 0)
    dplyr::filter(gap_df, method == base_method[1])
  else gap_df   # fall back to all rows only if no baseline estimator is present
  peak_gap <- gap_base %>%
    dplyr::semi_join(cell_map %>% dplyr::distinct(country), by = "country") %>%
    dplyr::filter(is.finite(gap_prop)) %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      peak_gap_pp      = max(gap_prop, na.rm = TRUE) * 100,
      peak_gap_antigen = vaccine[which.max(gap_prop)],
      peak_gap_year    = year[which.max(gap_prop)],
      .groups = "drop")
  
  out <- peak_yll %>%
    dplyr::left_join(peak_gap, by = "country") %>%
    dplyr::left_join(conflict_df %>%
                       dplyr::transmute(country,
                                        window_years = conflict_end - conflict_start + 1L),
                     by = "country") %>%
    dplyr::arrange(dplyr::desc(concentration_index)) %>%
    dplyr::mutate(
      profile = dplyr::case_when(
        is.na(concentration_index)        ~ NA_character_,
        concentration_index >= 0.50       ~ "acute (single-cohort)",
        concentration_index >= 0.25       ~ "intermediate",
        TRUE                              ~ "chronic (multi-cohort)"))
  
  .rev_save(out, "Table_S8_temporal_concentration.csv")
  message("  (R5) temporal concentration: ", nrow(out), " countries. ",
          "Most concentrated: ",
          paste(utils::head(out$country[!is.na(out$concentration_index)], 3),
                collapse = ", "), ".")
  invisible(out)
}

# ============================================================================
# (R1a) PANDEMIC-COHORT SHARE  (mean-level, no re-run)
# ============================================================================
# What fraction of each country's composite YLL accrues in pandemic cohort-years?
# Computed from the yearly summary on the SAME composite cells, at the mean level
# (exact for a mean, since the period mean is the sum of yearly means). This is
# the "how much of the gap survives dropping COVID years" figure, country by
# country, with no Monte Carlo re-run.
revision_pandemic_share <- function(yll_summary, yll_draws,
                                    pandemic_years = 2020:2021,
                                    pandemic_years_wide = 2020:2022,
                                    framing_focus = .rev_framing()) {
  if (is.null(yll_summary) || nrow(yll_summary) == 0) {
    message("  (R1a) pandemic share: empty summary -> skip."); return(invisible(NULL))
  }
  pick <- .rev_composite_pick(yll_draws, framing_focus = framing_focus)
  if (nrow(pick) == 0) { message("  (R1a): empty composite pick -> skip."); return(invisible(NULL)) }
  cell_map <- pick %>% dplyr::distinct(country, disease, method)
  
  yr <- yll_summary %>%
    dplyr::filter(!is.na(year), catchup == 0, framing == framing_focus) %>%
    dplyr::semi_join(cell_map, by = c("country", "disease", "method")) %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(yll_year = sum(yll_undisc_mean, na.rm = TRUE), .groups = "drop")
  
  by_country <- yr %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      total_yll          = sum(yll_year, na.rm = TRUE),
      pandemic_yll       = sum(yll_year[year %in% pandemic_years], na.rm = TRUE),
      pandemic_yll_wide  = sum(yll_year[year %in% pandemic_years_wide], na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(
      pct_in_pandemic       = ifelse(total_yll > 0, 100 * pandemic_yll / total_yll, NA_real_),
      pct_in_pandemic_wide  = ifelse(total_yll > 0, 100 * pandemic_yll_wide / total_yll, NA_real_),
      pct_surviving_drop    = ifelse(total_yll > 0, 100 * (total_yll - pandemic_yll) / total_yll, NA_real_)) %>%
    dplyr::arrange(dplyr::desc(pct_in_pandemic))
  
  overall <- tibble::tibble(
    country               = "ALL (composite)",
    total_yll             = sum(by_country$total_yll, na.rm = TRUE),
    pandemic_yll          = sum(by_country$pandemic_yll, na.rm = TRUE),
    pandemic_yll_wide     = sum(by_country$pandemic_yll_wide, na.rm = TRUE)) %>%
    dplyr::mutate(
      pct_in_pandemic      = 100 * pandemic_yll / total_yll,
      pct_in_pandemic_wide = 100 * pandemic_yll_wide / total_yll,
      pct_surviving_drop   = 100 * (total_yll - pandemic_yll) / total_yll)
  
  out <- dplyr::bind_rows(by_country, overall)
  .rev_save(out, "Table_S9_pandemic_cohort_share.csv")
  message(sprintf("  (R1a) pandemic-cohort share (mean-level): %.1f%% of composite YLL in %d-%d; %.1f%% survives dropping them.",
                  overall$pct_in_pandemic, min(pandemic_years), max(pandemic_years),
                  overall$pct_surviving_drop))
  invisible(out)
}

# ============================================================================
# (R2)  LEAVE-ONE-OUT COMPOSITE HEADLINE
# ============================================================================
# Composite total with each country removed. Exact and draw-level: the composite
# is the element-wise sum of the per-country picked draw vectors, so dropping a
# country just removes its rows. Nigeria-excluded is the row reviewers asked for.
revision_leave_one_out <- function(yll_draws, framing_focus = .rev_framing()) {
  pick <- .rev_composite_pick(yll_draws, framing_focus = framing_focus)
  if (nrow(pick) == 0) { message("  (R2) LOO: empty composite pick -> skip."); return(invisible(NULL)) }
  
  full <- .rev_sum_draws(pick$d_undisc)
  full_q <- .rev_q(full)
  countries <- sort(unique(pick$country))
  
  rows <- lapply(countries, function(cc) {
    keep_rows <- pick %>% dplyr::filter(country != cc)
    drop_rows <- pick %>% dplyr::filter(country == cc)
    loo  <- .rev_sum_draws(keep_rows$d_undisc)
    ctry <- .rev_sum_draws(drop_rows$d_undisc)
    loo_q  <- .rev_q(loo)
    ctry_q <- .rev_q(ctry)
    tibble::tibble(
      excluded_country       = cc,
      country_yll_mean       = ctry_q["mean"],
      country_yll_lo         = ctry_q["lo"],
      country_yll_hi         = ctry_q["hi"],
      pct_of_headline        = 100 * ctry_q["mean"] / full_q["mean"],
      headline_excl_mean     = loo_q["mean"],
      headline_excl_lo       = loo_q["lo"],
      headline_excl_hi       = loo_q["hi"])
  })
  
  out <- dplyr::bind_rows(
    tibble::tibble(
      excluded_country   = "(none: full composite)",
      country_yll_mean   = NA_real_, country_yll_lo = NA_real_, country_yll_hi = NA_real_,
      pct_of_headline    = NA_real_,
      headline_excl_mean = full_q["mean"],
      headline_excl_lo   = full_q["lo"],
      headline_excl_hi   = full_q["hi"]),
    dplyr::bind_rows(rows) %>% dplyr::arrange(dplyr::desc(pct_of_headline)))
  
  .rev_save(out, "Table_S10_leave_one_out.csv")
  ng <- out %>% dplyr::filter(excluded_country == "Nigeria")
  if (nrow(ng) == 1)
    message(sprintf("  (R2) Nigeria = %.0f%% of headline; composite excluding Nigeria = %s (95%% CI %s-%s).",
                    ng$pct_of_headline,
                    format(round(ng$headline_excl_mean), big.mark = ","),
                    format(round(ng$headline_excl_lo),   big.mark = ","),
                    format(round(ng$headline_excl_hi),   big.mark = ",")))
  invisible(out)
}

# ============================================================================
# (R10) ESTIMATOR-SELECTION SENSITIVITY OF THE COMPOSITE
# ============================================================================
# The 10-country composite recomputed with EACH estimator placed first in the
# priority (best-available fills cells that estimator cannot). Shows how far the
# headline moves under the selection rule. ITS is included for completeness but
# is a lower bound by construction (structural zeros in Ethiopia/Iraq).
revision_composite_by_anchor <- function(yll_draws, framing_focus = .rev_framing()) {
  d0 <- yll_draws %>% dplyr::filter(catchup == 0, framing == framing_focus)
  keep <- vapply(d0$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d0 <- d0[keep, , drop = FALSE]
  if (nrow(d0) == 0) { message("  (R10): empty draws -> skip."); return(invisible(NULL)) }
  present <- unique(d0$method)
  # A sensible global fallback ordering (donor -> DiD -> baseline -> ITS), used to
  # fill cells the first-placed estimator cannot supply.
  pref <- c("Synthetic Control", "SC + covariates", "Augmented SC",
            "DiD (CS)", "Baseline (3yr)", "Baseline (1yr)", "Baseline (5yr)", "ITS")
  pref <- intersect(pref, present)
  
  rows <- lapply(present, function(M) {
    priority <- unique(c(M, pref))
    pick <- .rev_composite_pick(yll_draws, anchor_priority = priority,
                                framing_focus = framing_focus)
    if (nrow(pick) == 0) return(NULL)
    qq <- .rev_q(.rev_sum_draws(pick$d_undisc))
    mix <- pick %>% dplyr::count(method, name = "n") %>%
      dplyr::arrange(dplyr::desc(n))
    tibble::tibble(
      anchor_first   = M,
      n_countries    = dplyr::n_distinct(pick$country),
      n_cells_anchor = sum(pick$method == M),
      yll_mean       = qq["mean"], yll_lo = qq["lo"], yll_hi = qq["hi"],
      cell_mix       = paste(sprintf("%s:%d", mix$method, mix$n), collapse = "; "))
  })
  out <- dplyr::bind_rows(rows) %>% dplyr::arrange(dplyr::desc(yll_mean))
  .rev_save(out, "Table_S11_composite_by_anchor.csv")
  if (nrow(out) > 0)
    message(sprintf("  (R10) composite spans %s - %s YLL across estimator-first rules.",
                    format(round(min(out$yll_mean)), big.mark = ","),
                    format(round(max(out$yll_mean)), big.mark = ",")))
  invisible(out)
}

# ============================================================================
# (R7)  CONVENTION CORNERS  (exact, draw-level)
# ============================================================================
# Four corners on the SAME composite cells:
#   undisc x reference LE   = d_undisc                       (the headline)
#   disc   x reference LE   = d_disc                         (exact)
#   undisc x national LE    = deaths * L_nat                 (exact; deaths=d_undisc/L_ref)
#   disc   x national LE    = deaths * L_nat * disc(L_nat)   (exact; analytic discount)
# L_ref is the per-cell reference ex at the disease's age at death (= ex0 for the
# infant cohort). L_nat is the mean national life expectancy over the country's
# conflict window (matching 09's life-table sensitivity convention).
revision_convention_corners <- function(yll_draws,
                                        covariates = covariates_panel,
                                        conflict_df = conflict_info,
                                        framing_focus = .rev_framing()) {
  pick <- .rev_composite_pick(yll_draws, framing_focus = framing_focus)
  if (nrow(pick) == 0) { message("  (R7): empty composite pick -> skip."); return(invisible(NULL)) }
  r <- .rev_rate()
  disc_factor <- function(L) ifelse(L > 0, (1 - exp(-r * L)) / (r * L), 1)
  
  # Per-cell reference L (robust to the under-5 age-band mode).
  L_ref_cell <- if (exists("age_at_death_for_disease") && exists("ref_ex_at_age"))
    vapply(pick$disease, function(dis) ref_ex_at_age(age_at_death_for_disease(dis)), numeric(1))
  else rep(.rev_ref_ex0(), nrow(pick))
  
  # Per-country national LE over the conflict window (keyed by ISO3).
  le_nat <- conflict_df %>%
    dplyr::select(ISO3, conflict_start, conflict_end) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(L_nat = mean(covariates$Life_Expectancy[
      covariates$ISO3 == ISO3 &
        covariates$year >= conflict_start &
        covariates$year <= conflict_end], na.rm = TRUE)) %>%
    dplyr::ungroup()
  le_map <- stats::setNames(le_nat$L_nat, le_nat$ISO3)
  
  S_uu <- S_ud <- S_nu <- S_nd <- NULL
  for (i in seq_len(nrow(pick))) {
    v  <- pick$d_undisc[[i]]
    vd <- pick$d_disc[[i]]
    if (is.null(v) || length(v) == 0) next
    Lr <- L_ref_cell[i]
    Ln <- le_map[[pick$ISO3[i]]]
    if (is.null(Ln) || !is.finite(Ln)) Ln <- Lr   # no rescale if LE missing
    deaths <- v / Lr
    nat_u  <- deaths * Ln
    nat_d  <- deaths * Ln * disc_factor(Ln)
    S_uu <- if (is.null(S_uu)) v      else S_uu + v
    S_ud <- if (is.null(S_ud)) vd     else S_ud + vd
    S_nu <- if (is.null(S_nu)) nat_u  else S_nu + nat_u
    S_nd <- if (is.null(S_nd)) nat_d  else S_nd + nat_d
  }
  mk <- function(label, v, conservatism) {
    q <- .rev_q(v)
    tibble::tibble(scenario = label, life_table = conservatism$lt,
                   discounting = conservatism$disc,
                   yll_mean = q["mean"], yll_lo = q["lo"], yll_hi = q["hi"])
  }
  out <- dplyr::bind_rows(
    mk("Maximal headline (undisc, reference LE)", S_uu, list(lt = "reference", disc = "none")),
    mk("Discounted, reference LE",                S_ud, list(lt = "reference", disc = "3%")),
    mk("Undiscounted, national LE",              S_nu, list(lt = "national",  disc = "none")),
    mk("Conservative (disc 3% + national LE)",   S_nd, list(lt = "national",  disc = "3%")))
  head_mean <- out$yll_mean[1]
  out <- out %>% dplyr::mutate(pct_vs_headline = 100 * (yll_mean - head_mean) / head_mean)
  
  .rev_save(out, "Table_S12_convention_corners.csv")
  message(sprintf("  (R7) corners: maximal = %s; conservative (disc + national LE) = %s (%.0f%% of headline).",
                  format(round(out$yll_mean[1]), big.mark = ","),
                  format(round(out$yll_mean[4]), big.mark = ","),
                  100 * out$yll_mean[4] / head_mean))
  invisible(out)
}

# ============================================================================
# RE-RUN HARNESS  (used by R1b and R3; opt-in)
# ============================================================================
# Reuses run_full_yll + assemble_headline_totals (canonical). gap_transform and
# cov_transform are identity by default; pass a function to modify the gap panel
# (drop cohort-years) or the covariate panel (inflate an incidence sigma). The
# seed is held at the headline value so the SHARED structural draws are
# CRN-matched to the main run.
rerun_yll_headline <- function(gap_df, covariates = covariates_panel,
                               conflict_df = conflict_info,
                               vaccine_disease = vaccine_disease_links,
                               gap_transform = function(x) x,
                               cov_transform = function(x) x,
                               n_sim = 3000L, seed = 42, methods = NULL,
                               anchor_priority = c("Synthetic Control",
                                                   "DiD (CS)", "Baseline (3yr)")) {
  gd <- gap_transform(gap_df)
  cv <- cov_transform(covariates)
  yr <- run_full_yll(gd, cv, conflict_df = conflict_df,
                     vaccine_disease = vaccine_disease,
                     n_sim = n_sim, seed = seed, methods = methods)
  list(headline = assemble_headline_totals(yr$draws, anchor_priority = anchor_priority),
       draws = yr$draws)
}

# ---- (R1b) draw-level COVID-excluded headline ------------------------------
revision_covid_excluded_headline <- function(gap_df, yll_draws,
                                             covariates = covariates_panel,
                                             drop_years = 2020:2021,
                                             n_sim = 3000L,
                                             framing_focus = .rev_framing()) {
  full <- assemble_headline_totals(yll_draws)
  if (is.null(full)) { message("  (R1b): no full headline -> skip."); return(invisible(NULL)) }
  drop_fun <- function(x) dplyr::filter(x, !(year %in% drop_years))
  rr <- tryCatch(
    rerun_yll_headline(gap_df, covariates = covariates,
                       gap_transform = drop_fun, n_sim = n_sim),
    error = function(e) { message("  (R1b) re-run failed: ", conditionMessage(e)); NULL })
  if (is.null(rr) || is.null(rr$headline)) return(invisible(NULL))
  f <- full$composite_10country; x <- rr$headline$composite_10country
  out <- tibble::tibble(
    scenario = c("Full composite (all cohort-years)",
                 sprintf("COVID-excluded (drop %d-%d cohort-years)",
                         min(drop_years), max(drop_years))),
    yll_mean = c(f["mean"], x["mean"]),
    yll_lo   = c(f["lo"],   x["lo"]),
    yll_hi   = c(f["hi"],   x["hi"]),
    n_sim    = c(NA_integer_, n_sim))
  out <- out %>% dplyr::mutate(
    pct_of_full = 100 * yll_mean / yll_mean[1])
  .rev_save(out, "Table_S9b_covid_excluded_headline.csv")
  message(sprintf("  (R1b) COVID-excluded composite = %s (%.0f%% of the %s headline).",
                  format(round(x["mean"]), big.mark = ","),
                  100 * x["mean"] / f["mean"],
                  format(round(f["mean"]), big.mark = ",")))
  invisible(out)
}

# ---- (R3) pertussis incidence-uncertainty stress ---------------------------
# Floors the pertussis log-sigma at `target_sigma` (default ~0.73, from the
# published under-5 pertussis death range 38k-670k: ln(670/38)/3.92). Incidence
# is already lognormally sampled from the GBD UI; this widens ONLY pertussis to
# the broader external modelling range and reports the change in the headline CI.
revision_pertussis_uncertainty <- function(gap_df, yll_draws,
                                           covariates = covariates_panel,
                                           target_sigma = log(670 / 38) / 3.92,
                                           n_sim = 3000L) {
  full <- assemble_headline_totals(yll_draws)
  if (is.null(full)) { message("  (R3): no full headline -> skip."); return(invisible(NULL)) }
  if (!"Pertussis_Inc_Sigma" %in% names(covariates)) {
    message("  (R3): Pertussis_Inc_Sigma absent -> skip."); return(invisible(NULL))
  }
  infl <- function(x) {
    s <- x$Pertussis_Inc_Sigma
    s[is.na(s)] <- target_sigma
    x$Pertussis_Inc_Sigma <- pmax(s, target_sigma)
    x
  }
  rr <- tryCatch(
    rerun_yll_headline(gap_df, covariates = covariates,
                       cov_transform = infl, n_sim = n_sim),
    error = function(e) { message("  (R3) re-run failed: ", conditionMessage(e)); NULL })
  if (is.null(rr) || is.null(rr$headline)) return(invisible(NULL))
  f <- full$composite_10country; x <- rr$headline$composite_10country
  out <- tibble::tibble(
    scenario = c("Headline (GBD-UI incidence sigma)",
                 sprintf("Pertussis sigma floored at %.2f (38k-670k range)", target_sigma)),
    yll_mean = c(f["mean"], x["mean"]),
    yll_lo   = c(f["lo"],   x["lo"]),
    yll_hi   = c(f["hi"],   x["hi"])) %>%
    dplyr::mutate(ci_width = yll_hi - yll_lo)
  .rev_save(out, "Table_S13_pertussis_uncertainty.csv")
  message(sprintf("  (R3) pertussis-sigma stress: CI %s-%s -> %s-%s (mean %s -> %s).",
                  format(round(f["lo"]), big.mark = ","), format(round(f["hi"]), big.mark = ","),
                  format(round(x["lo"]), big.mark = ","), format(round(x["hi"]), big.mark = ","),
                  format(round(f["mean"]), big.mark = ","), format(round(x["mean"]), big.mark = ",")))
  invisible(out)
}

# ============================================================================
# DRIVER
# ============================================================================
# No-re-run blocks (R5, R1a, R2, R10, R7) always run. The re-run blocks (R1b, R3)
# are opt-in via CVD_RUN_REVISION_RERUNS=TRUE (each re-runs run_full_yll once).
run_revision_analyses <- function(yll_raw, gap_df,
                                  covariates = covariates_panel,
                                  conflict_df = conflict_info,
                                  rerun_nsim = 3000L) {
  message("\n>>> Revision analyses (R5/R1/R2/R10/R7) <<<\n")
  res <- list()
  res$temporal   <- tryCatch(revision_temporal_concentration(yll_raw$summary, gap_df, yll_raw$draws,
                                                             conflict_df = conflict_df),
                             error = function(e) { message("  (R5) skipped: ", conditionMessage(e)); NULL })
  res$pandemic   <- tryCatch(revision_pandemic_share(yll_raw$summary, yll_raw$draws),
                             error = function(e) { message("  (R1a) skipped: ", conditionMessage(e)); NULL })
  res$loo        <- tryCatch(revision_leave_one_out(yll_raw$draws),
                             error = function(e) { message("  (R2) skipped: ", conditionMessage(e)); NULL })
  res$by_anchor  <- tryCatch(revision_composite_by_anchor(yll_raw$draws),
                             error = function(e) { message("  (R10) skipped: ", conditionMessage(e)); NULL })
  res$corners    <- tryCatch(revision_convention_corners(yll_raw$draws, covariates = covariates,
                                                         conflict_df = conflict_df),
                             error = function(e) { message("  (R7) skipped: ", conditionMessage(e)); NULL })
  
  run_reruns <- isTRUE(as.logical(Sys.getenv("CVD_RUN_REVISION_RERUNS", unset = "TRUE")))
  if (run_reruns) {
    message("\n  Re-run blocks ON (CVD_RUN_REVISION_RERUNS=TRUE; n_sim = ", rerun_nsim, ").")
    res$covid      <- tryCatch(revision_covid_excluded_headline(gap_df, yll_raw$draws,
                                                                covariates = covariates, n_sim = rerun_nsim),
                               error = function(e) { message("  (R1b) skipped: ", conditionMessage(e)); NULL })
    res$pertussis  <- tryCatch(revision_pertussis_uncertainty(gap_df, yll_raw$draws,
                                                              covariates = covariates, n_sim = rerun_nsim),
                               error = function(e) { message("  (R3) skipped: ", conditionMessage(e)); NULL })
  } else {
    message("\n  Re-run blocks OFF (set CVD_RUN_REVISION_RERUNS=TRUE for R1b COVID-excluded + R3 pertussis-sigma; each re-runs run_full_yll once).")
  }
  invisible(res)
}