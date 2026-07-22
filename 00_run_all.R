# ============================================================================
# 00_RUN_ALL.R
# ----------------------------------------------------------------------------
# Master driver. Sources every file in the pipeline in order and writes all
# outputs to OUT_DIR.
#
# USAGE:
#   1. Edit DATA_DIR / WHO_DIR in 01_setup_and_parameters.R if your data is
#      in a non-default location.
#   2. In R/RStudio: setwd() to the directory containing these scripts.
#   3. source("00_run_all.R")
# ============================================================================

script_dir <- '~/Documents/Conflict Vaccination DALYs/Output_New/'
# Fall back gracefully if that path isn't present on this machine (e.g. a
# co-author's checkout) without overriding it when it IS: prefer the env var,
# then the hard-coded path, then the working directory.
.cand_dirs <- c(Sys.getenv("CVD_SCRIPT_DIR", unset = NA),
                script_dir, getwd())
.cand_dirs <- .cand_dirs[!is.na(.cand_dirs)]
.hit <- .cand_dirs[file.exists(file.path(.cand_dirs, "01_setup_and_parameters.R"))]
if (length(.hit) > 0) script_dir <- .hit[1]
message("Pipeline directory: ", script_dir)
if (!file.exists(file.path(script_dir, "01_setup_and_parameters.R"))) {
  stop("01_setup_and_parameters.R not found in ", script_dir,
       "\n  Set CVD_SCRIPT_DIR or setwd() to the directory containing the pipeline scripts.")
}

# ---- Source pipeline ------------------------------------------------------
source(file.path(script_dir, "01_setup_and_parameters.R"))
source(file.path(script_dir, "02_load_data.R"))
source(file.path(script_dir, "02b_econ_covariates.R"))   # GDP/HDI predictors for SC+cov & DiD (after 02)
source(file.path(script_dir, "01b_life_expectancy.R"))   # L source layer (after 02, before 04)
source(file.path(script_dir, "03_coverage_gap_methods.R"))
source(file.path(script_dir, "03b_advanced_methods.R"))  # SC+cov, Augmented SC, DiD (CS), baseline-window band (after 03)
source(file.path(script_dir, "04_yll_monte_carlo.R"))
source(file.path(script_dir, "05_figures_and_tables.R"))
source(file.path(script_dir, "05b_fig4_triangulation.R")) # Nature-style triangulation Fig 4 (new headline; supersedes make_fig4_method_comparison)
source(file.path(script_dir, "05c_figures_v2.R"))        # redesigned main Figs 1-4 (v2); built in STEP 8
source(file.path(script_dir, "06_sensitivity_start.R"))  # onset-timing sensitivity (functions only; opt-in run below)
# ---- Supplementary robustness/validation modules (additive; PRIMARY spec
#      = manual-onset shift-0 stays the headline everywhere). Each self-skips if
#      an optional package (did/HonestDiD/fixest/ggplot2) is missing. -----------
source(file.path(script_dir, "07_robustness_honest_did.R")) # HonestDiD PT sensitivity on CS event study
source(file.path(script_dir, "08_yll_per_birth.R"))         # per-1,000-births intensity normalisation
source(file.path(script_dir, "09_structured_sensitivity.R"))# tornado over scenario axes
source(file.path(script_dir, "10_external_validation.R"))   # implied measles/pertussis deaths vs benchmark
source(file.path(script_dir, "11_placebo_in_time.R"))       # in-time placebo falsification (all methods)

# ---- Donor-pool scope metadata (F1) --------------------------------------
# ISO3 -> WB region + GDP-derived income tier, consulted by build_donor_pool /
# the AugSC & DiD donor sets ONLY when CVD_DONOR_SCOPE is region/income/
# region_income. Under the default "global" scope nothing changes. Setting
# CVD_DONOR_SCOPE="region_income" (or "no_neighbors") runs the comparable-donor
# / SUTVA-screened sensitivity the reviewers asked for, with no other edits.
.donor_meta <- tryCatch(build_donor_meta(coverage_long),
                        error = function(e) { message("  .donor_meta skipped: ",
                                                      conditionMessage(e)); NULL })
message("  Donor scope: ", DONOR_SCOPE,
        if (!identical(DONOR_SCOPE, "global"))
          "  (neighbours excluded; region/income matched per CVD_DONOR_SCOPE)" else
            "  (worldwide non-conflict donors)")

# ---- Method selection -----------------------------------------------------
# Which counterfactual estimators to run for the gap step. ALL_METHODS (from
# 03b) is the full triangulation set: Baseline, ITS, Synthetic Control,
# SC + covariates, Augmented SC, DiD (CS). Trim here to skip an estimator (e.g.
# drop "Augmented SC" if {augsynth} is not installed) without touching 03b.
GAP_METHODS <- ALL_METHODS
# Which methods to carry into the (heavy) YLL Monte Carlo. Default NULL = every
# method present in gap_df. Set to a subset to keep the MC to a headline panel.
YLL_METHODS <- NULL
# Onset-timing sensitivity (06): re-runs the gap -> YLL chain with conflict_start
# shifted -2/-1/0/+1/+2 yr. ON by default now (the onset year is a researcher
# degree of freedom and reviewers expect this robustness check). It re-runs the
# heaviest part of the pipeline once per shift; set FALSE to skip, or reduce
# START_SHIFT_NSIM. Produces Table_S2 + Figure_S_onset_sensitivity.
RUN_START_SENSITIVITY <- FALSE
START_SHIFT_NSIM      <- 2000L   # MC draws for the sensitivity (point-estimate stable; CRN-matched across shifts)

# ---- STEP 1: Coverage gap estimation -------------------------------------
message("\n>>> STEP 1: Coverage gap estimation (triangulation) <<<\n")
# estimate_all_gaps_plus runs every estimator in GAP_METHODS and stacks them
# under a common gap contract. The CS event-study aggregation rides along as an
# attribute ("event_study") for Figure 4 panel A, so keep gap_df in memory.
gap_df <- estimate_all_gaps_plus(coverage_long, conflict_info, methods = GAP_METHODS)
saveRDS(gap_df, file.path(OUT_DIR, "coverage_gaps.rds"))
utils::write.csv(gap_df, file.path(OUT_DIR, "coverage_gaps.csv"), row.names = FALSE)
message("  Methods present in gap_df: ",
        paste(unique(gap_df$method), collapse = ", "))
message("  Rows: ", nrow(gap_df))

# Reconcile requested vs produced. The advanced estimators self-skip (returning
# no rows) when their optional package is absent, so a method can silently drop
# out of the whole downstream (YLL, figures). Make that loud and actionable.
.missing_methods <- setdiff(GAP_METHODS, unique(gap_df$method))
if (length(.missing_methods) > 0) {
  .pkg_hint <- c(
    "Augmented SC" = "install 'augsynth' (GitHub: remotes::install_github('ebenmichael/augsynth'))",
    "DiD (CS)"     = "install 'did' (CRAN: install.packages('did'))",
    "SC + covariates" = "needs 'tidysynth' AND econ_panel from 02b (check GDP/HDI files loaded)",
    "Synthetic Control" = "install 'tidysynth'")
  warning("Requested methods missing from gap_df (excluded from YLL/figures): ",
          paste(.missing_methods, collapse = ", "), call. = FALSE)
  for (m in .missing_methods)
    message("    - ", m, ": ",
            if (!is.na(.pkg_hint[m])) .pkg_hint[m] else "see gap-step log above for the skip/error reason")
}

# ---- STEP 2: YLL Monte Carlo ---------------------------------------------
message("\n>>> STEP 2: YLL Monte Carlo <<<\n")
# run_full_yll now returns a list: $summary (yearly + country-period rows) and
# $draws (per-country period-total draw vectors for draw-level aggregation).
yll_raw <- run_full_yll(gap_df, covariates_panel, conflict_info, n_sim = N_SIM,
                        methods = YLL_METHODS)
saveRDS(yll_raw$summary, file.path(OUT_DIR, "yll_yearly_all.rds"))
message("  Summary rows: ", nrow(yll_raw$summary),
        " | draw rows: ", nrow(yll_raw$draws))

# Separate yearly trend rows from country-period totals
yll_yearly_only <- yearly_yll(yll_raw$summary)
yll_country     <- aggregate_country_yll(yll_raw$summary)
# Disease/global aggregation is now draw-level (recommendation #4), consuming
# the per-country draw vectors rather than recombining summary CIs.
yll_disease     <- aggregate_disease_yll(yll_raw$draws)

saveRDS(yll_yearly_only, file.path(OUT_DIR, "yll_yearly.rds"))
saveRDS(yll_country,     file.path(OUT_DIR, "yll_country.rds"))
saveRDS(yll_disease,     file.path(OUT_DIR, "yll_disease.rds"))
utils::write.csv(yll_yearly_only, file.path(OUT_DIR, "yll_yearly.csv"), row.names = FALSE)
utils::write.csv(yll_country,     file.path(OUT_DIR, "yll_country.csv"), row.names = FALSE)
utils::write.csv(yll_disease,     file.path(OUT_DIR, "yll_disease.csv"), row.names = FALSE)
message("  yll_country rows: ", nrow(yll_country))
message("  yll_disease rows: ", nrow(yll_disease))

# ---- STEP 3: Figures ------------------------------------------------------
message("\n>>> STEP 3: Figures <<<\n")
# The redesigned main figures (Figs 1-4) are built in STEP 8 by build_figures_v2()
# once every input exists (tornado_res, S4/S2c/S6c tables). The ORIGINAL figure
# calls below are retained as opt-in legacy: set CVD_LEGACY_FIGS=TRUE to emit
# them (distinct filenames, no collision with the v2 outputs).
LEGACY_FIGS <- isTRUE(as.logical(Sys.getenv("CVD_LEGACY_FIGS", unset = "FALSE")))
if (LEGACY_FIGS) {
  fig1 <- make_fig1(coverage_long, gap_df, conflict_info)
  save_fig(fig1, "Figure_1_coverage_decline_legacy", width = 11, height = 9)
  
  fig2 <- make_fig2(yll_country, lambda_focus = 0, framing_focus = FRAMING_HEADLINE)
  save_fig(fig2, "Figure_2_yll_forest_legacy", width = 11, height = 6.5)
}

# Catch-up sensitivity grid is not part of the redesigned Figs 1-4; keep it as a
# supplementary figure (always emitted).
fig3 <- make_fig3(yll_disease)
save_fig(fig3, "Figure_S_catchup_sensitivity", width = 11, height = 6.5)

# ---- Figure 4 (legacy): method triangulation -----------------------------
# Superseded by the v2 triangulation (STEP 8, fig3_v2). Opt-in via CVD_LEGACY_FIGS.
if (LEGACY_FIGS) {
  message("  Building legacy Figure 4 (method comparison)...")
  global_by_method <- aggregate_global_by_method_draws(
    yll_raw$draws, methods = NULL, framing_focus = FRAMING_HEADLINE)
  baseline_windows <- NULL
  sa_es <- tryCatch(
    did_event_study_sa(coverage_long, conflict_info, vaccine = "MCV1"),
    error = function(e) { message("    SA event study skipped: ", conditionMessage(e)); NULL })
  fig4 <- tryCatch(
    make_fig4_method_comparison(gap_df, global_by_method, baseline_windows,
                                sa_df = sa_es, vaccine_focus = "MCV1"),
    error = function(e) { message("    Figure 4 skipped: ", conditionMessage(e)); NULL })
  if (!is.null(fig4)) save_fig(fig4, "Figure_4_method_comparison_legacy", width = 11, height = 9)
}

# ---- STEP 4: Tables -------------------------------------------------------
message("\n>>> STEP 4: Tables <<<\n")
# F6: flag ITS cells floored to 0 by a rising pre-trend (cf >= observed), so they
# print "0\u2020" and aren't read as a true null. Footnote in the manuscript:
#   "\u2020 ITS counterfactual met or exceeded observed coverage (rising
#    pre-conflict trend); attributable YLL floored at 0 and not comparable to the
#    donor-comparison estimators - see Methods."
its_fl <- tryCatch(its_floored_flags(gap_df, vaccine_disease_links),
                   error = function(e) { message("  ITS floored-flag skipped: ",
                                                 conditionMessage(e)); NULL })
save_tab(make_table1(yll_country, its_floored = its_fl), "Table_1_country_disease_yll.csv")
save_tab(make_table2(yll_raw$draws), "Table_2_method_comparison.csv")
save_tab(make_table_s1(),           "Table_S1_parameters.csv")

# F5: headline totals. The SC anchor covers 7 countries only (no SC cell for
# Pakistan/Somalia/Sri Lanka); report the labelled SC complete-case total AND a
# 10-country composite (per-country best-available estimator) so no country is
# silently zero. Both with draw-based 95% CIs.
headline <- tryCatch(assemble_headline_totals(yll_raw$draws),
                     error = function(e) { message("  headline totals skipped: ",
                                                   conditionMessage(e)); NULL })
if (!is.null(headline)) {
  ht <- tibble::tibble(
    estimate        = c("SC (7-country complete-case)",
                        "Composite (10-country, best-available)"),
    n_countries     = c(headline$n_countries_sc, headline$n_countries_comp),
    yll_undisc_mean = c(headline$sc_7country["mean"], headline$composite_10country["mean"]),
    yll_undisc_lo   = c(headline$sc_7country["lo"],   headline$composite_10country["lo"]),
    yll_undisc_hi   = c(headline$sc_7country["hi"],   headline$composite_10country["hi"]))
  save_tab(ht, "Table_headline_totals.csv")
  save_tab(headline$anchor_mix, "Table_headline_anchor_mix.csv")
  message(sprintf("  Headline (undisc): SC 7-country = %s; composite 10-country = %s.",
                  formatC(round(headline$sc_7country["mean"]), big.mark = ",", format = "d"),
                  formatC(round(headline$composite_10country["mean"]), big.mark = ",", format = "d")))
}

# ---- STEP 5 (optional): onset-timing sensitivity -------------------------
if (isTRUE(RUN_START_SENSITIVITY)) {
  message("\n>>> STEP 5: Conflict-onset sensitivity (\u00b11 yr) <<<\n")
  onset_sens <- run_start_sensitivity(
    coverage_long, covariates_panel, conflict_info,
    shifts = c(-2L, -1L, 0L, 1L, 2L), methods = GAP_METHODS,
    n_sim = START_SHIFT_NSIM, seed = 42)
  saveRDS(onset_sens, file.path(OUT_DIR, "onset_sensitivity.rds"))
  save_tab(onset_sens$by_method,  "Table_S2_onset_sensitivity_by_method.csv")
  save_tab(onset_sens$by_country, "Table_S2b_onset_sensitivity_by_country.csv")
  # Composition-STABLE onset table (fixed cell set per method across shifts). This
  # is the version to read: the by_method table above can change its complete-case
  # cell set across shifts (cells_match_ref=FALSE), mixing composition with the
  # onset effect; by_method_fixed holds the cell set constant within each method.
  if (!is.null(onset_sens$by_method_fixed) && nrow(onset_sens$by_method_fixed) > 0)
    save_tab(onset_sens$by_method_fixed, "Table_S2c_onset_sensitivity_fixed_cells.csv")
  # Model-agnostic post-onset event-time coverage-gap profile (uses the PRIMARY
  # manual-onset gap_df; tau = year - conflict_start).
  et_profile <- tryCatch(event_time_gap_profile(gap_df, conflict_info),
                         error = function(e) { message("  event-time profile skipped: ",
                                                       conditionMessage(e)); NULL })
  if (!is.null(et_profile)) save_tab(et_profile, "Table_S2d_event_time_gap_profile.csv")
  fig_s <- tryCatch(make_fig_start_sensitivity(onset_sens),
                    error = function(e) { message("  onset fig skipped: ",
                                                  conditionMessage(e)); NULL })
  if (!is.null(fig_s)) save_fig(fig_s, "Figure_S_onset_sensitivity", width = 9, height = 6)
  # Headline robustness line for the log.
  rng <- onset_sens$by_method
  if (nrow(rng) > 0) {
    rng <- rng[is.finite(rng$pct_change_undisc), ]
    message("  Onset \u00b11yr: undiscounted global YLL moves by ",
            sprintf("%+.1f to %+.1f%%", min(rng$pct_change_undisc),
                    max(rng$pct_change_undisc)), " across methods.")
  }
}

# ---- STEP 6: Supplementary robustness / validation -----------------------
# All additive; the headline (Tables 1-2, Figures 1-4) above is unchanged. Each
# call self-skips if its optional dependency is absent.
message("\n>>> STEP 6: Supplementary robustness / validation <<<\n")

# (S3) HonestDiD parallel-trends sensitivity on the CS event study.
honest_did_res <- tryCatch(
  run_honest_did(gap_df),
  error = function(e) { message("  HonestDiD step skipped: ", conditionMessage(e)); NULL })

# (S4) Per-1,000-births intensity normalisation.
per_birth_res <- tryCatch(
  run_yll_per_birth(yll_country, covariates_panel, conflict_info),
  error = function(e) { message("  per-birth step skipped: ", conditionMessage(e)); NULL })

# (S5) Structured (tornado) sensitivity of the headline total. Anchored on the
# COMPOSITE best-available cells (the same 12.6M the headline reports); the SC
# 7-country anchor is emitted alongside as Table_S5b for the appendix.
tornado_res <- tryCatch(
  run_structured_sensitivity(yll_raw$draws, anchor_mode = "composite",
                             anchor_method = "Synthetic Control",
                             covariates = covariates_panel, conflict_df = conflict_info),
  error = function(e) { message("  tornado step skipped: ", conditionMessage(e)); NULL })

# (S6) External mortality validation (needs an optional benchmark CSV; see module).
validation_res <- tryCatch(
  run_external_validation(yll_country, conflict_info),
  error = function(e) { message("  validation step skipped: ", conditionMessage(e)); NULL })

# (S7) In-time placebo falsification. Heavy (re-runs SC/DiD on pre-onset windows),
# so opt-in. Set RUN_PLACEBO_IN_TIME <- TRUE to produce Table_S7 + figure.
RUN_PLACEBO_IN_TIME <- isTRUE(as.logical(Sys.getenv("CVD_RUN_PLACEBO", unset = "FALSE")))
if (RUN_PLACEBO_IN_TIME) {
  placebo_res <- tryCatch(
    run_placebo_in_time(coverage_long, conflict_info, gap_df_real = gap_df,
                        methods = GAP_METHODS),
    error = function(e) { message("  placebo step skipped: ", conditionMessage(e)); NULL })
} else {
  message("  Placebo-in-time falsification OFF (set CVD_RUN_PLACEBO=TRUE; re-runs SC/DiD).")
}

source(file.path(script_dir, "12_revision_analyses.R"))   # after 09/10
revision_res <- tryCatch(
  run_revision_analyses(yll_raw, gap_df,
                        covariates = covariates_panel,
                        conflict_df = conflict_info),
  error = function(e) { message("  revision analyses skipped: ", conditionMessage(e)); NULL })

# ---- STEP 7 (legacy): triangulation Fig 4 --------------------------------
# Superseded by fig3_v2 in STEP 8. Opt-in via CVD_LEGACY_FIGS.
if (LEGACY_FIGS) {
  message("\n>>> STEP 7 (legacy): Figure 4 (triangulation) <<<\n")
  fig4_tri <- tryCatch(
    build_fig4_triangulation(
      yll_draws   = yll_raw$draws,
      struct_sens = if (exists("tornado_res")) tornado_res else NULL,
      table2_csv  = file.path(TAB_DIR, "Table_2_method_comparison.csv"),
      s5_csv      = file.path(TAB_DIR, "Table_S5_structured_sensitivity.csv"),
      out_pdf     = file.path(FIG_DIR, "Figure_4_triangulation_legacy.pdf")),
    error = function(e) { message("  Fig 4 triangulation skipped: ",
                                  conditionMessage(e)); NULL })
}

# ---- STEP 8: Main figures (redesigned, Nature spec) ----------------------
# Figures 1-4 of the manuscript, built here at the end so every input exists:
#   Fig 1  coverage decline by antigen (coverage_long)
#   Fig 2  coverage gaps -> YLL, by country/method (gap_df, yll_raw$draws)
#   Fig 3  triangulation (yll_raw$draws + tornado_res / Table_S5)
#   Fig 4  per-capita (S4) + onset placebo (S2c, only if STEP 5 ran) + validation (S6c)
# In-memory objects are passed where available; S4/S2c/S6c fall back to the saved
# CSVs. Self-skips per panel if an optional package (patchwork/ggridges) or an
# input is missing.
message("\n>>> STEP 8: Main figures (v2) <<<\n")
figs_v2 <- tryCatch(
  build_figures_v2(
    coverage_long = coverage_long,
    gap_df        = gap_df,
    yll_draws     = yll_raw$draws,
    struct_sens   = if (exists("tornado_res")) tornado_res else NULL,
    drop_metric   = "trough"),
  error = function(e) { message("  STEP 8 figures skipped: ", conditionMessage(e)); NULL })

message("\n=========================================")
message(" Pipeline complete.")
message(" Outputs in: ", OUT_DIR)
message("=========================================\n")