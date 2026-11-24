t0 <- Sys.time()
PASS <- Sys.getenv("CVD_PASS", unset = "main")
stopifnot(PASS %in% c("main", "sens_donor"))
if (PASS == "sens_donor") Sys.setenv(CVD_DONOR_SCOPE = "region_income", CVD_OUT = "out/sens_donor")
MAIN <- PASS == "main"

for (f in c("01_setup.R", "01b_le.R", "02_data.R", "03_gaps.R", "04_yll.R", "05_figs.R", "06_onset.R",
            "07_honestdid.R", "08_perbirth.R", "09_tornado.R", "10_validation.R", "11_placebo.R",
            "12_revision.R", "13_phi_sim.R", "00_select.R"))
  source(file.path("R", f))
message("Pass: ", PASS, " | donor scope: ", DONOR_SCOPE, " | cores: ", N_CORES)

if (MAIN) select_countries()
.donor_meta <- build_donor_meta(coverage_long)

# ---- 1. coverage gaps ----
gap_df <- estimate_all_gaps_plus(coverage_long, conflict_info, methods = ALL_METHODS)
saveRDS(gap_df, file.path(OUT_DIR, "coverage_gaps.rds"))
utils::write.csv(gap_df, file.path(OUT_DIR, "coverage_gaps.csv"), row.names = FALSE)
stopifnot(setequal(unique(gap_df$method), ALL_METHODS))

# ---- 2. YLL Monte Carlo (headline: N_SIM, SEED_MC) ----
yll_raw <- run_full_yll(gap_df, covariates_panel, conflict_info, n_sim = N_SIM)
yll_yearly_only <- yearly_yll(yll_raw$summary)
yll_country     <- aggregate_country_yll(yll_raw$summary)
yll_disease     <- aggregate_disease_yll(yll_raw$draws)
utils::write.csv(yll_yearly_only, file.path(OUT_DIR, "yll_yearly.csv"), row.names = FALSE)
utils::write.csv(yll_country,     file.path(OUT_DIR, "yll_country.csv"), row.names = FALSE)
utils::write.csv(yll_disease,     file.path(OUT_DIR, "yll_disease.csv"), row.names = FALSE)

# ---- 3-4. catch-up figure, Tables 1, 2, S1, headline ----
save_fig(make_fig_catchup(yll_disease), "Figure_S_catchup_sensitivity", width = 11, height = 6.5)
save_tab(make_table1(yll_country, its_floored = its_floored_flags(gap_df, vaccine_disease_links)),
         "Table_1_country_disease_yll.csv")
save_tab(make_table2(yll_raw$draws), "Table_2_method_comparison.csv")
save_tab(make_table_s1(), "Table_S1_parameters.csv")
headline <- assemble_headline_totals(yll_raw$draws)
save_tab(tibble::tibble(
  estimate = c("SC (7-country complete-case)", "Composite (10-country, best-available)"),
  n_countries = c(headline$n_countries_sc, headline$n_countries_comp),
  yll_undisc_mean = c(headline$sc_7country["mean"], headline$composite_10country["mean"]),
  yll_undisc_lo   = c(headline$sc_7country["lo"],   headline$composite_10country["lo"]),
  yll_undisc_hi   = c(headline$sc_7country["hi"],   headline$composite_10country["hi"])),
  "Table_headline_totals.csv")
save_tab(headline$anchor_mix, "Table_headline_anchor_mix.csv")
message(sprintf("Headline (undisc): composite %s (%s-%s); SC 7-country %s",
                format(round(headline$composite_10country["mean"]), big.mark = ","),
                format(round(headline$composite_10country["lo"]), big.mark = ","),
                format(round(headline$composite_10country["hi"]), big.mark = ","),
                format(round(headline$sc_7country["mean"]), big.mark = ",")))

# ---- 5. (i) onset timing, shifts -2..+2 (n_sim 2000) ----
if (MAIN) {
  onset_sens <- run_start_sensitivity(coverage_long, covariates_panel, conflict_info, gap_df0 = gap_df)
  save_tab(onset_sens$by_method,       "Table_S2_onset_sensitivity_by_method.csv")
  save_tab(onset_sens$by_country,      "Table_S2b_onset_sensitivity_by_country.csv")
  save_tab(onset_sens$by_method_fixed, "Table_S2c_onset_sensitivity_fixed_cells.csv")
  save_tab(event_time_gap_profile(gap_df, conflict_info), "Table_S2d_event_time_gap_profile.csv")
  save_fig(make_fig_start_sensitivity(onset_sens), "Figure_S_onset_sensitivity", width = 9, height = 6)
}

# ---- 6. supplementary analyses ----
honest_did_res <- run_honest_did(gap_df)                                         # (ii)
per_birth_res  <- run_yll_per_birth(yll_country, covariates_panel, conflict_info) # (vi)
tornado_res    <- run_structured_sensitivity(yll_raw$draws, anchor_mode = "composite",  # (iii)
                                             covariates = covariates_panel, conflict_df = conflict_info)
validation_res <- run_external_validation(yll_country, conflict_info)            # (v)
if (MAIN) placebo_res <- run_placebo_in_time(coverage_long, conflict_info, gap_df_real = gap_df)  # (iv)
revision_res <- run_revision_analyses(yll_raw, gap_df, covariates = covariates_panel,   # (ix)-(xiv)
                                      conflict_df = conflict_info, reruns = MAIN)
if (MAIN) {                                                                      # (xv) S14
  inc_basis_res <- run_incidence_basis_sensitivity(gap_df, yll_raw$draws, covariates = covariates_panel,
                                                   conflict_df = conflict_info, tornado_tab = tornado_res)
  tornado_res <- inc_basis_res$tornado
}

# ---- 8. main figures + phi verification ----
if (MAIN) {
  build_figures_v2(coverage_long = coverage_long, gap_df = gap_df, yll_draws = yll_raw$draws,
                   struct_sens = tornado_res, drop_metric = "trough")
  run_phi_sim()
}

# ---- provenance ----
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
utils::write.csv(data.frame(file = MANIFEST$file,
                            sha256 = vapply(file.path(ROOT, "data", MANIFEST$file),
                                            function(f) digest::digest(file = f, algo = "sha256"), "")),
                 file.path(OUT_DIR, "input_sha256.csv"), row.names = FALSE)
rt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
writeLines(sprintf("pass=%s cores=%d runtime_min=%.1f finished=%s", PASS, N_CORES, rt, format(Sys.time())),
           file.path(OUT_DIR, "runtime.txt"))
message(sprintf("Pass %s done in %.1f min -> %s", PASS, rt, OUT_DIR))

# ---- (viii) donor-scope pass: separate process, same code, region_income donors ----
if (MAIN) {
  rc <- system2(file.path(R.home("bin"), "Rscript"), file.path("R", "00_run.R"),
                env = c("CVD_PASS=sens_donor", paste0("CVD_CORES=", N_CORES)))
  if (rc != 0) stop("donor-scope pass failed (exit ", rc, ")")
}
