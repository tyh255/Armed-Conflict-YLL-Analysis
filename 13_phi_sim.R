# 13_phi_sim.R | in: PHI_SCHEDULE, unvax_phi(), ve_incidence_draws() (01) | out: Table_phi_verification.csv | Methods: Years of life lost (conversion verified by individual-level simulation; Code availability)
# Infant cohort (N = 2e6) with known hazard h among unprotected infants; fraction c vaccinated at age a = 1 - phi,
# protected thereafter (all-or-nothing or leaky). Simulated all-infant incidence I plays GBD's role; each basis
# must recover h. PASS: timing-adjusted form within 1% of h in all 24 scenarios + helper checks.
#   legacy I | full-year I/(1 - c*VE_I) | timed I/(1 - phi*c*VE_I)

run_phi_sim <- function(N = 2e6, h = 0.08, CFR = 0.03, g = 0.10, seed = SEED_PHI_SIM, tol = 0.01) {
  basis0 <- INC_BASIS
  on.exit(assign("INC_BASIS", basis0, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  sim_incidence <- function(c, ve_case, a, mode) {
    vacc <- stats::runif(N) < c
    pre  <- ifelse(vacc, a, 1); post <- ifelse(vacc, 1 - a, 0)
    rate <- if (mode == "all_or_nothing") {
      responder <- vacc & (stats::runif(N) < ve_case)
      h * pre + h * post * ifelse(responder, 0, 1)
    } else h * pre + h * post * (1 - ve_case)
    mean(stats::rpois(N, rate))
  }
  grid <- expand.grid(vaccine = c("BCG", "DTP3", "MCV1", "MCV1_M12"), c = c(0.3, 0.6, 0.9),
                      mode = c("all_or_nothing", "leaky"), stringsAsFactors = FALSE)
  ve_death <- c(BCG = 0.71, DTP3 = 0.74, MCV1 = 0.83, MCV1_M12 = 0.92)   # numerator VE
  ve_case  <- ve_death; ve_case["BCG"] <- 0.51                           # divisor VE (BCG: disease)
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    v <- grid$vaccine[i]; cc <- grid$c[i]; md <- grid$mode[i]
    phi <- unname(PHI_SCHEDULE[v]); vc <- unname(ve_case[v]); vdth <- unname(ve_death[v])
    I <- sim_incidence(cc, vc, 1 - phi, md)
    hat <- c(legacy = I, full = I / (1 - cc * vc), timed = I / (1 - phi * cc * vc))
    D_unit <- N * g * h * CFR * vdth
    data.frame(vaccine = v, coverage = cc, vaccine_model = md, phi = round(phi, 3), I_sim = I, h_true = h,
               Iu_legacy = hat[["legacy"]], Iu_full = hat[["full"]], Iu_timed = hat[["timed"]],
               ratio_legacy = hat[["legacy"]] / h, ratio_full = hat[["full"]] / h, ratio_timed = hat[["timed"]] / h,
               D_unit_true = D_unit, D_infant_true = D_unit * phi,
               D_legacy = N * g * hat[["legacy"]] * CFR * vdth, D_full = N * g * hat[["full"]] * CFR * vdth,
               D_timed = N * g * hat[["timed"]] * CFR * vdth)
  }))
  # pipeline helpers
  assign("INC_BASIS", "unvaccinated_timed", envir = .GlobalEnv)
  phi_ok <- isTRUE(all.equal(unvax_phi("BCG", "Nigeria"), 1)) &&
    isTRUE(all.equal(unvax_phi("DTP3", "Nigeria"), 1 - 16/52)) &&
    isTRUE(all.equal(unvax_phi("MCV1", "Nigeria"), 1 - 41/52)) && isTRUE(all.equal(unvax_phi("MCV1", "Syria"), 0))
  assign("INC_BASIS", "unvaccinated", envir = .GlobalEnv)
  phi_ok <- phi_ok && unvax_phi("MCV1", "Syria") == 1 && unvax_phi("DTP3", "Yemen") == 1
  assign("INC_BASIS", "unvaccinated_timed", envir = .GlobalEnv)
  vd <- get_shared_draws(draw_shared_params(10000), "VE_BCG"); vi <- ve_incidence_draws(vd, "BCG")
  bcg_ok <- abs(mean(vi) - 0.51) < 0.01 && abs(mean(vd) - 0.71) < 0.01 &&
    stats::cor(vd, vi, method = "spearman") > 0.9999 && identical(ve_incidence_draws(vd, "DTP3"), vd)
  timed_ok <- all(abs(res$ratio_timed - 1) < tol)
  res$timed_within_1pct <- abs(res$ratio_timed - 1) < tol
  res$helpers_ok <- phi_ok && bcg_ok
  save_tab(res, "Table_phi_verification.csv")
  message(sprintf("  phi check: timed within %.0f%% in %d/%d scenarios (max error %.3f%%); helpers %s",
                  100 * tol, sum(res$timed_within_1pct), nrow(res), 100 * max(abs(res$ratio_timed - 1)),
                  phi_ok && bcg_ok))
  if (!(timed_ok && phi_ok && bcg_ok)) stop("incidence-basis verification FAILED")
  invisible(res)
}
