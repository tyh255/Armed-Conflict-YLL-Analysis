# ============================================================================
# 07_ROBUSTNESS_HONEST_DID.R   (supplementary; PRIMARY spec unchanged)
# ----------------------------------------------------------------------------
# Formal parallel-trends sensitivity for the Callaway & Sant'Anna (2021) DiD,
# using Rambachan & Roth (2023, RevEcon Studies) "honest" bounds. This does NOT
# change any headline estimate; it answers the question a Nature referee will
# ask of the DiD arm: "your pre-onset coefficients look ~flat, but eyeballing a
# pre-trend is underpowered and conditioning on it induces pre-test bias -- how
# large a violation of parallel trends would it take to overturn the conflict
# effect?" The answer is the BREAKDOWN VALUE (the smallest restriction parameter
# at which the robust CI first includes 0).
#
# Two restriction families are reported (Rambachan & Roth 2023):
#   * relative_magnitude (Mbar): post-onset PT violation is at most Mbar x the
#     MAX pre-onset violation. Mbar=1 means "the post-period deviation from
#     parallel trends is no larger than the worst pre-period deviation."
#   * smoothness (M): bounds the curvature (second difference) of the event-study
#     path; M=0 is exact linear extrapolation of the pre-trend. UNITS ARE THE
#     OUTCOME'S: coverage is in percentage points, so M is in pp of second
#     difference per period. The default grid is therefore on a pp scale
#     (0,0.5,1,1.5,2), NOT the 0.01-0.05 used previously (which is ~100x too
#     tight for pp coverage and made every robust CI collapse onto exact PT).
#     Sanity-check the top of the grid against the observed pre-trend curvature.
#
# This matters most because the whole comparison-group COVID argument is the
# parallel-trends assumption applied to the 2020-21 shock, and the onset
# sensitivity (06) shows pre-onset coverage degradation in several countries --
# both make the PT assumption worth bounding rather than asserting.
#
# IMPLEMENTATION NOTE: the bridge from a {did} dynamic aggregation (AGGTEobj) to
# {HonestDiD} is Pedro H. C. Sant'Anna's honest_did.AGGTEobj (github.com/
# pedrohcgs/CS_RR), not yet merged into either package, so it is vendored below
# verbatim-in-spirit with attribution. It REQUIRES base_period="universal" in
# att_gt (the pipeline already uses that, 03b) and the full att_gt object, which
# 03b attaches as attr(gap_df, "cs_attgt"). Roth's guidance: restrict the event
# study to e in [-5,5] because distant pre-period coefficients are noisy. That
# restriction is now applied INSIDE the bridge (on vectors already aligned to a
# single full aggregation) rather than via aggte(min_e/max_e), so the influence
# function and point estimates can never get trimmed out of step across {did}
# versions, and non-estimable (NA) event-times are dropped contiguously instead
# of poisoning the vcov.
#
# Source AFTER 03b (needs the cs_attgt attribute) and 05 (save_fig/save_tab,
# nm_palette/theme). Optional packages: did, HonestDiD, ggplot2. Self-skips if
# any is absent.
# ============================================================================

# ---- vendored bridge: honest_did for {did} AGGTEobj (Sant'Anna, CS_RR) ------
honest_did <- function(...) UseMethod("honest_did")

honest_did.AGGTEobj <- function(es, e = 0,
                                type = c("smoothness", "relative_magnitude"),
                                min_e = -Inf, max_e = Inf,
                                avg_window = NULL,
                                gridPoints = 100, ...) {
  type <- match.arg(type)
  if (es$type != "dynamic")
    stop("honest_did: need a dynamic aggregation (aggte(..., type='dynamic')).")
  if (is.null(es$DIDparams$base_period) || es$DIDparams$base_period != "universal")
    stop("honest_did: att_gt must use base_period='universal'.")
  
  # Influence function of the dynamic event-study estimates -> full vcov, built
  # from the SAME (untrimmed) aggregation so its columns align 1:1 with es$egt
  # and es$att.egt. (Trimming is done below, on these aligned vectors.)
  es_inf_func <- es$inf.function$dynamic.inf.func.e
  if (is.null(es_inf_func))
    stop("honest_did: es$inf.function$dynamic.inf.func.e is NULL; aggte did not ",
         "return an event-study influence function (check the {did} version).")
  n <- nrow(es_inf_func)
  V <- t(es_inf_func) %*% es_inf_func / n / n
  
  egt  <- es$egt
  beta <- es$att.egt
  if (length(beta) != ncol(V))
    stop("honest_did: att.egt length (", length(beta), ") != vcov dim (",
         ncol(V), "); influence function and point estimates are misaligned.")
  
  # ---- assemble the event-study window, dimension-safe -------------------
  # HonestDiD's convention: betahat = (pre-coefs ascending, post-coefs ascending)
  # with the universal-base reference (e = -1, att == 0, zero-variance row/col)
  # OMITTED; pre block is e <= -2, post block is e >= 0. Retaining e = -1 makes
  # sigma singular. We grow OUTWARD from the reference on each side, keeping only
  # estimable, contiguous integer event-times within [min_e, max_e]; this honours
  # Roth's [-5,5] guidance and silently survives any non-estimable interior cell
  # (which would otherwise feed NA into HonestDiD).
  estimable <- is.finite(beta) & is.finite(diag(V))
  est_at <- function(et) {
    idx <- match(et, egt)
    length(idx) == 1L && !is.na(idx) && estimable[idx]
  }
  post_e <- integer(0); k <- 0L
  while (k <= max_e && est_at(k)) { post_e <- c(post_e, k); k <- k + 1L }
  pre_e <- integer(0); k <- -2L                       # -1 is the omitted reference
  while (k >= min_e && est_at(k)) { pre_e <- c(pre_e, k); k <- k - 1L }
  pre_e <- rev(pre_e)                                  # ascending
  if (length(pre_e) < 1L || length(post_e) < 1L)
    stop("honest_did: need >=1 estimable pre (e<=-2) and >=1 post (e>=0) ",
         "coefficient contiguous with the reference (npre=", length(pre_e),
         ", npost=", length(post_e), "). Widen min_e/max_e or inspect the ",
         "event study for non-estimable cells near impact.")
  
  sel_e   <- c(pre_e, post_e)
  sel_idx <- match(sel_e, egt)
  beta_k  <- beta[sel_idx]
  V_k     <- V[sel_idx, sel_idx, drop = FALSE]
  npre    <- length(pre_e)
  npost   <- length(post_e)
  
  # ---- target l_vec: single horizon e, OR an average over an event-time window
  # avg_window = c(lo, hi). The averaged target is the estimand the BURDEN uses
  # (the gap accumulates over event-time; the on-impact coefficient understates
  # it for the lagged MCV1/DTP3 arms). l_vec puts equal weight on the estimable
  # post-periods inside the window (clipped to the estimable post block), which is
  # exactly HonestDiD's mechanism for a weighted post-period effect.
  if (!is.null(avg_window)) {
    lv <- .hd_avg_lvec(post_e, avg_window)        # equal-weight over window
    baseVec1 <- lv$l_vec
    target_e <- lv$used_e
  } else {
    if (max(post_e) < e)
      stop("honest_did: requested on-impact horizon e=", e, " exceeds the ",
           "estimable post-window (max e=", max(post_e), ").")
    # post_e is contiguous ascending from 0, so the position of event-time e is
    # match(e, post_e) -- robust even if post_e does not start at 0.
    pos_e <- match(e, post_e)
    if (is.na(pos_e))
      stop("honest_did: on-impact horizon e=", e, " is not an estimable post-period.")
    baseVec1 <- HonestDiD::basisVector(index = pos_e, size = npost)
    target_e <- e
  }
  orig_ci  <- HonestDiD::constructOriginalCS(
    betahat = beta_k, sigma = V_k, numPrePeriods = npre,
    numPostPeriods = npost, l_vec = baseVec1)
  
  robust_ci <- if (type == "relative_magnitude")
    HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat = beta_k, sigma = V_k, numPrePeriods = npre,
      numPostPeriods = npost, l_vec = baseVec1, gridPoints = gridPoints, ...)
  else
    HonestDiD::createSensitivityResults(
      betahat = beta_k, sigma = V_k, numPrePeriods = npre,
      numPostPeriods = npost, l_vec = baseVec1, ...)
  
  list(robust_ci = robust_ci, orig_ci = orig_ci, type = type,
       npre = npre, npost = npost, egt = sel_e, target_e = target_e)
}

# Equal-weight averaging l_vec over the estimable post-periods inside an event-
# time window. `post_e` is the ascending vector of estimable post event-times
# (length npost); the window is clipped to it. Pure base R (unit-tested).
.hd_avg_lvec <- function(post_e, window) {
  if (length(window) != 2L || any(!is.finite(window)))
    stop(".hd_avg_lvec: window must be c(lo, hi) finite.")
  lo <- min(window); hi <- max(window)
  sel <- which(post_e >= lo & post_e <= hi)
  if (length(sel) == 0L)
    stop(".hd_avg_lvec: no estimable post-period in [", lo, ",", hi, "] ",
         "(estimable post-times: ", paste(post_e, collapse = ","), ").")
  l <- numeric(length(post_e))
  l[sel] <- 1 / length(sel)
  list(l_vec = l, used_e = post_e[sel])
}

# ---- breakdown value -------------------------------------------------------
# Smallest restriction parameter (Mbar or M) at which the robust CI first
# includes 0 (i.e. the conflict effect is no longer distinguishable from 0).
# Returns Inf if every tested restriction still excludes 0 (very robust),
# the smallest grid value if even the tightest robust CI already includes 0
# (not identified at any tested restriction). ALWAYS returns a numeric scalar:
# returning c(value=<num>, par=<chr>) silently coerced the value to character,
# which then broke sprintf("%.2f", .) downstream and aborted the table write.
.breakdown_value <- function(robust_ci) {
  rc <- robust_ci
  par_col <- if ("Mbar" %in% names(rc)) "Mbar" else "M"
  rc <- rc[order(rc[[par_col]]), ]
  includes0 <- rc$lb <= 0 & rc$ub >= 0
  if (all(!includes0)) return(Inf)
  if (all(includes0))  return(min(rc[[par_col]]))
  rc[[par_col]][which(includes0)[1]]
}

# Local table writer: prefer the pipeline's save_tab (05); fall back to a plain
# CSV write so a standalone run of this module still persists Table_S3/S3b.
.hd_save_tab <- function(df, file, out_dir) {
  if (exists("save_tab", mode = "function")) return(save_tab(df, file))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(out_dir, file)
  utils::write.csv(df, path, row.names = FALSE)
  message("  [HonestDiD] wrote ", path, " (save_tab() not found; used write.csv).")
}

# ---- per-(vaccine,target) runner -------------------------------------------
# Runs both restriction families for ONE target (a single horizon e, or an
# event-time average window) and returns tagged robust-CI rows. Factored out so
# the on-impact and post-period-average horizons share identical error handling.
.hd_run_target <- function(dyn, vacc, target, min_e, max_e, Mbar_grid, M_grid) {
  hz <- target$label
  errs <- character(0)
  mk_args <- function(type, grid_arg, grid) {
    a <- list(dyn, type = type, min_e = min_e, max_e = max_e)
    if (!is.null(target$avg_window)) a$avg_window <- target$avg_window else a$e <- target$e
    a[[grid_arg]] <- grid
    a
  }
  hd_rm <- tryCatch(
    do.call(honest_did, mk_args("relative_magnitude", "Mbarvec", Mbar_grid)),
    error = function(err) {
      pos <- Mbar_grid[Mbar_grid > 0]
      if (length(pos) >= 1L && length(pos) < length(Mbar_grid)) {
        tryCatch(do.call(honest_did, mk_args("relative_magnitude", "Mbarvec", pos)),
                 error = function(e2) {
                   m <- paste0("RM (", vacc, "/", hz, "): ", conditionMessage(e2))
                   message("    ", m); errs[[length(errs) + 1L]] <<- m; NULL })
      } else {
        m <- paste0("RM (", vacc, "/", hz, "): ", conditionMessage(err))
        message("    ", m); errs[[length(errs) + 1L]] <<- m; NULL
      }
    })
  hd_sm <- tryCatch(
    do.call(honest_did, mk_args("smoothness", "Mvec", M_grid)),
    error = function(err) {
      m <- paste0("smoothness (", vacc, "/", hz, "): ", conditionMessage(err))
      message("    ", m); errs[[length(errs) + 1L]] <<- m; NULL })
  rm_row <- NULL; sm_row <- NULL
  if (!is.null(hd_rm)) {
    bd <- .breakdown_value(hd_rm$robust_ci)
    rm_row <- hd_rm$robust_ci %>%
      dplyr::mutate(vaccine = vacc, horizon = hz,
                    target_e = paste(hd_rm$target_e, collapse = ","),
                    orig_lb = hd_rm$orig_ci$lb, orig_ub = hd_rm$orig_ci$ub,
                    breakdown_Mbar = bd)
  }
  if (!is.null(hd_sm)) {
    bd <- .breakdown_value(hd_sm$robust_ci)
    sm_row <- hd_sm$robust_ci %>%
      dplyr::mutate(vaccine = vacc, horizon = hz,
                    target_e = paste(hd_sm$target_e, collapse = ","),
                    orig_lb = hd_sm$orig_ci$lb, orig_ub = hd_sm$orig_ci$ub,
                    breakdown_M = bd)
  }
  list(rm = rm_row, sm = sm_row, errs = errs)
}

# ---- driver ----------------------------------------------------------------
# Consumes attr(gap_df, "cs_attgt") (named list of att_gt objects by vaccine).
# For each vaccine: aggregate the FULL dynamic event study, then for EACH target
# horizon run both restriction families over a grid (trimming to [min_e,max_e]
# inside the bridge) and record the original CI + breakdown.
#
# TWO horizons by default (a review fix): the on-impact coefficient (e=0) AND the
# post-period AVERAGE over `avg_window` (default e in [0,4]). The averaged target
# is the headline robustness read, because the YLL burden is driven by the gap
# ACCUMULATED over event time (S2d), not the on-impact year -- bounding PT
# violations only at e=0 understates robustness for the lagged MCV1/DTP3 arms,
# whose on-impact CI already straddles 0. Set avg_window = NULL to recover the
# e=0-only behaviour. Tables gain a `horizon` column.
run_honest_did <- function(gap_df,
                           e = 0, min_e = -5L, max_e = 5L,
                           avg_window = c(0L, 4L),
                           Mbar_grid = c(0, 0.5, 1, 1.5, 2),
                           # smoothness M is in OUTCOME (percentage-point) units;
                           # see header. Was c(0,0.01,..,0.05) -> ~100x too tight.
                           M_grid    = c(0, 0.5, 1, 1.5, 2),
                           out_dir = if (exists("TAB_DIR")) TAB_DIR
                           else if (exists("OUT_DIR")) OUT_DIR else ".") {
  if (!requireNamespace("did", quietly = TRUE) ||
      !requireNamespace("HonestDiD", quietly = TRUE)) {
    message("  HonestDiD: need both {did} and {HonestDiD} -> skipping ",
            "(install.packages(c('did','HonestDiD'))).")
    return(invisible(NULL))
  }
  attgt <- attr(gap_df, "cs_attgt")
  if (is.null(attgt) || length(attgt) == 0) {
    message("  HonestDiD: attr(gap_df,'cs_attgt') is empty -> DiD not run, or ",
            "an older gap_df. Re-run estimate_all_gaps_plus with DiD enabled.")
    return(invisible(NULL))
  }
  
  # Build the target list: on-impact (e) plus the post-period average (avg_window).
  targets <- list(list(label = sprintf("on_impact_e%d", as.integer(e)),
                       e = as.integer(e), avg_window = NULL))
  if (!is.null(avg_window)) {
    wlab <- sprintf("avg_e%d_e%d", as.integer(min(avg_window)), as.integer(max(avg_window)))
    targets <- c(targets, list(list(label = wlab, e = NULL,
                                    avg_window = as.integer(avg_window))))
  }
  
  rows_rm <- list(); rows_sm <- list(); es_tab <- list(); errs <- character(0)
  for (vacc in names(attgt)) {
    att <- attgt[[vacc]]
    # FULL dynamic aggregation (no min_e/max_e here): keeps att.egt, se.egt and
    # the influence function mutually aligned. The [-5,5] restriction is applied
    # inside honest_did().
    dyn <- tryCatch(
      did::aggte(att, type = "dynamic", na.rm = TRUE),
      error = function(err) { message("    aggte failed (", vacc, "): ",
                                      conditionMessage(err)); NULL })
    if (is.null(dyn)) next
    
    # store the (trimmed) event study for the figure
    keep_es <- dyn$egt >= min_e & dyn$egt <= max_e
    es_tab[[vacc]] <- tibble::tibble(
      vaccine = vacc, event_time = dyn$egt[keep_es],
      att = dyn$att.egt[keep_es], se = dyn$se.egt[keep_es])
    
    for (tg in targets) {
      rt <- .hd_run_target(dyn, vacc, tg, min_e, max_e, Mbar_grid, M_grid)
      key <- paste(vacc, tg$label, sep = "|")
      if (!is.null(rt$rm)) rows_rm[[key]] <- rt$rm
      if (!is.null(rt$sm)) rows_sm[[key]] <- rt$sm
      if (length(rt$errs)) errs <- c(errs, rt$errs)
    }
  }
  
  rm_tab <- dplyr::bind_rows(rows_rm)
  sm_tab <- dplyr::bind_rows(rows_sm)
  es_all <- dplyr::bind_rows(es_tab)
  
  # WRITE TABLES FIRST. The deliverables are persisted before any (non-essential)
  # console diagnostics, so a formatting slip in a diagnostic can never again
  # abort the save. Both writes go through the save_tab/write.csv fallback.
  if (nrow(rm_tab) > 0) .hd_save_tab(rm_tab, "Table_S3_honestdid_relative_magnitude.csv", out_dir)
  if (nrow(sm_tab) > 0) .hd_save_tab(sm_tab, "Table_S3b_honestdid_smoothness.csv",        out_dir)
  
  if (nrow(rm_tab) == 0 && nrow(sm_tab) == 0) {
    message("  [HonestDiD] produced NO rows -> Table_S3/S3b not written. ",
            "Captured errors:")
    if (length(errs) == 0)
      message("    (none captured; every aggte() returned NULL -- check the DiD ",
              "event study has a universal base period and estimable pre/post cells.)")
    else for (m in errs) message("    - ", m)
  }
  
  # SELF-TEST / DIAGNOSTICS (wrapped: must never block the outputs above) ----
  # (1) robust CIs should (weakly) widen as the restriction relaxes;
  # (2) breakdown is a finite grid value or Inf.
  tryCatch({
    if (nrow(rm_tab) > 0) {
      chk <- rm_tab %>% dplyr::group_by(vaccine, horizon) %>%
        dplyr::summarise(widen_ok = all(diff(ub - lb) >= -1e-8), .groups = "drop")
      if (any(!chk$widen_ok))
        message("  [HonestDiD diag] non-monotone RM CI width for: ",
                paste(sprintf("%s/%s", chk$vaccine[!chk$widen_ok], chk$horizon[!chk$widen_ok]),
                      collapse = ", "),
                " (inspect; can happen with very noisy pre-periods).")
      bd_by_v <- rm_tab %>% dplyr::group_by(vaccine, horizon) %>%
        dplyr::summarise(bd = breakdown_Mbar[1], .groups = "drop")
      message("  [HonestDiD] relative-magnitude breakdown Mbar (vaccine/horizon): ",
              paste(sprintf("%s/%s=%.2f", bd_by_v$vaccine, bd_by_v$horizon, bd_by_v$bd),
                    collapse = ", "))
      # The interpretive read: compare on-impact vs averaged horizon per vaccine.
      av <- rm_tab %>% dplyr::filter(grepl("^avg_", horizon)) %>%
        dplyr::group_by(vaccine) %>% dplyr::summarise(bd = breakdown_Mbar[1], .groups = "drop")
      if (nrow(av) > 0)
        message("  [HonestDiD] post-period-average breakdown Mbar (headline robustness): ",
                paste(sprintf("%s=%.2f", av$vaccine, av$bd), collapse = ", "),
                " -- averaged target reflects the accumulated gap, not the on-impact year.")
    }
  }, error = function(err)
    message("  [HonestDiD diag] skipped: ", conditionMessage(err)))
  
  fig <- tryCatch(.fig_honest_did(rm_tab, es_all),
                  error = function(err) { message("  HonestDiD fig skipped: ",
                                                  conditionMessage(err)); NULL })
  if (!is.null(fig) && exists("save_fig"))
    save_fig(fig, "Figure_S_honestdid_sensitivity", width = 11, height = 5.5)
  
  invisible(list(relative_magnitude = rm_tab, smoothness = sm_tab,
                 event_study = es_all))
}

# Robust-CI-vs-restriction plot (one facet per vaccine), with the original
# (PT-assumed) CI as the leftmost interval and a 0 reference line.
.fig_honest_did <- function(rm_tab, es_all) {
  if (is.null(rm_tab) || nrow(rm_tab) == 0)
    stop(".fig_honest_did: empty relative-magnitude table.")
  if (!"horizon" %in% names(rm_tab)) rm_tab$horizon <- "on_impact_e0"
  orig <- rm_tab %>% dplyr::distinct(vaccine, horizon, orig_lb, orig_ub) %>%
    dplyr::mutate(Mbar = -0.25, lb = orig_lb, ub = orig_ub, kind = "Original (PT)")
  rob <- rm_tab %>% dplyr::transmute(vaccine, horizon, Mbar, lb, ub, kind = "Robust (R&R)")
  df <- dplyr::bind_rows(orig, rob)
  ggplot2::ggplot(df, ggplot2::aes(x = Mbar, ymin = lb, ymax = ub, colour = kind)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_linerange(linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = (lb + ub) / 2), size = 1.6) +
    ggplot2::facet_grid(horizon ~ vaccine, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("Original (PT)" = "#444444",
                                            "Robust (R&R)" = "#B5651D"),
                                 name = NULL) +
    ggplot2::labs(
      x = expression(paste("Relative-magnitude restriction  ", bar(M))),
      y = "Conflict effect on coverage (ATT, % pts)",
      title = "HonestDiD sensitivity of the DiD conflict effect to parallel-trends violations",
      subtitle = paste0("Rows: on-impact (e=0) vs post-period average. Robust 95% CI ",
                        "widens as Mbar relaxes; breakdown = smallest Mbar whose CI includes 0."),
      caption = "Rambachan & Roth (2023). Original = CI under exact parallel trends.") +
    { if (exists("theme_nm")) theme_nm() else ggplot2::theme_minimal() }
}