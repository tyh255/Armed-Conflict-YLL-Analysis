# ============================================================================
# 01b_LIFE_EXPECTANCY.R   (new module)
# ----------------------------------------------------------------------------
# Life-expectancy source layer for the YLL loss function (the "L" in
# YLL = deaths * L). SOURCE THIS AFTER 02_load_data.R AND BEFORE
# 04_yll_monte_carlo.R (it depends on harmonise_country() / the ISO3 lookup
# defined in 02, and is consumed by simulate_yll_one() in 04).
#
# WHY THIS MODULE EXISTS
# ----------------------
# The prior pipeline used World Bank NATIONAL period life expectancy as L. For a
# descriptive, GBD-comparable burden estimate that is the wrong input: GBD/WHO
# compute YLL against a single STANDARD REFERENCE life table (the Theoretical
# Minimum Risk Life Table, TMRLT), identical for every country, so that a death
# at a given age is valued equally regardless of where it occurs. Using local LE
# structurally downweights deaths in high-mortality settings - exactly the
# conflict countries in this study - and breaks cross-country comparability.
#
# IMPORTANT - WHICH FILE IS WHICH
# -------------------------------
# Files named IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_28_* are GBD 2019
# COUNTRY-YEAR period life tables (age group <1). They are LOCAL life expectancy
# (IHME-sourced instead of World Bank). They are NOT the reference table and do
# NOT, on their own, implement the reference-table recommendation. They are wired
# in here only as the LOCAL-LE sensitivity (mode "local_gbd").
#
# The HEADLINE reference table (TMRLT, ~88.87 yr at birth) is a SEPARATE GHDx
# download: "Global Burden of Disease Study 2019 (GBD 2019) Reference Life Table".
# Export it to a 2-column CSV (columns: age, ex) at REFERENCE_LT_FILE; otherwise
# the embedded fallback below is used. The infant cohort needs only ex at age 0
# (= 88.8719), which the fallback supplies; the older-age rows in the fallback
# are approximate and should be VERIFIED against the official CSV before being
# relied on for any age-band (under-5) sensitivity.
#
# THREE MODES (set via env var CVD_LIFE_TABLE_MODE, or edit LIFE_TABLE_MODE):
#   "reference" (DEFAULT, headline) - GBD 2019 TMRLT; constant across countries.
#   "local_gbd" (sensitivity)       - IHME GBD 2019 country LE at birth, NO-SHOCK.
#   "local_wb"  (prior behaviour)   - World Bank national LE at birth.
#
# NO-SHOCK vs WITH-SHOCK (for local_gbd): we use NO-SHOCK. Using with-shock LE
# (depressed by war/disaster deaths) to value conflict-attributable YLLs is
# circular - it lets the very conflict being studied shrink the life expectancy
# used to weight its own burden - and biases the headline downward.
# ============================================================================

stopifnot(requireNamespace("dplyr", quietly = TRUE),
          requireNamespace("tibble", quietly = TRUE))

# ---- Mode -----------------------------------------------------------------
LIFE_TABLE_MODE <- Sys.getenv("CVD_LIFE_TABLE_MODE", unset = "reference")
if (!LIFE_TABLE_MODE %in% c("reference", "local_gbd", "local_wb")) {
  stop("LIFE_TABLE_MODE must be one of 'reference', 'local_gbd', 'local_wb'; got '",
       LIFE_TABLE_MODE, "'.")
}

# ---- Reference life table (GBD 2019 TMRLT) --------------------------------
# Optional authoritative override: a 2-column CSV (age, ex). If present it is
# used; otherwise the embedded values below (the official GBD 2019 Reference
# Life Table, confirmed by the user) are used. Either way ex at age 0 = 88.872
# drives the infant-cohort headline; older-age rows support the under-5 / age-
# band sensitivity (recommendation #7).
REFERENCE_LT_FILE <- file.path(YLL_DATA_DIR, "GBD_2019_Reference_Life_Table.csv")

# Official GBD 2019 Reference Life Table (TMRLT), ages 0-95 at standard intervals.
.gbd2019_tmrlt_embedded <- tibble::tribble(
  ~age,  ~ex,
  0,     88.8718951,
  1,     88.00051053,
  5,     84.03008056,
  10,    79.04633476,
  15,    74.0665492,
  20,    69.10756792,
  25,    64.14930031,
  30,    59.1962771,
  35,    54.25261364,
  40,    49.31739311,
  45,    44.43332057,
  50,    39.63473787,
  55,    34.91488095,
  60,    30.25343822,
  65,    25.68089534,
  70,    21.28820012,
  75,    17.10351469,
  80,    13.23872477,
  85,    9.990181244,
  90,    7.617724915,
  95,    5.922359078
)

load_reference_lt <- function() {
  if (file.exists(REFERENCE_LT_FILE)) {
    df <- utils::read.csv(REFERENCE_LT_FILE, stringsAsFactors = FALSE)
    nm <- tolower(trimws(names(df)))
    age_idx <- which(nm %in% c("age", "age_group", "age_start", "age_group_years_start", "x"))
    ex_idx  <- which(nm %in% c("ex", "life_expectancy", "lifeexpectancy",
                               "standard_le", "le", "val", "value"))
    if (length(age_idx) == 0 || length(ex_idx) == 0) {
      stop("REFERENCE_LT_FILE found but could not identify 'age' and 'ex' columns ",
           "in ", REFERENCE_LT_FILE, ". Provide a 2-column CSV with headers ",
           "'age' and 'ex'.")
    }
    out <- tibble::tibble(age = suppressWarnings(as.numeric(df[[age_idx[1]]])),
                          ex  = suppressWarnings(as.numeric(df[[ex_idx[1]]]))) %>%
      dplyr::filter(!is.na(age), !is.na(ex)) %>%
      dplyr::arrange(age)
    if (!0 %in% out$age)
      warning("Reference life table has no age-0 row; infant-cohort L will be ",
              "interpolated/clamped.")
    message("Reference life table: loaded ", nrow(out), " rows from ",
            basename(REFERENCE_LT_FILE),
            " (ex at age 0 = ",
            round(stats::approx(out$age, out$ex, xout = 0, rule = 2)$y, 4), ").")
    return(out)
  }
  message("Reference life table: official file not present at\n  ", REFERENCE_LT_FILE,
          "\n  -> using EMBEDDED GBD 2019 TMRLT (official values; ex at age 0 = 88.872).")
  .gbd2019_tmrlt_embedded
}
REFERENCE_LT <- load_reference_lt()

# Reference ex at an arbitrary age (linear interpolation, clamped to table range).
ref_ex_at_age <- function(age) {
  lt <- REFERENCE_LT
  stats::approx(lt$age, lt$ex, xout = age, rule = 2)$y
}

# ---- Age at death per disease ---------------------------------------------
# The current pipeline is an INFANT (<1) cohort, so age at death is 0 for all
# diseases and the reference loss is ex(0) = 88.87. These are hooks for the
# under-5 / age-band sensitivity (recommendation #7): when incidence, the
# denominator N, and the CFR are jointly switched to the 0-4 band for TB and
# diphtheria, set their age here to the mean age at death within 0-4 (deaths in
# these diseases cluster at the youngest ages, so ~0.5-1 yr is typical) so that L
# stays consistent with the band used for incidence/CFR/N.
AGE_AT_DEATH <- list(
  Tuberculosis = 0,   # -> ~1 if/when TB moves to the 0-4 band
  Measles      = 0,
  Diphtheria   = 0,   # -> ~1 if/when diphtheria moves to the 0-4 band
  Pertussis    = 0,
  Tetanus      = 0
)
age_at_death_for_disease <- function(disease) {
  # Under the under-5 sensitivity (recommendation #7), TB and diphtheria deaths
  # are valued at ~age 1 (deaths cluster at the youngest ages within 0-4) so L
  # stays consistent with the 0-4 incidence/CFR/denominator. Otherwise age 0.
  if (exists("AGE_BAND_MODE") && AGE_BAND_MODE == "under5_tb_diph" &&
      disease %in% c("Tuberculosis", "Diphtheria")) {
    return(1)
  }
  a <- AGE_AT_DEATH[[disease]]
  if (is.null(a)) stop("No age-at-death mapping for disease: ", disease)
  a
}

# ---- Local GBD 2019 country life table (sensitivity, NO-SHOCK) ------------
# Default points at the NO-SHOCK ID_28 file. Override the filename via the env
# var CVD_GBD_LOCAL_LT_FILE if needed. measure_name for no-shock LE in the ID_28
# release is "Life expectancy no-shock with hiv" (codebook measure_id 31); the
# with-shock equivalent is "Life expectancy" (measure_id 26).
GBD_LOCAL_LT_FILE <- Sys.getenv(
  "CVD_GBD_LOCAL_LT_FILE",
  unset = file.path(
    YLL_DATA_DIR,
    "IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_28_NOSHOCK_Y2020M07D31.csv"))
GBD_LOCAL_LE_MEASURE <- "Life expectancy no-shock with hiv"

load_gbd_local_le <- function() {
  if (!file.exists(GBD_LOCAL_LT_FILE)) {
    message("local_gbd LE file not found at\n  ", GBD_LOCAL_LT_FILE,
            "\n  -> 'local_gbd' mode will return NA (years skipped).")
    return(NULL)
  }
  if (!exists("harmonise_country")) {
    stop("harmonise_country() not found; source 02_load_data.R before this module.")
  }
  df <- utils::read.csv(GBD_LOCAL_LT_FILE, stringsAsFactors = FALSE)
  need <- c("location_name", "sex_name", "age_group_id", "measure_name",
            "year_id", "val")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0)
    stop("GBD local LE file missing columns: ", paste(miss, collapse = ", "))
  out <- df %>%
    dplyr::filter(tolower(.data$sex_name) == "both",
                  .data$age_group_id == 28L,
                  .data$measure_name == GBD_LOCAL_LE_MEASURE) %>%
    dplyr::mutate(country      = harmonise_country(.data$location_name),
                  year         = as.integer(.data$year_id),
                  le_local_gbd = as.numeric(.data$val)) %>%
    dplyr::select(country, year, le_local_gbd)
  # Attach ISO3 via the lookup built in 02 (covariate-/coverage-derived).
  if (exists("country_iso_lookup")) {
    out <- out %>% dplyr::left_join(country_iso_lookup, by = "country")
  } else if (exists("iso3_lookup")) {
    out <- out %>% dplyr::left_join(iso3_lookup, by = "country")
  } else {
    warning("No ISO3 lookup available; local_gbd lookups will fail to match.")
    out$ISO3 <- NA_character_
  }
  out <- out %>% dplyr::filter(!is.na(le_local_gbd))
  message("local_gbd LE: loaded ", nrow(out), " country-year rows (no-shock) from ",
          basename(GBD_LOCAL_LT_FILE), ".")
  out
}
GBD_LOCAL_LE <- if (LIFE_TABLE_MODE == "local_gbd") load_gbd_local_le() else NULL

gbd_local_le_lookup <- function(iso3, year) {
  if (is.null(GBD_LOCAL_LE)) return(NA_real_)
  hit <- GBD_LOCAL_LE %>%
    dplyr::filter(.data$ISO3 == iso3, .data$year == !!year)
  if (nrow(hit) == 0) return(NA_real_)
  mean(hit$le_local_gbd, na.rm = TRUE)
}

# ---- Resolver -------------------------------------------------------------
# Returns L (years) for a death attributed to `disease` in `iso3`, `year`.
#   reference : reference ex at the disease's age at death (country-invariant)
#   local_gbd : GBD 2019 country LE at birth (no-shock) for that ISO3-year
#   local_wb  : the World Bank LE passed in from the covariate panel
# Callers in 04 pass g$Life_Expectancy[yi] as local_wb_le; it is ignored unless
# mode == "local_wb".
resolve_life_expectancy <- function(disease, iso3, year, local_wb_le = NA_real_) {
  if (LIFE_TABLE_MODE == "reference") {
    return(ref_ex_at_age(age_at_death_for_disease(disease)))
  }
  if (LIFE_TABLE_MODE == "local_gbd") {
    return(gbd_local_le_lookup(iso3, year))
  }
  local_wb_le  # local_wb
}

message("Life-expectancy layer ready: mode = '", LIFE_TABLE_MODE, "'",
        if (LIFE_TABLE_MODE == "reference")
          paste0(" (TMRLT ex at age 0 = ", round(ref_ex_at_age(0), 4), " yr).")
        else ".")