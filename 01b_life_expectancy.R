

stopifnot(requireNamespace("dplyr", quietly = TRUE),
          requireNamespace("tibble", quietly = TRUE))

# ---- Mode -----------------------------------------------------------------
LIFE_TABLE_MODE <- Sys.getenv("CVD_LIFE_TABLE_MODE", unset = "reference")
if (!LIFE_TABLE_MODE %in% c("reference", "local_gbd", "local_wb")) {
  stop("LIFE_TABLE_MODE must be one of 'reference', 'local_gbd', 'local_wb'; got '",
       LIFE_TABLE_MODE, "'.")
}

# ---- Reference life table (GBD 2019 TMRLT) --------------------------------
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
