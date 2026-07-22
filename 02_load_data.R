# ============================================================================
# 02_LOAD_DATA.R
# ----------------------------------------------------------------------------
# Reads all source data (coverage, IHME GBD incidence, World Bank covariates)
# and builds a long-format coverage panel + a wide covariates panel for use
# downstream.
#
# CHANGED in this version (vs prior): disease incidence now comes from IHME
# GBD (true incidence estimates with UIs) instead of WHO surveillance reported
# cases. WHO surveillance under-reports VPDs by 10-100x in LMIC and even more
# in conflict settings; using it as the incidence input was the dominant
# source of magnitude error in prior outputs.
#
# IHME files expected at <DATA_DIR>/IHME Data/IHME_<Disease>_GBD.csv with
# columns: measure_name, location_name, sex_name, age_name, metric_name,
# year, val, upper, lower.
#
# Pre-requisite: 01_setup_and_parameters.R must already be sourced.
# ============================================================================

# ---- 1. Read raw files ----------------------------------------------------
read_csv_safe <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# Vaccine coverage (WHO/UNICEF WUENIC)
bcg_data    <- read_csv_safe(file.path(COVERAGE_DIR, "BCG_WHO-UNICEF.csv"))
mcv1_data   <- read_csv_safe(file.path(COVERAGE_DIR, "MCV1_WHO-UNICEF.csv"))
dpt_data    <- read_csv_safe(file.path(COVERAGE_DIR, "DPT3_WHO-UNICEF.csv"))

# World Bank demographic covariates (births, population, life expectancy)
births_data <- read_csv_safe(file.path(YLL_DATA_DIR, "WB_Birth Rate.csv"))
pop_data    <- read_csv_safe(file.path(YLL_DATA_DIR, "WB_Population Data.csv"))
life_data   <- read_csv_safe(file.path(YLL_DATA_DIR, "WB_Life Expectancy.csv"))

# IHME GBD incidence (replaces all prior WHO surveillance / WB-TB rates)
IHME_DIR <- file.path(DATA_DIR, "IHME Data")

# ---- 2. Country name harmonisation ----------------------------------------
harmonise_country <- function(x) {
  dplyr::case_when(
    x == "American Samoa"                                                ~ "Samoa",
    x == "C\u00f4te d'Ivoire"                                            ~ "Cote d'Ivoire",
    x == "occupied Palestinian territory, including east Jerusalem"      ~ "State of Palestine",
    x == "T\u00fcrkiye"                                                  ~ "Turkiye",
    x == "United Kingdom of Great Britain and Northern Ireland"          ~ "United Kingdom",
    x == "United States of America"                                      ~ "United States",
    x == "Syrian Arab Republic"                                          ~ "Syria",
    x == "Iran (Islamic Republic of)"                                    ~ "Iran",
    x == "Lao People's Democratic Republic"                              ~ "Laos",
    x == "Viet Nam"                                                      ~ "Vietnam",
    x == "Republic of Korea"                                             ~ "South Korea",
    x == "Democratic People's Republic of Korea"                         ~ "North Korea",
    x == "Russian Federation"                                            ~ "Russia",
    x == "United Republic of Tanzania"                                   ~ "Tanzania",
    x == "Bolivia (Plurinational State of)"                              ~ "Bolivia",
    x == "Venezuela (Bolivarian Republic of)"                            ~ "Venezuela",
    x == "Republic of Moldova"                                           ~ "Moldova",
    x == "The former Yugoslav Republic of Macedonia"                     ~ "North Macedonia",
    x == "Czech Republic"                                                ~ "Czechia",
    x == "Eswatini (Swaziland)"                                          ~ "Eswatini",
    x == "Brunei Darussalam"                                             ~ "Brunei",
    x == "Micronesia (Federated States of)"                              ~ "Micronesia",
    TRUE                                                                 ~ x
  )
}

# ---- 3. Reshape World Bank wide-to-long -----------------------------------
to_long_wb <- function(df, value_name) {
  df %>%
    tidyr::pivot_longer(
      cols = starts_with("X"),
      names_to = "year",
      values_to = value_name
    ) %>%
    dplyr::mutate(year = as.integer(stringr::str_extract(year, "\\d{4}"))) %>%
    dplyr::rename(ISO3 = Country.Code) %>%
    dplyr::select(ISO3, year, dplyr::all_of(value_name))
}

WB_births_long <- to_long_wb(births_data, "Birth_Rate")
WB_pop_long    <- to_long_wb(pop_data, "Population") %>%
  dplyr::mutate(Population_per_1000 = Population / 1000)
WB_life_long   <- to_long_wb(life_data, "Life_Expectancy")

# ---- 4. Load IHME GBD incidence ------------------------------------------
# Filter: measure="Incidence", metric="Rate", sex="Both", age="<1 year".
# Rate is per 100,000 children under 1 -> divide by 1e5 for per-person.
# Compute per-row lognormal sigma from (lower, upper) for Monte Carlo.
#
# Age band note: prior pipeline convention is <1 year for all diseases
# (matches Birth_Cohort multiplier in the YLL formula). This is appropriate
# for pertussis and tetanus (mortality concentrated in infants) but
# under-counts measles/diphtheria/TB mortality from older children.
# Reviewers may want a sensitivity using under-5 age band; that requires
# also switching Birth_Cohort to under-5 population.
read_ihme_incidence <- function(filename, out_prefix,
                                age_band = "<1 year") {
  path <- file.path(IHME_DIR, filename)
  if (!file.exists(path)) stop("IHME file not found: ", path)
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  
  required_cols <- c("measure_name", "metric_name", "sex_name", "age_name",
                     "location_name", "year", "val", "upper", "lower")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("IHME file ", filename, " missing columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  out <- df %>%
    dplyr::filter(
      measure_name == "Incidence",
      metric_name  == "Rate",
      sex_name     == "Both",
      age_name     == age_band
    ) %>%
    dplyr::mutate(
      country = harmonise_country(location_name),
      year    = as.integer(year),
      # GBD rate is per 100,000; convert to per-person
      inc_rate  = val   / 1e5,
      inc_upper = upper / 1e5,
      inc_lower = lower / 1e5
    ) %>%
    dplyr::mutate(
      # Per-row lognormal sigma from 95% UI.
      # If lower is zero (some IHME rows for rare diseases), fall back to
      # NA and let the MC engine use the global INC_LOGNORM_SIGMA default.
      inc_sigma = dplyr::if_else(
        inc_upper > 0 & inc_lower > 0 & inc_upper > inc_lower,
        (log(inc_upper) - log(inc_lower)) / 3.92,
        NA_real_
      )
    ) %>%
    dplyr::select(country, year, inc_rate, inc_upper, inc_lower, inc_sigma) %>%
    # Dedupe (location/age/sex/metric filters should make this 1:1 already,
    # but guard against duplicate rows in source files)
    dplyr::group_by(country, year) %>%
    dplyr::summarise(
      inc_rate  = mean(inc_rate,  na.rm = TRUE),
      inc_upper = mean(inc_upper, na.rm = TRUE),
      inc_lower = mean(inc_lower, na.rm = TRUE),
      inc_sigma = mean(inc_sigma, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Rename with disease prefix for the join into covariates_panel
  names(out)[names(out) == "inc_rate"]  <- paste0(out_prefix, "_Inc_Rate")
  names(out)[names(out) == "inc_upper"] <- paste0(out_prefix, "_Inc_Upper")
  names(out)[names(out) == "inc_lower"] <- paste0(out_prefix, "_Inc_Lower")
  names(out)[names(out) == "inc_sigma"] <- paste0(out_prefix, "_Inc_Sigma")
  out
}

message("Loading IHME GBD incidence...")
# Age band per disease comes from disease_age_band() (01). In the default
# 'infant' mode all bands are '<1 year'; under 'under5_tb_diph' TB and
# diphtheria use 'Under 5'. (Verify the exact age_name strings in your IHME
# files; "Under 5" vs "<5 years" vs separate "<1 year"+"1 to 4" rows differ by
# extraction. Combining separate rows requires population-weighting, not
# averaging of rates.)
ihme_tb       <- read_ihme_incidence("IHME_Tuberculosis_GBD.csv", "TB",
                                     age_band = disease_age_band("Tuberculosis"))
ihme_measles  <- read_ihme_incidence("IHME_Measles_GBD.csv",      "Measles",
                                     age_band = disease_age_band("Measles"))
ihme_diphth   <- read_ihme_incidence("IHME_Diphtheria_GBD.csv",   "Diphtheria",
                                     age_band = disease_age_band("Diphtheria"))
ihme_pert     <- read_ihme_incidence("IHME_Pertussis_GBD.csv",    "Pertussis",
                                     age_band = disease_age_band("Pertussis"))
# TETANUS (recommendation #8): this file MUST be GBD cause "Tetanus" (which
# EXCLUDES neonatal tetanus), read at the <1 band -> post-neonatal infant
# tetanus, matching the post-neonatal CFR and the DTP3 lever. If it is
# "Neonatal tetanus" or a combined series the causal chain is INVALID; set
# TETANUS_ARM='drop' in 01, or supply the non-neonatal cause file.
message("  NOTE: IHME_Tetanus_GBD.csv must be GBD cause 'Tetanus' (excludes ",
        "neonatal tetanus). See TETANUS_ARM (01).")
ihme_tet      <- read_ihme_incidence("IHME_Tetanus_GBD.csv",      "Tetanus",
                                     age_band = disease_age_band("Tetanus"))
message("  IHME rows: TB=", nrow(ihme_tb),
        " Measles=", nrow(ihme_measles),
        " Diphtheria=", nrow(ihme_diphth),
        " Pertussis=", nrow(ihme_pert),
        " Tetanus=", nrow(ihme_tet))

# ---- 5. Reshape WHO-UNICEF coverage ---------------------------------------
reshape_coverage <- function(df, label) {
  iso_col <- if ("SpatialDimValueCode" %in% names(df)) df$SpatialDimValueCode
  else NA_character_
  df %>%
    dplyr::mutate(
      country = harmonise_country(Location),
      ISO3    = iso_col,
      year    = suppressWarnings(as.integer(Period)),
      coverage = suppressWarnings(as.numeric(Value)),
      vaccine = label
    ) %>%
    dplyr::filter(!is.na(year), !is.na(coverage)) %>%
    dplyr::select(country, ISO3, vaccine, year, coverage)
}

bcg_long  <- reshape_coverage(bcg_data,  "BCG")
mcv1_long <- reshape_coverage(mcv1_data, "MCV1")
dpt_long  <- reshape_coverage(dpt_data,  "DTP3")

# Backfill ISO3 in MCV1/DTP3 using BCG as the canonical source
iso3_lookup <- bcg_long %>%
  dplyr::filter(!is.na(ISO3), ISO3 != "") %>%
  dplyr::distinct(country, ISO3)

backfill_iso3 <- function(df) {
  df %>%
    dplyr::left_join(iso3_lookup, by = "country", suffix = c("", "_lookup")) %>%
    dplyr::mutate(ISO3 = ifelse(is.na(ISO3) | ISO3 == "", ISO3_lookup, ISO3)) %>%
    dplyr::select(-ISO3_lookup)
}
mcv1_long <- backfill_iso3(mcv1_long)
dpt_long  <- backfill_iso3(dpt_long)

coverage_long <- dplyr::bind_rows(bcg_long, mcv1_long, dpt_long) %>%
  dplyr::filter(!is.na(ISO3), ISO3 != "") %>%
  dplyr::arrange(country, vaccine, year)

# ---- 6. Build country-year covariates panel -------------------------------
covariates_panel <- WB_births_long %>%
  dplyr::full_join(WB_pop_long %>%
                     dplyr::select(ISO3, year, Population, Population_per_1000),
                   by = c("ISO3", "year")) %>%
  dplyr::full_join(WB_life_long, by = c("ISO3", "year")) %>%
  dplyr::arrange(ISO3, year) %>%
  dplyr::mutate(Birth_Cohort = Birth_Rate * Population_per_1000)

# Attach country names from the coverage panel (WB files only carry ISO3)
country_iso_lookup <- coverage_long %>% dplyr::distinct(country, ISO3)
covariates_panel <- covariates_panel %>%
  dplyr::left_join(country_iso_lookup, by = "ISO3")

# Attach IHME incidence (per-country-year, all 5 diseases with UIs)
covariates_panel <- covariates_panel %>%
  dplyr::left_join(ihme_tb,      by = c("country", "year")) %>%
  dplyr::left_join(ihme_measles, by = c("country", "year")) %>%
  dplyr::left_join(ihme_diphth,  by = c("country", "year")) %>%
  dplyr::left_join(ihme_pert,    by = c("country", "year")) %>%
  dplyr::left_join(ihme_tet,     by = c("country", "year"))

# ---- 6b. Under-5 population denominator (recommendation #7) ----------------
# Only needed when AGE_BAND_MODE == 'under5_tb_diph'. The base 'infant' mode is
# unaffected. We do NOT approximate under-5 population from births (that would
# bias the denominator ~5x); a real series must be supplied. Expected file:
# <YLL Data>/WB_Population_Under5.csv in the same wide format as WB_Population
# Data.csv (e.g. World Bank SP.POP.0004.MA/FE summed, or UN WPP under-5 pop).
if (exists("AGE_BAND_MODE") && AGE_BAND_MODE == "under5_tb_diph") {
  u5_path <- file.path(YLL_DATA_DIR, "WB_Population_Under5.csv")
  if (!file.exists(u5_path)) {
    stop("AGE_BAND_MODE='under5_tb_diph' needs an under-5 population series at\n  ",
         u5_path, "\n  (wide format with Country.Code + yearly columns). ",
         "Add it, or set CVD_AGE_BAND_MODE='infant'.")
  }
  WB_u5_long <- to_long_wb(read_csv_safe(u5_path), "Pop_Under5")
  covariates_panel <- covariates_panel %>%
    dplyr::left_join(WB_u5_long, by = c("ISO3", "year"))
  message("  Loaded under-5 population denominator (Pop_Under5) for ",
          "AGE_BAND_MODE='under5_tb_diph'.")
}

# ---- 7. Interpolation + LOCF for terminal years ---------------------------
# Use linear interpolation for interior gaps. For trailing NAs (e.g. GBD
# data ending in 2021 while a conflict runs through 2024), use rule = 2 to
# carry the last observation forward. This is necessary because the most
# recent GBD release typically lags 2-3 years behind the present.
interpolate_panel <- function(df, vars) {
  df %>%
    dplyr::group_by(ISO3) %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(vars),
                                ~ {
                                  x <- .x
                                  if (sum(!is.na(x)) < 2) return(x)
                                  zoo::na.approx(x, x = year,
                                                 na.rm = FALSE, rule = 2)
                                })) %>%
    dplyr::ungroup()
}

incidence_vars <- c(
  "TB_Inc_Rate", "TB_Inc_Upper", "TB_Inc_Lower", "TB_Inc_Sigma",
  "Measles_Inc_Rate", "Measles_Inc_Upper", "Measles_Inc_Lower", "Measles_Inc_Sigma",
  "Diphtheria_Inc_Rate", "Diphtheria_Inc_Upper", "Diphtheria_Inc_Lower", "Diphtheria_Inc_Sigma",
  "Pertussis_Inc_Rate", "Pertussis_Inc_Upper", "Pertussis_Inc_Lower", "Pertussis_Inc_Sigma",
  "Tetanus_Inc_Rate", "Tetanus_Inc_Upper", "Tetanus_Inc_Lower", "Tetanus_Inc_Sigma"
)
.interp_vars <- c("Birth_Rate", "Population_per_1000", "Birth_Cohort",
                  "Life_Expectancy", incidence_vars)
if ("Pop_Under5" %in% names(covariates_panel))
  .interp_vars <- c(.interp_vars, "Pop_Under5")
covariates_panel <- interpolate_panel(covariates_panel, .interp_vars)

# ---- 8. Validation --------------------------------------------------------
study_iso <- conflict_info$ISO3
missing_iso <- setdiff(study_iso, unique(covariates_panel$ISO3))
if (length(missing_iso) > 0) {
  warning("ISO3 codes missing from covariates_panel: ",
          paste(missing_iso, collapse = ", "))
}

# Coverage availability check for the pre-conflict 3-year window
coverage_check <- coverage_long %>%
  dplyr::inner_join(conflict_info %>% dplyr::select(ISO3, country, preconflict_year),
                    by = c("ISO3", "country")) %>%
  dplyr::filter(year >= preconflict_year - 2L, year <= preconflict_year) %>%
  dplyr::count(country, vaccine)
message("\nPre-conflict 3-year coverage availability:")
print(coverage_check %>% tidyr::pivot_wider(names_from = vaccine, values_from = n))

# Incidence availability check for the study countries during conflict years
inc_check <- covariates_panel %>%
  dplyr::filter(ISO3 %in% study_iso) %>%
  dplyr::inner_join(conflict_info %>%
                      dplyr::select(ISO3, conflict_start, conflict_end),
                    by = "ISO3") %>%
  dplyr::filter(year >= conflict_start, year <= conflict_end) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    n_yr           = dplyr::n(),
    TB_avail       = sum(!is.na(TB_Inc_Rate)),
    Measles_avail  = sum(!is.na(Measles_Inc_Rate)),
    Diphth_avail   = sum(!is.na(Diphtheria_Inc_Rate)),
    Pert_avail     = sum(!is.na(Pertussis_Inc_Rate)),
    Tet_avail      = sum(!is.na(Tetanus_Inc_Rate)),
    .groups = "drop"
  )
message("\nIHME incidence availability during conflict years (cols = years available / total conflict years):")
print(inc_check)

# Duplicate (ISO3, year) check
dup_check <- covariates_panel %>%
  dplyr::count(ISO3, year) %>%
  dplyr::filter(n > 1)
if (nrow(dup_check) > 0) {
  warning("Duplicate (ISO3, year) rows in covariates_panel: ", nrow(dup_check),
          " - first few:\n",
          paste(utils::capture.output(print(utils::head(dup_check))),
                collapse = "\n"))
}

message("\nData loading complete.")
message("  coverage_long rows:    ", nrow(coverage_long))
message("  covariates_panel rows: ", nrow(covariates_panel))