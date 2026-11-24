# 01b_le.R | in: none | out: REF_L0, resolve_life_expectancy() | Methods: Years of life lost (reference life table)

# GBD 2019 reference (theoretical-minimum-risk) life table, ex at age 0. Infant cohort ->
# every death is valued at L0. National-LE sensitivities rescale draws in 09/12.
REF_L0 <- 88.8718951

ref_ex_at_age <- function(age) {
  stopifnot(all(age == 0))   # infant cohort only
  rep(REF_L0, length(age))
}
age_at_death_for_disease <- function(disease) 0
resolve_life_expectancy <- function(disease, iso3, year) REF_L0
