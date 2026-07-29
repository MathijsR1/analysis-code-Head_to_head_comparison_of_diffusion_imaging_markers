
## Head_to_head_comparison_of_diffusion_imaging_markers
## Mathijs T. Rosbergen


# Replace the paths below with the corresponding locations of your datasets.

# Set working directory containing participant folders with PSMD files
setwd("path_to_participant_folders")

# Get list of participant folders
participant_folders <- list.files()

# Read PSMD values for each participant
psmd_list <- list()

for (folder in participant_folders) {
  file_path <- file.path(folder, "psmd")
  
  # Read PSMD file and convert value(s) to numeric
  psmd_list[[basename(folder)]] <- as.numeric(readLines(file_path))
}

# Combine all participant PSMD values into a single dataframe
df_psmd <- do.call(
  rbind,
  lapply(names(psmd_list), function(participant_name) {
    data.frame(
      Participant = participant_name,
      PSMD = psmd_list[[participant_name]]
    )
  })
)

# Save PSMD dataframe
df_psmd_save <- df_psmd
save(
  df_psmd_save,
  file = "PSMD.RData"
)


# Set working directory containing input data files
setwd("path_input_data")

library(haven)





# ------------------------------------------------------------
# Load eligibility data and select participants with available scans
# ------------------------------------------------------------

eligib <- read_sav("path_eligibility_file")

eligib <- eligib[eligib$scan_present == 1, ]

# Number of unique participants with available scan data
length(unique(eligib$ergoid))


# ------------------------------------------------------------
# Merge eligibility data with global DTI measures
# ------------------------------------------------------------

dti_global <- read_sav("path_global_dti_file")

df_eligib_dti <- merge(
  eligib,
  dti_global,
  by.x = "name_x",
  by.y = "name_y",
  all = TRUE
)

length(unique(df_eligib_dti$ergoid))


# ------------------------------------------------------------
# Add tract-specific DTI measures (IFO and PTR)
# ------------------------------------------------------------

dti_tract <- read_sav("path_dti_tracts_file")

# Keep only relevant tract measures
dti_tract <- dti_tract[c(
  "case",
  "ifo_FA",
  "ifo_MD",
  "ptr_FA",
  "ptr_MD"
)]

df_eligib_dti_compl <- merge(
  df_eligib_dti,
  dti_tract,
  by.x = "name_x",
  by.y = "name_y",
  all = TRUE
)


# ------------------------------------------------------------
# Add PSMD values
# ------------------------------------------------------------

df_eligib_dti_psmd <- merge(
  df_eligib_dti_compl,
  df_psmd,
  by.x = "name_x",
  by.y = "name_x",
  all = TRUE
)


# ------------------------------------------------------------
# Remove participants without DTI/eligibility identifier
# ------------------------------------------------------------

length(unique(df_eligib_dti_psmd$ergoid))

df_eligib_dti_psmd <- df_eligib_dti_psmd[
  !is.na(df_eligib_dti_psmd$ergoid),
]


# ------------------------------------------------------------
# Check availability of DTI and PSMD measurements
# ------------------------------------------------------------

length(unique(df_eligib_dti_psmd$ergoid))

measurement_check <- c(
  FA_DTI  = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$FA_dti)])),
  IFO_FA  = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$ifo_FA)])),
  PTR_FA  = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$ptr_FA)])),
  MD_DTI  = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$MD_dti)])),
  IFO_MD  = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$ifo_MD)])),
  PTR_MD  = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$ptr_MD)])),
  PSMD    = length(unique(df_eligib_dti_psmd$ergoid[!is.na(df_eligib_dti_psmd$PSMD)]))
)

measurement_check


# Number of participants with complete MD/PSMD information
length(unique(df_eligib_dti_psmd$ergoid[
  !is.na(df_eligib_dti_psmd$PSMD) &
    !is.na(df_eligib_dti_psmd$MD_dti)
]))


# Number of participants with complete MD information across all tracts
length(unique(df_eligib_dti_psmd$ergoid[
  !is.na(df_eligib_dti_psmd$PSMD) &
    !is.na(df_eligib_dti_psmd$MD_dti) &
    !is.na(df_eligib_dti_psmd$ifo_MD) &
    !is.na(df_eligib_dti_psmd$ptr_MD)
]))


# Participants with both PSMD and global MD available
unique(df_eligib_dti_psmd$ergoid[
  !is.na(df_eligib_dti_psmd$PSMD) &
    !is.na(df_eligib_dti_psmd$MD_dti)
])


# ------------------------------------------------------------
# Inspect PSMD distribution and remove extreme values
# ------------------------------------------------------------

hist(df_eligib_dti_psmd$PSMD, breaks = 500)
boxplot(df_eligib_dti_psmd$PSMD)

median(df_eligib_dti_psmd$PSMD, na.rm = TRUE)
min(df_eligib_dti_psmd$PSMD, na.rm = TRUE)
max(df_eligib_dti_psmd$PSMD, na.rm = TRUE)


# Remove implausibly high PSMD values (based on literature)
df_eligib_dti_psmd <- df_eligib_dti_psmd[
  df_eligib_dti_psmd$PSMD < 0.0007,
]

length(unique(df_eligib_dti_psmd$ergoid))

hist(df_eligib_dti_psmd$PSMD, breaks = 500)


# Re-check number of participants after PSMD filtering

length(unique(df_eligib_dti_psmd$ergoid[
  !is.na(df_eligib_dti_psmd$PSMD) &
    !is.na(df_eligib_dti_psmd$MD_dti)
]))

length(unique(df_eligib_dti_psmd$ergoid[
  !is.na(df_eligib_dti_psmd$PSMD) &
    !is.na(df_eligib_dti_psmd$MD_dti) &
    !is.na(df_eligib_dti_psmd$ifo_MD) &
    !is.na(df_eligib_dti_psmd$ptr_MD)
]))








# ============================================================
# Create analysis dataset
# ============================================================

library(haven)
library(data.table)
library(lubridate)
library(dplyr)


# ------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------

setwd("path_for_working_directory")


# ------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------

# Vascular markers
infarct_data <- read_sav("path_infarct_data")

# Covariates
load("path_covariates_data")

# Clinical datasets
mci_data   <- read_sav("path_mci_data")
dem_data   <- read_sav("path_dementia_data")
diab_data  <- read_sav("path_diabetes_data")
apoe_data  <- read_sav("path_apoe_data")

# White matter hyperintensities
wmh_data <- read_sav("path_whitematterhyperintensities_data")

# Stroke and mortality data
stroke_data <- read_sav("path_stroke_data")
mort_data   <- read_sav("path_mortality_data")


# ------------------------------------------------------------
# Load visit dates
# ------------------------------------------------------------

visit_path <- "path_visit_data"

visit_data_e4 <- read_sav(
  paste0(visit_path, "/path_for_cohort_wave_e4")
) %>%
  select(ergoid, e4_3494)

visit_data_ej <- read_sav(
  paste0(visit_path, "/path_for_cohort_wave_ej")
) %>%
  select(ergoid, ej_3494) %>%
  rename(e4_3494 = ej_3494)

# Combine E4 and early follow-up visits
visit_data_e4 <- rbind(
  visit_data_e4,
  visit_data_ej
)

visit_data_e5 <- read_sav(
  paste0(visit_path, "/path_for_cohort_wave_e5")
)

visit_data_e6 <- read_sav(
  paste0(visit_path, "/path_for_cohort_wave_e6")
)

# Merge all visit dates
visit_data <- merge(
  visit_data_e4,
  visit_data_e5,
  by = "ID",
  all = TRUE
)

visit_data <- merge(
  visit_data,
  visit_data_e6,
  by = "ID",
  all = TRUE
)

# Remove duplicate cohort variables created by merge
visit_data <- visit_data %>%
  select(-rs_cohort.y) %>%
  rename(rs_cohort = rs_cohort.x)



# ============================================================
# Merge all datasets
# ============================================================

# Add visit dates
df_psmd_visit <- merge(
  df_eligib_dti_psmd,
  visit_data,
  by = "ID",
  all.x = TRUE
)


# Add baseline covariates
df_psmd_cov <- merge(
  df_psmd_visit,
  cov_data %>% select(-rs_cohort, -rs_cohort.x, -rs_cohort.y),
  by = "ID",
  all.x = TRUE
)


# Add APOE genotype
df_psmd_cov_apoe <- merge(
  df_psmd_cov,
  apoe_data,
  by = "ID",
  all.x = TRUE
)


# Add diabetes information
df_psmd_cov_apoe_dm <- merge(
  df_psmd_cov_apoe,
  diab_data %>% select(-rs_cohort),
  by = "ID",
  all.x = TRUE
)


# Add MCI information
df_psmd_cov_apoe_dm_mci <- merge(
  df_psmd_cov_apoe_dm,
  mci_data %>% select(-rs_cohort),
  by = "ID",
  all.x = TRUE
)


# Add dementia information
df_psmd_cov_apoe_dm_mci_dem <- merge(
  df_psmd_cov_apoe_dm_mci,
  dem_data %>% select(-rs_cohort),
  by = "ID",
  all.x = TRUE
)


# Add stroke and mortality information
df_psmd_cov_apoe_dm_mci_dem_str <- merge(
  df_psmd_cov_apoe_dm_mci_dem,
  stroke_data %>% select(-rs_cohort),
  by = "ID",
  all.x = TRUE
)

df_psmd_cov_apoe_dm_mci_dem_str_mort <- merge(
  df_psmd_cov_apoe_dm_mci_dem_str,
  mort_data %>% select(-rs_cohort),
  by = "ID",
  all.x = TRUE
)


# Add WMH measures
df_psmd_cov_apoe_dm_mci_dem_wmh <- merge(
  df_psmd_cov_apoe_dm_mci_dem_str_mort,
  wmh_data,
  by = "ID",
  all.x = TRUE
)


# Add infarct information
df_wide <- merge(
  df_psmd_cov_apoe_dm_mci_dem_wmh,
  infarct_data,
  by = "ID",
  all.x = TRUE
)



# ============================================================
# Add cognitive data
# ============================================================

load("path_to_cognition_data")

df_wide <- merge(
  df_wide,
  cogn_data %>% select(-rs_cohort),
  by.x = "ID.x",
  by.y = "ID",
  all.x = TRUE
)


# Rename visit date variables
df_wide <- df_wide %>%
  rename(
    visitdate_e4 = e4_3494,
    visitdate_e5 = e5_3494,
    visitdate_e6 = e6_3494
  )


# Convert visit dates to Date format
df_wide <- df_wide %>%
  mutate(
    visitdate_e4 = as.Date(visitdate_e4),
    visitdate_e5 = as.Date(visitdate_e5),
    visitdate_e6 = as.Date(visitdate_e6)
  )



# ============================================================
# Convert wide dataset to longitudinal format
# ============================================================

variables <- c(
  "visitdate_e",
  "bmi",
  "smoke",
  "oh",
  "sbp",
  "dbp",
  "ht",
  "htdrug",
  "chol",
  "hdl",
  "lip_e",
  "MCI",
  "LDST_e",
  "PPBsum_e",
  "STR3_adjusted_e",
  "WFT_e",
  "WLTdel_e"
)


assign_cross_sectional <- function(data_function, variables) {
  
  scanround4 <- c("RS1_4X", "RS2_2", "RS2_2X", "RS3_1")
  scanround5 <- c("RS3_1X", "RS3_2", "RS1_5", "RS2_3")
  scanround6 <- c("RS1_6", "RS2_4")
  
  for (variable in variables) {
    
    data_function[[variable]] <- NA
    
    for (row in 1:nrow(data_function)) {
      
      if (data_function$scanround[row] %in% scanround4) {
        data_function[row, variable] <-
          data_function[row, paste0(variable, "4")]
      }
      
      if (data_function$scanround[row] %in% scanround5) {
        data_function[row, variable] <-
          data_function[row, paste0(variable, "5")]
      }
      
      if (data_function$scanround[row] %in% scanround6 &
          !variable %in% c("chol", "hdl")) {
        data_function[row, variable] <-
          data_function[row, paste0(variable, "6")]
      }
    }
  }
  
  data_function
}


long_data <- assign_cross_sectional(df_wide, variables)

names(long_data)[names(long_data) == "ID.x"] <- "ID"

long_data$visitdate_e <- as.Date(long_data$visitdate_e)






# ============================================================
# Define DTI availability and create baseline population
# ============================================================

# Participant considered to have complete DTI data if all required
# global and tract-specific measures are available
long_data$DTI_present <- ifelse(
  !is.na(long_data$MD_dti) &
    !is.na(long_data$PSMD) &
    !is.na(long_data$FA_dti) &
    !is.na(long_data$ifo_MD) &
    !is.na(long_data$ptr_MD),
  1,
  0
)


# Select participants with complete DTI data
# and retain first available visit per participant

df_allvisits <- long_data %>%
  filter(DTI_present == 1) %>%
  arrange(ergoid, visitdate_e)

df_first <- df_allvisits %>%
  group_by(ergoid) %>%
  slice(1)



# ============================================================
# Standardize DTI measures
# ============================================================

# Standardization allows interpretation of hazard ratios per
# standard deviation increase

dti_variables <- c(
  "PSMD",
  "FA_dti",
  "MD_dti",
  "ifo_FA",
  "ptr_FA",
  "ifo_MD",
  "ptr_MD"
)

for (var in dti_variables) {
  df_first[[paste0(var, "_std")]] <-
    (df_first[[var]] - mean(df_first[[var]])) /
    sd(df_first[[var]])
}


# Reverse FA measures so that higher values represent worse WM integrity
df_first$FA_std     <- 0 - df_first$FA_dti_std
df_first$ifo_FA_std <- 0 - df_first$ifo_FA_std
df_first$ptr_FA_std <- 0 - df_first$ptr_FA_std



# ============================================================
# Calculate baseline characteristics
# ============================================================

# Age at time of MRI scan
df_first$age_scan <- difftime(
  df_first$scandate.x,
  df_first$date_of_birth,
  units = "days"
)

df_first$age_scan <- as.numeric(df_first$age_scan) / 365.25



# ------------------------------------------------------------
# Diabetes status at time of scan
# ------------------------------------------------------------

table(df_first$prevalent_DM)

df_first$DM_scan <- df_first$prevalent_DM

table(df_first$DM_scan)

# Add incident diabetes if diagnosis occurred before scan date
for (i in seq_len(nrow(df_first))) {
  
  if (is.na(df_first$incident_DM[i])) {
    print(df_first$ergoid[i])
    
  } else if (
    df_first$incident_DM[i] == 1 &
    df_first$incident_DM_date[i] < df_first$scandate.x[i]
  ) {
    df_first$DM_scan[i] <- 1
  }
}

table(df_first$DM_scan)



# ------------------------------------------------------------
# Recode APOE genotype
# ------------------------------------------------------------

table(df_first$apoe)

df_first$apoe[df_first$apoe %in% c(24, 34)] <- 1
df_first$apoe[df_first$apoe == 44] <- 2
df_first$apoe[df_first$apoe %in% c(22, 23)] <- 3
df_first$apoe[df_first$apoe == 33] <- 4

table(df_first$apoe)



# ------------------------------------------------------------
# Check missingness of cardiovascular variables
# ------------------------------------------------------------

table(df_first$ht)
table(df_first$htdrug)

sum(is.na(df_first$ht))
sum(is.na(df_first$sbp))
sum(is.na(df_first$dbp))
sum(is.na(df_first$htdrug))
sum(is.na(df_first$lip_e))



# ============================================================
# Multiple imputation of covariates
# ============================================================

library(mice)

# Select variables used for imputation diagnostics
df_c_v <- df_first[c(
  "age_scan",
  "sex",
  "bmi",
  "ht",
  "sbp",
  "dbp",
  "htdrug",
  "smoke",
  "oh",
  "chol",
  "hdl",
  "lip_e",
  "education",
  "DM_scan",
  "apoe",
  "microbleeds_present",
  "lacunar_infarcts_present"
)]

# Inspect missing data pattern
md.pattern(df_c_v)



# Lipid lowering medication:
# value 9 represents missing/not available and is recoded as NA

sum(is.na(df_first$lip_e))

df_first$lip_e[df_first$lip_e == 9] <- NA

sum(is.na(df_first$lip_e))



# Prepare categorical variables for imputation

library(sjmisc)

df_first$education <- as.factor(df_first$education)
df_first$apoe <- as.factor(df_first$apoe)



# Dataset used for multiple imputation

cox_data_imp <- df_first[c(
  "ergoid",
  "sex",
  "education",
  "bmi",
  "apoe",
  "smoke",
  "oh",
  "ht",
  "sbp",
  "dbp",
  "htdrug",
  "chol",
  "hdl",
  "lip_e",
  "DM_scan",
  "microbleeds_present",
  "lacunar_infarcts_present"
)]


# Predictor matrix:
# participant identifier is excluded as predictor

pred_matrix <- make.predictorMatrix(cox_data_imp)
pred_matrix[, "ergoid"] <- 0


# Perform multiple imputation

imp <- mice(
  cox_data_imp,
  maxit = 5,
  predictorMatrix = pred_matrix
)


# Combine imputed variables with original dataset

imp_single <- merge_imputations(
  cox_data_imp,
  imp = imp,
  ori = df_first
)



# ============================================================
# Baseline characteristics table
# ============================================================

library(tableone)

CreateTableOne(
  vars = c(
    "age_scan",
    "sex.x",
    "bmi_imp",
    "ht_imp",
    "sbp_imp",
    "dbp_imp",
    "htdrug_imp",
    "smoke_imp",
    "oh_imp",
    "chol_imp",
    "hdl_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  factorVars = c(
    "sex.x",
    "ht_imp",
    "htdrug_imp",
    "smoke_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  data = imp_single
)



# ============================================================
# Exploratory plots
# ============================================================

plot(imp_single$age_scan, df_first$PSMD)
plot(imp_single$MD_dti, df_first$PSMD)












# ============================================================
# Cognitive analysis population
# ============================================================

# Start with imputed baseline dataset
df_cogn <- imp_single


# ------------------------------------------------------------
# Exclude participants with prevalent dementia
# ------------------------------------------------------------

df_cogn <- df_cogn %>%
  arrange(ergoid)

# Check cohort membership
table(df_cogn$IC.x)

# Keep only IC = 2 participants
df_IC <- df_cogn$ergoid[df_cogn$IC.x == 1]

df_cogn <- df_cogn[df_cogn$IC.x == 2, ]


# Remove participants with dementia at baseline MRI

table(df_cogn$dementia_prevalent)

df_cogn <- df_cogn[
  df_cogn$dementia_prevalent == 0,
]


# Exclude dementia diagnoses within 3 months after MRI scan
# (possible dementia already present at time of scanning)

df_cogn$date_scan3 <- df_cogn$scandate.x %m+% months(3)

df_cogn$dem_excl <- 0

for (i in seq_len(nrow(df_cogn))) {
  
  if (
    df_cogn$dementia_incident[i] == 1 &
    df_cogn$dementia_date[i] < df_cogn$date_scan3[i]
  ) {
    
    df_cogn$dem_excl[i] <- 1
    print(df_cogn$ergoid[i])
  }
}

table(df_cogn$dem_excl)

df_cogn <- df_cogn[df_cogn$dem_excl == 0, ]

length(unique(df_cogn$ergoid))



# ============================================================
# Check availability of baseline cognitive tests
# ============================================================

cognitive_tests <- c(
  "PPBsum_e",
  "WLTdel_e",
  "STR3_adjusted_e",
  "WFT_e",
  "LDST_e"
)


# Participants with at least one cognitive test available

has_one_test <- rowSums(
  !is.na(df_cogn[, cognitive_tests])
) > 0


# Participants with complete cognitive data

has_complete_tests <- rowSums(
  !is.na(df_cogn[, cognitive_tests])
) == length(cognitive_tests)


length(df_cogn$ergoid[has_one_test])
length(df_cogn$ergoid[has_complete_tests])


# Compare participants with incomplete versus complete testing

df_check_cogn_tests_cross <- df_cogn

df_check_completetest <- df_check_cogn_tests_cross[has_complete_tests, ]

df_check_1test <- df_check_cogn_tests_cross[has_one_test, ]

df_check_1test <- df_check_1test[
  !df_check_1test$ergoid %in% df_check_completetest$ergoid,
]


# Compare baseline characteristics

CreateTableOne(
  vars = c(
    "age_scan",
    "sex.x",
    "bmi_imp",
    "ht_imp",
    "sbp_imp",
    "dbp_imp",
    "htdrug_imp",
    "smoke_imp",
    "oh_imp",
    "chol_imp",
    "hdl_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  factorVars = c(
    "sex.x",
    "ht_imp",
    "htdrug_imp",
    "smoke_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  data = df_check_1test
)


CreateTableOne(
  vars = c(
    "age_scan",
    "sex.x",
    "bmi_imp",
    "ht_imp",
    "sbp_imp",
    "dbp_imp",
    "htdrug_imp",
    "smoke_imp",
    "oh_imp",
    "chol_imp",
    "hdl_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  factorVars = c(
    "sex.x",
    "ht_imp",
    "htdrug_imp",
    "smoke_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  data = df_check_completetest
)


median(df_check_1test$oh_imp)
quantile(df_check_1test$oh_imp)



# ============================================================
# Longitudinal cognitive follow-up availability
# ============================================================

library(tidyr)


# Variables for each follow-up visit

variables_e4 <- c(
  "visitdate_e4",
  "PPBsum_e4",
  "WLTdel_e4",
  "STR3_adjusted_e4",
  "WFT_e4",
  "LDST_e4"
)

variables_e5 <- c(
  "visitdate_e5",
  "PPBsum_e5",
  "WLTdel_e5",
  "STR3_adjusted_e5",
  "WFT_e5",
  "LDST_e5"
)

variables_e6 <- c(
  "visitdate_e6",
  "PPBsum_e6",
  "WLTdel_e6",
  "STR3_adjusted_e6",
  "WFT_e6",
  "LDST_e6"
)


# Convert wide cognitive data to long format

df_check_cogn_tests_long <- df_cogn %>%
  pivot_longer(
    cols = all_of(c(
      variables_e4,
      variables_e5,
      variables_e6
    )),
    names_to = c(".value", "visit"),
    names_pattern = "(.*)(e[456])"
  )


# Keep visits where participant attended centre

df_check_cogn_tests_long <- df_check_cogn_tests_long[
  !is.na(df_check_cogn_tests_long$visitdate_),
]


df_check_cogn_tests_long <- df_check_cogn_tests_long %>%
  arrange(ergoid, visitdate_)



# Keep only visits after baseline MRI

df_check_cogn_tests_long <- df_check_cogn_tests_long[
  df_check_cogn_tests_long$visitdate_ >= 
    df_check_cogn_tests_long$visitdate_e,
]



# Select first two available visits

df_check_cogn_tests_long_12visit <- df_check_cogn_tests_long %>%
  group_by(ergoid) %>%
  filter(row_number() <= 2)



# Identify participants with follow-up

followup_ids <- df_check_cogn_tests_long_12visit$ergoid[
  duplicated(df_check_cogn_tests_long_12visit$ergoid)
]


df_check_cogn_tests_long_2visits <-
  df_check_cogn_tests_long_12visit[
    df_check_cogn_tests_long_12visit$ergoid %in% followup_ids,
  ]


df_check_cogn_tests_long_1visit <-
  df_check_cogn_tests_long_12visit[
    !df_check_cogn_tests_long_12visit$ergoid %in% followup_ids,
  ]



# Select follow-up visit only

df_check_cogn_tests_long <- df_check_cogn_tests_long_2visits %>%
  group_by(ergoid) %>%
  filter(row_number() == 2)



# ============================================================
# Create final longitudinal cognition dataset
# ============================================================

df_cogn <- df_cogn[
  rowSums(!is.na(df_cogn[, cognitive_tests])) == length(cognitive_tests),
]

length(unique(df_cogn$ergoid))


df_long_cogn <- df_check_cogn_tests_long_12visit


# Keep only visits with complete cognitive tests

cognitive_tests_long <- c(
  "PPBsum_e",
  "WLTdel_e",
  "STR3_adjusted_e",
  "WFT_e",
  "LDST_e"
)

df_long_cogn <- df_long_cogn[
  rowSums(!is.na(df_long_cogn[, cognitive_tests_long])) ==
    length(cognitive_tests_long),
]


df_long_cogn <- df_long_cogn %>%
  arrange(ergoid, visitdate_)



# ============================================================
# Calculate g-factor using PCA
# ============================================================

tests <- c(
  "LDST_",
  "PPBsum_",
  "STR3_adjusted_",
  "WFT_",
  "WLTdel_"
)


# Ensure numeric format

df_long_cogn[tests] <- lapply(
  df_long_cogn[tests],
  as.numeric
)


# Principal component analysis

pca_result <- prcomp(
  df_long_cogn[, tests],
  center = TRUE,
  scale = TRUE
)


# First principal component represents general cognition

df_long_cogn$g_factor <- scale(
  pca_result$x[, 1]
)


summary(pca_result)

loadings_g_factor <- pca_result$rotation[, 1]

print(loadings_g_factor)



# ============================================================
# Standardize individual cognitive tests
# ============================================================

for (test in tests) {
  
  df_long_cogn[[paste0(test, "std")]] <-
    (
      df_long_cogn[[test]] -
        mean(df_long_cogn[[test]])
    ) /
    sd(df_long_cogn[[test]])
}












############################################################
#### Cross-sectional analyses
############################################################

library(splines)
library(car)
library(ggplot2)
library(tableone)


############################################################
# Prepare WMH variable
############################################################

# Correct WMH volume for intracranial volume
df_long_cogn$wmh_icv <- df_long_cogn$Total_wml /
  df_long_cogn$ICV_from_mask

# Log-transform and standardize WMH
df_long_cogn$wmh_log <- log(df_long_cogn$wmh_icv)

df_long_cogn$wmh_std <- (
  df_long_cogn$wmh_log - mean(df_long_cogn$wmh_log)
) / sd(df_long_cogn$wmh_log)



############################################################
# Select first visit per participant
############################################################

df_long_cogn_f <- df_long_cogn %>%
  group_by(ergoid) %>%
  slice(1)



############################################################
# Baseline characteristics
############################################################

CreateTableOne(
  c(
    "age_scan",
    "sex.x",
    "bmi_imp",
    "ht_imp",
    "sbp_imp",
    "dbp_imp",
    "htdrug_imp",
    "smoke_imp",
    "oh_imp",
    "chol_imp",
    "hdl_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  factorVars = c(
    "sex.x",
    "ht_imp",
    "htdrug_imp",
    "smoke_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  data = df_long_cogn_f
)

median(df_long_cogn_f$oh_imp)
quantile(df_long_cogn_f$oh_imp)



############################################################
# Age versus cognition plots
############################################################

plot_age <- function(y_variable, y_label, title, filename) {
  
  pdf(
    paste0(
      "path_output",
      filename
    ),
    width = 5,
    height = 5
  )
  
  print(
    ggplot(
      df_long_cogn_f,
      aes(
        x = age_scan,
        y = .data[[y_variable]]
      )
    ) +
      geom_point(size = 0.5) +
      geom_smooth(
        method = "lm",
        formula = y ~ ns(x, df = 2),
        se = FALSE
      ) +
      labs(
        x = "Age [years]",
        y = y_label
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line()
      ) +
      ggtitle(title)
  )
  
  dev.off()
}


plot_age(
  "g_factor",
  "G-factor",
  "G-factor vs. age",
  "Gfactor_age.pdf"
)

plot_age(
  "LDST_std",
  "LDST",
  "LDST vs. age",
  "LDST_age.pdf"
)

plot_age(
  "PPBsum_std",
  "PPBsum",
  "PPBsum vs. age",
  "PPBsum_age.pdf"
)

plot_age(
  "STR3_adjusted_std",
  "Stroop",
  "Stroop vs. age",
  "STR3_adjusted_age.pdf"
)

plot_age(
  "WFT_std",
  "WFT",
  "WFT vs. age",
  "WFT_age.pdf"
)

plot_age(
  "WLTdel_std",
  "WLT",
  "WLT vs. age",
  "WLT_del_age.pdf"
)



############################################################
# WMH analysis
############################################################

wmh_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

wmh_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(wmh_age, wmh_age_2)


lm_WMH_G_1 <- lm(
  g_factor ~ wmh_std,
  data = df_long_cogn_f
)

lm_WMH_G_2 <- lm(
  g_factor ~ wmh_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_WMH_G_3 <- lm(
  g_factor ~ wmh_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_wmh_G_1 <- AIC(lm_WMH_G_1)
AIC_wmh_G_2 <- AIC(lm_WMH_G_2)
AIC_wmh_G_3 <- AIC(lm_WMH_G_3)

summary(lm_WMH_G_3)
confint(lm_WMH_G_3)



############################################################
# PSMD analysis
############################################################

hist(df_long_cogn_f$PSMD_std, prob = TRUE)

lines(
  density(df_long_cogn_f$PSMD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$PSMD_std
)


# Log transformed PSMD

df_long_cogn_f$PSMD_log <- log(df_long_cogn_f$PSMD)

df_long_cogn_f$PSMD_log_std <- (
  df_long_cogn_f$PSMD_log -
    mean(df_long_cogn_f$PSMD_log)
) / sd(df_long_cogn_f$PSMD_log)



psmd_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

psmd_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(psmd_age, psmd_age_2)



lm_PSMD_G_1 <- lm(
  g_factor ~ PSMD_std + pedir,
  data = df_long_cogn_f
)

lm_PSMD_G_2 <- lm(
  g_factor ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_PSMD_G_3 <- lm(
  g_factor ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_psmd_G_1 <- AIC(lm_PSMD_G_1)
AIC_psmd_G_2 <- AIC(lm_PSMD_G_2)
AIC_psmd_G_3 <- AIC(lm_PSMD_G_3)



lm_PSMD_log_1 <- lm(
  g_factor ~ PSMD_log_std,
  data = df_long_cogn_f
)

lm_PSMD_log_2 <- lm(
  g_factor ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_PSMD_log_3 <- lm(
  g_factor ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_psmd_log_1 <- AIC(lm_PSMD_log_1)
AIC_psmd_log_2 <- AIC(lm_PSMD_log_2)
AIC_psmd_log_3 <- AIC(lm_PSMD_log_3)


summary(lm_PSMD_G_3)
confint(lm_PSMD_G_3)
AIC_psmd_G_3

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)
AIC_psmd_log_3


hist(df_long_cogn_f$PSMD)
hist(df_long_cogn_f$PSMD_log_std)





############################################################
# MD analysis
############################################################

hist(
  df_long_cogn_f$MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$MD_std
)


# Test linear versus nonlinear age association

md_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(md_age, md_age_2)



# MD models

lm_MD_G_1 <- lm(
  g_factor ~ MD_std + pedir,
  data = df_long_cogn_f
)

lm_MD_G_2 <- lm(
  g_factor ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_MD_G_3 <- lm(
  g_factor ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_md_G_1 <- AIC(lm_MD_G_1)
AIC_md_G_2 <- AIC(lm_MD_G_2)
AIC_md_G_3 <- AIC(lm_MD_G_3)


summary(lm_MD_G_3)
confint(lm_MD_G_3)



############################################################
# FA analysis
############################################################

hist(
  df_long_cogn_f$FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$FA_std
)


fa_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(fa_age, fa_age_2)



# FA models

lm_FA_G_1 <- lm(
  g_factor ~ FA_std + pedir,
  data = df_long_cogn_f
)

lm_FA_G_2 <- lm(
  g_factor ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_FA_G_3 <- lm(
  g_factor ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_fa_G_1 <- AIC(lm_FA_G_1)
AIC_fa_G_2 <- AIC(lm_FA_G_2)
AIC_fa_G_3 <- AIC(lm_FA_G_3)


summary(lm_FA_G_3)
confint(lm_FA_G_3)



############################################################
# IFO MD analysis
############################################################

hist(
  df_long_cogn_f$ifo_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ifo_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_MD_std
)


ifo_md_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

ifo_md_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(ifo_md_age, ifo_md_age_2)



lm_ifo_MD_G_1 <- lm(
  g_factor ~ ifo_MD_std + pedir,
  data = df_long_cogn_f
)

lm_ifo_MD_G_2 <- lm(
  g_factor ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_ifo_MD_G_3 <- lm(
  g_factor ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ifo_md_G_1 <- AIC(lm_ifo_MD_G_1)
AIC_ifo_md_G_2 <- AIC(lm_ifo_MD_G_2)
AIC_ifo_md_G_3 <- AIC(lm_ifo_MD_G_3)


summary(lm_ifo_MD_G_3)
confint(lm_ifo_MD_G_3)



############################################################
# IFO FA analysis
############################################################

hist(
  df_long_cogn_f$ifo_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_FA_std
)


ifo_fa_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

ifo_fa_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(ifo_fa_age, ifo_fa_age_2)



lm_ifo_FA_G_1 <- lm(
  g_factor ~ ifo_FA_std + pedir,
  data = df_long_cogn_f
)

lm_ifo_FA_G_2 <- lm(
  g_factor ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_ifo_FA_G_3 <- lm(
  g_factor ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ifo_fa_G_1 <- AIC(lm_ifo_FA_G_1)
AIC_ifo_fa_G_2 <- AIC(lm_ifo_FA_G_2)
AIC_ifo_fa_G_3 <- AIC(lm_ifo_FA_G_3)


summary(lm_ifo_FA_G_3)
confint(lm_ifo_FA_G_3)
AIC_ifo_fa_G_3



############################################################
# PTR MD analysis
############################################################

hist(
  df_long_cogn_f$ptr_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ptr_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_MD_std
)


ptr_md_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

ptr_md_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(ptr_md_age, ptr_md_age_2)



lm_ptr_MD_G_3 <- lm(
  g_factor ~ ptr_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_md_G_3 <- AIC(lm_ptr_MD_G_3)

AIC_ptr_md_G_3



############################################################
# PTR FA analysis
############################################################

hist(
  df_long_cogn_f$ptr_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_FA_std
)


ptr_fa_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_f
)

ptr_fa_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(ptr_fa_age, ptr_fa_age_2)



lm_ptr_FA_G_3 <- lm(
  g_factor ~ ptr_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_fa_G_3 <- AIC(lm_ptr_FA_G_3)

AIC_ptr_fa_G_3



############################################################
# Model assumption checks
############################################################

check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  # Residuals versus fitted values
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(h = 0)
  
  
  # Independence
  cat("Durbin-Watson Test:\n")
  print(durbinWatsonTest(model))
  
  
  # Homoscedasticity
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(h = 0)
  
  
  # Normality residuals
  qqnorm(residuals(model))
  qqline(residuals(model))
  
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro.test(residuals(model)))
  
  
  # Multicollinearity
  cat("Variance Inflation Factor:\n")
  print(vif(model))
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_G_1)
check_assumptions(lm_FA_G_2)
check_assumptions(lm_FA_G_3)

check_assumptions(lm_MD_G_1)
check_assumptions(lm_MD_G_2)
check_assumptions(lm_MD_G_3)

check_assumptions(lm_PSMD_G_1)
check_assumptions(lm_PSMD_G_2)
check_assumptions(lm_PSMD_G_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Final model summaries
############################################################

summary(lm_FA_G_1)
confint(lm_FA_G_1)
AIC(lm_FA_G_1)

summary(lm_FA_G_2)
confint(lm_FA_G_2)
AIC(lm_FA_G_2)

summary(lm_FA_G_3)
confint(lm_FA_G_3)
AIC(lm_FA_G_3)


summary(lm_MD_G_1)
confint(lm_MD_G_1)
AIC(lm_MD_G_1)

summary(lm_MD_G_2)
confint(lm_MD_G_2)
AIC(lm_MD_G_2)

summary(lm_MD_G_3)
confint(lm_MD_G_3)
AIC(lm_MD_G_3)


summary(lm_PSMD_G_1)
confint(lm_PSMD_G_1)
AIC(lm_PSMD_G_1)

summary(lm_PSMD_G_2)
confint(lm_PSMD_G_2)
AIC(lm_PSMD_G_2)

summary(lm_PSMD_G_3)
confint(lm_PSMD_G_3)
AIC(lm_PSMD_G_3)


summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)















############################################################
#### Longitudinal analyses
############################################################


############################################################
# Create longitudinal cognitive dataset
############################################################

# Select participants with at least two cognitive measurements

df_long_cogn_long_ergoid <- df_long_cogn$ergoid[
  duplicated(df_long_cogn$ergoid)
]

length(unique(df_long_cogn$ergoid))


df_long_cogn_long <- df_long_cogn[
  df_long_cogn$ergoid %in% df_long_cogn_long_ergoid,
]

length(unique(df_long_cogn_long$ergoid))


# Order visits chronologically

df_long_cogn_long <- df_long_cogn_long %>%
  arrange(ergoid, visitdate_)



############################################################
# Select baseline and follow-up visits
############################################################

# First available cognitive visit

df_long_cogn_long_f <- df_long_cogn_long %>%
  group_by(ergoid) %>%
  filter(row_number() == 1)


# Second available cognitive visit

df_long_cogn_long_l <- df_long_cogn_long %>%
  group_by(ergoid) %>%
  filter(row_number() == 2)


# Use follow-up visit as analysis dataset

df_long_cogn_long <- df_long_cogn_long_l



############################################################
# Calculate follow-up time
############################################################

# Add baseline visit date

df_long_cogn_long$baseline_visitdate_ <-
  df_long_cogn_long$visitdate_e


# Follow-up duration in years

df_long_cogn_long$fup_cogn_time <-
  as.numeric(
    difftime(
      df_long_cogn_long$visitdate_,
      df_long_cogn_long$baseline_visitdate_,
      units = "days"
    )
  ) / 365.25


hist(
  df_long_cogn_long$fup_cogn_time,
  breaks = 30
)



############################################################
# Add baseline cognition measures
############################################################

# Correct follow-up cognition for baseline cognition

df_long_cogn_long$g_factor_baseline <-
  as.numeric(df_long_cogn_long_f$g_factor)

df_long_cogn_long$LDST_std_baseline <-
  as.numeric(df_long_cogn_long_f$LDST_std)

df_long_cogn_long$PPBsum_std_baseline <-
  as.numeric(df_long_cogn_long_f$PPBsum_std)

df_long_cogn_long$STR3_adjusted_std_baseline <-
  as.numeric(df_long_cogn_long_f$STR3_adjusted_std)

df_long_cogn_long$WFT_std_baseline <-
  as.numeric(df_long_cogn_long_f$WFT_std)

df_long_cogn_long$WLTdel_std_baseline <-
  as.numeric(df_long_cogn_long_f$WLTdel_std)



############################################################
# Identify participants without follow-up
############################################################

sum(df_long_cogn_long$fup_cogn_time == 0)

df_long_cogn_long$nofup <-
  ifelse(
    df_long_cogn_long$fup_cogn_time == 0,
    1,
    0
  )


# g-factor as continuous variable

df_long_cogn_long$g_factor <-
  as.numeric(df_long_cogn_long$g_factor)



############################################################
# Baseline characteristics
############################################################

library(tableone)

CreateTableOne(
  c(
    "age_scan",
    "sex.x",
    "bmi_imp",
    "ht_imp",
    "sbp_imp",
    "dbp_imp",
    "htdrug_imp",
    "smoke_imp",
    "oh_imp",
    "chol_imp",
    "hdl_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present",
    "g_factor",
    "LDST_",
    "PPBsum_",
    "STR3_adjusted_",
    "WFT_",
    "WLTdel_",
    "PSMD_std",
    "MD_std",
    "FA_std"
  ),
  factorVars = c(
    "sex.x",
    "ht_imp",
    "htdrug_imp",
    "smoke_imp",
    "lip_e_imp",
    "education_imp",
    "DM_scan_imp",
    "apoe_imp",
    "microbleeds_present",
    "lacunar_infarcts_present"
  ),
  data = df_long_cogn_long
)


median(df_long_cogn_long$oh_imp)
quantile(df_long_cogn_long$oh_imp)



############################################################
# Exclude participants without follow-up
############################################################

df_long_cogn_long <-
  df_long_cogn_long[
    df_long_cogn_long$nofup == 0,
  ]

hist(
  df_long_cogn_long$fup_cogn_time
)

mean(df_long_cogn_long$fup_cogn_time)



############################################################
# Calculate cognitive decline
############################################################

df_long_cogn_long$g_factor_diff <- NA


for (i in df_long_cogn_long$ergoid) {
  
  df_long_cogn_long$g_factor_diff[
    which(df_long_cogn_long$ergoid == i)
  ] <-
    df_long_cogn_long$g_factor[
      which(df_long_cogn_long$ergoid == i)
    ] -
    df_long_cogn_f$g_factor[
      which(df_long_cogn$ergoid == i)
    ]
}



############################################################
# WMH longitudinal model
############################################################

wmh_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

wmh_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  wmh_age,
  wmh_age_2
)



lm_PSMD_G_long_1 <- lm(
  g_factor ~ wmh_std +
    fup_cogn_time +
    g_factor_baseline,
  data = df_long_cogn_long
)


lm_PSMD_G_long_2 <- lm(
  g_factor ~ wmh_std +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_PSMD_G_long_3 <- lm(
  g_factor ~ wmh_std +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacun​ar_infarcts_present),
  data = df_long_cogn_long
)


AIC_psmd_G_long_1 <- AIC(lm_PSMD_G_long_1)
AIC_psmd_G_long_2 <- AIC(lm_PSMD_G_long_2)
AIC_psmd_G_long_3 <- AIC(lm_PSMD_G_long_3)


summary(lm_PSMD_G_long_3)
confint(lm_PSMD_G_long_3)
AIC_psmd_G_long_3



############################################################
# PSMD longitudinal model
############################################################

psmd_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

psmd_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  psmd_age,
  psmd_age_2
)


lm_PSMD_G_long_1 <- lm(
  g_factor ~ PSMD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline,
  data = df_long_cogn_long
)


lm_PSMD_G_long_2 <- lm(
  g_factor ~ PSMD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_PSMD_G_long_3 <- lm(
  g_factor ~ PSMD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


lm_PSMD_log_1 <- lm(
  g_factor ~ PSMD_log_std,
  data = df_long_cogn_long
)


lm_PSMD_log_2 <- lm(
  g_factor ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_PSMD_log_3 <- lm(
  g_factor ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_psmd_G_long_1 <- AIC(lm_PSMD_G_long_1)
AIC_psmd_G_long_2 <- AIC(lm_PSMD_G_long_2)
AIC_psmd_G_long_3 <- AIC(lm_PSMD_G_long_3)


summary(lm_PSMD_G_long_3)
confint(lm_PSMD_G_long_3)
AIC_psmd_G_long_3



############################################################
# MD longitudinal model
############################################################

md_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

md_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  md_age,
  md_age_2
)


lm_MD_G_long_1 <- lm(
  g_factor ~ MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline,
  data = df_long_cogn_long
)


lm_MD_G_long_2 <- lm(
  g_factor ~ MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_MD_G_long_3 <- lm(
  g_factor ~ MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_md_G_long_1 <- AIC(lm_MD_G_long_1)
AIC_md_G_long_2 <- AIC(lm_MD_G_long_2)
AIC_md_G_long_3 <- AIC(lm_MD_G_long_3)


summary(lm_MD_G_long_3)
confint(lm_MD_G_long_3)
AIC_md_G_long_3








############################################################
# FA longitudinal model
############################################################

fa_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

fa_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  fa_age,
  fa_age_2
)


lm_FA_G_long_1 <- lm(
  g_factor ~ FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline,
  data = df_long_cogn_long
)


lm_FA_G_long_2 <- lm(
  g_factor ~ FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_FA_G_long_3 <- lm(
  g_factor ~ FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_fa_G_long_1 <- AIC(lm_FA_G_long_1)
AIC_fa_G_long_2 <- AIC(lm_FA_G_long_2)
AIC_fa_G_long_3 <- AIC(lm_FA_G_long_3)


summary(lm_FA_G_long_3)
confint(lm_FA_G_long_3)
AIC_fa_G_long_3



############################################################
# IFO MD longitudinal model
############################################################

ifo_md_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

ifo_md_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  ifo_md_age,
  ifo_md_age_2
)


lm_ifo_MD_G_long_1 <- lm(
  g_factor ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline,
  data = df_long_cogn_long
)


lm_ifo_MD_G_long_2 <- lm(
  g_factor ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_ifo_MD_G_long_3 <- lm(
  g_factor ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_ifo_md_G_long_1 <- AIC(lm_ifo_MD_G_long_1)
AIC_ifo_md_G_long_2 <- AIC(lm_ifo_MD_G_long_2)
AIC_ifo_md_G_long_3 <- AIC(lm_ifo_MD_G_long_3)



summary(lm_ifo_MD_G_long_3)
confint(lm_ifo_MD_G_long_3)
AIC_ifo_md_G_long_3



############################################################
# IFO FA longitudinal model
############################################################

ifo_fa_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

ifo_fa_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  ifo_fa_age,
  ifo_fa_age_2
)


lm_ifo_FA_G_long_1 <- lm(
  g_factor ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline,
  data = df_long_cogn_long
)


lm_ifo_FA_G_long_2 <- lm(
  g_factor ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_long
)


lm_ifo_FA_G_long_3 <- lm(
  g_factor ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_ifo_fa_G_long_1 <- AIC(lm_ifo_FA_G_long_1)
AIC_ifo_fa_G_long_2 <- AIC(lm_ifo_FA_G_long_2)
AIC_ifo_fa_G_long_3 <- AIC(lm_ifo_FA_G_long_3)



summary(lm_ifo_FA_G_long_3)
confint(lm_ifo_FA_G_long_3)
AIC_ifo_fa_G_long_3



############################################################
# PTR MD longitudinal model
############################################################

ptr_md_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

ptr_md_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  ptr_md_age,
  ptr_md_age_2
)


lm_ptr_MD_G_long_3 <- lm(
  g_factor ~ ptr_MD_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_ptr_md_G_long_3 <- AIC(lm_ptr_MD_G_long_3)

AIC_ptr_md_G_long_3



############################################################
# PTR FA longitudinal model
############################################################

ptr_fa_age <- lm(
  g_factor ~ age_scan,
  data = df_long_cogn_long
)

ptr_fa_age_2 <- lm(
  g_factor ~ ns(age_scan, 2),
  data = df_long_cogn_long
)

anova(
  ptr_fa_age,
  ptr_fa_age_2
)


lm_ptr_FA_G_long_3 <- lm(
  g_factor ~ ptr_FA_std +
    pedir +
    fup_cogn_time +
    g_factor_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_long
)


AIC_ptr_fa_G_long_3 <- AIC(lm_ptr_FA_G_long_3)

AIC_ptr_fa_G_long_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  # Residuals versus fitted values
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(h = 0, col = "red")
  
  
  # Independence of residuals
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # Homoscedasticity
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(h = 0, col = "red")
  
  
  # Normality residuals
  qqnorm(residuals(model))
  qqline(residuals(model), col = "red")
  
  
  # Shapiro-Wilk test
  shapiro_test <- shapiro.test(residuals(model))
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # Multicollinearity
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_G_long_1)
check_assumptions(lm_FA_G_long_2)
check_assumptions(lm_FA_G_long_3)

check_assumptions(lm_MD_G_long_1)
check_assumptions(lm_MD_G_long_2)
check_assumptions(lm_MD_G_long_3)

check_assumptions(lm_PSMD_G_long_1)
check_assumptions(lm_PSMD_G_long_2)
check_assumptions(lm_PSMD_G_long_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

summary(lm_FA_G_long_1)
confint(lm_FA_G_long_1)
AIC(lm_FA_G_long_1)

summary(lm_FA_G_long_2)
confint(lm_FA_G_long_2)
AIC(lm_FA_G_long_2)

summary(lm_FA_G_long_3)
confint(lm_FA_G_long_3)
AIC(lm_FA_G_long_3)


summary(lm_MD_G_long_1)
confint(lm_MD_G_long_1)
AIC(lm_MD_G_long_1)

summary(lm_MD_G_long_2)
confint(lm_MD_G_long_2)
AIC(lm_MD_G_long_2)

summary(lm_MD_G_long_3)
confint(lm_MD_G_long_3)
AIC(lm_MD_G_long_3)


summary(lm_PSMD_G_long_1)
confint(lm_PSMD_G_long_1)
AIC(lm_PSMD_G_long_1)

summary(lm_PSMD_G_long_2)
confint(lm_PSMD_G_long_2)
AIC(lm_PSMD_G_long_2)

summary(lm_PSMD_G_long_3)
confint(lm_PSMD_G_long_3)
AIC(lm_PSMD_G_long_3)


summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)














############################################################
#### Separate cognitive tests
#### Cross-sectional analysis
#### Outcome: LDST
############################################################


############################################################
# PSMD analysis
############################################################

# Distribution of PSMD

hist(
  df_long_cogn_f$PSMD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$PSMD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$PSMD_std
)


# Log-transform PSMD

df_long_cogn_f$PSMD_log <- log(df_long_cogn_f$PSMD)

df_long_cogn_f$PSMD_log_std <- (
  df_long_cogn_f$PSMD_log - mean(df_long_cogn_f$PSMD_log)
) / sd(df_long_cogn_f$PSMD_log)



# Test linear versus nonlinear age association

psmd_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

psmd_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  psmd_age,
  psmd_age_2
)



# PSMD models

lm_PSMD_LDST_1 <- lm(
  LDST_std ~ PSMD_std + pedir,
  data = df_long_cogn_f
)


lm_PSMD_LDST_2 <- lm(
  LDST_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_LDST_3 <- lm(
  LDST_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# Log-transformed PSMD models

lm_PSMD_log_1 <- lm(
  LDST_std ~ PSMD_log_std,
  data = df_long_cogn_f
)


lm_PSMD_log_2 <- lm(
  LDST_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_log_3 <- lm(
  LDST_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn_f
)



# AIC comparison

AIC_psmd_LDST_1 <- AIC(lm_PSMD_LDST_1)
AIC_psmd_LDST_2 <- AIC(lm_PSMD_LDST_2)
AIC_psmd_LDST_3 <- AIC(lm_PSMD_LDST_3)



############################################################
# MD analysis
############################################################

# Distribution of MD

hist(
  df_long_cogn_f$MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# MD models

lm_MD_LDST_1 <- lm(
  LDST_std ~ MD_std + pedir,
  data = df_long_cogn_f
)


lm_MD_LDST_2 <- lm(
  LDST_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_MD_LDST_3 <- lm(
  LDST_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


# AIC comparison

AIC_md_LDST_1 <- AIC(lm_MD_LDST_1)
AIC_md_LDST_2 <- AIC(lm_MD_LDST_2)
AIC_md_LDST_3 <- AIC(lm_MD_LDST_3)



############################################################
# FA analysis
############################################################

# Distribution of FA

hist(
  df_long_cogn_f$FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# FA models

lm_FA_LDST_1 <- lm(
  LDST_std ~ FA_std + pedir,
  data = df_long_cogn_f
)


lm_FA_LDST_2 <- lm(
  LDST_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_FA_LDST_3 <- lm(
  LDST_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_fa_LDST_1 <- AIC(lm_FA_LDST_1)
AIC_fa_LDST_2 <- AIC(lm_FA_LDST_2)
AIC_fa_LDST_3 <- AIC(lm_FA_LDST_3)



############################################################
# IFO MD analysis
############################################################

# Distribution of IFO MD

hist(
  df_long_cogn_f$ifo_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ifo_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# IFO MD models

lm_ifo_MD_LDST_1 <- lm(
  LDST_std ~ ifo_MD_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_MD_LDST_2 <- lm(
  LDST_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_MD_LDST_3 <- lm(
  LDST_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


# AIC comparison

AIC_ifo_md_LDST_1 <- AIC(lm_ifo_MD_LDST_1)
AIC_ifo_md_LDST_2 <- AIC(lm_ifo_MD_LDST_2)
AIC_ifo_md_LDST_3 <- AIC(lm_ifo_MD_LDST_3)




############################################################
# IFO FA analysis
############################################################

# Distribution of IFO FA

hist(
  df_long_cogn_f$ifo_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# IFO FA models

lm_ifo_FA_LDST_1 <- lm(
  LDST_std ~ ifo_FA_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_FA_LDST_2 <- lm(
  LDST_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_FA_LDST_3 <- lm(
  LDST_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


# AIC comparison

AIC_ifo_fa_LDST_1 <- AIC(lm_ifo_FA_LDST_1)
AIC_ifo_fa_LDST_2 <- AIC(lm_ifo_FA_LDST_2)
AIC_ifo_fa_LDST_3 <- AIC(lm_ifo_FA_LDST_3)



############################################################
# PTR MD analysis
############################################################

# Distribution of PTR MD

hist(
  df_long_cogn_f$ptr_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ptr_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# PTR MD model

lm_ptr_MD_LDST_3 <- lm(
  LDST_std ~ ptr_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_md_LDST_3 <- AIC(lm_ptr_MD_LDST_3)

AIC_ptr_md_LDST_3



############################################################
# PTR FA analysis
############################################################

# Distribution of PTR FA

hist(
  df_long_cogn_f$ptr_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  LDST_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  LDST_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# PTR FA model

lm_ptr_FA_LDST_3 <- lm(
  LDST_std ~ ptr_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_fa_LDST_3 <- AIC(lm_ptr_FA_LDST_3)

AIC_ptr_fa_LDST_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # Residuals versus fitted values
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # Independence of residuals
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  
  # Homoscedasticity
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  
  # Normality of residuals
  qqnorm(
    residuals(model)
  )
  
  qqline(
    residuals(model),
    col = "red"
  )
  
  
  
  # Shapiro-Wilk test
  shapiro_test <- shapiro.test(
    residuals(model)
  )
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  
  # Multicollinearity
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_LDST_1)
check_assumptions(lm_FA_LDST_2)
check_assumptions(lm_FA_LDST_3)

check_assumptions(lm_MD_LDST_1)
check_assumptions(lm_MD_LDST_2)
check_assumptions(lm_MD_LDST_3)

check_assumptions(lm_PSMD_LDST_1)
check_assumptions(lm_PSMD_LDST_2)
check_assumptions(lm_PSMD_LDST_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

# NOTE:
# In the original script these refer to lm_FA_1, lm_MD_1,
# lm_PSMD_1, etc. These objects are not created in this block.
# The lines are kept unchanged to preserve the original workflow.

summary(lm_FA_1)
confint(lm_FA_1)
AIC(lm_FA_1)

summary(lm_FA_2)
confint(lm_FA_2)
AIC(lm_FA_2)

summary(lm_FA_3)
confint(lm_FA_3)
AIC(lm_FA_3)


summary(lm_MD_1)
confint(lm_MD_1)
AIC(lm_MD_1)

summary(lm_MD_2)
confint(lm_MD_2)
AIC(lm_MD_2)

summary(lm_MD_3)
confint(lm_MD_3)
AIC(lm_MD_3)


summary(lm_PSMD_1)
confint(lm_PSMD_1)
AIC(lm_PSMD_1)

summary(lm_PSMD_2)
confint(lm_PSMD_2)
AIC(lm_PSMD_2)

summary(lm_PSMD_3)
confint(lm_PSMD_3)
AIC(lm_PSMD_3)


summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)





############################################################
#### Separate cognitive tests
#### Cross-sectional analysis
#### Outcome: PPBsum
############################################################


############################################################
# PSMD analysis
############################################################

# Inspect PSMD distribution

hist(
  df_long_cogn_f$PSMD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$PSMD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$PSMD_std
)



# Log-transform PSMD

df_long_cogn_f$PSMD_log <- log(df_long_cogn_f$PSMD)

df_long_cogn_f$PSMD_log_std <- (
  df_long_cogn_f$PSMD_log -
    mean(df_long_cogn_f$PSMD_log)
) /
  sd(df_long_cogn_f$PSMD_log)



# Test linear versus nonlinear age association

psmd_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

psmd_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  psmd_age,
  psmd_age_2
)



# PSMD models

lm_PSMD_PPB_1 <- lm(
  PPBsum_std ~ PSMD_std + pedir,
  data = df_long_cogn_f
)


lm_PSMD_PPB_2 <- lm(
  PPBsum_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_PPB_3 <- lm(
  PPBsum_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# Log-transformed PSMD models

lm_PSMD_log_1 <- lm(
  PPBsum_std ~ PSMD_log_std,
  data = df_long_cogn_f
)


lm_PSMD_log_2 <- lm(
  PPBsum_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_log_3 <- lm(
  PPBsum_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn_f
)



# AIC comparison

AIC_psmd_PPB_1 <- AIC(lm_PSMD_PPB_1)
AIC_psmd_PPB_2 <- AIC(lm_PSMD_PPB_2)
AIC_psmd_PPB_3 <- AIC(lm_PSMD_PPB_3)



############################################################
# MD analysis
############################################################

# Inspect MD distribution

hist(
  df_long_cogn_f$MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# MD models

lm_MD_PPB_1 <- lm(
  PPBsum_std ~ MD_std + pedir,
  data = df_long_cogn_f
)


lm_MD_PPB_2 <- lm(
  PPBsum_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_MD_PPB_3 <- lm(
  PPBsum_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_md_PPB_1 <- AIC(lm_MD_PPB_1)
AIC_md_PPB_2 <- AIC(lm_MD_PPB_2)
AIC_md_PPB_3 <- AIC(lm_MD_PPB_3)



############################################################
# FA analysis
############################################################

# Inspect FA distribution

hist(
  df_long_cogn_f$FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# FA models

lm_FA_PPB_1 <- lm(
  PPBsum_std ~ FA_std + pedir,
  data = df_long_cogn_f
)


lm_FA_PPB_2 <- lm(
  PPBsum_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_FA_PPB_3 <- lm(
  PPBsum_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_fa_PPB_1 <- AIC(lm_FA_PPB_1)
AIC_fa_PPB_2 <- AIC(lm_FA_PPB_2)
AIC_fa_PPB_3 <- AIC(lm_FA_PPB_3)



############################################################
# IFO MD analysis
############################################################

# Inspect IFO MD distribution

hist(
  df_long_cogn_f$ifo_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ifo_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# IFO MD models

lm_ifo_MD_PPB_1 <- lm(
  PPBsum_std ~ ifo_MD_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_MD_PPB_2 <- lm(
  PPBsum_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_MD_PPB_3 <- lm(
  PPBsum_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_md_PPB_1 <- AIC(lm_ifo_MD_PPB_1)
AIC_ifo_md_PPB_2 <- AIC(lm_ifo_MD_PPB_2)
AIC_ifo_md_PPB_3 <- AIC(lm_ifo_MD_PPB_3)



############################################################
# IFO FA analysis
############################################################

# Inspect IFO FA distribution

hist(
  df_long_cogn_f$ifo_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# IFO FA models

lm_ifo_FA_PPB_1 <- lm(
  PPBsum_std ~ ifo_FA_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_FA_PPB_2 <- lm(
  PPBsum_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_FA_PPB_3 <- lm(
  PPBsum_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_fa_PPB_1 <- AIC(lm_ifo_FA_PPB_1)
AIC_ifo_fa_PPB_2 <- AIC(lm_ifo_FA_PPB_2)
AIC_ifo_fa_PPB_3 <- AIC(lm_ifo_FA_PPB_3)



############################################################
# PTR MD analysis
############################################################

# Inspect PTR MD distribution

hist(
  df_long_cogn_f$ptr_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ptr_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# PTR MD model

lm_ptr_MD_PPB_3 <- lm(
  PPBsum_std ~ ptr_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_md_PPB_3 <- AIC(lm_ptr_MD_PPB_3)

AIC_ptr_md_PPB_3



############################################################
# PTR FA analysis
############################################################

# Inspect PTR FA distribution

hist(
  df_long_cogn_f$ptr_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  PPBsum_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  PPBsum_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# PTR FA model

lm_ptr_FA_PPB_3 <- lm(
  PPBsum_std ~ ptr_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunars_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_fa_PPB_3 <- AIC(lm_ptr_FA_PPB_3)

AIC_ptr_fa_PPB_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # Linearity: residuals versus fitted values
  
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  
  # Independence: Durbin-Watson test
  
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  
  # Homoscedasticity: Scale-location plot
  
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  
  # Normality of residuals
  
  qqnorm(
    residuals(model)
  )
  
  qqline(
    residuals(model),
    col = "red"
  )
  
  
  
  # Shapiro-Wilk test
  
  shapiro_test <- shapiro.test(
    residuals(model)
  )
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  
  # Multicollinearity: VIF
  
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_PPB_1)
check_assumptions(lm_FA_PPB_2)
check_assumptions(lm_FA_PPB_3)

check_assumptions(lm_MD_PPB_1)
check_assumptions(lm_MD_PPB_2)
check_assumptions(lm_MD_PPB_3)

check_assumptions(lm_PSMD_PPB_1)
check_assumptions(lm_PSMD_PPB_2)
check_assumptions(lm_PSMD_PPB_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
#### Separate cognitive tests
#### Cross-sectional analysis
#### Outcome: Stroop (STR3_adjusted_std)
############################################################


############################################################
# PSMD analysis
############################################################

# Inspect PSMD distribution

hist(
  df_long_cogn_f$PSMD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$PSMD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$PSMD_std
)



# Log-transform PSMD

df_long_cogn_f$PSMD_log <- log(df_long_cogn_f$PSMD)

df_long_cogn_f$PSMD_log_std <- (
  df_long_cogn_f$PSMD_log -
    mean(df_long_cogn_f$PSMD_log)
) /
  sd(df_long_cogn_f$PSMD_log)



# Test linear versus nonlinear age association

psmd_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

psmd_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  psmd_age,
  psmd_age_2
)



# PSMD models

lm_PSMD_S_1 <- lm(
  STR3_adjusted_std ~ PSMD_std + pedir,
  data = df_long_cogn_f
)


lm_PSMD_S_2 <- lm(
  STR3_adjusted_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_S_3 <- lm(
  STR3_adjusted_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# Log-transformed PSMD models

lm_PSMD_log_1 <- lm(
  STR3_adjusted_std ~ PSMD_log_std,
  data = df_long_cogn_f
)


lm_PSMD_log_2 <- lm(
  STR3_adjusted_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_log_3 <- lm(
  STR3_adjusted_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn_f
)



# AIC comparison

AIC_psmd_S_1 <- AIC(lm_PSMD_S_1)
AIC_psmd_S_2 <- AIC(lm_PSMD_S_2)
AIC_psmd_S_3 <- AIC(lm_PSMD_S_3)



############################################################
# MD analysis
############################################################

# Inspect MD distribution

hist(
  df_long_cogn_f$MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# MD models

lm_MD_S_1 <- lm(
  STR3_adjusted_std ~ MD_std + pedir,
  data = df_long_cogn_f
)


lm_MD_S_2 <- lm(
  STR3_adjusted_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_MD_S_3 <- lm(
  STR3_adjusted_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_md_S_1 <- AIC(lm_MD_S_1)
AIC_md_S_2 <- AIC(lm_MD_S_2)
AIC_md_S_3 <- AIC(lm_MD_S_3)



############################################################
# FA analysis
############################################################

# Inspect FA distribution

hist(
  df_long_cogn_f$FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# FA models

lm_FA_S_1 <- lm(
  STR3_adjusted_std ~ FA_std + pedir,
  data = df_long_cogn_f
)


lm_FA_S_2 <- lm(
  STR3_adjusted_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_FA_S_3 <- lm(
  STR3_adjusted_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_fa_S_1 <- AIC(lm_FA_S_1)
AIC_fa_S_2 <- AIC(lm_FA_S_2)
AIC_fa_S_3 <- AIC(lm_FA_S_3)



############################################################
# IFO MD analysis
############################################################

# Inspect IFO MD distribution

hist(
  df_long_cogn_f$ifo_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ifo_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# IFO MD models

lm_ifo_MD_S_1 <- lm(
  STR3_adjusted_std ~ ifo_MD_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_MD_S_2 <- lm(
  STR3_adjusted_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_MD_S_3 <- lm(
  STR3_adjusted_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_md_S_1 <- AIC(lm_ifo_MD_S_1)
AIC_ifo_md_S_2 <- AIC(lm_ifo_MD_S_2)
AIC_ifo_md_S_3 <- AIC(lm_ifo_MD_S_3)



############################################################
# IFO FA analysis
############################################################

# Inspect IFO FA distribution

hist(
  df_long_cogn_f$ifo_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# IFO FA models

lm_ifo_FA_S_1 <- lm(
  STR3_adjusted_std ~ ifo_FA_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_FA_S_2 <- lm(
  STR3_adjusted_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_FA_S_3 <- lm(
  STR3_adjusted_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_fa_S_1 <- AIC(lm_ifo_FA_S_1)
AIC_ifo_fa_S_2 <- AIC(lm_ifo_FA_S_2)
AIC_ifo_fa_S_3 <- AIC(lm_ifo_FA_S_3)



############################################################
# PTR MD analysis
############################################################

# Inspect PTR MD distribution

hist(
  df_long_cogn_f$ptr_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ptr_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)



# PTR MD model

lm_ptr_MD_S_3 <- lm(
  STR3_adjusted_std ~ ptr_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_md_S_3 <- AIC(lm_ptr_MD_S_3)

AIC_ptr_md_S_3



############################################################
# PTR FA analysis
############################################################

# Inspect PTR FA distribution

hist(
  df_long_cogn_f$ptr_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  STR3_adjusted_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  STR3_adjusted_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)



# PTR FA model

lm_ptr_FA_S_3 <- lm(
  STR3_adjusted_std ~ ptr_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_fa_S_3 <- AIC(lm_ptr_FA_S_3)

AIC_ptr_fa_S_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # 1. Linearity:
  # residuals versus fitted values
  
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # 2. Independence:
  # Durbin-Watson test
  
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # 3. Homoscedasticity:
  # Scale-location plot
  
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # 4. Normality of residuals:
  # Q-Q plot
  
  qqnorm(
    residuals(model)
  )
  
  qqline(
    residuals(model),
    col = "red"
  )
  
  
  # Shapiro-Wilk test
  
  shapiro_test <- shapiro.test(
    residuals(model)
  )
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # 5. Multicollinearity:
  # Variance Inflation Factor
  
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_S_1)
check_assumptions(lm_FA_S_2)
check_assumptions(lm_FA_S_3)

check_assumptions(lm_MD_S_1)
check_assumptions(lm_MD_S_2)
check_assumptions(lm_MD_S_3)

check_assumptions(lm_PSMD_S_1)
check_assumptions(lm_PSMD_S_2)
check_assumptions(lm_PSMD_S_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

# Note:
# Original script refers to lm_FA_1, lm_MD_1 and lm_PSMD_1.
# These objects are not created in this Stroop block.
# They are therefore kept unchanged to preserve original workflow.

summary(lm_FA_1)
confint(lm_FA_1)
AIC(lm_FA_1)

summary(lm_FA_2)
confint(lm_FA_2)
AIC(lm_FA_2)

summary(lm_FA_3)
confint(lm_FA_3)
AIC(lm_FA_3)


summary(lm_MD_1)
confint(lm_MD_1)
AIC(lm_MD_1)

summary(lm_MD_2)
confint(lm_MD_2)
AIC(lm_MD_2)

summary(lm_MD_3)
confint(lm_MD_3)
AIC(lm_MD_3)


summary(lm_PSMD_1)
confint(lm_PSMD_1)
AIC(lm_PSMD_1)

summary(lm_PSMD_2)
confint(lm_PSMD_2)
AIC(lm_PSMD_2)

summary(lm_PSMD_3)
confint(lm_PSMD_3)
AIC(lm_PSMD_3)


summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)




############################################################
#### Separate cognitive tests
#### Cross-sectional analysis
#### Outcome: Word Fluency Test (WFT)
############################################################


############################################################
# PSMD analysis
############################################################

# Inspect PSMD distribution

hist(
  df_long_cogn_f$PSMD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$PSMD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$PSMD_std
)



# Log-transform PSMD

df_long_cogn_f$PSMD_log <- log(df_long_cogn_f$PSMD)

df_long_cogn_f$PSMD_log_std <- (
  df_long_cogn_f$PSMD_log -
    mean(df_long_cogn_f$PSMD_log)
) /
  sd(df_long_cogn_f$PSMD_log)



# Test linear versus nonlinear age association

psmd_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

psmd_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  psmd_age,
  psmd_age_2
)
# Geen verschil: GEEN splines gebruiken in verdere modellen



# PSMD models

lm_PSMD_WFT_1 <- lm(
  WFT_std ~ PSMD_std + pedir,
  data = df_long_cogn_f
)


lm_PSMD_WFT_2 <- lm(
  WFT_std ~ PSMD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_WFT_3 <- lm(
  WFT_std ~ PSMD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# Log-transformed PSMD models

lm_PSMD_log_1 <- lm(
  WFT_std ~ PSMD_log_std,
  data = df_long_cogn_f
)


lm_PSMD_log_2 <- lm(
  WFT_std ~ PSMD_log_std +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_log_3 <- lm(
  WFT_std ~ PSMD_log_std +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunarinfarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn_f
)



# AIC comparison

AIC_psmd_WFT_1 <- AIC(lm_PSMD_WFT_1)
AIC_psmd_WFT_2 <- AIC(lm_PSMD_WFT_2)
AIC_psmd_WFT_3 <- AIC(lm_PSMD_WFT_3)



############################################################
# MD analysis
############################################################

# Inspect MD distribution

hist(
  df_long_cogn_f$MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)
# Geen verschil: GEEN splines gebruiken



# MD models

lm_MD_WFT_1 <- lm(
  WFT_std ~ MD_std + pedir,
  data = df_long_cogn_f
)


lm_MD_WFT_2 <- lm(
  WFT_std ~ MD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_MD_WFT_3 <- lm(
  WFT_std ~ MD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_md_WFT_1 <- AIC(lm_MD_WFT_1)
AIC_md_WFT_2 <- AIC(lm_MD_WFT_2)
AIC_md_WFT_3 <- AIC(lm_MD_WFT_3)



############################################################
# FA analysis
############################################################

# Inspect FA distribution

hist(
  df_long_cogn_f$FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)
# Geen verschil: GEEN splines gebruiken



# FA models

lm_FA_WFT_1 <- lm(
  WFT_std ~ FA_std + pedir,
  data = df_long_cogn_f
)


lm_FA_WFT_2 <- lm(
  WFT_std ~ FA_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_FA_WFT_3 <- lm(
  WFT_std ~ FA_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_fa_WFT_1 <- AIC(lm_FA_WFT_1)
AIC_fa_WFT_2 <- AIC(lm_FA_WFT_2)
AIC_fa_WFT_3 <- AIC(lm_FA_WFT_3)


############################################################
# IFO MD analysis
############################################################

# Inspect IFO MD distribution

hist(
  df_long_cogn_f$ifo_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ifo_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)
# Geen verschil: GEEN splines gebruiken



# IFO MD models

lm_ifo_MD_WFT_1 <- lm(
  WFT_std ~ ifo_MD_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_MD_WFT_2 <- lm(
  WFT_std ~ ifo_MD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_MD_WFT_3 <- lm(
  WFT_std ~ ifo_MD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_md_WFT_1 <- AIC(lm_ifo_MD_WFT_1)
AIC_ifo_md_WFT_2 <- AIC(lm_ifo_MD_WFT_2)
AIC_ifo_md_WFT_3 <- AIC(lm_ifo_MD_WFT_3)



############################################################
# IFO FA analysis
############################################################

# Inspect IFO FA distribution

hist(
  df_long_cogn_f$ifo_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)
# Geen verschil: GEEN splines gebruiken



# IFO FA models

lm_ifo_FA_WFT_1 <- lm(
  WFT_std ~ ifo_FA_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_FA_WFT_2 <- lm(
  WFT_std ~ ifo_FA_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_FA_WFT_3 <- lm(
  WFT_std ~ ifo_FA_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_fa_WFT_1 <- AIC(lm_ifo_FA_WFT_1)
AIC_ifo_fa_WFT_2 <- AIC(lm_ifo_FA_WFT_2)
AIC_ifo_fa_WFT_3 <- AIC(lm_ifo_FA_WFT_3)



############################################################
# PTR MD analysis
############################################################

# Inspect PTR MD distribution

hist(
  df_long_cogn_f$ptr_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ptr_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)
# Geen verschil: GEEN splines gebruiken



# PTR MD model

lm_ptr_MD_WFT_3 <- lm(
  WFT_std ~ ptr_MD_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_md_WFT_3 <- AIC(lm_ptr_MD_WFT_3)

AIC_ptr_md_WFT_3



############################################################
# PTR FA analysis
############################################################

# Inspect PTR FA distribution

hist(
  df_long_cogn_f$ptr_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  WFT_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  WFT_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)
# Geen verschil: GEEN splines gebruiken



# PTR FA model

lm_ptr_FA_WFT_3 <- lm(
  WFT_std ~ ptr_FA_std +
    pedir +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunarinfarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_fa_WFT_3 <- AIC(lm_ptr_FA_WFT_3)

AIC_ptr_fa_WFT_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # Linearity: residuals versus fitted values
  
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # Independence: Durbin-Watson test
  
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # Homoscedasticity: Scale-location plot
  
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # Normality of residuals: Q-Q plot
  
  qqnorm(
    residuals(model)
  )
  
  qqline(
    residuals(model),
    col = "red"
  )
  
  
  # Shapiro-Wilk test
  
  shapiro_test <- shapiro.test(
    residuals(model)
  )
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # Multicollinearity: VIF
  
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_WFT_1)
check_assumptions(lm_FA_WFT_2)
check_assumptions(lm_FA_WFT_3)

check_assumptions(lm_MD_WFT_1)
check_assumptions(lm_MD_WFT_2)
check_assumptions(lm_MD_WFT_3)

check_assumptions(lm_PSMD_WFT_1)
check_assumptions(lm_PSMD_WFT_2)
check_assumptions(lm_PSMD_WFT_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

summary(lm_FA_1)
confint(lm_FA_1)
AIC(lm_FA_1)

summary(lm_FA_2)
confint(lm_FA_2)
AIC(lm_FA_2)

summary(lm_FA_3)
confint(lm_FA_3)
AIC(lm_FA_3)



summary(lm_MD_1)
confint(lm_MD_1)
AIC(lm_MD_1)

summary(lm_MD_2)
confint(lm_MD_2)
AIC(lm_MD_2)

summary(lm_MD_3)
confint(lm_MD_3)
AIC(lm_MD_3)



summary(lm_PSMD_1)
confint(lm_PSMD_1)
AIC(lm_PSMD_1)

summary(lm_PSMD_2)
confint(lm_PSMD_2)
AIC(lm_PSMD_2)

summary(lm_PSMD_3)
confint(lm_PSMD_3)
AIC(lm_PSMD_3)



summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)



############################################################
#### Separate cognitive tests
#### Cross-sectional analysis
#### Outcome: Word Learning Test delayed recall (WLTdel)
############################################################


############################################################
# PSMD analysis
############################################################

# Inspect PSMD distribution

hist(
  df_long_cogn_f$PSMD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$PSMD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$PSMD_std
)



# Log-transform PSMD

df_long_cogn_f$PSMD_log <- log(df_long_cogn_f$PSMD)

df_long_cogn_f$PSMD_log_std <- (
  df_long_cogn_f$PSMD_log -
    mean(df_long_cogn_f$PSMD_log)
) /
  sd(df_long_cogn_f$PSMD_log)



# Test linear versus nonlinear age association

psmd_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

psmd_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  psmd_age,
  psmd_age_2
)



# PSMD models

lm_PSMD_WLT_1 <- lm(
  WLTdel_std ~ PSMD_std + pedir,
  data = df_long_cogn_f
)


lm_PSMD_WLT_2 <- lm(
  WLTdel_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_WLT_3 <- lm(
  WLTdel_std ~ PSMD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# Log-transformed PSMD models

lm_PSMD_log_1 <- lm(
  WLTdel_std ~ PSMD_log_std,
  data = df_long_cogn_f
)


lm_PSMD_log_2 <- lm(
  WLTdel_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_PSMD_log_3 <- lm(
  WLTdel_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn_f
)



# AIC comparison

AIC_psmd_WLT_1 <- AIC(lm_PSMD_WLT_1)
AIC_psmd_WLT_2 <- AIC(lm_PSMD_WLT_2)
AIC_psmd_WLT_3 <- AIC(lm_PSMD_WLT_3)



############################################################
# MD analysis
############################################################

# Inspect MD distribution

hist(
  df_long_cogn_f$MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)
# Geen verschil: GEEN splines gebruiken



# MD models

lm_MD_WLT_1 <- lm(
  WLTdel_std ~ MD_std + pedir,
  data = df_long_cogn_f
)


lm_MD_WLT_2 <- lm(
  WLTdel_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_MD_WLT_3 <- lm(
  WLTdel_std ~ MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_md_WLT_1 <- AIC(lm_MD_WLT_1)
AIC_md_WLT_2 <- AIC(lm_MD_WLT_2)
AIC_md_WLT_3 <- AIC(lm_MD_WLT_3)



############################################################
# FA analysis
############################################################

# Inspect FA distribution

hist(
  df_long_cogn_f$FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)
# Geen verschil: GEEN splines gebruiken



# FA models

lm_FA_WLT_1 <- lm(
  WLTdel_std ~ FA_std + pedir,
  data = df_long_cogn_f
)


lm_FA_WLT_2 <- lm(
  WLTdel_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_FA_WLT_3 <- lm(
  WLTdel_std ~ FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_fa_WLT_1 <- AIC(lm_FA_WLT_1)
AIC_fa_WLT_2 <- AIC(lm_FA_WLT_2)
AIC_fa_WLT_3 <- AIC(lm_FA_WLT_3)



############################################################
# IFO MD analysis
############################################################

# Inspect IFO MD distribution

hist(
  df_long_cogn_f$ifo_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ifo_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)
# Geen verschil: GEEN splines gebruiken



# IFO MD models

lm_ifo_MD_WLT_1 <- lm(
  WLTdel_std ~ ifo_MD_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_MD_WLT_2 <- lm(
  WLTdel_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_MD_WLT_3 <- lm(
  WLTdel_std ~ ifo_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_md_WLT_1 <- AIC(lm_ifo_MD_WLT_1)
AIC_ifo_md_WLT_2 <- AIC(lm_ifo_MD_WLT_2)
AIC_ifo_md_WLT_3 <- AIC(lm_ifo_MD_WLT_3)



############################################################
# IFO FA analysis
############################################################

# Inspect IFO FA distribution

hist(
  df_long_cogn_f$ifo_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ifo_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)
# Geen verschil: GEEN splines gebruiken



# IFO FA models

lm_ifo_FA_WLT_1 <- lm(
  WLTdel_std ~ ifo_FA_std + pedir,
  data = df_long_cogn_f
)


lm_ifo_FA_WLT_2 <- lm(
  WLTdel_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)


lm_ifo_FA_WLT_3 <- lm(
  WLTdel_std ~ ifo_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)



# AIC comparison

AIC_ifo_fa_WLT_1 <- AIC(lm_ifo_FA_WLT_1)
AIC_ifo_fa_WLT_2 <- AIC(lm_ifo_FA_WLT_2)
AIC_ifo_fa_WLT_3 <- AIC(lm_ifo_FA_WLT_3)



############################################################
# PTR MD analysis
############################################################

# Inspect PTR MD distribution

hist(
  df_long_cogn_f$ptr_MD_std,
  prob = TRUE
)

lines(
  density(df_long_cogn_f$ptr_MD_std),
  col = "red",
  lwd = 2
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_MD_std
)



# Test linear versus nonlinear age association

md_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

md_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  md_age,
  md_age_2
)
# Geen verschil: GEEN splines gebruiken



# PTR MD model

lm_ptr_MD_WLT_3 <- lm(
  WLTdel_std ~ ptr_MD_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_md_WLT_3 <- AIC(lm_ptr_MD_WLT_3)

AIC_ptr_md_WLT_3



############################################################
# PTR FA analysis
############################################################

# Inspect PTR FA distribution

hist(
  df_long_cogn_f$ptr_FA_std,
  breaks = 50
)

plot(
  df_long_cogn_f$age_scan,
  df_long_cogn_f$ptr_FA_std
)



# Test linear versus nonlinear age association

fa_age <- lm(
  WLTdel_std ~ age_scan,
  data = df_long_cogn_f
)

fa_age_2 <- lm(
  WLTdel_std ~ ns(age_scan, 2),
  data = df_long_cogn_f
)

anova(
  fa_age,
  fa_age_2
)
# Geen verschil: GEEN splines gebruiken



# PTR FA model

lm_ptr_FA_WLT_3 <- lm(
  WLTdel_std ~ ptr_FA_std +
    pedir +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


AIC_ptr_fa_WLT_3 <- AIC(lm_ptr_FA_WLT_3)



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # Linearity: residuals versus fitted values
  
  plot(
    model$fitted.values,
    residuals(model),
    main = "Residuals vs Fitted",
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # Independence: Durbin-Watson test
  
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # Homoscedasticity: Scale-location plot
  
  plot(
    model$fitted.values,
    sqrt(abs(residuals(model))),
    main = "Scale-Location",
    xlab = "Fitted values",
    ylab = "Square Root of |Residuals|"
  )
  
  abline(
    h = 0,
    col = "red"
  )
  
  
  # Normality of residuals
  
  qqnorm(
    residuals(model)
  )
  
  qqline(
    residuals(model),
    col = "red"
  )
  
  
  # Shapiro-Wilk test
  
  shapiro_test <- shapiro.test(
    residuals(model)
  )
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # Multicollinearity (VIF)
  
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Run assumption checks
############################################################

check_assumptions(lm_FA_WLT_1)
check_assumptions(lm_FA_WLT_2)
check_assumptions(lm_FA_WLT_3)

check_assumptions(lm_MD_WLT_1)
check_assumptions(lm_MD_WLT_2)
check_assumptions(lm_MD_WLT_3)

check_assumptions(lm_PSMD_WLT_1)
check_assumptions(lm_PSMD_WLT_2)
check_assumptions(lm_PSMD_WLT_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

summary(lm_FA_1)
confint(lm_FA_1)
AIC(lm_FA_1)

summary(lm_FA_2)
confint(lm_FA_2)
AIC(lm_FA_2)

summary(lm_FA_3)
confint(lm_FA_3)
AIC(lm_FA_3)



summary(lm_MD_1)
confint(lm_MD_1)
AIC(lm_MD_1)

summary(lm_MD_2)
confint(lm_MD_2)
AIC(lm_MD_2)

summary(lm_MD_3)
confint(lm_MD_3)
AIC(lm_MD_3)



summary(lm_PSMD_1)
confint(lm_PSMD_1)
AIC(lm_PSMD_1)

summary(lm_PSMD_2)
confint(lm_PSMD_2)
AIC(lm_PSMD_2)

summary(lm_PSMD_3)
confint(lm_PSMD_3)
AIC(lm_PSMD_3)



summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)


















############################################################
#### Longitudinal analyses - LDST
############################################################

# Dataset for longitudinal analyses
df_long_cogn <- df_long_cogn_long


############################################################
# PSMD
############################################################

# Distribution and relation with age
hist(df_long_cogn$PSMD_std, prob = TRUE,
     main = "Distribution of standardized PSMD",
     xlab = "PSMD (standardized)")
lines(density(df_long_cogn$PSMD_std), col = "red", lwd = 2)

plot(df_long_cogn$age_scan, df_long_cogn$PSMD_std,
     xlab = "Age at scan",
     ylab = "PSMD (standardized)")


# Log transformation of PSMD
df_long_cogn$PSMD_log <- log(df_long_cogn$PSMD)

df_long_cogn$PSMD_log_std <- 
  (df_long_cogn$PSMD_log - mean(df_long_cogn$PSMD_log)) /
  sd(df_long_cogn$PSMD_log)


# Test linear vs nonlinear age effect
psmd_age <- lm(LDST_std ~ age_scan,
               data = df_long_cogn)

psmd_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
                 data = df_long_cogn)

anova(psmd_age, psmd_age_2)


# Linear models PSMD - LDST
lm_PSMD_LDST_long_1 <- lm(
  LDST_std ~ PSMD_std + pedir +
    fup_cogn_time +
    LDST_std_baseline,
  data = df_long_cogn
)

lm_PSMD_LDST_long_2 <- lm(
  LDST_std ~ PSMD_std + pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)

lm_PSMD_LDST_long_3 <- lm(
  LDST_std ~ PSMD_std + pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


# Log-transformed PSMD models
lm_PSMD_log_1 <- lm(
  LDST_std ~ PSMD_log_std,
  data = df_long_cogn_f
)

lm_PSMD_log_2 <- lm(
  LDST_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_PSMD_log_3 <- lm(
  LDST_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunarinfarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn_f
)


# Model comparison using AIC
AIC_psmd_LDST_long_1 <- AIC(lm_PSMD_LDST_long_1)
AIC_psmd_LDST_long_2 <- AIC(lm_PSMD_LDST_long_2)
AIC_psmd_LDST_long_3 <- AIC(lm_PSMD_LDST_long_3)



############################################################
# MD
############################################################

hist(df_long_cogn$MD_std, prob = TRUE,
     main = "Distribution of standardized MD",
     xlab = "MD (standardized)")
lines(density(df_long_cogn$MD_std), col = "red", lwd = 2)

plot(df_long_cogn$age_scan, df_long_cogn$MD_std,
     xlab = "Age at scan",
     ylab = "MD (standardized)")


md_age <- lm(LDST_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)


lm_MD_LDST_long_1 <- lm(
  LDST_std ~ MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline,
  data = df_long_cogn
)

lm_MD_LDST_long_2 <- lm(
  LDST_std ~ MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)

lm_MD_LDST_long_3 <- lm(
  LDST_std ~ MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_md_LDST_long_1 <- AIC(lm_MD_LDST_long_1)
AIC_md_LDST_long_2 <- AIC(lm_MD_LDST_long_2)
AIC_md_LDST_long_3 <- AIC(lm_MD_LDST_long_3)



############################################################
# FA
############################################################

hist(df_long_cogn$FA_std, breaks = 50,
     main = "Distribution of standardized FA",
     xlab = "FA (standardized)")

plot(df_long_cogn$age_scan, df_long_cogn$FA_std,
     xlab = "Age at scan",
     ylab = "FA (standardized)")


fa_age <- lm(LDST_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)


lm_FA_LDST_long_1 <- lm(
  LDST_std ~ FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline,
  data = df_long_cogn
)

lm_FA_LDST_long_2 <- lm(
  LDST_std ~ FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)

lm_FA_LDST_long_3 <- lm(
  LDST_std ~ FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_fa_LDST_long_1 <- AIC(lm_FA_LDST_long_1)
AIC_fa_LDST_long_2 <- AIC(lm_FA_LDST_long_2)
AIC_fa_LDST_long_3 <- AIC(lm_FA_LDST_long_3)



############################################################
# ifo MD
############################################################

hist(df_long_cogn$ifo_MD_std, prob = TRUE,
     main = "Distribution of standardized ifo MD",
     xlab = "ifo MD (standardized)")
lines(density(df_long_cogn$ifo_MD_std), col = "red", lwd = 2)

plot(df_long_cogn$age_scan, df_long_cogn$ifo_MD_std,
     xlab = "Age at scan",
     ylab = "ifo MD (standardized)")


md_age <- lm(LDST_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)


lm_ifo_MD_LDST_long_1 <- lm(
  LDST_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline,
  data = df_long_cogn
)

lm_ifo_MD_LDST_long_2 <- lm(
  LDST_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)

lm_ifo_MD_LDST_long_3 <- lm(
  LDST_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ifo_md_LDST_long_1 <- AIC(lm_ifo_MD_LDST_long_1)
AIC_ifo_md_LDST_long_2 <- AIC(lm_ifo_MD_LDST_long_2)
AIC_ifo_md_LDST_long_3 <- AIC(lm_ifo_MD_LDST_long_3)



############################################################
# ifo FA
############################################################

hist(df_long_cogn$ifo_FA_std, breaks = 50,
     main = "Distribution of standardized ifo FA",
     xlab = "ifo FA (standardized)")

plot(df_long_cogn$age_scan, df_long_cogn$ifo_FA_std,
     xlab = "Age at scan",
     ylab = "ifo FA (standardized)")


fa_age <- lm(LDST_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)


lm_ifo_FA_LDST_long_1 <- lm(
  LDST_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline,
  data = df_long_cogn
)

lm_ifo_FA_LDST_long_2 <- lm(
  LDST_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)

lm_ifo_FA_LDST_long_3 <- lm(
  LDST_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ifo_fa_LDST_long_1 <- AIC(lm_ifo_FA_LDST_long_1)
AIC_ifo_fa_LDST_long_2 <- AIC(lm_ifo_FA_LDST_long_2)
AIC_ifo_fa_LDST_long_3 <- AIC(lm_ifo_FA_LDST_long_3)



############################################################
# ptr MD
############################################################

hist(df_long_cogn$ptr_MD_std, prob = TRUE,
     main = "Distribution of standardized ptr MD",
     xlab = "ptr MD (standardized)")
lines(density(df_long_cogn$ptr_MD_std), col = "red", lwd = 2)

plot(df_long_cogn$age_scan, df_long_cogn$ptr_MD_std,
     xlab = "Age at scan",
     ylab = "ptr MD (standardized)")


md_age <- lm(LDST_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)


lm_ptr_MD_LDST_long_3 <- lm(
  LDST_std ~ ptr_MD_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ptr_md_LDST_long_3 <- AIC(lm_ptr_MD_LDST_long_3)
AIC_ptr_md_LDST_long_3



############################################################
# ptr FA
############################################################

hist(df_long_cogn$ptr_FA_std, breaks = 50,
     main = "Distribution of standardized ptr FA",
     xlab = "ptr FA (standardized)")

plot(df_long_cogn$age_scan, df_long_cogn$ptr_FA_std,
     xlab = "Age at scan",
     ylab = "ptr FA (standardized)")


fa_age <- lm(LDST_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(LDST_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)


lm_ptr_FA_LDST_long_3 <- lm(
  LDST_std ~ ptr_FA_std +
    pedir +
    fup_cogn_time +
    LDST_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ptr_fa_LDST_long_3 <- AIC(lm_ptr_FA_LDST_long_3)
AIC_ptr_fa_LDST_long_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)

check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  # Residuals vs fitted values: linearity
  plot(model$fitted.values,
       residuals(model),
       main = "Residuals vs Fitted",
       xlab = "Fitted values",
       ylab = "Residuals")
  abline(h = 0, col = "red")
  
  
  # Durbin-Watson test: independence
  dw_test <- durbinWatsonTest(model)
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # Scale-location plot: homoscedasticity
  plot(model$fitted.values,
       sqrt(abs(residuals(model))),
       main = "Scale-Location",
       xlab = "Fitted values",
       ylab = "Square Root of |Residuals|")
  abline(h = 0, col = "red")
  
  
  # Q-Q plot: normality of residuals
  qqnorm(residuals(model))
  qqline(residuals(model), col = "red")
  
  
  # Shapiro-Wilk test
  shapiro_test <- shapiro.test(residuals(model))
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # Variance Inflation Factor
  vif_values <- vif(model)
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Check assumptions for models
############################################################

check_assumptions(lm_FA_1)
check_assumptions(lm_FA_2)
check_assumptions(lm_FA_3)

check_assumptions(lm_MD_1)
check_assumptions(lm_MD_2)
check_assumptions(lm_MD_3)

check_assumptions(lm_PSMD_1)
check_assumptions(lm_PSMD_2)
check_assumptions(lm_PSMD_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

summary(lm_FA_1)
confint(lm_FA_1)
AIC(lm_FA_1)

summary(lm_FA_2)
confint(lm_FA_2)
AIC(lm_FA_2)

summary(lm_FA_3)
confint(lm_FA_3)
AIC(lm_FA_3)


summary(lm_MD_1)
confint(lm_MD_1)
AIC(lm_MD_1)

summary(lm_MD_2)
confint(lm_MD_2)
AIC(lm_MD_2)

summary(lm_MD_3)
confint(lm_MD_3)
AIC(lm_MD_3)


summary(lm_PSMD_1)
confint(lm_PSMD_1)
AIC(lm_PSMD_1)

summary(lm_PSMD_2)
confint(lm_PSMD_2)
AIC(lm_PSMD_2)

summary(lm_PSMD_3)
confint(lm_PSMD_3)
AIC(lm_PSMD_3)


summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)


############################################################
#### Longitudinal analyses - PPBsum
############################################################


############################################################
# PSMD
############################################################

# Distribution and relation with age
hist(df_long_cogn$PSMD_std, prob = TRUE,
     main = "Distribution of standardized PSMD",
     xlab = "PSMD (standardized)")

lines(density(df_long_cogn$PSMD_std),
      col = "red",
      lwd = 2)

plot(df_long_cogn$age_scan,
     df_long_cogn$PSMD_std,
     xlab = "Age at scan",
     ylab = "PSMD (standardized)")


# Log transformation of PSMD
df_long_cogn$PSMD_log <- log(df_long_cogn$PSMD)

df_long_cogn$PSMD_log_std <-
  (df_long_cogn$PSMD_log - mean(df_long_cogn$PSMD_log)) /
  sd(df_long_cogn$PSMD_log)


# Test linearity of age effect
psmd_age <- lm(PPBsum_std ~ age_scan,
               data = df_long_cogn)

psmd_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
                 data = df_long_cogn)

anova(psmd_age, psmd_age_2)


# PSMD models
lm_PSMD_PPB_long_1 <- lm(
  PPBsum_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline,
  data = df_long_cogn
)


lm_PSMD_PPB_long_2 <- lm(
  PPBsum_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_PPB_long_3 <- lm(
  PPBsum_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


# Log-transformed PSMD models
lm_PSMD_log_1 <- lm(
  PPBsum_std ~ PSMD_log_std,
  data = df_long_cogn_f
)

lm_PSMD_log_2 <- lm(
  PPBsum_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn_f
)

lm_PSMD_log_3 <- lm(
  PPBsum_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn_f
)


# AIC comparison
AIC_psmd_PPB_long_1 <- AIC(lm_PSMD_PPB_long_1)
AIC_psmd_PPB_long_2 <- AIC(lm_PSMD_PPB_long_2)
AIC_psmd_PPB_long_3 <- AIC(lm_PSMD_PPB_long_3)



############################################################
# MD
############################################################

hist(df_long_cogn$MD_std,
     prob = TRUE,
     main = "Distribution of standardized MD",
     xlab = "MD (standardized)")

lines(density(df_long_cogn$MD_std),
      col = "red",
      lwd = 2)

plot(df_long_cogn$age_scan,
     df_long_cogn$MD_std,
     xlab = "Age at scan",
     ylab = "MD (standardized)")


md_age <- lm(PPBsum_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)


lm_MD_PPB_long_1 <- lm(
  PPBsum_std ~ MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline,
  data = df_long_cogn
)


lm_MD_PPB_long_2 <- lm(
  PPBsum_std ~ MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_MD_PPB_long_3 <- lm(
  PPBsum_std ~ MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_md_PPB_long_1 <- AIC(lm_MD_PPB_long_1)
AIC_md_PPB_long_2 <- AIC(lm_MD_PPB_long_2)
AIC_md_PPB_long_3 <- AIC(lm_MD_PPB_long_3)



############################################################
# FA
############################################################

hist(df_long_cogn$FA_std,
     breaks = 50,
     main = "Distribution of standardized FA",
     xlab = "FA (standardized)")

plot(df_long_cogn$age_scan,
     df_long_cogn$FA_std,
     xlab = "Age at scan",
     ylab = "FA (standardized)")


fa_age <- lm(PPBsum_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)


lm_FA_PPB_long_1 <- lm(
  PPBsum_std ~ FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline,
  data = df_long_cogn
)


lm_FA_PPB_long_2 <- lm(
  PPBsum_std ~ FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_FA_PPB_long_3 <- lm(
  PPBsum_std ~ FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_fa_PPB_long_1 <- AIC(lm_FA_PPB_long_1)
AIC_fa_PPB_long_2 <- AIC(lm_FA_PPB_long_2)
AIC_fa_PPB_long_3 <- AIC(lm_FA_PPB_long_3)


############################################################
# ifo MD
############################################################

hist(df_long_cogn$ifo_MD_std,
     prob = TRUE,
     main = "Distribution of standardized ifo MD",
     xlab = "ifo MD (standardized)")

lines(density(df_long_cogn$ifo_MD_std),
      col = "red",
      lwd = 2)

plot(df_long_cogn$age_scan,
     df_long_cogn$ifo_MD_std,
     xlab = "Age at scan",
     ylab = "ifo MD (standardized)")


md_age <- lm(PPBsum_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)


lm_ifo_MD_PPB_long_1 <- lm(
  PPBsum_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline,
  data = df_long_cogn
)


lm_ifo_MD_PPB_long_2 <- lm(
  PPBsum_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_ifo_MD_PPB_long_3 <- lm(
  PPBsum_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ifo_md_PPB_long_1 <- AIC(lm_ifo_MD_PPB_long_1)
AIC_ifo_md_PPB_long_2 <- AIC(lm_ifo_MD_PPB_long_2)
AIC_ifo_md_PPB_long_3 <- AIC(lm_ifo_MD_PPB_long_3)



############################################################
# ifo FA
############################################################

hist(df_long_cogn$ifo_FA_std,
     breaks = 50,
     main = "Distribution of standardized ifo FA",
     xlab = "ifo FA (standardized)")

plot(df_long_cogn$age_scan,
     df_long_cogn$ifo_FA_std,
     xlab = "Age at scan",
     ylab = "ifo FA (standardized)")


fa_age <- lm(PPBsum_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)


lm_ifo_FA_PPB_long_1 <- lm(
  PPBsum_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline,
  data = df_long_cogn
)


lm_ifo_FA_PPB_long_2 <- lm(
  PPBsum_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_ifo_FA_PPB_long_3 <- lm(
  PPBsum_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ifo_fa_PPB_long_1 <- AIC(lm_ifo_FA_PPB_long_1)
AIC_ifo_fa_PPB_long_2 <- AIC(lm_ifo_FA_PPB_long_2)
AIC_ifo_fa_PPB_long_3 <- AIC(lm_ifo_FA_PPB_long_3)



############################################################
# ptr MD
############################################################

hist(df_long_cogn$ptr_MD_std,
     prob = TRUE,
     main = "Distribution of standardized ptr MD",
     xlab = "ptr MD (standardized)")

lines(density(df_long_cogn$ptr_MD_std),
      col = "red",
      lwd = 2)

plot(df_long_cogn$age_scan,
     df_long_cogn$ptr_MD_std,
     xlab = "Age at scan",
     ylab = "ptr MD (standardized)")


md_age <- lm(PPBsum_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)


lm_ptr_MD_PPB_long_3 <- lm(
  PPBsum_std ~ ptr_MD_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)


AIC_ptr_md_PPB_long_3 <- AIC(lm_ptr_MD_PPB_long_3)
AIC_ptr_md_PPB_long_3



############################################################
# ptr FA
############################################################

hist(df_long_cogn$ptr_FA_std,
     breaks = 50,
     main = "Distribution of standardized ptr FA",
     xlab = "ptr FA (standardized)")

plot(df_long_cogn$age_scan,
     df_long_cogn$ptr_FA_std,
     xlab = "Age at scan",
     ylab = "ptr FA (standardized)")


fa_age <- lm(PPBsum_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(PPBsum_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)


lm_ptr_FA_PPB_long_3 <- lm(
  PPBsum_std ~ ptr_FA_std +
    pedir +
    fup_cogn_time +
    PPBsum_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunarinfarcts_present),
  data = df_long_cogn
)


AIC_ptr_fa_PPB_long_3 <- AIC(lm_ptr_FA_PPB_long_3)
AIC_ptr_fa_PPB_long_3



############################################################
# Model assumption checks
############################################################

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # Residuals vs fitted values: linearity
  plot(model$fitted.values,
       residuals(model),
       main = "Residuals vs Fitted",
       xlab = "Fitted values",
       ylab = "Residuals")
  
  abline(h = 0, col = "red")
  
  
  # Durbin-Watson test: independence
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # Scale-location plot: homoscedasticity
  plot(model$fitted.values,
       sqrt(abs(residuals(model))),
       main = "Scale-Location",
       xlab = "Fitted values",
       ylab = "Square Root of |Residuals|")
  
  abline(h = 0, col = "red")
  
  
  # Q-Q plot: normality of residuals
  qqnorm(residuals(model))
  qqline(residuals(model), col = "red")
  
  
  # Shapiro-Wilk test
  shapiro_test <- shapiro.test(residuals(model))
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # Multicollinearity
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1, 1))
}



############################################################
# Check assumptions
############################################################

check_assumptions(lm_FA_1)
check_assumptions(lm_FA_2)
check_assumptions(lm_FA_3)

check_assumptions(lm_MD_1)
check_assumptions(lm_MD_2)
check_assumptions(lm_MD_3)

check_assumptions(lm_PSMD_1)
check_assumptions(lm_PSMD_2)
check_assumptions(lm_PSMD_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)



############################################################
# Model summaries
############################################################

summary(lm_FA_1)
confint(lm_FA_1)
AIC(lm_FA_1)

summary(lm_FA_2)
confint(lm_FA_2)
AIC(lm_FA_2)

summary(lm_FA_3)
confint(lm_FA_3)
AIC(lm_FA_3)


summary(lm_MD_1)
confint(lm_MD_1)
AIC(lm_MD_1)

summary(lm_MD_2)
confint(lm_MD_2)
AIC(lm_MD_2)

summary(lm_MD_3)
confint(lm_MD_3)
AIC(lm_MD_3)


summary(lm_PSMD_1)
confint(lm_PSMD_1)
AIC(lm_PSMD_1)

summary(lm_PSMD_2)
confint(lm_PSMD_2)
AIC(lm_PSMD_2)

summary(lm_PSMD_3)
confint(lm_PSMD_3)
AIC(lm_PSMD_3)


summary(lm_PSMD_log_1)
confint(lm_PSMD_log_1)

summary(lm_PSMD_log_2)
confint(lm_PSMD_log_2)

summary(lm_PSMD_log_3)
confint(lm_PSMD_log_3)




########## Stroop - longitudinal analyses


## PSMD

hist(df_long_cogn$PSMD_std, prob = TRUE)
lines(density(df_long_cogn$PSMD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$PSMD_std)


# log-transform PSMD
df_long_cogn$PSMD_log <- log(df_long_cogn$PSMD)

df_long_cogn$PSMD_log_std <- 
  (df_long_cogn$PSMD_log - mean(df_long_cogn$PSMD_log)) /
  sd(df_long_cogn$PSMD_log)


# age spline check
psmd_age <- lm(STR3_adjusted_std ~ age_scan, data = df_long_cogn)

psmd_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
                 data = df_long_cogn)

anova(psmd_age, psmd_age_2)



# Main models

lm_PSMD_S_long_1 <- lm(
  STR3_adjusted_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline,
  data = df_long_cogn
)


lm_PSMD_S_long_2 <- lm(
  STR3_adjusted_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_S_long_3 <- lm(
  STR3_adjusted_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



# log PSMD models

lm_PSMD_log_1 <- lm(
  STR3_adjusted_std ~ PSMD_log_std,
  data = df_long_cogn
)


lm_PSMD_log_2 <- lm(
  STR3_adjusted_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_log_3 <- lm(
  STR3_adjusted_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present) +
    sbp_imp +
    dbp_imp +
    as.factor(htdrug) +
    bmi_imp +
    as.factor(smoke_imp) +
    chol_imp +
    hdl_imp +
    as.factor(lip_e_imp) +
    oh_imp +
    as.factor(apoe_imp) +
    as.factor(DM_scan_imp),
  data = df_long_cogn
)



# AIC

AIC_psmd_S_long_1 <- AIC(lm_PSMD_S_long_1)
AIC_psmd_S_long_2 <- AIC(lm_PSMD_S_long_2)
AIC_psmd_S_long_3 <- AIC(lm_PSMD_S_long_3)



############################################################


## MD

hist(df_long_cogn$MD_std, prob = TRUE)
lines(density(df_long_cogn$MD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$MD_std)


md_age <- lm(STR3_adjusted_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)



lm_MD_S_long_1 <- lm(
  STR3_adjusted_std ~ MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline,
  data = df_long_cogn
)


lm_MD_S_long_2 <- lm(
  STR3_adjusted_std ~ MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_MD_S_long_3 <- lm(
  STR3_adjusted_std ~ MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_md_S_long_1 <- AIC(lm_MD_S_long_1)
AIC_md_S_long_2 <- AIC(lm_MD_S_long_2)
AIC_md_S_long_3 <- AIC(lm_MD_S_long_3)


############################################################


## FA

hist(df_long_cogn$FA_std, breaks = 50)
plot(df_long_cogn$age_scan, df_long_cogn$FA_std)


fa_age <- lm(STR3_adjusted_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)



lm_FA_S_long_1 <- lm(
  STR3_adjusted_std ~ FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline,
  data = df_long_cogn
)


lm_FA_S_long_2 <- lm(
  STR3_adjusted_std ~ FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_FA_S_long_3 <- lm(
  STR3_adjusted_std ~ FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_fa_S_long_1 <- AIC(lm_FA_S_long_1)
AIC_fa_S_long_2 <- AIC(lm_FA_S_long_2)
AIC_fa_S_long_3 <- AIC(lm_FA_S_long_3)



############################################################


## ifo MD

hist(df_long_cogn$ifo_MD_std, prob = TRUE)
lines(density(df_long_cogn$ifo_MD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$ifo_MD_std)


ifo_md_age <- lm(STR3_adjusted_std ~ age_scan,
                 data = df_long_cogn)

ifo_md_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
                   data = df_long_cogn)

anova(ifo_md_age, ifo_md_age_2)



lm_ifo_MD_S_long_1 <- lm(
  STR3_adjusted_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline,
  data = df_long_cogn
)


lm_ifo_MD_S_long_2 <- lm(
  STR3_adjusted_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_ifo_MD_S_long_3 <- lm(
  STR3_adjusted_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ifo_md_S_long_1 <- AIC(lm_ifo_MD_S_long_1)
AIC_ifo_md_S_long_2 <- AIC(lm_ifo_MD_S_long_2)
AIC_ifo_md_S_long_3 <- AIC(lm_ifo_MD_S_long_3)



############################################################


## ifo FA

hist(df_long_cogn$ifo_FA_std, breaks = 50)
plot(df_long_cogn$age_scan, df_long_cogn$ifo_FA_std)


ifo_fa_age <- lm(STR3_adjusted_std ~ age_scan,
                 data = df_long_cogn)

ifo_fa_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
                   data = df_long_cogn)

anova(ifo_fa_age, ifo_fa_age_2)



lm_ifo_FA_S_long_1 <- lm(
  STR3_adjusted_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline,
  data = df_long_cogn
)


lm_ifo_FA_S_long_2 <- lm(
  STR3_adjusted_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_ifo_FA_S_long_3 <- lm(
  STR3_adjusted_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ifo_fa_S_long_1 <- AIC(lm_ifo_FA_S_long_1)
AIC_ifo_fa_S_long_2 <- AIC(lm_ifo_FA_S_long_2)
AIC_ifo_fa_S_long_3 <- AIC(lm_ifo_FA_S_long_3)



############################################################


## ptr MD

hist(df_long_cogn$ptr_MD_std, prob = TRUE)
lines(density(df_long_cogn$ptr_MD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$ptr_MD_std)


ptr_md_age <- lm(STR3_adjusted_std ~ age_scan,
                 data = df_long_cogn)

ptr_md_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
                   data = df_long_cogn)

anova(ptr_md_age, ptr_md_age_2)



lm_ptr_MD_S_long_3 <- lm(
  STR3_adjusted_std ~ ptr_MD_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ptr_md_S_long_3 <- AIC(lm_ptr_MD_S_long_3)



############################################################


## ptr FA

hist(df_long_cogn$ptr_FA_std, breaks = 50)
plot(df_long_cogn$age_scan, df_long_cogn$ptr_FA_std)


ptr_fa_age <- lm(STR3_adjusted_std ~ age_scan,
                 data = df_long_cogn)

ptr_fa_age_2 <- lm(STR3_adjusted_std ~ ns(age_scan, 2),
                   data = df_long_cogn)

anova(ptr_fa_age, ptr_fa_age_2)



lm_ptr_FA_S_long_3 <- lm(
  STR3_adjusted_std ~ ptr_FA_std +
    pedir +
    fup_cogn_time +
    STR3_adjusted_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ptr_fa_S_long_3 <- AIC(lm_ptr_FA_S_long_3)



############################################################


# Assumptie functie

library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2,2))
  
  plot(model$fitted.values,
       residuals(model),
       main = "Residuals vs Fitted",
       xlab = "Fitted values",
       ylab = "Residuals")
  
  abline(h = 0, col = "red")
  
  
  dw_test <- durbinWatsonTest(model)
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  plot(model$fitted.values,
       sqrt(abs(residuals(model))),
       main = "Scale-Location",
       xlab = "Fitted values",
       ylab = "Square Root of |Residuals|")
  
  abline(h = 0, col = "red")
  
  
  qqnorm(residuals(model))
  qqline(residuals(model), col = "red")
  
  
  shapiro_test <- shapiro.test(residuals(model))
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  vif_values <- vif(model)
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1,1))
}



############################################################


# Check assumptions

check_assumptions(lm_FA_S_long_1)
check_assumptions(lm_FA_S_long_2)
check_assumptions(lm_FA_S_long_3)

check_assumptions(lm_MD_S_long_1)
check_assumptions(lm_MD_S_long_2)
check_assumptions(lm_MD_S_long_3)

check_assumptions(lm_PSMD_S_long_1)
check_assumptions(lm_PSMD_S_long_2)
check_assumptions(lm_PSMD_S_long_3)

check_assumptions(lm_PSMD_log_1)
check_assumptions(lm_PSMD_log_2)
check_assumptions(lm_PSMD_log_3)

check_assumptions(lm_ifo_MD_S_long_1)
check_assumptions(lm_ifo_MD_S_long_2)
check_assumptions(lm_ifo_MD_S_long_3)

check_assumptions(lm_ifo_FA_S_long_1)
check_assumptions(lm_ifo_FA_S_long_2)
check_assumptions(lm_ifo_FA_S_long_3)

check_assumptions(lm_ptr_MD_S_long_3)
check_assumptions(lm_ptr_FA_S_long_3)



############################################################


# Model summaries

summary(lm_PSMD_S_long_1)
confint(lm_PSMD_S_long_1)
AIC(lm_PSMD_S_long_1)

summary(lm_PSMD_S_long_2)
confint(lm_PSMD_S_long_2)
AIC(lm_PSMD_S_long_2)

summary(lm_PSMD_S_long_3)
confint(lm_PSMD_S_long_3)
AIC(lm_PSMD_S_long_3)



summary(lm_MD_S_long_1)
confint(lm_MD_S_long_1)
AIC(lm_MD_S_long_1)

summary(lm_MD_S_long_2)
confint(lm_MD_S_long_2)
AIC(lm_MD_S_long_2)

summary(lm_MD_S_long_3)
confint(lm_MD_S_long_3)
AIC(lm_MD_S_long_3)



summary(lm_FA_S_long_1)
confint(lm_FA_S_long_1)
AIC(lm_FA_S_long_1)

summary(lm_FA_S_long_2)
confint(lm_FA_S_long_2)
AIC(lm_FA_S_long_2)

summary(lm_FA_S_long_3)
confint(lm_FA_S_long_3)
AIC(lm_FA_S_long_3)




############# WFT - longitudinal analyses


## PSMD

hist(df_long_cogn$PSMD_std, prob = TRUE)
lines(density(df_long_cogn$PSMD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$PSMD_std)


# log-transform
df_long_cogn$PSMD_log <- log(df_long_cogn$PSMD)

df_long_cogn$PSMD_log_std <- 
  (df_long_cogn$PSMD_log - mean(df_long_cogn$PSMD_log)) /
  sd(df_long_cogn$PSMD_log)



# leeftijdsmodel check
psmd_age <- lm(WFT_std ~ age_scan,
               data = df_long_cogn)

psmd_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
                 data = df_long_cogn)

anova(psmd_age, psmd_age_2)
# geen verschil, dus geen splines gebruiken



## PSMD models

lm_PSMD_WFT_long_1 <- lm(
  WFT_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline,
  data = df_long_cogn
)


lm_PSMD_WFT_long_2 <- lm(
  WFT_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_WFT_long_3 <- lm(
  WFT_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



## PSMD log models

lm_PSMD_log_1 <- lm(
  WFT_std ~ PSMD_log_std,
  data = df_long_cogn
)


lm_PSMD_log_2 <- lm(
  WFT_std ~ PSMD_log_std +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_log_3 <- lm(
  WFT_std ~ PSMD_log_std +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



# AIC

AIC_psmd_WFT_long_1 <- AIC(lm_PSMD_WFT_long_1)
AIC_psmd_WFT_long_2 <- AIC(lm_PSMD_WFT_long_2)
AIC_psmd_WFT_long_3 <- AIC(lm_PSMD_WFT_long_3)





############################################################


## MD


hist(df_long_cogn$MD_std, prob = TRUE)
lines(density(df_long_cogn$MD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$MD_std)


md_age <- lm(WFT_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)
# geen verschil, dus geen splines gebruiken



lm_MD_WFT_long_1 <- lm(
  WFT_std ~ MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline,
  data = df_long_cogn
)


lm_MD_WFT_long_2 <- lm(
  WFT_std ~ MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_MD_WFT_long_3 <- lm(
  WFT_std ~ MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_md_WFT_long_1 <- AIC(lm_MD_WFT_long_1)
AIC_md_WFT_long_2 <- AIC(lm_MD_WFT_long_2)
AIC_md_WFT_long_3 <- AIC(lm_MD_WFT_long_3)





############################################################


## FA


hist(df_long_cogn$FA_std, breaks = 50)
plot(df_long_cogn$age_scan, df_long_cogn$FA_std)


fa_age <- lm(WFT_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)
# geen verschil, dus geen splines gebruiken



lm_FA_WFT_long_1 <- lm(
  WFT_std ~ FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline,
  data = df_long_cogn
)


lm_FA_WFT_long_2 <- lm(
  WFT_std ~ FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_FA_WFT_long_3 <- lm(
  WFT_std ~ FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_fa_WFT_long_1 <- AIC(lm_FA_WFT_long_1)
AIC_fa_WFT_long_2 <- AIC(lm_FA_WFT_long_2)
AIC_fa_WFT_long_3 <- AIC(lm_FA_WFT_long_3)



############################################################


## ifo MD


hist(df_long_cogn$ifo_MD_std, prob = TRUE)
lines(density(df_long_cogn$ifo_MD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$ifo_MD_std)


md_age <- lm(WFT_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)
# geen verschil, dus geen splines gebruiken



lm_ifo_MD_WFT_long_1 <- lm(
  WFT_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline,
  data = df_long_cogn
)


lm_ifo_MD_WFT_long_2 <- lm(
  WFT_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_ifo_MD_WFT_long_3 <- lm(
  WFT_std ~ ifo_MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ifo_md_WFT_long_1 <- AIC(lm_ifo_MD_WFT_long_1)
AIC_ifo_md_WFT_long_2 <- AIC(lm_ifo_MD_WFT_long_2)
AIC_ifo_md_WFT_long_3 <- AIC(lm_ifo_MD_WFT_long_3)





############################################################


## ifo FA


hist(df_long_cogn$ifo_FA_std, breaks = 50)
plot(df_long_cogn$age_scan, df_long_cogn$ifo_FA_std)


fa_age <- lm(WFT_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)
# geen verschil, dus geen splines gebruiken



lm_ifo_FA_WFT_long_1 <- lm(
  WFT_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline,
  data = df_long_cogn
)


lm_ifo_FA_WFT_long_2 <- lm(
  WFT_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_ifo_FA_WFT_long_3 <- lm(
  WFT_std ~ ifo_FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ifo_fa_WFT_long_1 <- AIC(lm_ifo_FA_WFT_long_1)
AIC_ifo_fa_WFT_long_2 <- AIC(lm_ifo_FA_WFT_long_2)
AIC_ifo_fa_WFT_long_3 <- AIC(lm_ifo_FA_WFT_long_3)





############################################################


## ptr MD


hist(df_long_cogn$ptr_MD_std, prob = TRUE)
lines(density(df_long_cogn$ptr_MD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$ptr_MD_std)


md_age <- lm(WFT_std ~ age_scan,
             data = df_long_cogn)

md_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(md_age, md_age_2)
# geen verschil, dus geen splines gebruiken



lm_ptr_MD_WFT_long_3 <- lm(
  WFT_std ~ ptr_MD_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ptr_md_WFT_long_3 <- AIC(lm_ptr_MD_WFT_long_3)
AIC_ptr_md_WFT_long_3





############################################################


## ptr FA


hist(df_long_cogn$ptr_FA_std, breaks = 50)
plot(df_long_cogn$age_scan, df_long_cogn$ptr_FA_std)


fa_age <- lm(WFT_std ~ age_scan,
             data = df_long_cogn)

fa_age_2 <- lm(WFT_std ~ ns(age_scan, 2),
               data = df_long_cogn)

anova(fa_age, fa_age_2)
# geen verschil, dus geen splines gebruiken



lm_ptr_FA_WFT_long_3 <- lm(
  WFT_std ~ ptr_FA_std +
    pedir +
    fup_cogn_time +
    WFT_std_baseline +
    age_scan +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_ptr_fa_WFT_long_3 <- AIC(lm_ptr_FA_WFT_long_3)
AIC_ptr_fa_WFT_long_3





############################################################
# Assumpties controleren


library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2,2))
  
  
  # Lineariteit
  plot(model$fitted.values,
       residuals(model),
       main = "Residuals vs Fitted",
       xlab = "Fitted values",
       ylab = "Residuals")
  
  abline(h = 0, col = "red")
  
  
  # Onafhankelijkheid
  dw_test <- durbinWatsonTest(model)
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # Homoscedasticiteit
  plot(model$fitted.values,
       sqrt(abs(residuals(model))),
       main = "Scale-Location",
       xlab = "Fitted values",
       ylab = "Square Root of |Residuals|")
  
  abline(h = 0, col = "red")
  
  
  # Normaliteit
  qqnorm(residuals(model))
  qqline(residuals(model), col = "red")
  
  
  shapiro_test <- shapiro.test(residuals(model))
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # Multicollineariteit
  vif_values <- vif(model)
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1,1))
}




############################################################
# check assumptions


check_assumptions(lm_FA_WFT_long_1)
check_assumptions(lm_FA_WFT_long_2)
check_assumptions(lm_FA_WFT_long_3)

check_assumptions(lm_MD_WFT_long_1)
check_assumptions(lm_MD_WFT_long_2)
check_assumptions(lm_MD_WFT_long_3)

check_assumptions(lm_PSMD_WFT_long_1)
check_assumptions(lm_PSMD_WFT_long_2)
check_assumptions(lm_PSMD_WFT_long_3)

check_assumptions(lm_ifo_MD_WFT_long_1)
check_assumptions(lm_ifo_MD_WFT_long_2)
check_assumptions(lm_ifo_MD_WFT_long_3)

check_assumptions(lm_ifo_FA_WFT_long_1)
check_assumptions(lm_ifo_FA_WFT_long_2)
check_assumptions(lm_ifo_FA_WFT_long_3)

check_assumptions(lm_ptr_MD_WFT_long_3)
check_assumptions(lm_ptr_FA_WFT_long_3)





############################################################
# Model summaries


summary(lm_FA_WFT_long_1)
confint(lm_FA_WFT_long_1)
AIC(lm_FA_WFT_long_1)

summary(lm_FA_WFT_long_2)
confint(lm_FA_WFT_long_2)
AIC(lm_FA_WFT_long_2)

summary(lm_FA_WFT_long_3)
confint(lm_FA_WFT_long_3)
AIC(lm_FA_WFT_long_3)



summary(lm_MD_WFT_long_1)
confint(lm_MD_WFT_long_1)
AIC(lm_MD_WFT_long_1)

summary(lm_MD_WFT_long_2)
confint(lm_MD_WFT_long_2)
AIC(lm_MD_WFT_long_2)

summary(lm_MD_WFT_long_3)
confint(lm_MD_WFT_long_3)
AIC(lm_MD_WFT_long_3)



summary(lm_PSMD_WFT_long_1)
confint(lm_PSMD_WFT_long_1)
AIC(lm_PSMD_WFT_long_1)

summary(lm_PSMD_WFT_long_2)
confint(lm_PSMD_WFT_long_2)
AIC(lm_PSMD_WFT_long_2)

summary(lm_PSMD_WFT_long_3)
confint(lm_PSMD_WFT_long_3)
AIC(lm_PSMD_WFT_long_3)



summary(lm_ifo_MD_WFT_long_3)
confint(lm_ifo_MD_WFT_long_3)
AIC(lm_ifo_MD_WFT_long_3)


summary(lm_ifo_FA_WFT_long_3)
confint(lm_ifo_FA_WFT_long_3)
AIC(lm_ifo_FA_WFT_long_3)


summary(lm_ptr_MD_WFT_long_3)
confint(lm_ptr_MD_WFT_long_3)
AIC(lm_ptr_MD_WFT_long_3)


summary(lm_ptr_FA_WFT_long_3)
confint(lm_ptr_FA_WFT_long_3)
AIC(lm_ptr_FA_WFT_long_3)




############# WLTdel


## PSMD

hist(df_long_cogn$PSMD_std, prob = TRUE)
lines(density(df_long_cogn$PSMD_std), col = "red", lwd = 2)
plot(df_long_cogn$age_scan, df_long_cogn$PSMD_std)


# log-transform
df_long_cogn$PSMD_log = log(df_long_cogn$PSMD)

df_long_cogn$PSMD_log_std = (
  df_long_cogn$PSMD_log - mean(df_long_cogn$PSMD_log)
) / sd(df_long_cogn$PSMD_log)



psmd_age = lm(WLTdel_std ~ age_scan,
              data = df_long_cogn)

psmd_age_2 = lm(WLTdel_std ~ ns(age_scan, 2),
                data = df_long_cogn)

anova(psmd_age, psmd_age_2)



lm_PSMD_WLT_long_1 = lm(
  WLTdel_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    WLTdel_std_baseline,
  data = df_long_cogn
)


lm_PSMD_WLT_long_2 = lm(
  WLTdel_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    WLTdel_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_WLT_long_3 = lm(
  WLTdel_std ~ PSMD_std +
    pedir +
    fup_cogn_time +
    WLTdel_std_baseline +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



lm_PSMD_log_WLT_long_1 = lm(
  WLTdel_std ~ PSMD_log_std,
  data = df_long_cogn
)


lm_PSMD_log_WLT_long_2 = lm(
  WLTdel_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp),
  data = df_long_cogn
)


lm_PSMD_log_WLT_long_3 = lm(
  WLTdel_std ~ PSMD_log_std +
    ns(age_scan, 2) +
    as.factor(sex.x) +
    as.factor(education_imp) +
    as.factor(microbleeds_present) +
    as.factor(lacunar_infarcts_present),
  data = df_long_cogn
)



AIC_psmd_WLT_long_1 = AIC(lm_PSMD_WLT_long_1)
AIC_psmd_WLT_long_2 = AIC(lm_PSMD_WLT_long_2)
AIC_psmd_WLT_long_3 = AIC(lm_PSMD_WLT_long_3)

############################################################
## Assumpties controleren


library(ggplot2)


check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  
  
  # 1. Lineariteit
  plot(model$fitted.values,
       residuals(model),
       main = "Residuals vs Fitted",
       xlab = "Fitted values",
       ylab = "Residuals")
  
  abline(h = 0, col = "red")
  
  
  # 2. Onafhankelijkheid
  dw_test <- durbinWatsonTest(model)
  
  cat("Durbin-Watson Test:\n")
  print(dw_test)
  
  
  # 3. Homoscedasticiteit
  plot(model$fitted.values,
       sqrt(abs(residuals(model))),
       main = "Scale-Location",
       xlab = "Fitted values",
       ylab = "Square Root of |Residuals|")
  
  abline(h = 0, col = "red")
  
  
  # 4. Normaliteit
  qqnorm(residuals(model))
  qqline(residuals(model), col = "red")
  
  
  shapiro_test <- shapiro.test(residuals(model))
  
  cat("Shapiro-Wilk Test:\n")
  print(shapiro_test)
  
  
  # 5. Multicollineariteit
  vif_values <- vif(model)
  
  cat("Variance Inflation Factor (VIF):\n")
  print(vif_values)
  
  
  par(mfrow = c(1,1))
}





############################################################
# Check assumptions


check_assumptions(lm_FA_WLT_long_1)
check_assumptions(lm_FA_WLT_long_2)
check_assumptions(lm_FA_WLT_long_3)


check_assumptions(lm_MD_WLT_long_1)
check_assumptions(lm_MD_WLT_long_2)
check_assumptions(lm_MD_WLT_long_3)


check_assumptions(lm_PSMD_WLT_long_1)
check_assumptions(lm_PSMD_WLT_long_2)
check_assumptions(lm_PSMD_WLT_long_3)


check_assumptions(lm_PSMD_log_WLT_long_1)
check_assumptions(lm_PSMD_log_WLT_long_2)
check_assumptions(lm_PSMD_log_WLT_long_3)


check_assumptions(lm_ifo_MD_WLT_long_1)
check_assumptions(lm_ifo_MD_WLT_long_2)
check_assumptions(lm_ifo_MD_WLT_long_3)


check_assumptions(lm_ifo_FA_WLT_long_1)
check_assumptions(lm_ifo_FA_WLT_long_2)
check_assumptions(lm_ifo_FA_WLT_long_3)


check_assumptions(lm_ptr_MD_WLT_long_3)
check_assumptions(lm_ptr_FA_WLT_long_3)






############################################################
# Model summaries


summary(lm_FA_WLT_long_1)
confint(lm_FA_WLT_long_1)
AIC(lm_FA_WLT_long_1)


summary(lm_FA_WLT_long_2)
confint(lm_FA_WLT_long_2)
AIC(lm_FA_WLT_long_2)


summary(lm_FA_WLT_long_3)
confint(lm_FA_WLT_long_3)
AIC(lm_FA_WLT_long_3)




summary(lm_MD_WLT_long_1)
confint(lm_MD_WLT_long_1)
AIC(lm_MD_WLT_long_1)


summary(lm_MD_WLT_long_2)
confint(lm_MD_WLT_long_2)
AIC(lm_MD_WLT_long_2)


summary(lm_MD_WLT_long_3)
confint(lm_MD_WLT_long_3)
AIC(lm_MD_WLT_long_3)





summary(lm_PSMD_WLT_long_1)
confint(lm_PSMD_WLT_long_1)
AIC(lm_PSMD_WLT_long_1)


summary(lm_PSMD_WLT_long_2)
confint(lm_PSMD_WLT_long_2)
AIC(lm_PSMD_WLT_long_2)


summary(lm_PSMD_WLT_long_3)
confint(lm_PSMD_WLT_long_3)
AIC(lm_PSMD_WLT_long_3)





summary(lm_PSMD_log_WLT_long_1)
confint(lm_PSMD_log_WLT_long_1)


summary(lm_PSMD_log_WLT_long_2)
confint(lm_PSMD_log_WLT_long_2)


summary(lm_PSMD_log_WLT_long_3)
confint(lm_PSMD_log_WLT_long_3)





summary(lm_ifo_MD_WLT_long_1)
confint(lm_ifo_MD_WLT_long_1)
AIC(lm_ifo_MD_WLT_long_1)


summary(lm_ifo_MD_WLT_long_2)
confint(lm_ifo_MD_WLT_long_2)
AIC(lm_ifo_MD_WLT_long_2)


summary(lm_ifo_MD_WLT_long_3)
confint(lm_ifo_MD_WLT_long_3)
AIC(lm_ifo_MD_WLT_long_3)





summary(lm_ifo_FA_WLT_long_1)
confint(lm_ifo_FA_WLT_long_1)
AIC(lm_ifo_FA_WLT_long_1)


summary(lm_ifo_FA_WLT_long_2)
confint(lm_ifo_FA_WLT_long_2)
AIC(lm_ifo_FA_WLT_long_2)


summary(lm_ifo_FA_WLT_long_3)
confint(lm_ifo_FA_WLT_long_3)
AIC(lm_ifo_FA_WLT_long_3)





summary(lm_ptr_MD_WLT_long_3)
confint(lm_ptr_MD_WLT_long_3)
AIC(lm_ptr_MD_WLT_long_3)


summary(lm_ptr_FA_WLT_long_3)
confint(lm_ptr_FA_WLT_long_3)
AIC(lm_ptr_FA_WLT_long_3)