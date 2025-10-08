

# -- Load packages
library(dplyr)
library(purrr)
library(stringr)
library(openxlsx)
library(progressr)
library(lme4)

# Formatter: if p < 2.2e-16, report as "<2.2e-16"
p_fmt <- function(x, machine = 2.2e-16) {
  out <- ifelse(x < machine,
                paste0("<", formatC(machine, format = "e", digits = 1)),
                formatC(signif(x, 2), format = "e"))
  out[is.na(x)] <- NA_character_
  out
}

# -- Set working directory
wkdir <- "/QRISdata/Q7280/pharmacogenomics"
output_dir <- "/QRISdata/Q7280/pharmacogenomics/pharma_summaries"

# -- Load data
pheno <- read.csv(file.path(wkdir, "phenotypes/survey_phenotypes.csv"))
pheno$timesincefirstpreg = pheno$AGE - pheno$PREG1AGE
pheno <- pheno %>%
  mutate(
    earliest_last_end = PREG1AGE + 2 *pmax(CHILDREN, 0),
    latest_last_end = AGE,
    preg_last5_deterministic = case_when(
      CHILDREN == 0 ~ 0,
      earliest_last_end >= AGE - 5 ~ 1,
      TRUE ~ NA
    ),
    BFEEDANY = case_when(
      BFEEDANY == 1 ~ 0,
      BFEEDANY == 2 | BFEEDANY == 3 ~ 1,
      TRUE ~ NA
    )
  )
dat <- read.csv(file.path(wkdir, "data/AGDSAcceptabilityTreatmentGroups_25082025.csv"))

# -- Drug mapping as tibble
source("/QRISdata/Q7280/pharmacogenomics/Drug_Reference/Drug_Reference_Table.R")
drug_mapping <- drug_ref %>%
  select(ATCCode, DrugName, DrugClass)

# -- Map ATC codes
ad_mapped <- dat %>%
  left_join(drug_mapping, by = "ATCCode")

# -- Create participant summary
pharma <- ad_mapped %>%
  group_by(ParticipantID) %>%
  summarize(
    num_ATC = n_distinct(ATCCode),
    num_class = n_distinct(DrugClass),
    Total_TEs = sum(NumberOfPrescriptionEpisodes, na.rm = TRUE),
    TotalPrescriptionDays = sum(PrescriptionDays, na.rm = TRUE)
  )

# Join data
pheno <- pheno %>% mutate(ParticipantID = STUDYID)
pharma_full <- ad_mapped %>%
  left_join(pheno, by = "ParticipantID") %>%
  left_join(pharma, by = "ParticipantID") %>%
  filter(MDD==1)

# -- Create recurrent depression variables
pharma_full <- pharma_full %>%
  mutate(
    Recurrent = as.integer(TIMES2WK > 1),
    Recurrent_TEs = as.integer(Total_TEs > 1)
  )

# -- Creating a mapping for drug names
drug_col_map <- drug_ref %>%
  select(DrugName, Drug) %>%
  filter(DrugName %in% pharma_full$DrugName)


# -- Create derived columns function
create_derived_column <- function(df, prefix) {
  col_name <- paste0(prefix, "AD")
  
  # Initialize column with NA
  result <- df %>% mutate(!!sym(col_name) := NA_real_)
  
  # Update for each drug
  for (i in 1:nrow(drug_col_map)) {
    drug <- drug_col_map$DrugName[i]
    suffix <- drug_col_map$Drug[i]
    
    # Set values for matching drugs
    result <- result %>%
      mutate(!!sym(col_name) := case_when(
        DrugName == drug & !!sym(suffix) == 1 ~ 1,
        DrugName == drug & !!sym(suffix) == 0 ~ 0,
        TRUE ~ !!sym(col_name)
      ))
  }
  
  return(result)
}

#-- Add WELLAD and STOPAD columns
pharma_full <- pharma_full %>%
  create_derived_column("WELL") %>%
  create_derived_column("STOP")

#-- Define independent variables
#-- Dementia related medications was excluded as 0% were dispensing
independent <- c("AGE", "SEX", "AGE2WKF", "Recurrent", "circadian",
                 "BMI", "SUICIDEA", "SELFHARM", "EDU", "PHYSHLTH", "REGSMK", "DRK3FRQ", "preg_last5_deterministic", "BFEEDANY",
                 "TYPE11B", "TYPE33C","DXPDMD", "DXANX", "DXPERSD", "DXBPD2", "DXSUD", "DXADHD", "DXOCD", "DXSAD", "DXANOR", "DXSCZ",
                 "DXPHYS3", "DXPHYS6", "DXPHYS12", "DXPHYS34", "MIGEVR", "DXENDO", "DXFIBRO", "DXPCOS",
                 "LOWINT2W", "DEP2WK", "FATIGUED.x", "GUILTY", "NOFOCUS", "DEATHTHK", "APWTCHANGE", "SLEEP", "MOVEMENT", "atypical",
                 "Augmentation", "num_conditions", "ADHD", "Analgesics", "Anxiety", "Asthma_and_COPD", "Cancer", "Cardiovascular",
                 "Diabetes", "Dyslipidemia", "Hepatic", "Immunosuppressants", "Siezures", "Sleep", "Thyroid")

#-- Vector of non-binary variables
non_binary <- c("AGE", "AGE2WKF", "BMI", "EDU", "PHYSHLTH", "DRK3FRQ",  "GA_score", "num_conditions")

# -- Convert SEX
pharma_full <- pharma_full %>%
  mutate(SEX = case_when(
    SEX == "Female" ~ 0,
    SEX == "Male" ~ 1,
    TRUE ~ NA_real_
  ))

# -- Linear model function
run_lm <- function(var, data, dependent_var, adjust_vars = NULL) {
  
  # Filter only for this specific variable
  model_data <- data %>% filter(!is.na(!!sym(var)))
  
  #-- Number of cases and controls if variable is binary
  if (!(var %in% non_binary)) {
    Total_Cases <- model_data %>%
      filter(!!sym(var) == 1) %>%
      nrow()
    
    Total_Controls <- model_data %>%
      filter(!!sym(var) == 0) %>%
      nrow()
  } else {
    Total_Cases <- NA
    Total_Controls <- NA
  }
  
  if (is.null(adjust_vars)) {
    # No adjustment variables
    formula_str <- as.formula(paste0(dependent_var, " ~ ", var))
  } else {
    # With adjustment variables
    formula_str <- as.formula(paste0(dependent_var, " ~ ", paste(adjust_vars, collapse = " + "), " + ", var))
  }
  
  # Run model
  lm_model <- lm(formula_str, data = model_data)
  
  # Extract results
  model_summary <- summary(lm_model)
  coefficients <- tibble::rownames_to_column(as.data.frame(coef(model_summary)), "Term") %>%
    as_tibble()
  
  # Add metadata
  coefficients <- coefficients %>%
    mutate(
      Independent = var,
      Total_N = nrow(model_data),
      Total_Cases = Total_Cases,
      Total_Controls = Total_Controls
    )
  
  return(coefficients)
}

# -- Prepare base data for modeling
tab <- pharma_full %>%
  select(ParticipantID, SEX, AGE, num_ATC, num_class, TotalPrescriptionDays, all_of(independent)) %>%
  distinct()

#=========== Analysis 1: Total Dispense (Cumulative Total Prescription Days) =============================================

lm1_base_data <- tab %>% 
  filter(!is.na(AGE) & !is.na(SEX) & !is.na(TotalPrescriptionDays))

# Apply run_lm to each independent variable
with_progress({
  p <- progressor(steps = length(independent))
  
  lm1_results <- map_dfr(independent, function(var) {
    result <- if (var %in% c("AGE", "SEX")) {
      run_lm(var, lm1_base_data, "TotalPrescriptionDays")
    } else {
      run_lm(var, lm1_base_data, "TotalPrescriptionDays", c("AGE", "SEX"))
    }
    p(sprintf("Completed %s", var))
    result
  })
})

# -- Format results
results1 <- lm1_results %>%
  mutate(Dependent = "Cumulative Prescription Dispense (days)") %>%
  arrange(Independent) %>%
  select(Dependent, Independent, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, Total_N, Total_Cases, Total_Controls)

# Calculate FDR and Bonferroni p-values
terms <- results1 %>%
  filter(!(Term %in% c("(Intercept)", "AGE", "SEX")) | (Independent %in% c("AGE", "SEX") & Term != "(Intercept)")) %>%
  mutate(
    FDR_P = p.adjust(`Pr(>|t|)`, method = "fdr"),
    Bonf_P = p.adjust(`Pr(>|t|)`, method = "bonferroni"),
    Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
    Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
  ) %>%
  select(Dependent, Independent, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)

# Join with main results
results1 <- results1 %>%
  left_join(terms, by = c("Dependent", "Independent", "Term")) %>%
  arrange(Dependent, Independent) %>%
  group_by(Dependent, Independent) %>%
  mutate(
    Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
    Independent = if_else(Independent != lag(Independent, default = ""), Independent, NA_character_),
    Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_),
    Total_Cases = if_else(Total_Cases != lag(Total_Cases, default = 0), Total_Cases, NA_integer_),
    Total_Controls = if_else(Total_Controls != lag(Total_Controls, default = 0), Total_Controls, NA_integer_)
  ) %>%
  mutate(
    across(c(`Pr(>|t|)`, FDR_P, Bonf_P), p_fmt),
    across(c(`Std. Error`,`t value`), ~ round(.x, 2)),
    across(c(Estimate), ~ round(.x, 1))
  )

#================= Medication-specific effects analysis =====================================================
filtered_pheno <- pharma_full %>%
  filter(!is.na(AGE) & !is.na(SEX)) %>%
  select(ParticipantID, SEX, AGE, PrescriptionDays, DrugName, DrugClass) %>%
  distinct() %>%
  filter(ParticipantID %in% pheno$ParticipantID) # Filter for those with phenotypes

# Count participants per medication
Medication_N <- filtered_pheno %>%
  count(DrugName, name = "Medication_N")

# Run medication effects model (set reference level to Escitalopram)
filtered_pheno <- filtered_pheno %>%
  mutate(DrugName = factor(DrugName, levels = c("SSRI:Escitalopram", 
                                                setdiff(unique(drug_ref$DrugName), "SSRI:Escitalopram"))))

#lm_model <- lm(PrescriptionDays ~ AGE + SEX + DrugName, data = filtered_pheno)
#lm_output <- tibble::rownames_to_column(as.data.frame(coef(summary(lm_model))), "Term") %>%
#  as_tibble()

mixed_model <- lmer(PrescriptionDays ~ AGE + SEX + DrugName + (1|ParticipantID), 
                    data = filtered_pheno)
mixed_output <- tibble::rownames_to_column(
  as.data.frame(coef(summary(mixed_model))), 
  "Term"
) %>%
  as_tibble()

# Clean up variable names
lm_output <- mixed_output %>%
  mutate(
    Term = str_replace(Term, "DrugName", ""),
    Term = if_else(Term == "(Intercept)", "SSRI:Escitalopram", Term),
    Total_N = length(unique(filtered_pheno$ParticipantID))
  )

# Join with participant counts
med_results <- lm_output %>%
  left_join(Medication_N, by = c("Term" = "DrugName")) %>%
  mutate(
    #`Pr(>|t|)` = format(signif(`Pr(>|t|)`, 2),scientific = TRUE),
    across(c(`Std. Error`,`t value`), ~ round(.x, 2)),
    across(c(Estimate), ~ round(.x, 1)),
    Dependent = "Cumulative Prescription Dispense (days)"
  ) %>%
  arrange(Dependent, Term) %>%
  select(Dependent, Term, Estimate, `Std. Error`, `t value`, Total_N, Medication_N)

# Format for reporting
med_results <- med_results %>%
  group_by(Dependent) %>%
  mutate(
    Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
    Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_)
  ) %>%
  mutate(
    Term = if_else(Term == "SSRI:Escitalopram", "Reference:SSRI:Escitalopram",
                   if_else(Term == "AGE", "AGE (yrs)", 
                           if_else(Term == "SEX", "SEX(Female=1, Male=2)", Term)
                   )
    )
  )
    
# Save the results
#write.csv(med_results, file.path(output_dir, "LM_PrescriptionDays_MedicationEffects_Escitalopram_Reference.csv"), quote = FALSE, row.names = FALSE)

#============ Class-specific effects ==============================================================
# Set reference level for drug class
filtered_pheno <- filtered_pheno %>%
  mutate(DrugClass = factor(DrugClass, levels = c("SSRI", 
                                                setdiff(unique(drug_ref$DrugClass), "SSRI")))) %>%
  filter(ParticipantID %in% pheno$ParticipantID) # Filter for those with phenotypes

# Count by class
Class_N <- filtered_pheno %>%
  count(DrugClass, name = "Class_N")

# Run class effects model
#lm_model <- lm(PrescriptionDays ~ AGE + SEX + DrugClass, data = filtered_pheno)
#lm_output <- tibble::rownames_to_column(as.data.frame(coef(summary(lm_model))), "Term") %>%
#  as_tibble()

mixed_model <- lmer(PrescriptionDays ~ AGE + SEX + DrugClass + (1|ParticipantID), 
                    data = filtered_pheno)
mixed_output <- tibble::rownames_to_column(
  as.data.frame(coef(summary(mixed_model))), 
  "Term"
) %>%
  as_tibble()


# Clean up variable names
lm_output <- mixed_output %>%
  mutate(
    Term = str_replace(Term, "DrugClass", ""),
    Term = if_else(Term == "(Intercept)", "SSRI", Term),
    Total_N = length(unique(filtered_pheno$ParticipantID))
  )

# Join with class counts
class_results <- lm_output %>%
  left_join(Class_N, by = c("Term" = "DrugClass")) %>%
  mutate(
    #`Pr(>|t|)` = format(signif(`Pr(>|t|)`, 2),scientific = TRUE),
    across(c(`Std. Error`,`t value`), ~ round(.x, 2)),
    across(c(Estimate), ~ round(.x, 1)),
    Dependent = "Cumulative Prescription Dispense (days)"
  ) %>%
  arrange(Dependent, Term) %>%
  select(Dependent, Term, Estimate, `Std. Error`, `t value`, Total_N, Class_N)

# Format for reporting
class_results <- class_results %>%
  group_by(Dependent) %>%
  mutate(
    Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
    Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_)
  )  %>%
  mutate(
    Term = if_else(Term == "SSRI", "Reference:SSRI",
                   if_else(Term == "AGE", "AGE (yrs)", 
                           if_else(Term == "SEX", "SEX(Female=1, Male=2)", Term)
                   )
    )
  )
# Save the results
#write.csv(class_results, file.path(output_dir, "LM_PrescriptionDays_ClassEffects_SSRI_Reference.csv"), quote = FALSE, row.names = FALSE)

#============== Analysis 2: Medication Diversity =============================================

run_analysis <- function(data, dependent_var, dependent_label) {
  # Filter data
  base_data <- data %>% 
    filter(!is.na(AGE) & !is.na(SEX) & !is.na(!!sym(dependent_var)))
  
  # Run models for each independent variable
  with_progress({
    p <- progressor(steps = length(independent))
    
    results <- map_dfr(independent, function(var) {
      result <- if (var %in% c("AGE", "SEX")) {
        run_lm(var, base_data, dependent_var)
      } else {
        run_lm(var, base_data, dependent_var, c("AGE", "SEX"))
      }
      p()  # Just increment the counter
      result
    })
  })
  
  # Format results
  results <- results %>%
    mutate(Dependent = dependent_label) %>%
    arrange(Independent) %>%
    select(Dependent, Independent, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, Total_N, Total_Cases, Total_Controls)
  
  # Calculate FDR and Bonferroni p-values
  terms <- results %>%
    filter(!(Term %in% c("(Intercept)", "AGE", "SEX")) | (Independent %in% c("AGE", "SEX") & Term != "(Intercept)")) %>%
    mutate(
      FDR_P = p.adjust(`Pr(>|t|)`, method = "fdr"),
      Bonf_P = p.adjust(`Pr(>|t|)`, method = "bonferroni"),
      Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
      Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
    ) %>%
    select(Dependent, Independent, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)
  
  # Join and format
  results <- results %>%
    left_join(terms, by = c("Dependent", "Independent", "Term")) %>%
    arrange(Dependent, Independent) %>%
    group_by(Dependent, Independent) %>%
    mutate(
      Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
      Independent = if_else(Independent != lag(Independent, default = ""), Independent, NA_character_),
      Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_),
      Total_Cases = if_else(Total_Cases != lag(Total_Cases, default = 0), Total_Cases, NA_integer_),
      Total_Controls = if_else(Total_Controls != lag(Total_Controls, default = 0), Total_Controls, NA_integer_)
    ) %>%
    mutate(
      across(c(`Pr(>|t|)`, FDR_P, Bonf_P), p_fmt),
      across(c(`Std. Error`, Estimate), ~ round(.x, 4)),
      across(c(`t value`), ~ round(.x, 2))
    )
}

# Run Analysis 2 and 3 with the function
results2 <- run_analysis(tab, "num_ATC", "Medication Diversity") 
results3 <- run_analysis(tab, "num_class", "Class Diversity")

#========= Rename variables for readability ==========================
# Remapping of self-report variables
rename_mapping <- c(
  "AGE" = "Age", 
  "SEX" = "Sex", 
  "AGE2WKF" = "Age of MDD Onset",
  "Recurrent" = "Recurrent 2-Weeks of MDD",
  "BMI" = "Body Mass Index",
  "SUICIDEA" = "Suicidal Ideation",
  "SELFHARM" = "Self harm",
  "circadian" = "Circadian Subtype",
  "EDU" = "Education level",
  "PHYSHLTH" = "Physical health",
  "REGSMK" = "Regular Smoker",
  "DRK3FRQ" = "Drinks over 3 Months",
  "preg_last5_deterministic" = "Likely Pregnant",
  "BFEEDANY" = "Breastfeed Ever",
  "TYPE11B" = "Type 2 Diabetes", 
  "TYPE33C" = "Stomach Ulcers",
  "DXANX" = "Anxiety Disorder",
  "DXPERSD" = "Personality Disorder",
  "DXBPD2" = "Bipolar Disorder",
  "DXSCZ" = "Schizophrenia",
  "DXANOR" = "Anorexia Nervosa",
  "DXSUD" = "Substance Use Disorder",
  "DXADHD" = "ADHD",
  "DXOCD" = "Obsessive-compulsive Disorder",
  "DXSAD" = "Seasonal Affective Disorder",
  "DXPHYS3" = "Back pain",
  "DXPHYS6" = "Chronic Fatigue Syndrome",
  "DXPHYS12" = "Epilepsy",
  "DXPHYS34" = "Chronic pain",
  "MIGEVR" = "Migraines or Headaches",
  "DXPDMD" = "Premenstrual Dysphoric Disorder",
  "DXENDO" = "Endometriosis",
  "DXFIBRO" = "Fibroids (uterus)",
  "DXPCOS" = "PCOS",
  "LOWINT2W" = "Low Interest 2-Weeks",
  "DEP2WK" = "Depressed 2-Weeks",
  "FATIGUED.x" = "Fatigued",
  "GUILTY" = "Guilty Feelings",
  "NOFOCUS" = "No Focus",
  "DEATHTHK" = "Death Thoughts",
  "APWTCHANGE" = "Appetite/Weight Change",
  "SLEEP" = "Sleep Disturbances",
  "MOVEMENT" = "Movement Changes",
  "atypical" = "Atypical Subtype",
  "Augmentation" = "Augmentation",
  "num_conditions" = "Number Co-occurring Conditions",
  "ADHD" = "ADHD-related Medications",
  "Analgesics" = "Analgesics-related Medications",
  "Anxiety" = "Anxiety-related Medications",
  "Asthma_and_COPD" = "Asthma/COPD-related Medications",
  "Cancer" = "Cancer-related Medications",
  "Cardiovascular" = "Cardiovascular-related Medications",
  "Diabetes" = "Diabetes-related Medications",
  "Dyslipidemia" = "Dyslipidemia-related Medications",
  "Hepatic" = "Hepatic-related Medications",
  "Immunosuppressants" = "Immunosuppressants",
  "Siezures" = "Siezure-related Medications",
  "Sleep" = "Sleep-related Medications",
  "Thyroid" = "Thyroid-related Medications",
  "(Intercept)" = "Intercept"
)

# 59 phenotypes
apply_renaming <- function(df) {
  df %>%
    mutate(
      Independent = recode(Independent, !!!rename_mapping),
      Term = recode(Term, !!!rename_mapping)
    )
}

results1_renamed <- apply_renaming(results1)
results2_renamed <- apply_renaming(results2)
results3_renamed <- apply_renaming(results3)

results_all <- bind_rows(results1_renamed, results2_renamed, results3_renamed)

#================= Write Results =========================================================
# Create a new workbook
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")

# Medication and class effects
drug_results <- bind_rows(med_results, class_results)
removeWorksheet(wb, "Table2")
addWorksheet(wb, "Table2")
writeData(wb, "Table2", drug_results)

# Main results
removeWorksheet(wb, "Table3")
addWorksheet(wb, "Table3")
writeData(wb, "Table3", results_all)

# Auto-adjust column widths for better readability
for(sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:20, widths = "auto")
}

# Save the workbook
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)