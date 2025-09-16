# -- Load packages
library(dplyr)
library(purrr)
library(stringr)
library(openxlsx)
library(progressr)

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

#-- Define independent variables (removed SEX since we're stratifying by it)
independent_base <- c("AGE", "AGE2WKF", "Recurrent", "circadian",
                      "BMI", "SUICIDEA", "SELFHARM", "EDU", "PHYSHLTH", "REGSMK", "DRK3FRQ",
                      "TYPE11B", "TYPE33C","DXANX", "DXPERSD", "DXBPD2", "DXSUD", "DXADHD", "DXOCD", "DXSAD", "DXSCZ",
                      "DXPHYS3", "DXPHYS6", "DXPHYS12", "DXPHYS34", "MIGEVR", 
                      "LOWINT2W", "DEP2WK", "FATIGUED.x", "GUILTY", "NOFOCUS", "DEATHTHK", "APWTCHANGE", "SLEEP", "MOVEMENT", "atypical",
                      "Augmentation", "num_conditions", "ADHD", "Analgesics", "Anxiety", "Asthma_and_COPD", "Cancer", "Cardiovascular",
                      "Diabetes", "Dyslipidemia", "Hepatic", "Immunosuppressants", "Siezures", "Sleep", "Thyroid")

# Sex-specific variables (only for females)
female_specific <- c("preg_last5_deterministic", "BFEEDANY",  "DXENDO", "DXFIBRO", "DXPCOS", "DXPDMD", "DXANOR")

#-- Vector of non-binary variables
non_binary <- c("AGE", "AGE2WKF", "BMI", "EDU", "PHYSHLTH", "DRK3FRQ", "GA_score", "num_conditions")

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

# -- Function to run analysis for a specific sex
run_sex_stratified_analysis <- function(sex_value, sex_label) {
  cat(paste("Running analysis for", sex_label, "\n"))
  
  # Filter data by sex
  pharma_sex <- pharma_full %>% filter(SEX == sex_value)
  
  # Select appropriate variables based on sex
  if (sex_value == 0) { # Female
    independent <- c(independent_base, female_specific)
  } else { # Male
    independent <- independent_base
  }
  
  # Prepare base data for modeling
  tab_sex <- pharma_sex %>%
    select(ParticipantID, AGE, num_ATC, num_class, TotalPrescriptionDays, all_of(independent)) %>%
    distinct()
  
  #=========== Analysis 1: Total Dispense (Cumulative Total Prescription Days) =============================================
  
  lm1_base_data <- tab_sex %>% 
    filter(!is.na(AGE) & !is.na(TotalPrescriptionDays))
  
  # Apply run_lm to each independent variable
  with_progress({
    p <- progressor(steps = length(independent))
    
    lm1_results <- map_dfr(independent, function(var) {
      result <- if (var == "AGE") {
        run_lm(var, lm1_base_data, "TotalPrescriptionDays")
      } else {
        run_lm(var, lm1_base_data, "TotalPrescriptionDays", c("AGE"))
      }
      p(sprintf("Completed %s for %s", var, sex_label))
      result
    })
  })
  
  # -- Format results
  results1 <- lm1_results %>%
    mutate(Dependent = "Cumulative Prescription Dispense (days)",
           Sex_Stratum = sex_label) %>%
    arrange(Independent) %>%
    select(Dependent, Sex_Stratum, Independent, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, Total_N, Total_Cases, Total_Controls)
  
  # Calculate FDR and Bonferroni p-values
  terms <- results1 %>%
    filter(!(Term %in% c("(Intercept)", "AGE")) | (Independent == "AGE" & Term != "(Intercept)")) %>%
    mutate(
      FDR_P = p.adjust(`Pr(>|t|)`, method = "fdr"),
      Bonf_P = p.adjust(`Pr(>|t|)`, method = "bonferroni"),
      Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
      Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
    ) %>%
    select(Dependent, Sex_Stratum, Independent, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)
  
  # Join with main results
  results1 <- results1 %>%
    left_join(terms, by = c("Dependent", "Sex_Stratum", "Independent", "Term")) %>%
    arrange(Dependent, Independent) %>%
    distinct() %>%  # Remove any duplicate rows
    group_by(Dependent, Independent) %>%
    mutate(
      Dependent = if_else(row_number() == 1, Dependent, NA_character_),
      Sex_Stratum = if_else(row_number() == 1, Sex_Stratum, NA_character_),
      Independent = if_else(row_number() == 1, Independent, NA_character_),
      Total_N = if_else(row_number() == 1, Total_N, NA_integer_),
      Total_Cases = if_else(row_number() == 1, Total_Cases, NA_integer_),
      Total_Controls = if_else(row_number() == 1, Total_Controls, NA_integer_)
    ) %>%
    mutate(
      across(c(`Pr(>|t|)`, FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE)),
      across(c(`Std. Error`,`t value`), ~ round(.x, 2)),
      across(c(Estimate), ~ round(.x, 1))
    )
  
  #================= Medication-specific effects analysis =====================================================
  filtered_pheno_sex <- pharma_sex %>%
    filter(!is.na(AGE)) %>%
    select(ParticipantID, AGE, PrescriptionDays, DrugName, DrugClass) %>%
    distinct() %>%
    filter(ParticipantID %in% pheno$ParticipantID) # Filter for those with phenotypes
  
  # Count participants per medication
  Medication_N <- filtered_pheno_sex %>%
    count(DrugName, name = "Medication_N")
  
  # Run medication effects model (set reference level to Escitalopram)
  filtered_pheno_sex <- filtered_pheno_sex %>%
    mutate(DrugName = factor(DrugName, levels = c("SSRI:Escitalopram", 
                                                  setdiff(unique(drug_ref$DrugName), "SSRI:Escitalopram"))))
  
  # Check if we have enough data for medication analysis
  if (nrow(filtered_pheno_sex) > 0 && length(unique(filtered_pheno_sex$DrugName)) > 1) {
    lm_model <- lm(PrescriptionDays ~ AGE + DrugName, data = filtered_pheno_sex)
    lm_output <- tibble::rownames_to_column(as.data.frame(coef(summary(lm_model))), "Term") %>%
      as_tibble()
    
    # Clean up variable names
    lm_output <- lm_output %>%
      mutate(
        Term = str_replace(Term, "DrugName", ""),
        Term = if_else(Term == "(Intercept)", "SSRI:Escitalopram", Term),
        Total_N = length(unique(filtered_pheno_sex$ParticipantID))
      )
    
    # Join with participant counts
    med_results <- lm_output %>%
      left_join(Medication_N, by = c("Term" = "DrugName")) %>%
      mutate(
        `Pr(>|t|)` = format(signif(`Pr(>|t|)`, 2), scientific = TRUE),
        across(c(`Std. Error`,`t value`), ~ round(.x, 2)),
        across(c(Estimate), ~ round(.x, 1)),
        Dependent = "Cumulative Prescription Dispense (days)",
        Sex_Stratum = sex_label
      ) %>%
      arrange(Dependent, Term) %>%
      select(Dependent, Sex_Stratum, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, Total_N, Medication_N)
    
    # Format for reporting
    med_results <- med_results %>%
      group_by(Dependent) %>%
      mutate(
        Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
        Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_)
      ) %>%
      mutate(
        Term = if_else(Term == "SSRI:Escitalopram", "Reference:SSRI:Escitalopram",
                       if_else(Term == "AGE", "AGE (yrs)", Term))
      )
  } else {
    med_results <- tibble()
  }
  
  #============ Class-specific effects ==============================================================
  # Set reference level for drug class
  filtered_pheno_sex <- filtered_pheno_sex %>%
    mutate(DrugClass = factor(DrugClass, levels = c("SSRI", 
                                                    setdiff(unique(drug_ref$DrugClass), "SSRI")))) %>%
    filter(ParticipantID %in% pheno$ParticipantID) # Filter for those with phenotypes
  
  # Count by class
  Class_N <- filtered_pheno_sex %>%
    count(DrugClass, name = "Class_N")
  
  # Run class effects model
  if (nrow(filtered_pheno_sex) > 0 && length(unique(filtered_pheno_sex$DrugClass)) > 1) {
    lm_model <- lm(PrescriptionDays ~ AGE + DrugClass, data = filtered_pheno_sex)
    lm_output <- tibble::rownames_to_column(as.data.frame(coef(summary(lm_model))), "Term") %>%
      as_tibble()
    
    # Clean up variable names
    lm_output <- lm_output %>%
      mutate(
        Term = str_replace(Term, "DrugClass", ""),
        Term = if_else(Term == "(Intercept)", "SSRI", Term),
        Total_N = length(unique(filtered_pheno_sex$ParticipantID))
      )
    
    # Join with class counts
    class_results <- lm_output %>%
      left_join(Class_N, by = c("Term" = "DrugClass")) %>%
      mutate(
        `Pr(>|t|)` = format(signif(`Pr(>|t|)`, 2), scientific = TRUE),
        across(c(`Std. Error`,`t value`), ~ round(.x, 2)),
        across(c(Estimate), ~ round(.x, 1)),
        Dependent = "Cumulative Prescription Dispense (days)",
        Sex_Stratum = sex_label
      ) %>%
      arrange(Dependent, Term) %>%
      select(Dependent, Sex_Stratum, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, Total_N, Class_N)
    
    # Format for reporting
    class_results <- class_results %>%
      group_by(Dependent) %>%
      mutate(
        Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
        Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_)
      ) %>%
      mutate(
        Term = if_else(Term == "SSRI", "Reference:SSRI",
                       if_else(Term == "AGE", "AGE (yrs)", Term))
      )
  } else {
    class_results <- tibble()
  }
  
  #============== Analysis 2 & 3: Medication Diversity =============================================
  
  run_analysis_stratified <- function(data, dependent_var, dependent_label) {
    # Filter data
    base_data <- data %>% 
      filter(!is.na(AGE) & !is.na(!!sym(dependent_var)))
    
    # Run models for each independent variable
    with_progress({
      p <- progressor(steps = length(independent))
      
      results <- map_dfr(independent, function(var) {
        result <- if (var == "AGE") {
          run_lm(var, base_data, dependent_var)
        } else {
          run_lm(var, base_data, dependent_var, c("AGE"))
        }
        p()
        result
      })
    })
    
    # Format results
    results <- results %>%
      mutate(Dependent = dependent_label,
             Sex_Stratum = sex_label) %>%
      arrange(Independent) %>%
      select(Dependent, Sex_Stratum, Independent, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, Total_N, Total_Cases, Total_Controls)
    
    # Calculate FDR and Bonferroni p-values
    terms <- results %>%
      filter(!(Term %in% c("(Intercept)", "AGE")) | (Independent == "AGE" & Term != "(Intercept)")) %>%
      distinct(Dependent, Sex_Stratum, Independent, Term, .keep_all = TRUE) %>%
      mutate(
        FDR_P = p.adjust(`Pr(>|t|)`, method = "fdr"),
        Bonf_P = p.adjust(`Pr(>|t|)`, method = "bonferroni"),
        Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
        Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
      ) %>%
      select(Dependent, Sex_Stratum, Independent, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)
    
    # Join and format
    results <- results %>%
      left_join(terms, by = c("Dependent", "Sex_Stratum", "Independent", "Term"), relationship = "many-to-one") %>%
      arrange(Dependent, Independent) %>%
      distinct() %>%  # Remove any duplicate rows
      group_by(Dependent, Independent) %>%
      mutate(
        Dependent = if_else(row_number() == 1, Dependent, NA_character_),
        Sex_Stratum = if_else(row_number() == 1, Sex_Stratum, NA_character_),
        Independent = if_else(row_number() == 1, Independent, NA_character_),
        Total_N = if_else(row_number() == 1, Total_N, NA_integer_),
        Total_Cases = if_else(row_number() == 1, Total_Cases, NA_integer_),
        Total_Controls = if_else(row_number() == 1, Total_Controls, NA_integer_)
      ) %>%
      ungroup() %>%
      mutate(
        across(c(`Pr(>|t|)`, FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE)),
        across(c(`Std. Error`, Estimate), ~ round(.x, 4)),
        across(c(`t value`), ~ round(.x, 2))
      )
    
    return(results)
  }
  
  # Run Analysis 2 and 3 with the function
  results2 <- run_analysis_stratified(tab_sex, "num_ATC", "Medication Diversity")
  results3 <- run_analysis_stratified(tab_sex, "num_class", "Class Diversity")
  
  # Combine results
  results_combined <- bind_rows(results1, results2, results3)
  drug_results_combined <- bind_rows(med_results, class_results)
  
  return(list(
    main_results = results_combined,
    drug_results = drug_results_combined
  ))
}

#========= Rename variables for readability ==========================
# Remapping of self-report variables (removed SEX since we're stratifying by it)
rename_mapping <- c(
  "AGE" = "Age", 
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

apply_renaming <- function(df) {
  df %>%
    mutate(
      Independent = recode(Independent, !!!rename_mapping),
      Term = recode(Term, !!!rename_mapping)
    )
}

# Run stratified analyses
cat("Running sex-stratified analyses...\n")

# Run analysis for females (SEX = 0)
female_results <- run_sex_stratified_analysis(0, "Female")

# Run analysis for males (SEX = 1)
male_results <- run_sex_stratified_analysis(1, "Male")

# Apply renaming to results
female_results$main_results <- apply_renaming(female_results$main_results)
male_results$main_results <- apply_renaming(male_results$main_results)

# Combine all results
all_main_results <- bind_rows(female_results$main_results, male_results$main_results)
all_drug_results <- bind_rows(female_results$drug_results, male_results$drug_results)

#================= Write Results =========================================================
# Create output file paths
output_main <- file.path(output_dir, "Sex_Stratified_Main_Results.csv")
output_drug <- file.path(output_dir, "Sex_Stratified_Drug_Results.csv")

# Write CSV files
write.csv(all_main_results, output_main, quote = FALSE, row.names = FALSE)
write.csv(all_drug_results, output_drug, quote = FALSE, row.names = FALSE)

# Create Excel workbook with stratified results
wb_stratified <- createWorkbook()

# Add worksheets for each sex
addWorksheet(wb_stratified, "Female_Main_Results")
writeData(wb_stratified, "Female_Main_Results", female_results$main_results)

addWorksheet(wb_stratified, "Male_Main_Results")
writeData(wb_stratified, "Male_Main_Results", male_results$main_results)

addWorksheet(wb_stratified, "Female_Drug_Results")
writeData(wb_stratified, "Female_Drug_Results", female_results$drug_results)

addWorksheet(wb_stratified, "Male_Drug_Results")
writeData(wb_stratified, "Male_Drug_Results", male_results$drug_results)

addWorksheet(wb_stratified, "Combined_Main_Results")
writeData(wb_stratified, "Combined_Main_Results", all_main_results)

addWorksheet(wb_stratified, "Combined_Drug_Results")
writeData(wb_stratified, "Combined_Drug_Results", all_drug_results)

# Auto-adjust column widths for better readability
for(sheet in names(wb_stratified)) {
  setColWidths(wb_stratified, sheet, cols = 1:20, widths = "auto")
}

# Save the workbook
saveWorkbook(wb_stratified, file.path(output_dir, "Sex_Stratified_All_Results.xlsx"), overwrite = TRUE)

cat("Analysis complete! Results saved to:\n")
cat("- Excel file:", file.path(output_dir, "Sex_Stratified_All_Results.xlsx"), "\n")
cat("- Main results CSV:", output_main, "\n")
cat("- Drug results CSV:", output_dir, "\n")