# -- Load R packages (optimized selection)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
  library(broom)
  library(purrr)
})

# -- Set path and duration parameters
wkdir <- "/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes"
duration_list <- c(360, 600)
adherence_thresholds <- c("All", "Adherent")  # Now using both thresholds
output_dir <- "/scratch/user/uqawal15"

# -- Pre-define constants for efficiency
drug_classes <- c("SSRI", "SNRI", "TCA", "TeCA", "BIP+L", "BIP-L", "Various")
pharma_classes <- c("SSRI", "SNRI", "TCA", "TeCA")
pharma_drugs <- c("SSRI:Sertraline", "SSRI:Escitalopram", "SSRI:Citalopram", 
                  "SNRI:Venlafaxine", "SNRI:Duloxetine", "TCA:Amitriptyline", 
                  "SSRI:Paroxetine", "TeCA:Mirtazapine", "SNRI:Desvenlafaxine", 
                  "SSRI:Fluoxetine")

all_quantitative_outcomes <- c(
  "AGE", "AGE2WKF_scaled", "TIMES2WK_scaled", "BMI", "EDU_scaled", "PHYSHLTH_scaled",
  "num_ATC_scaled", "num_class_scaled", "NumberOfPrescriptionEpisodes_scaled", 
  "AveragePrescriptionEpisodesDays_scaled", "PrescriptionDays_scaled", "GA_score_scaled", "num_conditions_scaled"
)

# ATC mapping (pre-defined for efficiency)
atc_drug_ref <- tibble(
  DrugName = pharma_drugs,
  ATCCodes = c('N06AB06', 'N06AB10', 'N06AB04', 'N06AX16', 'N06AX21', 
               'N06AA09', 'N06AB05', 'N06AX11', 'N06AX23', 'N06AB03'),
  DrugClass = c('SSRI', 'SSRI', 'SSRI', 'SNRI', 'SNRI', 'TCA', 'SSRI', 'TeCA', 'SNRI', 'SSRI')
)

# -- Optimized function to process phenotypes (vectorized operations)
process_phenotypes <- function(pheno_data, exposome_vars) {
  # Efficient SEX conversion using case_when
  pheno_data$SEX <- case_when(
    pheno_data$SEX == "Female" ~ 0,
    pheno_data$SEX == "Male" ~ 1,
    TRUE ~ NA_real_
  )
  
  # Vectorized scaling for all _scaled variables at once
  scale_vars <- grep("_scaled$", exposome_vars, value = TRUE)
  base_vars <- gsub("_scaled$", "", scale_vars)
  
  # Check which base variables exist in the data
  existing_base_vars <- intersect(base_vars, names(pheno_data))
  existing_scale_vars <- paste0(existing_base_vars, "_scaled")
  
  # Apply scaling vectorized
  if (length(existing_base_vars) > 0) {
    scaled_data <- pheno_data[existing_base_vars] %>%
      mutate(across(everything(), ~ as.numeric(scale(.x)[,1])))
    names(scaled_data) <- existing_scale_vars
    
    # Bind with original data
    pheno_data <- bind_cols(pheno_data, scaled_data)
  }
  
  # Select required columns efficiently
  required_cols <- c("ParticipantID", "SEX", "AGE", "DrugName", "DrugClass")
  available_vars <- intersect(exposome_vars, names(pheno_data))
  final_cols <- unique(c(required_cols, available_vars))
  
  return(pheno_data[final_cols])
}

# -- Optimized model fitting function with better error handling
fit_models_optimized <- function(data, outcome_var, reference_groups, group_col = "DrugName", age_adj = TRUE) {
  
  # Pre-filter data once
  data_filtered <- data %>%
    filter(
      !is.na(AGE), 
      !is.na(SEX), 
      !is.na(.data[[outcome_var]]),
      !is.na(.data[[group_col]])
    )
  
  if (nrow(data_filtered) == 0) {
    return(NULL)
  }
  
  # Build formula once
  if (age_adj && outcome_var != "AGE") {
    formula_str <- paste0(outcome_var, " ~ AGE + SEX + ", group_col)
  } else {
    formula_str <- paste0(outcome_var, " ~ ", group_col)
  }
  formula_obj <- as.formula(formula_str)
  
  # Vectorized processing of all reference groups
  results_list <- map_dfr(reference_groups, function(ref_group) {
    tryCatch({
      # Set reference level efficiently
      temp_data <- data_filtered
      temp_data[[group_col]] <- factor(temp_data[[group_col]], 
                                       levels = c(ref_group, setdiff(unique(temp_data[[group_col]]), ref_group)))
      
      # Fit model
      model <- lm(formula_obj, data = temp_data)
      
      # Extract results efficiently
      coef_summary <- broom::tidy(model)
      coef_summary$term <- gsub(paste0("^", group_col), "", coef_summary$term)
      coef_summary$reference_term <- ref_group
      coef_summary$term[1] <- ref_group
      
      # Get group counts efficiently
      group_counts <- temp_data %>%
        count(.data[[group_col]], name = "group_n")
      names(group_counts)[1] <- "term"
      
      # Combine results
      result <- coef_summary %>%
        left_join(group_counts, by = "term") %>%
        mutate(outcome = outcome_var)
      
      return(result)
      
    }, error = function(e) {
      cat(sprintf("Error fitting model for %s with reference %s: %s\n", 
                  outcome_var, ref_group, e$message))
      return(NULL)
    })
  })
  
  return(results_list)
}

# -- Load pharma data once (avoid reloading in loop)
cat("Loading pharmaceutical data...\n")
ad <- read.csv("/QRISdata/Q7280/pharmacogenomics/data/AGDSAcceptabilityTreatmentGroups_25082025.csv",
               stringsAsFactors = FALSE)

# Convert date columns manually
ad$EarliestPrescription <- as.Date(ad$EarliestPrescription, format = "%d/%m/%Y")
ad$LatestPrescription <- as.Date(ad$LatestPrescription, format = "%d/%m/%Y")

# Pre-process pharma data once
ad_mapped <- ad %>%
  left_join(atc_drug_ref, by = c("ATCCode" = "ATCCodes"))

pharma_summary <- ad_mapped %>%
  group_by(ParticipantID) %>%
  summarise(
    num_ATC = n_distinct(ATCCode),
    num_class = n_distinct(DrugClass),
    .groups = "drop"
  )

# -- Main processing loop (NOW INCLUDES BOTH ADHERENCE GROUPS)
all_results <- map_dfr(duration_list, function(duration) {
  cat(sprintf("Processing duration: %d days\n", duration))
  
  # Process both adherence thresholds
  map_dfr(adherence_thresholds, function(adherence) {
    cat(sprintf("  Processing adherence group: %s\n", adherence))
    
    # Load appropriate phenotype data based on adherence
    if (adherence == "All") {
      pheno <- read.csv(file.path(wkdir, paste0("inner_data_", duration, "days.csv")),
                        stringsAsFactors = FALSE)
    } else {
      pheno <- read.csv(file.path(wkdir, paste0("inner_data_", duration, "days_Adherent.csv")),
                        stringsAsFactors = FALSE)
    }
    
    # Get unique drugs for this duration/adherence combination
    unique_drugs <- unique(pheno$DrugName)
    
    # Combine with pharma data efficiently
    pheno_combined <- left_join(pheno, pharma_summary, by = "ParticipantID")
    
    # Process phenotypes
    combined_tab <- process_phenotypes(pheno_combined, all_quantitative_outcomes)
    
    # Process all outcomes
    adherence_results <- map_dfr(all_quantitative_outcomes, function(outcome) {
      if (!outcome %in% names(combined_tab)) {
        return(NULL)
      }
      
      cat(sprintf("    Processing outcome: %s\n", outcome))
      
      # Determine if pharma outcome
      is_pharma <- grepl("AveragePrescriptionEpisodesDays|NumberOfPrescriptionEpisodes", outcome)
      
      # Drug-level analysis
      if (is_pharma) {
        drug_data <- combined_tab %>% 
          filter(!DrugName %in% c("BIP+L", "BIP-L", "Various", "Combination"))
        drug_refs <- pharma_drugs
      } else {
        drug_data <- combined_tab
        drug_refs <- unique_drugs
      }
      
      drug_results <- fit_models_optimized(drug_data, outcome, drug_refs, "DrugName", outcome != "AGE")
      if (!is.null(drug_results)) {
        drug_results$analysis_level <- "Drug"
      }
      
      # Class-level analysis
      class_data <- combined_tab %>% filter(!is.na(DrugClass))
      if (is_pharma) {
        class_data <- class_data %>%
          filter(!DrugClass %in% c("BIP+L", "BIP-L", "Various"))
        class_refs <- pharma_classes
      } else {
        class_refs <- drug_classes
      }
      
      class_results <- fit_models_optimized(class_data, outcome, class_refs, "DrugClass", outcome != "AGE")
      if (!is.null(class_results)) {
        class_results$analysis_level <- "Class"
      }
      
      # Combine drug and class results
      combined_results <- bind_rows(drug_results, class_results)
      if (!is.null(combined_results)) {
        combined_results$threshold <- paste0(duration, " days")
        combined_results$adherence <- adherence  # Add adherence info
      }
      
      return(combined_results)
    })
    
    return(adherence_results)
  })
})

# -- Separate drug and class results efficiently
drug_results <- all_results %>% 
  filter(analysis_level == "Drug") %>%
  select(-analysis_level)

class_results <- all_results %>% 
  filter(analysis_level == "Class") %>%
  select(-analysis_level)

# -- Updated formatting function to handle adherence
format_results_optimized <- function(results, ref_term) {
  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }
  
  # Calculate p-value adjustments efficiently (now grouped by adherence too)
  results_with_adj <- results %>%
    filter(!(term %in% c("AGE", "SEX")) | (outcome %in% c("AGE", "SEX"))) %>%
    group_by(adherence, threshold, reference_term, term) %>%
    mutate(
      fdr_p = p.adjust(p.value, method = "fdr"),
      bonf_p = p.adjust(p.value, method = "bonferroni"),
      sig_fdr = ifelse(fdr_p < 0.05, "*", ""),
      sig_bonf = ifelse(bonf_p < 0.05, "*", "")
    ) %>%
    ungroup() %>%
    select(adherence, threshold, reference_term, term, outcome, fdr_p, bonf_p, sig_fdr, sig_bonf)
  
  # Calculate total N efficiently (now grouped by adherence too)
  total_n <- results %>%
    filter(reference_term == ref_term) %>%
    group_by(adherence, threshold, outcome) %>%
    summarise(total_n = sum(group_n, na.rm = TRUE), .groups = "drop")
  
  # Final formatting
  formatted <- results %>%
    left_join(results_with_adj, by = c("adherence", "threshold", "reference_term", "term", "outcome")) %>%
    left_join(total_n, by = c("adherence", "threshold", "outcome")) %>%
    mutate(
      #across(c(p.value, fdr_p, bonf_p), ~ format(signif(.x, 2), scientific = TRUE)),
      across(c(std.error), ~ signif(.x, 2)),
      across(c(estimate, statistic), ~ round(.x, 2))
    ) %>%
    select(adherence, threshold, reference_term, outcome, term, estimate, std.error, 
           statistic, p.value, fdr_p, bonf_p, sig_fdr, sig_bonf, total_n, group_n) %>%
    arrange(adherence, threshold, reference_term, outcome, term) %>%
    group_by(adherence, threshold, reference_term, outcome) %>%
    mutate(
      adherence = ifelse(adherence != lag(adherence, default = ""), adherence, NA_character_),
      threshold = ifelse(threshold != lag(threshold, default = ""), threshold, NA_character_),
      reference_term = ifelse(reference_term != lag(reference_term, default = ""), reference_term, NA_character_),
      outcome = ifelse(outcome != lag(outcome, default = ""), outcome, NA_character_),
      total_n = ifelse(total_n != lag(total_n, default = 0), total_n, NA_integer_)
    ) %>%
    ungroup()
  
  return(formatted)
}

# -- Format results
cat("Formatting results...\n")
drug_formatted <- format_results_optimized(drug_results, "SSRI:Sertraline")
class_formatted <- format_results_optimized(class_results, "SSRI")

# -- Renaming mapping
rename_mapping <- c(
  "GA_score_scaled" = "General Adherence Score", 
  "num_conditions_scaled" = "Number of Co-occuring Conditions",
  "num_ATC_scaled" = "AD Diversity",
  "num_class_scaled" = "Class Diversity",
  "NumberOfPrescriptionEpisodes_scaled" = "Number of Prescription Episodes",
  "AveragePrescriptionEpisodesDays_scaled" = "Average Prescription Episode Length",
  "PrescriptionDays_scaled" = "Cumulative AD Dispense",
  "AGE2WKF_scaled" = "Age of MDD Onset",
  "TIMES2WK_scaled" = "Times 2-Weeks of MDD",
  "EDU_scaled" = "Education Level",
  "PHYSHLTH_scaled" = "Physical health", 
  "AGE" = "Age", 
  "BMI" = "BMI"
)

# -- Save to Excel efficiently
cat("Saving to Excel...\n")
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")

# Process and save class results
if (!is.null(class_formatted)) {
  class_renamed <- class_formatted %>%
    mutate(outcome = recode(outcome, !!!rename_mapping))
  removeWorksheet(wb, "Table9")
  addWorksheet(wb, "Table9")
  writeData(wb, "Table9", class_renamed)
}

# Process and save drug results
if (!is.null(drug_formatted)) {
  drug_renamed <- drug_formatted %>%
    mutate(outcome = recode(outcome, !!!rename_mapping))
  removeWorksheet(wb, "Table10")
  addWorksheet(wb, "Table10")
  writeData(wb, "Table10", drug_renamed)
}

saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)
cat("Analysis complete! Results saved to All_Results.xlsx\n")






