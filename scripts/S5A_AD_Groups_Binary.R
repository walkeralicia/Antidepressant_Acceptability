# -- Load R packages (removed data.table)
suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(tidyr)
  library(forcats)
  library(tibble)
  library(parallel)
  library(openxlsx)
})

# -- Set path and duration parameters
wkdir <- "/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes"
output_dir <- "/scratch/user/uqawal15"

#=============== Binary Survey Traits - Optimized Processing =============================

# -- Define independent variables
# remove bfeedany due to only one level???
exposome_binary <- c("SEX","WELLAD", "STOPAD", "circadian", # clinical characteristics
                     "SUICIDEA", "REGSMK", "preg_last5_deterministic",  # risk factors
                     "TYPE11B", "TYPE33C", "DXPDMD", "DXANX", "DXPERSD", "DXSUD", "DXADHD", "DXOCD", "DXSAD","DXPHYS3", "DXANOR", "DXSCZ",
                     "DXPHYS6", "DXPHYS12", "DXPHYS34", "MIGEVR", "DXFIBRO", "DXPCOS", "DXENDO", # comorbidities
                     "LOWINT2W", "DEP2WK", "FATIGUED.x", "GUILTY", "NOFOCUS", "DEATHTHK", "APWTCHANGE", "SLEEP", "MOVEMENT", "atypical", # mdd symptoms
                     "Augmentation_ADspecific", "Augmentation", "ADHD", "Analgesics", "Anxiety", "Asthma_and_COPD", "Cancer", "Cardiovascular",
                     "Diabetes", "Dyslipidemia", "Hepatic", "Immunosuppressants", "Siezures", "Sleep", "Thyroid") # drug dispenses

# -- Optimized function to run GLM and extract results
run_glm_analysis <- function(data, column_name, group_var, group_name, reference_group) {
  tryCatch({
    # Create formula
    if (column_name != "SEX") {
      formula_str <- paste0(column_name, " ~ AGE + SEX + ", group_var)
    } else {
      formula_str <- paste0(column_name, " ~ ", group_var)
    }
    
    # Run GLM
    model <- glm(as.formula(formula_str), family = binomial(link = 'logit'), data = data)
    coefs <- coef(summary(model))
    
    # Get confidence intervals more efficiently
    conf_ints <- tryCatch({
      exp(confint.default(model))  # Use default method which is faster
    }, error = function(e) {
      matrix(NA, nrow = length(coef(model)), ncol = 2,
             dimnames = list(names(coef(model)), c("2.5 %", "97.5 %")))
    })
    
    # Create summary dataframe
    result_df <- data.frame(
      Group = rownames(coefs),
      OR = exp(coefs[, "Estimate"]),
      LCI = conf_ints[, "2.5 %"],
      UCI = conf_ints[, "97.5 %"],
      pValue = coefs[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
    
    # Clean group names
    result_df$Group <- gsub(paste0("^", group_var), "", result_df$Group)
    result_df$Reference <- reference_group
    result_df[1, "Group"] <- reference_group
    
    # Calculate group counts more efficiently - avoid duplicate column names
    count_prefix <- if(group_name == "Drug") "Group_N_" else "Class_N_"
    
    # Create count data using a safer approach with temporary column name
    temp_data <- data
    temp_outcome <- paste0("temp_", column_name)
    temp_data[[temp_outcome]] <- temp_data[[column_name]]
    
    count_data <- temp_data %>%
      count(.data[[group_var]], .data[[temp_outcome]], name = "Count_N") %>%
      pivot_wider(names_from = .data[[temp_outcome]], 
                  values_from = Count_N, 
                  names_prefix = count_prefix,
                  values_fill = 0)
    
    # Merge results
    names(count_data)[1] <- "Group"
    result_df <- left_join(result_df, count_data, by = "Group")
    
    return(result_df)
    
  }, error = function(e) {
    cat(sprintf("    Error in model for %s with reference %s: %s\n", 
                column_name, reference_group, e$message))
    return(NULL)
  })
}

# -- Optimized data preprocessing function
preprocess_data <- function(pheno, column_name) {
  # Convert SEX to factor more efficiently
  pheno$SEX <- factor(ifelse(pheno$SEX == "Female", 0, 
                             ifelse(pheno$SEX == "Male", 1, NA)))
  
  # Select required columns - avoid duplicates
  base_cols <- c("ParticipantID", "SEX", "AGE", "DrugName", "DrugClass")
  if (column_name %in% base_cols) {
    required_cols <- base_cols
  } else {
    required_cols <- c(base_cols, column_name)
  }
  
  tab <- pheno[, required_cols, drop = FALSE]
  
  # Filter out missing values
  filter_cols <- if (column_name %in% c("AGE", "SEX")) {
    c("AGE", "SEX", column_name)
  } else {
    c("AGE", "SEX", column_name)
  }
  filter_cols <- unique(filter_cols)  # Remove duplicates
  
  complete_cases <- complete.cases(tab[, filter_cols, drop = FALSE])
  tab <- tab[complete_cases, ]
  
  return(tab)
}

# -- Main processing loop with optimizations
duration_list <- c(360, 600)
adherence_thresholds <- c("All", "Adherent")
all_drug_results <- list()
all_class_results <- list()

for (duration in duration_list) {
  cat(sprintf("\nProcessing duration: %d days\n", duration))
  
  for (adherence in adherence_thresholds) {
    # -- Load treatment groups with readr (faster than base R)
    if (adherence == "All") {
      pheno <- read.csv(file.path(wkdir, paste0("inner_data_", duration, "days.csv")))
    } else {
      pheno <- read.csv(file.path(wkdir, paste0("inner_data_", duration, "days_Adherent.csv")))
    }
    # -- Prepare phenotypes
    pheno$timesincefirstpreg = pheno$AGE - pheno$PREG1AGE
    pheno <- pheno %>%
      mutate(
        earliest_last_end = PREG1AGE + 2 *pmax(CHILDREN, 0),
        latest_last_end = AGE,
        preg_last5_deterministic = case_when(
          CHILDREN == 0 ~ 0,
          earliest_last_end >= AGE - 5 ~ 1,
          TRUE ~ NA
        )
      )
    
    total_traits <- length(exposome_binary)
    
    # Process each trait
    for (i in seq_along(exposome_binary)) {
      column_name <- exposome_binary[i]
      cat(sprintf("[%d/%d] Processing trait: %s\n", i, total_traits, column_name))
      
      # Preprocess data
      tab <- preprocess_data(pheno, column_name)
      class_dat <- tab
      
      # Filter for sustained AD groups for specific traits
      if (column_name %in% c("WELLAD", "STOPAD", "Augmentation_ADspecific")) {
        excluded_drugs <- c("BIP-L", "BIP+L", "Various", "Combination")
        tab <- tab[!tab$DrugName %in% excluded_drugs, ]
        
        excluded_classes <- c("BIP+L", "BIP-L", "Various")
        class_dat <- class_dat[!class_dat$DrugClass %in% excluded_classes, ]
      }
      
      # Get unique groups
      drugs <- unique(tab$DrugName)
      classes <- unique(class_dat$DrugClass)
      
      # Process drug-level models using vectorized approach
      cat("Processing drug-level models:\n")
      drug_results <- lapply(drugs, function(drug) {
        cat(sprintf("  Processing drug: %s\n", drug))
        
        # Set reference level
        tab_copy <- tab
        tab_copy$DrugName <- factor(tab_copy$DrugName, levels = c(drug, setdiff(drugs, drug)))
        
        result <- run_glm_analysis(tab_copy, column_name, "DrugName", "Drug", drug)
        if (!is.null(result)) {
          result$Dependent <- column_name
          result$Threshold <- duration
          result$Adherence <- adherence
        }
        return(result)
      })
      
      # Process class-level models
      cat("Processing class-level models:\n")
      class_results <- lapply(classes, function(class) {
        cat(sprintf("  Processing class: %s\n", class))
        
        # Set reference level
        class_copy <- class_dat
        class_copy$DrugClass <- factor(class_copy$DrugClass, levels = c(class, setdiff(classes, class)))
        
        result <- run_glm_analysis(class_copy, column_name, "DrugClass", "Class", class)
        if (!is.null(result)) {
          result$Dependent <- column_name
          result$Threshold <- duration
          result$Adherence <- adherence
        }
        return(result)
      })
      
      # Store results
      all_drug_results <- c(all_drug_results, drug_results)
      all_class_results <- c(all_class_results, class_results)
    }
  }
}

# Combine all results efficiently
cat("\nCombining results...\n")
EX_glm <- bind_rows(all_drug_results)
class_EX_glm <- bind_rows(all_class_results)

#########################################################
# Optimized result processing and Excel export

# -- Efficient renaming mapping
rename_mapping <- c(
  "AGE" = "Age", 
  "SEX" = "Sex", 
  "AGE2WKF" = "Age of MDD Onset",
  "WELLAD" = "Self-report Responders",
  "STOPAD" = "Self-report Discontinuation", 
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
  "Augmentation_ADspecific" = "AD-specific Augmentation",
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
  "Thyroid" = "Thyroid-related Medications"
)

# Optimized function to process results for Excel
process_results_for_excel <- function(df, data_type) {
  cat(sprintf("Processing %s results for Excel...\n", data_type))
  
  # Rename columns based on data type
  if (data_type == "class") {
    df <- df %>% rename(Reference = Reference, Term = Group)
    n_cols <- paste0("Class_N_", 0:1)
    ref_filter_val <- "SSRI"
  } else {
    df <- df %>% rename(Reference = Reference, Term = Group)
    n_cols <- paste0("Group_N_", 0:1)  # Updated to match new column names
    ref_filter_val <- "SSRI:Sertraline"
  }
  
  # Calculate p-value adjustments efficiently
  terms <- df %>%
    filter(!(Term %in% c("AGE", "SEX1")) | (Dependent %in% c("AGE", "SEX"))) %>%
    group_by(Adherence, Threshold, Reference, Term) %>%
    mutate(
      FDR_P = p.adjust(pValue, method = "fdr"),
      Bonf_P = p.adjust(pValue, method = "bonferroni"),
      Sig_FDR = ifelse(FDR_P < 0.05, "*", ""),
      Sig_Bonf = ifelse(Bonf_P < 0.05, "*", "")
    ) %>%
    ungroup() %>%
    select(Adherence, Threshold, Reference, Dependent, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)
  
  # Calculate total N more efficiently - handle dynamic column names
  available_n_cols <- intersect(n_cols, names(df))
  
  total_n <- df %>%
    filter(Reference == ref_filter_val, Threshold == 360, Adherence == "All") %>%
    rowwise() %>%
    mutate(Total_Participants = if(length(available_n_cols) > 0) {
      sum(c_across(all_of(available_n_cols)), na.rm = TRUE)
    } else {
      0
    }) %>%
    group_by(Dependent) %>%
    summarise(Total_N = sum(Total_Participants, na.rm = TRUE), .groups = "drop")
  
  # Final result processing - handle dynamic column selection
  available_n_cols <- intersect(n_cols, names(df))
  
  results <- df %>%
    left_join(terms, by = c("Adherence", "Threshold", "Reference", "Dependent", "Term")) %>%
    left_join(total_n, by = "Dependent") %>%
    mutate(across(c(OR, LCI, UCI), ~ signif(.x, 2))#,
           #across(c(pValue, FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE))
           ) %>%
    select(Adherence, Threshold, Reference, Dependent, Term, OR, LCI, UCI, pValue, 
           FDR_P, Bonf_P, Sig_FDR, Sig_Bonf, Total_N, all_of(available_n_cols)) %>%
    arrange(Adherence, Threshold, Reference, Dependent, Term) %>%
    group_by(Adherence, Threshold, Reference, Dependent) %>%
    mutate(
      Adherence = ifelse(Adherence != lag(Adherence, default = ""), Adherence, NA_character_),
      Threshold = ifelse(Threshold != lag(Threshold, default = 0), Threshold, NA_integer_),
      Reference = ifelse(Reference != lag(Reference, default = ""), Reference, NA_character_),
      Dependent = ifelse(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
      Total_N = ifelse(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_)
    ) %>%
    ungroup() %>%
    mutate(Dependent = recode(Dependent, !!!rename_mapping),
           Term = ifelse(Term == "SEX1", "SEX", Term))
  
  return(results)
}

# Process and save results
cat("\nloading Excel workbook...\n")
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")

# Process class results
class_results <- process_results_for_excel(class_EX_glm, "class")
removeWorksheet(wb, "Table11")
addWorksheet(wb, "Table11")
writeData(wb, "Table11", class_results)

# Process drug results  
drug_results <- process_results_for_excel(EX_glm, "drug")
removeWorksheet(wb, "Table12")
addWorksheet(wb, "Table12")
writeData(wb, "Table12", drug_results)

# Save workbook
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)
cat("Analysis complete! Results saved to All_Results.xlsx\n")