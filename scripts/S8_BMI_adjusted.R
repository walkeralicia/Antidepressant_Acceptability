# -- Load R packages
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
output_dir <- "/scratch/user/uqawal15"
duration_list <- c(360, 600)
adherence_thresholds <- c("All", "Adherent")  # Add adherence groups

# -- Load and process genetic data once
bmi_pgs <- read.table("/QRISdata/Q7280/pharmacogenomics/pgs/scores/BMI_LOO/BMI_LOO_agds_sbrc_gctb_plink2.sscore")

eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id")

pcs <- read.table("/QRISdata/Q5338/Ancestry_analysis/AGDS_R11_TOPMedr2_pruned.05.common_pca3.proj.eigenvec") %>%
  select(-V1) %>%
  rename(
    IID = V2, 
    PC1 = V3,
    PC2 = V4, 
    PC3 = V5
  )

# Filter and standardize PGS data
bmi_pgs_processed <- bmi_pgs %>%
  filter(V2 %in% eur$V2) %>%
  mutate(
    std_pgs = as.numeric(scale(V6)),
    sd_pgs = V6 / sd(V6, na.rm = TRUE)
  ) %>%
  select(V2, std_pgs, sd_pgs) %>%
  left_join(pcs, by = c("V2" = "IID"))

# -- Function to run BMI analysis
run_bmi_analysis <- function(data, group_var, reference_groups, duration, adherence) {
  filtered_data <- data %>%
    filter(!is.na(AGE), !is.na(SEX), !is.na(sd_pgs), !is.na(BMI), !is.na(!!sym(group_var)))
  
  if (nrow(filtered_data) == 0) return(NULL)
  
  formula_obj <- as.formula(paste0("BMI ~ AGE + SEX + PC1 + PC2 + PC3 + sd_pgs + ", group_var))
  
  map_dfr(reference_groups, function(ref_group) {
    tryCatch({
      temp_data <- filtered_data
      temp_data[[group_var]] <- factor(temp_data[[group_var]], 
                                       levels = c(ref_group, setdiff(unique(temp_data[[group_var]]), ref_group)))
      
      model <- lm(formula_obj, data = temp_data)
      
      coef_results <- tidy(model) %>%
        mutate(
          term = gsub(paste0("^", group_var), "", term),
          reference_term = ref_group,
          outcome = "BMI"
        )
      
      coef_results$term[1] <- ref_group
      
      group_counts <- temp_data %>%
        count(!!sym(group_var), name = "group_n") %>%
        rename(term = !!sym(group_var))
      
      coef_results %>%
        left_join(group_counts, by = "term") %>%
        mutate(
          threshold = paste0(duration, " days"),
          adherence = adherence
        ) %>%
        select(adherence, threshold, term, estimate, std.error, statistic, p.value, reference_term, outcome, group_n)
      
    }, error = function(e) NULL)
  })
}

# -- Process each duration and adherence combination
process_duration_adherence <- function(duration, adherence) {
  cat("Processing duration:", duration, "days, adherence:", adherence, "\n")
  
  # Load appropriate phenotype data based on adherence
  if (adherence == "All") {
    pheno_path <- file.path(wkdir, paste0("inner_data_", duration, "days.csv"))
  } else {
    pheno_path <- file.path(wkdir, paste0("inner_data_", duration, "days_Adherent.csv"))
  }
  
  pheno <- read.csv(pheno_path)
  
  pheno_processed <- pheno %>%
    transmute(
      IID,
      ParticipantID,
      SEX = case_when(SEX == "Female" ~ 0, SEX == "Male" ~ 1, TRUE ~ NA_real_),
      AGE, DrugName, DrugClass, BMI
    )
  
  combined_data <- pheno_processed %>%
    inner_join(bmi_pgs_processed, by = c("IID" = "V2"))
  
  unique_drugs <- unique(combined_data$DrugName[!is.na(combined_data$DrugName)])
  unique_classes <- unique(combined_data$DrugClass[!is.na(combined_data$DrugClass)])
  
  drug_results <- run_bmi_analysis(combined_data, "DrugName", unique_drugs, duration, adherence)
  if (!is.null(drug_results)) {
    drug_results$analysis_level <- "Drug"
  }
  
  class_results <- run_bmi_analysis(combined_data, "DrugClass", unique_classes, duration, adherence)
  if (!is.null(class_results)) {
    class_results$analysis_level <- "Class"
  }
  
  bind_rows(drug_results, class_results)
}

# -- Process all combinations using nested map
all_results <- map_dfr(duration_list, function(duration) {
  map_dfr(adherence_thresholds, function(adherence) {
    process_duration_adherence(duration, adherence)
  })
})

# -- Separate results
drug_results <- all_results %>%
  filter(analysis_level == "Drug") %>%
  select(-analysis_level)

class_results <- all_results %>%
  filter(analysis_level == "Class") %>%
  select(-analysis_level)

# -- Process results function (updated for adherence)
process_results_optimized <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df <- df %>%
    rename(
      Term = term, ReferenceTerm = reference_term, Outcome = outcome,
      Threshold = threshold, Group_N = group_n, `Std..Error` = std.error,
      `t.value` = statistic, `Pr...t..` = p.value, Estimate = estimate,
      Adherence = adherence
    )
  
  df %>% 
    mutate(
      across(c(`Pr...t..`), ~ format(signif(.x, 2), scientific = TRUE)),
      across(c(`Std..Error`), ~ signif(.x, 2)),
      across(c(Estimate, `t.value`), ~ round(.x, 2))
    ) %>%
    select(Adherence, Threshold, ReferenceTerm, Outcome, Term, Estimate, `Std..Error`, 
           `t.value`, `Pr...t..`, Group_N) %>%
    filter((Adherence == "All" & Threshold == "360 days") | 
             (Adherence == "Adherent" & Threshold == "360 days" & ReferenceTerm == "SSRI:Sertraline") |
             (Adherence == "All" & Threshold == "600 days" & ReferenceTerm == "SSRI:Sertraline")) %>%
    mutate(Adherence = factor(Adherence, levels = c("All", "Adherent"))) %>%
    arrange(Adherence, Threshold, ReferenceTerm, Outcome, Term) %>%
    group_by(Adherence, Threshold, ReferenceTerm, Outcome) %>%
    mutate(
      Adherence = ifelse(as.character(Adherence) != lag(as.character(Adherence), default = ""), 
                         as.character(Adherence), 
                         NA_character_),
      Threshold = ifelse(Threshold != lag(Threshold, default = ""), Threshold, NA_character_),
      ReferenceTerm = ifelse(ReferenceTerm != lag(ReferenceTerm, default = ""), ReferenceTerm, NA_character_),
      Outcome = ifelse(Outcome != lag(Outcome, default = ""), Outcome, NA_character_)
    ) %>%
    ungroup()
}

# -- Process final results
drug_results_final <- process_results_optimized(drug_results)
class_results_final <- process_results_optimized(class_results)

# -- Save to Excel
wb_path <- file.path(output_dir, "All_Results.xlsx")
wb <- if (file.exists(wb_path)) loadWorkbook(wb_path) else createWorkbook()

safe_add_worksheet <- function(wb, sheet_name, data_to_write) {
  if (is.null(data_to_write) || nrow(data_to_write) == 0) return(FALSE)
  
  if (sheet_name %in% names(wb)) removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data_to_write)
  TRUE
}

safe_add_worksheet(wb, "Table15", class_results_final)
safe_add_worksheet(wb, "Table16", drug_results_final)

saveWorkbook(wb, wb_path, overwrite = TRUE)

cat("BMI PGS analysis complete! Results saved for both All and Adherent groups.\n")