# -- Load R packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(openxlsx)
  library(broom)
  library(purrr)
})

# -- Set path and duration parameters
wkdir <- "/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes"
output_dir <- "/scratch/user/uqawal15"

# -- Load and process genetic data once
bmi_pgs <- read_table("/QRISdata/Q7280/pharmacogenomics/pgs/scores/BMI_LOO/BMI_LOO_agds_sbrc_gctb_plink2.sscore",
                      show_col_types = FALSE)

eur <- read_table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id", 
                  col_names = c("V1", "V2"), show_col_types = FALSE)

# Filter and standardize PGS data
bmi_pgs_processed <- bmi_pgs %>%
  filter(IID %in% eur$V2) %>%
  mutate(
    std_pgs = as.numeric(scale(SCORE1_SUM)),
    sd_pgs = SCORE1_SUM / sd(SCORE1_SUM, na.rm = TRUE)
  ) %>%
  select(IID, std_pgs, sd_pgs)

# -- Function to run BMI analysis
run_bmi_analysis <- function(data, group_var, reference_groups) {
  filtered_data <- data %>%
    filter(!is.na(AGE), !is.na(SEX), !is.na(sd_pgs), !is.na(BMI), !is.na(.data[[group_var]]))
  
  if (nrow(filtered_data) == 0) return(NULL)
  
  formula_obj <- as.formula(paste0("BMI ~ AGE + SEX + sd_pgs + ", group_var))
  
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
        count(.data[[group_var]], name = "group_n") %>%
        rename(term = .data[[group_var]])
      
      coef_results %>%
        left_join(group_counts, by = "term") %>%
        select(term, estimate, std.error, statistic, p.value, reference_term, outcome, group_n)
      
    }, error = function(e) NULL)
  })
}

# -- Process each duration
process_duration <- function(threshold) {
  pheno_path <- file.path(wkdir, paste0("inner_data_", threshold, "days.csv"))
  pheno <- read_csv(pheno_path, show_col_types = FALSE)
  
  pheno_processed <- pheno %>%
    transmute(
      ParticipantID,
      SEX = case_when(SEX == "Female" ~ 0, SEX == "Male" ~ 1, TRUE ~ NA_real_),
      AGE, DrugName, DrugClass, BMI
    )
  
  combined_data <- pheno_processed %>%
    inner_join(bmi_pgs_processed, by = c("ParticipantID" = "IID"))
  
  unique_drugs <- unique(combined_data$DrugName[!is.na(combined_data$DrugName)])
  unique_classes <- unique(combined_data$DrugClass[!is.na(combined_data$DrugClass)])
  
  drug_results <- run_bmi_analysis(combined_data, "DrugName", unique_drugs)
  if (!is.null(drug_results)) {
    drug_results$threshold <- paste0(threshold, " days")
    drug_results$analysis_level <- "Drug"
  }
  
  class_results <- run_bmi_analysis(combined_data, "DrugClass", unique_classes)
  if (!is.null(class_results)) {
    class_results$threshold <- paste0(threshold, " days")
    class_results$analysis_level <- "Class"
  }
  
  bind_rows(drug_results, class_results)
}

# -- Process all durations
durations <- c(360, 600)
all_duration_results <- map_dfr(durations, process_duration)

# -- Separate results
drug_results <- all_duration_results %>%
  filter(analysis_level == "Drug") %>%
  select(-analysis_level)

class_results <- all_duration_results %>%
  filter(analysis_level == "Class") %>%
  select(-analysis_level)

# -- Process results function
process_results_optimized <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    rename(
      Term = term, ReferenceTerm = reference_term, Outcome = outcome,
      Threshold = threshold, Group_N = group_n, `Std..Error` = std.error,
      `t.value` = statistic, `Pr...t..` = p.value, Estimate = estimate
    ) %>%
    group_by(Threshold, ReferenceTerm, Outcome, Term) %>%
    mutate(
      FDR_P = p.adjust(`Pr...t..`, method = "fdr"),
      Bonf_P = p.adjust(`Pr...t..`, method = "bonferroni"),
      Sig_FDR = ifelse(FDR_P < 0.05, "*", ""),
      Sig_Bonf = ifelse(Bonf_P < 0.05, "*", "")
    ) %>%
    ungroup() %>%
    mutate(
      across(c(`Pr...t..`, FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE)),
      across(c(`Std..Error`), ~ signif(.x, 2)),
      across(c(Estimate, `t.value`), ~ round(.x, 2))
    ) %>%
    select(Threshold, ReferenceTerm, Outcome, Term, Estimate, `Std..Error`, 
           `t.value`, `Pr...t..`, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf, Group_N) %>%
    arrange(Threshold, ReferenceTerm, Outcome, Term) %>%
    group_by(Threshold, ReferenceTerm, Outcome) %>%
    mutate(
      Threshold = ifelse(Threshold != lag(Threshold, default = ""), Threshold, NA_character_),
      ReferenceTerm = ifelse(ReferenceTerm != lag(ReferenceTerm, default = ""), ReferenceTerm, NA_character_),
      Outcome = ifelse(Outcome != lag(Outcome, default = ""), Outcome, NA_character_)
    ) %>%
    ungroup()
}

# -- Process final results
drug_results_final <- process_results_optimized(drug_results)
class_results_final <- process_results_optimized(class_results)
ls

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