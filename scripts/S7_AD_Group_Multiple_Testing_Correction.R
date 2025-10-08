# ===============================================================================
# MEMORY-EFFICIENT COMBINED CORRECTION ACROSS QUANTITATIVE AND BINARY FEATURES
# ===============================================================================

library(dplyr)
library(openxlsx)

# Set paths
output_dir <- "C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability"
excel_file <- file.path(output_dir, "All_Results.xlsx")

# ===============================================================================
# STEP 1: READ AND PROCESS DATA IN CHUNKS
# ===============================================================================

cat("Reading data from Excel sheets...\n")

# Function to read and immediately clean data
read_and_clean <- function(file, sheet) {
  df <- read.xlsx(file, sheet = sheet)
  # Remove empty rows and clean data
  df %>%
    filter(if_any(everything(), ~ !is.na(.x) & .x != "")) %>%
    mutate(across(everything(), ~ ifelse(.x == "", NA, .x)))
}

# Read data
quant_class_data <- read_and_clean(excel_file, "Table9") %>%
  mutate(row_id = row_number())
quant_drug_data <- read_and_clean(excel_file, "Table10") %>%
  mutate(row_id = row_number())
binary_class_data <- read_and_clean(excel_file, "Table11") %>%
  mutate(row_id = row_number())
binary_drug_data <- read_and_clean(excel_file, "Table12") %>%
  mutate(row_id = row_number())

# Debug data types
cat("=== DATA TYPE DEBUG INFO ===\n")
cat("Quant class - adherence:", class(quant_class_data$adherence), "threshold:", class(quant_class_data$threshold), "\n")
cat("Quant drug - adherence:", class(quant_drug_data$adherence), "threshold:", class(quant_drug_data$threshold), "\n")
cat("Binary class - Adherence:", class(binary_class_data$Adherence), "Threshold:", class(binary_class_data$Threshold), "\n")
cat("Binary drug - Adherence:", class(binary_drug_data$Adherence), "Threshold:", class(binary_drug_data$Threshold), "\n")
cat("===============================\n")

# ===============================================================================
# STEP 2: MEMORY-EFFICIENT PROCESSING FUNCTION
# ===============================================================================

process_corrections_efficiently <- function(quant_data, binary_data, data_type) {
  cat(sprintf("Processing %s data...\n", data_type))
  
  # Extract only needed columns for correction
  if (data_type == "class") {
    quant_subset <- quant_data %>%
      filter(!is.na(p.value)) %>%
      select(adherence, threshold, reference_term, term, outcome, p.value, row_id) %>%
      filter(!(term %in% c("AGE", "SEX")) | (outcome %in% c("Age", "AGE"))) %>%
      mutate(
        analysis_type = "quantitative",
        pvalue_num = as.numeric(p.value),
        adherence = as.character(adherence),
        threshold = as.character(threshold),
        reference_term = as.character(reference_term),
        term = as.character(term)
      )
    
    binary_subset <- binary_data %>%
      filter(!is.na(pValue)) %>%
      select(Adherence, Threshold, Reference, Term, Dependent, pValue, row_id) %>%
      filter(!(Term %in% c("AGE", "SEX1", "SEX")) | (Dependent %in% c("Age", "Sex", "AGE", "SEX"))) %>%
      mutate(
        analysis_type = "binary",
        pvalue_num = as.numeric(pValue),
        adherence = as.character(Adherence),
        threshold = as.character(Threshold),
        reference_term = as.character(Reference),
        term = as.character(Term),
        outcome = as.character(Dependent)
      ) %>%
      select(analysis_type, adherence, threshold, reference_term, term, outcome, pvalue_num, row_id)
    
  } else {  # drug
    quant_subset <- quant_data %>%
      filter(!is.na(p.value)) %>%
      select(adherence, threshold, reference_term, term, outcome, p.value, row_id) %>%
      filter(!(term %in% c("AGE", "SEX")) | (outcome %in% c("Age", "AGE"))) %>%
      mutate(
        analysis_type = "quantitative",
        pvalue_num = as.numeric(p.value),
        adherence = as.character(adherence),
        threshold = as.character(threshold),
        reference_term = as.character(reference_term),
        term = as.character(term)
      )
    
    binary_subset <- binary_data %>%
      filter(!is.na(pValue)) %>%
      select(Adherence, Threshold, Reference, Term, Dependent, pValue, row_id) %>%
      filter(!(Term %in% c("AGE", "SEX1", "SEX")) | (Dependent %in% c("Age", "Sex", "AGE", "SEX"))) %>%
      mutate(
        analysis_type = "binary",
        pvalue_num = as.numeric(pValue),
        adherence = as.character(Adherence),
        threshold = as.character(Threshold),
        reference_term = as.character(Reference),
        term = as.character(Term),
        outcome = as.character(Dependent)
      ) %>%
      select(analysis_type, adherence, threshold, reference_term, term, outcome, pvalue_num, row_id)
  }
  
  # Combine for correction
  combined_subset <- bind_rows(quant_subset, binary_subset)
  
  # Apply correction
  corrected_subset <- combined_subset %>%
    group_by(adherence, threshold, reference_term, term) %>%
    mutate(
      unified_fdr_p = p.adjust(pvalue_num, method = "fdr"),
      unified_bonf_p = p.adjust(pvalue_num, method = "bonferroni"),
      unified_sig_fdr = ifelse(unified_fdr_p < 0.05, "*", ""),
      unified_sig_bonf = ifelse(unified_bonf_p < 0.05, "*", "")
    ) %>%
    ungroup()
  
  return(corrected_subset)
}

# Process each data type separately
class_corrections <- process_corrections_efficiently(quant_class_data, binary_class_data, "class")
drug_corrections <- process_corrections_efficiently(quant_drug_data, binary_drug_data, "drug")

# ===============================================================================
# STEP 3: FORMAT CORRECTIONS 
# ===============================================================================

#-- Class quantitative corrections
quant_class_corrections <- class_corrections %>%
  filter(analysis_type == "quantitative") %>%
  select(unified_fdr_p, unified_bonf_p, unified_sig_fdr, unified_sig_bonf, row_id)

quant_class_temp <- quant_class_data %>%
  left_join(quant_class_corrections, 
            by = c("row_id")) %>%
  arrange(row_id) %>%
  select(-row_id)

quant_class_combined <- quant_class_temp %>%
  select(-fdr_p, - bonf_p, -sig_fdr, -sig_bonf) %>%
  rename(
    fdr_p = unified_fdr_p, 
    bonf_p = unified_bonf_p, 
    sig_fdr = unified_sig_fdr, 
    sig_bonf = unified_sig_bonf
  ) %>%
  mutate(across(c(fdr_p, bonf_p), ~ format(signif(.x, 2), scientific = TRUE)))

quant_class_combined_filtered <- quant_class_combined %>%
  fill(adherence, threshold, reference_term, outcome) %>%
  filter((adherence == "All" & threshold == "360 days") | 
           (adherence == "Adherent" & threshold == "360 days" & reference_term == "SSRI") |
           (adherence == "All" & threshold == "600 days" & reference_term == "SSRI")) %>%
  mutate(adherence = factor(adherence, levels = c("All", "Adherent"))) %>%
  arrange(adherence) %>%
  group_by(adherence, threshold, reference_term, outcome) %>%
  mutate(adherence = ifelse(as.character(adherence) != lag(as.character(adherence), default = ""), 
                                    as.character(adherence), 
                                    NA_character_),
         threshold = ifelse(threshold != lag(threshold, default = ""), threshold, NA_character_),
         reference_term = ifelse(reference_term != lag(reference_term, default = ""), reference_term, NA_character_),
         outcome = ifelse(outcome != lag(outcome, default = ""), outcome, NA_character_))


wb <- loadWorkbook(file.path(output_dir, "All_Results.xlsx"))
removeWorksheet(wb, "Table9Unified")
addWorksheet(wb, "Table9Unified")
writeData(wb, "Table9Unified", quant_class_combined_filtered)
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)

#=========================================

#-- Class binary corrections
binary_class_corrections <- class_corrections %>%
  filter(analysis_type != "quantitative") %>%
  select(unified_fdr_p, unified_bonf_p, unified_sig_fdr, unified_sig_bonf, row_id)

binary_class_temp <- binary_class_data %>%
  left_join(binary_class_corrections, 
            by = c("row_id")) %>%
  arrange(row_id) %>%
  select(-row_id)

binary_class_combined <- binary_class_temp %>%
  select(-FDR_P, -Bonf_P, -Sig_FDR, -Sig_Bonf) %>%
  rename(
    FDR_P = unified_fdr_p, 
    Bonf_P = unified_bonf_p, 
    Sig_FDR = unified_sig_fdr, 
    Sig_Bonf = unified_sig_bonf
  ) %>%
  mutate(across(c(FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE)))

binary_class_combined_filtered <- binary_class_combined %>%
  fill(Adherence, Threshold, Reference, Dependent) %>%
  filter((Adherence == "All" & Threshold == "360") | 
           (Adherence == "Adherent" & Threshold == "360" & Reference == "SSRI") |
           (Adherence == "All" & Threshold == "600" & Reference == "SSRI")) %>%
  mutate(Adherence = factor(Adherence, levels = c("All", "Adherent"))) %>%
  arrange(Adherence) %>%
  group_by(Adherence, Threshold, Reference, Dependent) %>%
  mutate(Adherence = ifelse(as.character(Adherence) != lag(as.character(Adherence), default = ""), 
                            as.character(Adherence), 
                            NA_character_),
         Threshold = ifelse(Threshold != lag(Threshold, default = 0), Threshold, NA_integer_),
         Reference = ifelse(Reference != lag(Reference, default = ""), Reference, NA_character_),
         Dependent = ifelse(Dependent != lag(Dependent, default = ""), Dependent, NA_character_))

wb <- loadWorkbook(file.path(output_dir, "All_Results.xlsx"))
removeWorksheet(wb, "Table11Unified")
addWorksheet(wb, "Table11Unified")
writeData(wb, "Table11Unified", binary_class_combined_filtered)
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)

#=======================================

#-- drug quantitative corrections
quant_drug_corrections <- drug_corrections %>%
  filter(analysis_type == "quantitative") %>%
  select(unified_fdr_p, unified_bonf_p, unified_sig_fdr, unified_sig_bonf, row_id)

quant_drug_temp <- quant_drug_data %>%
  left_join(quant_drug_corrections, 
            by = c("row_id")) %>%
  arrange(row_id) %>%
  select(-row_id)
quant_drug_temp$p.value <- as.numeric(quant_drug_temp$p.value)
quant_drug_temp$unified_fdr_p <- as.numeric(quant_drug_temp$unified_fdr_p)
quant_drug_temp$unified_bonf_p <- as.numeric(quant_drug_temp$unified_bonf_p)


quant_drug_combined <- quant_drug_temp %>%
  select(-fdr_p, - bonf_p, -sig_fdr, -sig_bonf) %>%
  rename(
    fdr_p = unified_fdr_p, 
    bonf_p = unified_bonf_p, 
    sig_fdr = unified_sig_fdr, 
    sig_bonf = unified_sig_bonf
  ) %>%
  mutate(across(c(p.value, fdr_p, bonf_p), ~ sprintf("%.1e", .x)))

quant_drug_combined_filtered <- quant_drug_combined %>%
  fill(adherence, threshold, reference_term, outcome) %>%
  filter((adherence == "All" & threshold == "360 days") | 
           (adherence == "Adherent" & threshold == "360 days" & reference_term == "SSRI:Sertraline") |
           (adherence == "All" & threshold == "600 days" & reference_term == "SSRI:Sertraline")) %>%
  mutate(adherence = factor(adherence, levels = c("All", "Adherent"))) %>%
  arrange(adherence) %>%
  group_by(adherence, threshold, reference_term, outcome) %>%
  mutate(adherence = ifelse(as.character(adherence) != lag(as.character(adherence), default = ""), 
                            as.character(adherence), 
                            NA_character_),
         threshold = ifelse(threshold != lag(threshold, default = ""), threshold, NA_character_),
         reference_term = ifelse(reference_term != lag(reference_term, default = ""), reference_term, NA_character_),
         outcome = ifelse(outcome != lag(outcome, default = ""), outcome, NA_character_))



wb <- loadWorkbook(file.path(output_dir, "All_Results.xlsx"))
removeWorksheet(wb, "Table10Unified")
addWorksheet(wb, "Table10Unified")
writeData(wb, "Table10Unified", quant_drug_combined_filtered)
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)

#=========================================

#-- Drug binary corrections
binary_drug_corrections <- drug_corrections %>%
  filter(analysis_type != "quantitative") %>%
  select(unified_fdr_p, unified_bonf_p, unified_sig_fdr, unified_sig_bonf, row_id)

binary_drug_temp <- binary_drug_data %>%
  left_join(binary_drug_corrections, 
            by = c("row_id")) %>%
  arrange(row_id) %>%
  select(-row_id)

binary_drug_combined <- binary_drug_temp %>%
  select(-FDR_P, -Bonf_P, -Sig_FDR, -Sig_Bonf) %>%
  rename(
    FDR_P = unified_fdr_p, 
    Bonf_P = unified_bonf_p, 
    Sig_FDR = unified_sig_fdr, 
    Sig_Bonf = unified_sig_bonf
  ) %>%
  mutate(across(c(FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE)))

binary_drug_combined_filtered <- binary_drug_combined %>%
  fill(Adherence, Threshold, Reference, Dependent) %>%
  filter((Adherence == "All" & Threshold == "360") | 
           (Adherence == "Adherent" & Threshold == "360" & Reference == "SSRI:Sertraline") |
           (Adherence == "All" & Threshold == "600" & Reference == "SSRI:Sertraline")) %>%
  mutate(Adherence = factor(Adherence, levels = c("All", "Adherent"))) %>%
  arrange(Adherence) %>%
  group_by(Adherence, Threshold, Reference, Dependent) %>%
  mutate(Adherence = ifelse(as.character(Adherence) != lag(as.character(Adherence), default = ""), 
                            as.character(Adherence), 
                            NA_character_),
         Threshold = ifelse(Threshold != lag(Threshold, default = 0), Threshold, NA_integer_),
         Reference = ifelse(Reference != lag(Reference, default = ""), Reference, NA_character_),
         Dependent = ifelse(Dependent != lag(Dependent, default = ""), Dependent, NA_character_))



wb <- loadWorkbook(file.path(output_dir, "All_Results.xlsx"))
removeWorksheet(wb, "Table12Unified")
addWorksheet(wb, "Table12Unified")
writeData(wb, "Table12Unified", binary_drug_combined_filtered)
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)

