# -- Load libraries
library(emmeans)
library(openxlsx)
library(broom)
library(dplyr)
library(stringr)
library(purrr)

# -- Set paths and parameters
wkdir <- "/QRISdata/Q7280/pharmacogenomics"
duration_list <- c(360, 600)
adherence_thresholds <- c("All", "Adherent")  # Add adherence groups

# -- List of Europeans
eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id", header = FALSE)
eur_ids <- as.character(eur$V2)

#-- PC file
pcs <- read.table("/QRISdata/Q5338/Ancestry_analysis/AGDS_R11_TOPMedr2_pruned.05.common_pca3.proj.eigenvec") %>%
  select(-V1) %>%
  rename(
    IID = V2, 
    PC1 = V3,
    PC2 = V4, 
    PC3 = V5
  )

# -- Find PGS files
pgs_codes <- c("PUD_01", "T2D_03", "CNT_03", "ADHD_01", "BIP_LOO", "BMI_LOO", "MDD_LOO", 
               "SCZ_03", "Migraine_01", "SBP_01", "ANO_LOO", 
               "UKB_35BM_2021_LDL_direct_adjstatins", "CRP_01", "Neuroticism_01", "LRA_01")
pgslist <- system('find /QRISdata/Q7280/pharmacogenomics/pgs -name "*agds_sbrc_gctb_plink2.sscore" -type f -exec ls {} \\;', intern = TRUE)
selected_pgs <- pgslist[grep(paste(pgs_codes, collapse = "|"), pgslist)]

# -- Functions for processing PGS files
read_pgs_file <- function(file_path, j) {
  # Standard PGS file format
  pgs <- read.table(file_path, header = FALSE)
  
  # Extract trait information
  pgs_filename <- basename(file_path)
  trait_parts <- str_split(str_split(pgs_filename, "\\.")[[1]][1], "_")[[1]]
  trait <- paste0(trait_parts[1], "_", trait_parts[2])
  
  # Filter for Europeans and standardize
  pgs_filtered <- pgs %>%
    filter(V2 %in% eur_ids) %>%
    mutate(V6 = as.numeric(V6),
           std_pgs = scale(V6)[,1])
  
  # Add PCs to PGS data
  pgs_filtered_pcs <- pgs_filtered %>%
    left_join(pcs, by = c("V2" = "IID"))
  
  return(list(
    data = pgs_filtered_pcs,
    trait = trait,
    id_col = "V2"
  ))
}

# Function to run models for drugs and classes
run_models <- function(pgs_data, trait, ad_data, duration, adherence) {
  # Identify join column
  id_col <- pgs_data$id_col
  
  # Join PGS with treatment data
  if (id_col == "V2") {
    pgs_atc <- ad_data %>%
      inner_join(pgs_data$data, by = c("IID" = "V2"))
  } else {
    pgs_atc <- ad_data %>%
      inner_join(pgs_data$data, by = "IID")
  }
  
  # Get unique drugs and classes
  drugs <- unique(pgs_atc$DrugName)
  classes <- c(unique(pgs_atc$DrugClass))
  
  # Function to create model results
  create_model_results <- function(model_data, reference_col, reference_val, grouping_col) {
    # Set reference level
    model_data[[reference_col]] <- factor(model_data[[reference_col]])
    model_data[[reference_col]] <- relevel(model_data[[reference_col]], ref = reference_val)
    
    # Run linear model
    formula <- as.formula(paste("std_pgs ~ PC1 + PC2 + PC3 + ", reference_col))
    model <- lm(formula, data = model_data)
    
    # Extract coefficients
    coefs <- tidy(model) %>%
      rename(!!grouping_col := term) %>%
      mutate(
        !!grouping_col := str_replace(!!sym(grouping_col), paste0(reference_col), ""),
        !!grouping_col := if_else(!!sym(grouping_col) == "(Intercept)", reference_val, !!sym(grouping_col))
      )
    
    # Get participant counts per group
    counts <- model_data %>%
      count(!!sym(reference_col), name = "Group_N") %>%
      rename(!!grouping_col := !!sym(reference_col))
    
    # Combine results
    results <- coefs %>%
      full_join(counts, by = grouping_col) %>%
      mutate(
        Reference = reference_val,
        PGS = trait,
        Threshold = paste0(duration, " days"),
        Adherence = adherence  # Add adherence info
      )
    
    return(results)
  }
  
  # Run models for each drug
  drug_results <- map_dfr(drugs, ~{
    create_model_results(pgs_atc, "DrugName", .x, "DrugName")
  })
  
  # Run models for each class
  class_results <- map_dfr(classes, ~{
    create_model_results(pgs_atc, "DrugClass", .x, "DrugClass")
  })
  
  return(list(drug_results = drug_results, class_results = class_results))
}

# -- Initialize lists for result storage
all_results <- tibble()
all_class_results <- tibble()

# -- Main processing loop (NOW INCLUDES BOTH ADHERENCE GROUPS)
for (duration in duration_list) {
  cat("Processing duration:", duration, "days\n")
  
  # Process both adherence thresholds
  for (adherence in adherence_thresholds) {
    cat("  Processing adherence group:", adherence, "\n")
    
    # Load appropriate treatment data based on adherence
    if (adherence == "All") {
      ad_data <- read.csv(paste0("/QRISdata/Q7280/pharmacogenomics/treatment_groups/Final_Treatmentgroups_", duration, "days.csv"))
    } else {
      ad_data <- read.csv(paste0("/QRISdata/Q7280/pharmacogenomics/treatment_groups/Final_Treatmentgroups_", duration, "days_Adherent.csv"))
    }
    
    # Process each PGS trait
    for (j in 1:length(selected_pgs)) {
      # Read and prepare PGS data (same for both adherence groups)
      pgs_file_path <- selected_pgs[j]
      pgs_data <- read_pgs_file(pgs_file_path, j)
      cat("    Processing trait:", pgs_data$trait, "(", j, ")\n")
      
      # Run models and collect results
      model_results <- run_models(pgs_data, pgs_data$trait, ad_data, duration, adherence)
      
      # Add to result collection
      all_results <- bind_rows(all_results, model_results$drug_results)
      all_class_results <- bind_rows(all_class_results, model_results$class_results)
    }
  }
}

# -- Process drug-level results (now grouped by adherence too)
results_lm <- all_results %>%
  arrange(Adherence, Threshold, PGS, DrugName)

# Calculate participant counts per threshold/PGS/adherence combination
num_dependent <- results_lm %>%
  filter(Reference == "SSRI:Escitalopram") %>%
  group_by(Adherence, Threshold, PGS) %>%
  summarise(Total_N = sum(Group_N, na.rm = TRUE), .groups = "drop")

# Join with main results
EX_labeled <- results_lm %>%
  left_join(num_dependent, by = c("Adherence", "Threshold", "PGS"))

# -- Process class-level results (now grouped by adherence too)
results_class_lm <- all_class_results %>%
  arrange(Adherence, Threshold, PGS, DrugClass)

# Calculate participant counts per threshold/PGS/adherence combination
class_num_dependent <- results_class_lm %>%
  filter(Reference == "SSRI") %>%
  group_by(Adherence, Threshold, PGS) %>%
  summarise(Total_N = sum(Group_N, na.rm = TRUE), .groups = "drop")

# Join with main results
class_EX_labeled <- results_class_lm %>%
  left_join(class_num_dependent, by = c("Adherence", "Threshold", "PGS"))

# -- Updated function to add results to Excel workbook
process_and_add_to_workbook <- function(data_type) {
  # Read appropriate file
  if (data_type == "class") {
    df <- class_EX_labeled %>%
      rename(Term = DrugClass)
    sheet_name <- "Table13"
  } else {
    df <- EX_labeled %>%
      rename(Term = DrugName)
    sheet_name <- "Table14" 
  }
  
  # Calculate FDR and Bonferroni p-values (now grouped by adherence too)
  terms <- df %>%
    filter(!Term %in% c("PC1", "PC2", "PC3")) %>%
    group_by(Adherence, Threshold, Reference, Term) %>%
    mutate(
      FDR_P = p.adjust(p.value, method = "fdr"),
      Bonf_P = p.adjust(p.value, method = "bonferroni"),
      Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
      Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
    ) %>%
    select(Adherence, Threshold, Reference, PGS, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf) %>%
    ungroup()
  
  # Format results
  results <- df %>%
    left_join(terms, by = c("Adherence", "Threshold", "Reference", "PGS", "Term")) %>%
    arrange(Adherence, Threshold, Reference, PGS, Term) %>%
    mutate(across(c(p.value, FDR_P, Bonf_P), ~format(signif(.x, 2), scientific = TRUE)),
           across(c(estimate, statistic), ~round(., 3)),
           across(c(std.error), ~signif(., 2))) %>%
    select(Adherence, Threshold, Reference, PGS, Term, estimate, std.error, statistic, 
           p.value, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf, Total_N, Group_N)
  
  # Format for display (handling consecutive duplicates)
  results_formatted <- results %>%
    filter((Adherence == "All" & Threshold == "360 days") | 
             (Adherence == "Adherent" & Threshold == "360 days" & Reference %in%  c("SSRI", "SSRI:Sertraline")) |
             (Adherence == "All" & Threshold == "600 days" & Reference %in% c("SSRI", "SSRI:Sertraline"))) %>%
    mutate(Adherence = factor(Adherence, levels = c("All", "Adherent"))) %>%
    arrange(Adherence) %>%
    group_by(Adherence, Threshold, Reference, PGS) %>%
    mutate(
      Adherence = ifelse(as.character(Adherence) != lag(as.character(Adherence), default = ""), 
                         as.character(Adherence), 
                         NA_character_),
      Threshold = if_else(Threshold != lag(Threshold, default = ""), Threshold, NA_character_),
      Reference = if_else(Reference != lag(Reference, default = ""), Reference, NA_character_),
      PGS = if_else(PGS != lag(PGS, default = ""), PGS, NA_character_),
      Total_N = if_else(Total_N != lag(Total_N, default = 0), Total_N, NA_integer_)
    )
  
  # Rename PGS for display
  results_renamed <- results_formatted %>%
    mutate(PGS = str_extract(PGS, "^[^_]+")) %>%
    mutate(PGS = ifelse(PGS == "UKB", "LDL-c", PGS))
  
  # Add to workbook
  if (sheet_name %in% names(wb)) {
    removeWorksheet(wb, sheet_name)
  }
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, results_renamed)
  
  return(wb)
}

# Process both data types and save workbook
wb <- loadWorkbook("All_Results.xlsx")
wb <- process_and_add_to_workbook("class")
wb <- process_and_add_to_workbook("drug")
saveWorkbook(wb, "/scratch/user/uqawal15/All_Results.xlsx", overwrite = TRUE)

cat("Analysis complete! Results saved for both All and Adherent groups.\n")