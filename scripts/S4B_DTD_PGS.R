# -- Load packages
library(dplyr)
library(purrr)
library(stringr)
library(openxlsx)
library(progressr)
library(emmeans)
library(broom)

# -- Set working directory
wkdir <- "/QRISdata/Q7280/pharmacogenomics"
output_dir <- "/scratch/user/uqawal15"

# -- Load treatment group data and extract AGE, SNPSEX and PC information
pcs <- read.table("/QRISdata/Q5338/Ancestry_analysis/AGDS_R11_TOPMedr2_pruned.05.common_pca3.proj.eigenvec") %>%
  select(-V1) %>%
  rename(
    IID = V2, 
    PC1 = V3,
    PC2 = V4, 
    PC3 = V5
  )
dat <- read.csv(file.path(wkdir, "data/AGDSAcceptabilityTreatmentGroups_25082025.csv"))
pheno <- read.csv(file.path(wkdir, "phenotypes/survey_phenotypes.csv")) %>%
  select(STUDYID, AGE, SEX, BMI, DXBPD2, MDD) %>% 
  rename(ParticipantID = STUDYID)

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

# -- Genetic data setup
pgslist <- system('find /QRISdata/Q7280/pharmacogenomics/pgs -name "*agds_sbrc_gctb_plink2.sscore" -type f -exec ls {} \\;', intern = TRUE)
pgs_codes <- c("PUD_01", "T2D_03", "CNT_03", "ADHD_01", "BIP_LOO", "BMI_LOO", "MDD_LOO", "SCZ_03", "Migraine_01", "SBP_01",  "ANO_LOO", "UKB_35BM_2021_LDL_direct_adjstatins", "CRP_01", "Neuroticism_01", "LRA_01")
selected_pgs <- pgslist[grep(paste(pgs_codes, collapse = "|"), pgslist)]

# -- List of Europeans
eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id")

# -- Linkage file for those with genetic data
link <- read.table(file.path("/QRISdata/Q7280/pharmacogenomics", "data/StudyID_GenoID_link.txt"), header = TRUE) %>%
  rename("ParticipantID" = "StudyCodeID") %>%
  distinct()

#=========== FUNCTION TO RUN PGS ANALYSIS =================================================================

run_pgs_analysis <- function(exclude_bip = FALSE) {
  
  # Start with the base dataset
  pharma_full <- pharma %>%
    left_join(pheno, by = "ParticipantID") %>%
    filter(MDD==1)
  
  # Apply filters based on parameters
  if (exclude_bip) {
    pharma_full <- pharma_full %>%
      filter(DXBPD2 != 1) %>%
      select(-DXBPD2)
    analysis_label <- "BIP_Excluded"
  } else {
    pharma_full <- pharma_full %>%
      select(-DXBPD2)
    analysis_label <- "Full_Sample"
  }
  
  # -- Convert SEX
  pharma_full <- pharma_full %>%
    mutate(SEX = case_when(
      SEX == "Female" ~ 0,
      SEX == "Male" ~ 1,
      TRUE ~ NA_real_
    ))
  
  # -- Set up results list
  total_lm_list <- list()
  numAD_lm_list <- list()
  numClass_lm_list <- list()
  
  # -- For each PGS trait fit a linear model
  for (j in 1:length(selected_pgs)) {
    
    # -- Read in PGS file for the selected trait
    pgs_file_path <- selected_pgs[j]
    pgs <- read.table(selected_pgs[j], header = FALSE)
    
    # -- Extract PGS trait name
    pgs_filename <- basename(pgs_file_path)
    trait_parts <- unlist(strsplit(unlist(strsplit(pgs_filename, "\\."))[1], "_"))
    trait <- paste0(trait_parts[1], "_", trait_parts[2])
    cat(trait, j, analysis_label, '\n')
    
    # -- Filter PGS for Europeans
    ids_to_keep <- as.character(eur$V2)
    pgs_e <- pgs %>% 
      filter(V2 %in% ids_to_keep)
    
    # -- Standardize PGS
    pgs_e <- pgs_e %>%
      mutate(SCORE1_SUM = as.numeric(V6),
             std_pgs = scale(V6)[,1])
    
    # -- Join PGS file with ad treatment group data
    pharma_link <- inner_join(link, pharma_full, by = "ParticipantID")
    pgs_atc <- inner_join(pharma_link, pgs_e, by = c("IID" = "V2"))
    
    #-- Join with PCs
    pgs_atc <- left_join(pgs_atc, pcs, by = "IID")
    
    #======= Analysis 1: Total Prescription Days
    total_lm <- lm(TotalPrescriptionDays ~ AGE + SEX + PC1 + PC2 + PC3 + std_pgs, data = pgs_atc) 
    lm_summary <- tidy(total_lm) %>% rename(Term = term)
    
    lm_dat <- tibble(
      Dependent = "TotalPrescriptions",
      PGS = trait,
      Analysis = analysis_label
    ) %>%
      bind_cols(lm_summary)
    
    total_lm_list[[length(total_lm_list) + 1]] <- lm_dat
    
    #====== Analysis 2: Number of unique Antidepressants
    numAD_lm <- lm(num_ATC ~ AGE + SEX + PC1 + PC2 + PC3 + std_pgs, data = pgs_atc) 
    lm_summary <- tidy(numAD_lm) %>% rename(Term = term)
    
    lm_dat <- tibble(
      Dependent = "NumADs",
      PGS = trait,
      Analysis = analysis_label
    ) %>%
      bind_cols(lm_summary)
    
    numAD_lm_list[[length(numAD_lm_list) + 1]] <- lm_dat
    
    #====== Analysis 3: Number of unique Antidepressant Classes
    numClass_lm <- lm(num_class ~ AGE + SEX + PC1 + PC2 + PC3 + std_pgs, data = pgs_atc) 
    lm_summary <- tidy(numClass_lm) %>% rename(Term = term)
    
    lm_dat <- tibble(
      Dependent = "NumClasses",
      PGS = trait,
      Analysis = analysis_label
    ) %>%
      bind_cols(lm_summary)
    
    numClass_lm_list[[length(numClass_lm_list) + 1]] <- lm_dat
  }
  
  # -- Combine results
  results <- bind_rows(total_lm_list, numAD_lm_list, numClass_lm_list) %>%
    rename(P.value = p.value, t.value = statistic)
  
  return(results)
}

#=========== RUN BOTH ANALYSES =================================================================

# Run analysis with BIP included (full sample)
cat("Running analysis with BIP included...\n")
results_bip_included <- run_pgs_analysis(exclude_bip = FALSE)

# Run analysis with BIP excluded
cat("Running analysis with BIP excluded...\n")
results_bip_excluded <- run_pgs_analysis(exclude_bip = TRUE)

# Combine both results
results_combined <- bind_rows(results_bip_included, results_bip_excluded)

#=========== PROCESS COMBINED RESULTS =================================================================

# Calculate adjusted p-values only for std_pgs terms
pgs_terms <- results_combined %>%
  filter(Term == "std_pgs") %>%
  group_by(Dependent, Analysis) %>%
  mutate(
    FDR_P = p.adjust(P.value, method = "fdr"),
    Bonf_P = p.adjust(P.value, method = "bonferroni"),
    Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
    Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
  ) %>%
  ungroup() %>%
  select(Dependent, PGS, Term, Analysis, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)

# Join with full data
results_annotated <- results_combined %>%
  left_join(pgs_terms, by = c("Dependent", "PGS", "Term", "Analysis"))

# -- Remapping of variables
pgs_rename_mapping <- c(
  "ADHD_01" = "ADHD",
  "ANO_LOO" = "Anorexia",
  "BIP_LOO" = "Bipolar Disorder",
  "BMI_LOO" = "BMI",
  "CNT_03" = "Chronotype",
  "CRP_01" = "C-reactive protein",
  "MDD_LOO" = "Depression",
  "Migraine_01" = "Migraine",
  "Neuroticism_01" = "Neuroticism",
  "PUD_01" = "Peptic Ulcer Disease",
  "SBP_01" = "Systolic Blood Pressure",
  "SCZ_03" = "Schizophrenia",
  "T2D_03" = "T2D",
  "UKB_35BM" = "LDL-c", 
  "LRA_01" = "LRA"
)

dependent_mapping <- c(
  "NumADs" = "Medication Diversity",
  "NumClasses" = "Class Diversity",
  "TotalPrescriptions" = "Cumulative Prescription Dispense (days)"
)

results_final <- results_annotated %>%
  mutate(
    PGS = recode(PGS, !!!pgs_rename_mapping),
    Dependent = recode(Dependent, !!!dependent_mapping)
  ) %>%
  mutate(
    across(c(P.value, FDR_P, Bonf_P), ~ format(signif(.x, 2), scientific = TRUE)),
    std.error = if_else(Dependent %in% c("Medication Diversity", "Class Diversity"),
                        round(std.error, 4), round(std.error, 2)),
    estimate = if_else(Dependent %in% c("Medication Diversity", "Class Diversity"),
                       round(estimate, 4), round(estimate, 1)),
    t.value = round(t.value, 2)
  )

# Reorder columns for clarity
results_final <- results_final %>%
  select(Analysis,Dependent, PGS, Term, estimate, std.error, t.value, P.value, 
         FDR_P, Bonf_P, Sig_FDR, Sig_Bonf, everything()) %>%
  group_by(Analysis, Dependent, PGS) %>%
  mutate(
    Analysis = if_else(Analysis != lag(Analysis, default = ""), Analysis, NA_character_),
    Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
    PGS = if_else(PGS != lag(PGS, default = ""), PGS, NA_character_)
  )

# Save to Excel
wb <- loadWorkbook(file.path("/scratch/user/uqawal15", "All_Results.xlsx"))
removeWorksheet(wb, "Table4")
addWorksheet(wb, "Table4")
writeData(wb, "Table4", results_final)
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)


