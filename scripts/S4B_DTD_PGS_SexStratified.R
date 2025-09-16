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

#=========== FUNCTION TO RUN SEX-STRATIFIED PGS ANALYSIS =================================================================

run_sex_stratified_pgs_analysis <- function(sex_value, sex_label) {
  cat(paste("Running PGS analysis for", sex_label, "\n"))
  
  # Start with the base dataset
  pharma_full <- pharma %>%
    left_join(pheno, by = "ParticipantID") %>%
    filter(MDD==1) %>%
    select(-DXBPD2)  # Remove BIP variable entirely
  
  # -- Convert SEX and filter by sex
  pharma_full <- pharma_full %>%
    mutate(SEX = case_when(
      SEX == "Female" ~ 0,
      SEX == "Male" ~ 1,
      TRUE ~ NA_real_
    )) %>%
    filter(SEX == sex_value)
  
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
    cat(trait, j, sex_label, '\n')
    
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
    
    # Check if we have sufficient data
    if (nrow(pgs_atc) < 10) {
      cat(sprintf("Insufficient data for %s in %s: only %d observations\n", trait, sex_label, nrow(pgs_atc)))
      next
    }
    
    #======= Analysis 1: Total Prescription Days (without SEX since we're stratified)
    total_lm <- lm(TotalPrescriptionDays ~ AGE + PC1 + PC2 + PC3 + std_pgs, data = pgs_atc) 
    lm_summary <- tidy(total_lm) %>% rename(Term = term)
    
    lm_dat <- tibble(
      Dependent = "TotalPrescriptions",
      PGS = trait,
      Sex_Stratum = sex_label
    ) %>%
      bind_cols(lm_summary)
    
    total_lm_list[[length(total_lm_list) + 1]] <- lm_dat
    
    #====== Analysis 2: Number of unique Antidepressants
    numAD_lm <- lm(num_ATC ~ AGE + PC1 + PC2 + PC3 + std_pgs, data = pgs_atc) 
    lm_summary <- tidy(numAD_lm) %>% rename(Term = term)
    
    lm_dat <- tibble(
      Dependent = "NumADs",
      PGS = trait,
      Sex_Stratum = sex_label
    ) %>%
      bind_cols(lm_summary)
    
    numAD_lm_list[[length(numAD_lm_list) + 1]] <- lm_dat
    
    #====== Analysis 3: Number of unique Antidepressant Classes
    numClass_lm <- lm(num_class ~ AGE + PC1 + PC2 + PC3 + std_pgs, data = pgs_atc) 
    lm_summary <- tidy(numClass_lm) %>% rename(Term = term)
    
    lm_dat <- tibble(
      Dependent = "NumClasses",
      PGS = trait,
      Sex_Stratum = sex_label
    ) %>%
      bind_cols(lm_summary)
    
    numClass_lm_list[[length(numClass_lm_list) + 1]] <- lm_dat
  }
  
  # -- Combine results for this sex
  if (length(total_lm_list) > 0 || length(numAD_lm_list) > 0 || length(numClass_lm_list) > 0) {
    results <- bind_rows(total_lm_list, numAD_lm_list, numClass_lm_list) %>%
      rename(P.value = p.value, t.value = statistic)
  } else {
    results <- tibble()
  }
  
  return(results)
}

#=========== RUN SEX-STRATIFIED ANALYSES =================================================================

cat("Running sex-stratified PGS analyses...\n")

# Run analysis for females (SEX = 0)
results_female <- run_sex_stratified_pgs_analysis(0, "Female")

# Run analysis for males (SEX = 1)  
results_male <- run_sex_stratified_pgs_analysis(1, "Male")

# Combine both results
results_combined <- bind_rows(results_female, results_male)

#=========== PROCESS COMBINED RESULTS =================================================================

# Calculate adjusted p-values only for std_pgs terms, stratified by sex
pgs_terms <- results_combined %>%
  filter(Term == "std_pgs") %>%
  group_by(Dependent, Sex_Stratum) %>%
  mutate(
    FDR_P = p.adjust(P.value, method = "fdr"),
    Bonf_P = p.adjust(P.value, method = "bonferroni"),
    Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
    Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
  ) %>%
  ungroup() %>%
  select(Dependent, PGS, Term, Sex_Stratum, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)

# Join with full data
results_annotated <- results_combined %>%
  left_join(pgs_terms, by = c("Dependent", "PGS", "Term", "Sex_Stratum"))

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
    across(c(P.value, FDR_P, Bonf_P), ~ ifelse(is.na(.x), NA_character_, format(signif(.x, 2), scientific = TRUE))),
    std.error = if_else(Dependent %in% c("Medication Diversity", "Class Diversity"),
                        round(std.error, 4), round(std.error, 2)),
    estimate = if_else(Dependent %in% c("Medication Diversity", "Class Diversity"),
                       round(estimate, 4), round(estimate, 1)),
    t.value = round(t.value, 2)
  )

# Reorder columns for clarity
results_final <- results_final %>%
  select(Sex_Stratum, Dependent, PGS, Term, estimate, std.error, t.value, P.value, 
         FDR_P, Bonf_P, Sig_FDR, Sig_Bonf, everything()) %>%
  arrange(Sex_Stratum, Dependent, PGS) %>%
  group_by(Sex_Stratum, Dependent, PGS) %>%
  mutate(
    Sex_Stratum = if_else(Sex_Stratum != lag(Sex_Stratum, default = ""), Sex_Stratum, NA_character_),
    Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
    PGS = if_else(PGS != lag(PGS, default = ""), PGS, NA_character_)
  ) %>%
  ungroup()

# Save to Excel with separate sheets for each sex
wb <- createWorkbook()

# Add new sheets for sex-stratified results
addWorksheet(wb, "Table4_Female_PGS")
writeData(wb, "Table4_Female_PGS", results_final %>% filter(!is.na(Sex_Stratum) & Sex_Stratum == "Female" | is.na(Sex_Stratum)))

addWorksheet(wb, "Table4_Male_PGS")
writeData(wb, "Table4_Male_PGS", results_final %>% filter(!is.na(Sex_Stratum) & Sex_Stratum == "Male" | is.na(Sex_Stratum)))

addWorksheet(wb, "Table4_Combined_PGS")
writeData(wb, "Table4_Combined_PGS", results_final)

# Auto-adjust column widths
for(sheet in c("Table4_Female_PGS", "Table4_Male_PGS", "Table4_Combined_PGS")) {
  setColWidths(wb, sheet, cols = 1:15, widths = "auto")
}

saveWorkbook(wb, file.path(output_dir, "All_Results_PGS_SexStratified.xlsx"), overwrite = TRUE)

cat("Sex-stratified PGS analysis complete! Results saved to:\n")
cat("- Excel file:", file.path(output_dir, "All_Results.xlsx"), "\n")
cat("- Sheets: Table4_Female_PGS, Table4_Male_PGS, Table4_Combined_PGS\n")