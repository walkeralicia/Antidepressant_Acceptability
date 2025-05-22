# -- Load packages
library(tidyverse)
library(emmeans)
library(openxlsx)
library(broom)

# -- Set working directory
wkdir <- "/QRISdata/Q7280/pharmacogenomics"
output_dir <- "/scratch/user/uqawal15"

# -- Load treatment group data and extract AGE and SNPSEX information
dat <- read_csv(file.path(wkdir, "data/AGDSAcceptabilityTreatmentGroups_14122024.csv"))
pheno <- read_csv(file.path(wkdir, "phenotypes/survey_phenotypes.csv")) %>%
  select(STUDYID, AGE, SEX, BMI)

# -- Drug mapping as tibble
drug_mapping <- tibble(
  ATCCode = c('N06AB06', 'N06AB10', 'N06AB04', 'N06AX16', 'N06AX21', 
              'N06AA09', 'N06AB05', 'N06AX11', 'N06AX23', 'N06AB03'),
  DrugName = c("SSRI:Sertraline", "SSRI:Escitalopram", "SSRI:Citalopram", "SNRI:Venlafaxine", "SNRI:Duloxetine", 
               "TCA:Amitriptyline", "SSRI:Paroxetine", "TeCA:Mirtazapine", "SNRI:Desvenlafaxine", "SSRI:Fluoxetine"),
  DrugClass = c('SSRI', 'SSRI', 'SSRI', 'SNRI', 'SNRI', 'TCA', 'SSRI', 'TeCA', 'SNRI', 'SSRI')
)

# -- Map ATC codes
ad_mapped <- dat %>%
  left_join(drug_mapping, by = "ATCCode")

# -- Create participant summary
pharma <- ad_mapped %>%
  group_by(ParticipantID) %>%
  summarize(
    num_ATC = n_distinct(ATCCode),
    num_class = n_distinct(DrugClass),
    Total_TEs = sum(NumberOfTreatmentPeriods, na.rm = TRUE),
    TotalPrescriptionDays = sum(PrescriptionDays, na.rm = TRUE)
  )

# Join data
pheno <- pheno %>% 
  mutate(ParticipantID = STUDYID)

pharma_full <- ad_mapped %>%
  left_join(pheno, by = "ParticipantID") %>%
  left_join(pharma, by = "ParticipantID")

# -- Convert SEX
pharma_full <- pharma_full %>%
  mutate(SEX = case_when(
    SEX == "Female" ~ 0,
    SEX == "Male" ~ 1,
    TRUE ~ NA_real_
  ))

# -- Genetic data
pgslist <- system('find /QRISdata/Q7280/pharmacogenomics/pgs -name "*agds_sbrc_gctb_plink2.sscore" -type f -exec ls {} \\;', intern = TRUE)
pgs_codes <- c("PUD_01", "T2D_03", "CNT_03", "ADHD_01", "BIP_LOO", "BMI_LOO", "MDD_LOO", "MDD_07_hsa04726", "SCZ_02", "Migraine_01", "SBP_01", "ANX_LOO", "ANO_LOO", "UKB_35BM_2021_LDL_direct_adjstatins", "CRP_01", "OCD_2024", "Neuroticism_01", "LRA_01")
num <- length(pgs_codes)
selected_pgs <- pgslist[grep(paste(pgs_codes, collapse = "|"), pgslist)]

# -- List of Europeans
eur <- read_table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id", col_names = FALSE)

# -- List of all antidepressant medications
drugs <- unique(pharma_full$DrugName)
classes <- c("SSRI", "SNRI", "TCA", "TeCA", "Lithium")

# -- Linkage file for those with genetic data
link <- read_table(file.path("/QRISdata/Q7280/pharmacogenomics", "data/StudyID_GenoID_link.txt")) %>%
  rename("ParticipantID" = "StudyCodeID") %>%
  distinct()

#=========== LM MODELS =================================================================

# -- Set up results list
total_lm_list <- list()
numAD_lm_list <- list()
numClass_lm_list <- list()

# -- For each PGS trait fit a linear model with treatment group as the independent variable
for (j in 1:length(selected_pgs)) {
  
  if (j != 16) {
    # -- Read in PGS file for the selected trait
    pgs_file_path <- selected_pgs[j]
    pgs <- read_table(selected_pgs[j], col_names = FALSE)
    
    # -- Extract PGS trait name
    pgs_filename <- basename(pgs_file_path)
    trait_parts <- unlist(strsplit(unlist(strsplit(pgs_filename, "\\."))[1], "_"))
    trait <- paste0(trait_parts[1], "_", trait_parts[2])
    cat(trait, j, '\n')
    
    # -- Filter PGS for Europeans (loss of 806 individuals)
    ids_to_keep <- as.character(eur$X2)
    pgs_e <- pgs %>% 
      filter(X2 %in% ids_to_keep)
    
    # -- Standardize PGS
    pgs_e <- pgs_e %>%
      mutate(X6 = as.numeric(X6),
             std_pgs = scale(X6)[,1])
    
    # -- Join PGS file with ad treatment group data
    pharma_link <- inner_join(link, pharma_full, by = "ParticipantID")
    pgs_atc <- inner_join(pharma_link, pgs_e, by = c("IID" = "X2"))
  } else {
    
    # -- Read in PGS file for the selected trait
    pgs_file_path <- selected_pgs[j]
    pgs <- read_table(selected_pgs[j])
    
    # -- Extract PGS trait name
    pgs_filename <- basename(pgs_file_path)
    trait_parts <- unlist(strsplit(unlist(strsplit(pgs_filename, "\\."))[1], "_"))
    trait <- paste0(trait_parts[1], "_", trait_parts[2], "_", trait_parts[3])
    cat(trait, j, '\n')
    
    # -- Filter PGS for Europeans (loss of 806 individuals)
    ids_to_keep <- as.character(eur$X2)
    pgs_e <- pgs %>% 
      filter(IID %in% ids_to_keep)
    
    # -- Standardize PGS
    pgs_e <- pgs_e %>%
      mutate(std_pgs = scale(SCORESUM)[,1])
    
    # -- Join PGS file with ad treatment group data
    pharma_link <- inner_join(link, pharma_full, by = "ParticipantID")
    pgs_atc <- inner_join(pharma_link, pgs_e, by = "IID")
  }
  
  #======= Analysis 1: Total Prescription Days
  
  # -- LM: Total Prescription Duration ~ AGE + SEX + PGS
  total_lm <- lm(TotalPrescriptionDays ~ AGE + SEX + std_pgs, data = pgs_atc) 
  lm_summary <- tidy(total_lm) %>% 
    rename(Term = term)
  
  # -- Create data frame
  lm_dat <- tibble(
    Dependent = "TotalPrescriptions",
    PGS = trait
  ) %>%
    bind_cols(lm_summary)
  
  # -- Add to list
  total_lm_list[[length(total_lm_list) + 1]] <- lm_dat
  
  #====== Analysis 2: Number of unique Antidepressants
  
  # --LM: Number of unique antidepressants ~ AGE + SEX + PGS
  numAD_lm <- lm(num_ATC ~ AGE + SEX + std_pgs, data = pgs_atc) 
  lm_summary <- tidy(numAD_lm) %>% 
    rename(Term = term)
  
  # -- Create data frame
  lm_dat <- tibble(
    Dependent = "NumADs",
    PGS = trait
  ) %>%
    bind_cols(lm_summary)
  
  # -- Add to list
  numAD_lm_list[[length(numAD_lm_list) + 1]] <- lm_dat
  
  #====== Analysis 3: Number of unique Antidepressant Classes
  
  # --LM: Number of unique antidepressants ~ AGE + SEX +PGS
  numClass_lm <- lm(num_class ~ AGE + SEX + std_pgs, data = pgs_atc) 
  lm_summary <- tidy(numClass_lm) %>% 
    rename(Term = term)
  
  # -- Create data frame
  lm_dat <- tibble(
    Dependent = "NumClasses",
    PGS = trait
  ) %>%
    bind_cols(lm_summary)
  
  # -- Add to list
  numClass_lm_list[[length(numClass_lm_list) + 1]] <- lm_dat
}

# -- Combine results from every model
results_lm <- bind_rows(total_lm_list, numAD_lm_list, numClass_lm_list) %>%
  rename(P.value = p.value, 
         t.value = statistic)

# Step 1: Calculate adjusted p-values only for std_pgs terms
pgs_terms <- results_lm %>%
  filter(Term == "std_pgs") %>%
  group_by(Dependent) %>%
  mutate(
    FDR_P = p.adjust(P.value, method = "fdr"),
    Bonf_P = p.adjust(P.value, method = "bonferroni"),
    Sig_FDR = if_else(FDR_P < 0.05, "*", ""),
    Sig_Bonf = if_else(Bonf_P < 0.05, "*", "")
  ) %>%
  ungroup() %>%
  select(Dependent, PGS, Term, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf)

# Step 2: Join with full data — but include Term in the join key
results_lm_annotated <- results_lm %>%
  left_join(pgs_terms, by = c("Dependent", "PGS", "Term"))

# -- Remapping of variables
pgs_rename_mapping <- c(
  "ADHD_01" = "ADHD",
  "ANO_LOO" = "Anorexia",
  "ANX_LOO" = "Anxiety",
  "BIP_LOO" = "Bipolar Disorder",
  "BMI_LOO" = "BMI",
  "CNT_03" = "Chronotype",
  "CRP_01" = "C-reactive protein",
  "MDD_07_hsa04726" = "Depression KEGG:hsa04726",
  "MDD_LOO" = "Depression",
  "Migraine_01" = "Migraine",
  "Neuroticism_01" = "Neuroticism",
  "OCD_2024" = "OCD",
  "PUD_01" = "Peptic Ulcer Disease",
  "SBP_01" = "Systolic Blood Pressure",
  "SCZ_02" = "Schizophrenia",
  "T2D_03" = "T2D",
  "UKB_35BM" = "LDL-c", 
  "LRA_01" = "LRA"
)

results_lm_renamed <- results_lm_annotated %>%
  mutate(PGS = recode(PGS, !!!pgs_rename_mapping))

dependent_mapping <- c(
  "NumADs" = "Medication Diversity",
  "NumClasses" = "Class Diversity",
  "TotalPrescriptions" = "Cumulative Prescription Dispense (days)"
)

results_lm_renamed <- results_lm_renamed %>%
  mutate(Dependent = recode(Dependent, !!!dependent_mapping))

# Format for Excel display
results_lm_renamed <- results_lm_renamed %>%
  group_by(Dependent, PGS) %>%
  mutate(
    Dependent = if_else(Dependent != lag(Dependent, default = ""), Dependent, NA_character_),
    PGS = if_else(PGS != lag(PGS, default = ""), PGS, NA_character_)
  )

results_lm_renamed <- results_lm_renamed %>%
  mutate(
    across(c(P.value, FDR_P, Bonf_P), ~ signif(.x, 3)),
    across(c(std.error, estimate, t.value), ~ round(.x, 2))
  )

# Add to existing workbook
wb <- loadWorkbook(file.path("/scratch/user/uqawal15", "All_Results.xlsx"))
addWorksheet(wb, "Table3")
writeData(wb, "Table3", results_lm_renamed)
saveWorkbook(wb, file.path(output_dir, "All_Results.xlsx"), overwrite = TRUE)