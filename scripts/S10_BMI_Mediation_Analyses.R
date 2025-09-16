# -- Load R packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
  library(broom)
  library(purrr)
})

#BMI PGS ──→ BMI ←── Treatment Group
#\           /
#  └─────────┘
#(Collider)


# If both BMI PGS and treatment group influence BMI, then:
# controlling for BMI PGS creates an artificial pathway
# treatment group appears more strongly associated with BMI
# conditioned on a collider (BMI PGS), creating spurious associations


# Load in phenotype data
pheno <- read.csv("/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes/inner_data_360days.csv")
pheno_processed <- pheno %>%
  transmute(
    IID,
    ParticipantID,
    SEX = case_when(SEX == "Female" ~ 0, SEX == "Male" ~ 1, TRUE ~ NA_real_),
    AGE, DrugName, DrugClass, BMI
  )

# Load in genetic data
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

# Combine data
combined_data <- pheno_processed %>%
  inner_join(bmi_pgs_processed, by = c("IID" = "V2")) %>%
  filter(DrugClass %in% c("SSRI", "SNRI")) %>%
  mutate(
    treatment_group = case_when(
      DrugClass == "SSRI" ~ 0,
      DrugClass == "SNRI" ~ 1,
      TRUE ~ NA
    )
  )
analysis_vars <- c("treatment_group", "std_pgs", "BMI", "AGE", "SEX", "PC1", "PC2", "PC3")
analysis_data <- combined_data[complete.cases(combined_data[analysis_vars]), ]

# Mediation analysis
# BMI PGS → Treatment Group → BMI


# Using mediation package in R
library(mediation)

# Step 1: Test BMI PGS → SNRI vs SSRI
mediator_model <- lm(treatment_group ~ std_pgs + AGE + SEX + PC1 + PC2 + PC3, 
                      data = analysis_data)

# Step 2: BMI outcome model
outcome_model <- lm(BMI ~ treatment_group + std_pgs + AGE + SEX + PC1 + PC2 + PC3, 
                    data = analysis_data)

# Step 3: Mediation analysis
mediation_result <- mediate(mediator_model, outcome_model, 
                            treat = "std_pgs", 
                            mediator = "treatment_group")

summary(mediation_result)

Causal Mediation Analysis

Quasi-Bayesian Confidence Intervals

Estimate 95% CI Lower 95% CI Upper p-value
ACME           0.00659479   0.00030880   0.01772769   0.024 *
  ADE            2.82595852   2.65146110   2.99586168  <2e-16 ***
  Total Effect   2.83255332   2.65879319   3.00441695  <2e-16 ***
  Prop. Mediated 0.00209130   0.00010666   0.00627004   0.024 *
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Sample Size Used: 4920


Simulations: 1000
