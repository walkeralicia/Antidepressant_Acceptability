#---- Create a covariate and phenotype file
R
library(dplyr)
library(openxlsx)

# Define file paths
PC_FILE <- "/QRISdata/Q5338/Ancestry_analysis/AGDS_R11_TOPMedr2_pruned.05.common_pca3.proj.eigenvec"
DEMO_FILE <- "/QRISdata/Q7280/pharmacogenomics/phenotypes/BMI_age_sex.csv"
TREATMENT_FILE <- "/QRISdata/Q7280/pharmacogenomics/treatment_groups/Final_Treatmentgroups_360days.csv"
PRESCRIPTION_FILE <- "/QRISdata/Q7280/pharmacogenomics/data/AGDSAcceptabilityTreatmentGroups_14122024.csv"
PHENO_FILE <- "/QRISdata/Q7280/pharmacogenomics/phenotypes/survey_phenotypes.csv"
OUTPUT_DIR <- "/scratch/user/uqawal15"

# Load and format principal components
pcs <- read.table(PC_FILE) %>%
  setNames(c("FID", "IID", "PC1", "PC2", "PC3"))

# Load demographic data
pheno <- read.csv(DEMO_FILE) %>%
  select(STUDYID, AGE, SEX)

# Load treatment group data
groups <- read.csv(TREATMENT_FILE) %>%
  select(IID, ParticipantID, DrugName, DrugClass)
  
# Load prescription file
drugs <- read.csv(PRESCRIPTION_FILE) %>%
  mutate(
    EarliestPrescription = as.Date(EarliestPrescription, format = "%d/%m/%Y"),
    LatestPrescription = as.Date(LatestPrescription, format = "%d/%m/%Y")
  )
  
#======== GWAS Outcomes =======================

#------- First Option (SSRI acceptability) --------------------

# Identify participants not in SSRI, BIP-L, or BIP+L groups (i.e., controls, n = 5205)
controls <- groups %>% 
  filter(!(DrugClass %in% c("SSRI") | DrugName %in% c("BIP-L", "BIP+L")))

# Filter prescription data to include only those from controls
drugs_controls <- drugs %>% 
  filter(ParticipantID %in% controls$ParticipantID)

# Identify controls with long-term SSRI use (>360 days, n = 476)
long_term_controls <- drugs_controls %>%
  filter(ATCCode %in% c("N06AB03", "N06AB04", "N06AB05", "N06AB06", "N06AB10") & PrescriptionDays > 360) %>%
  distinct(ParticipantID)

# Remove long-term users from the controls dataset (n = 4727)
controls_cleaned <- controls %>%
  filter(!ParticipantID %in% long_term_controls$ParticipantID)

# Create SSRI responder phenotype
# 1 = SSRI responder, 0 = non-responder, NA = not applicable
groups$SSRI_Acceptability <- ifelse(
  groups$ParticipantID %in% controls_cleaned$ParticipantID, 0,
  ifelse(groups$DrugClass == "SSRI", 1, NA)
)


#------ Second option (SSRI+SNRI acceptability) -------------------------

# Identify participants not in SSRI, SNRI, BIP-L, or BIP+L groups (i.e., controls)
controls <- groups %>% 
  filter(!(DrugClass %in% c("SSRI", "SNRI") | DrugName %in% c("BIP-L", "BIP+L")))

# Filter prescriptions data to include only those from the controls
drugs_controls <- drugs %>% 
  filter(ParticipantID %in% controls$ParticipantID)

# Identify controls with long-term SSRI or SNRI use (>360 days)
long_term_controls <- drugs_controls %>%
  filter(ATCCode %in% c("N06AB03", "N06AB04", "N06AB05", "N06AB06", "N06AB10", "N06AX16", "N06AX21", "N06AX23") & PrescriptionDays > 360) %>%
  distinct(ParticipantID)

# Remove long-term users from the controls dataset
controls_cleaned <- drugs_controls %>%
  filter(!ParticipantID %in% long_term_controls$ParticipantID)

# Create SSRI+SNRI responder phenotype
# 1 = SSRI+SNRI responder, 0 = non-responder, NA = not applicable
groups$SSRI_SNRI_Acceptability <- ifelse(
  groups$ParticipantID %in% controls_cleaned$ParticipantID, 0,
  ifelse(groups$DrugClass %in% c("SSRI", "SNRI"), 1, NA)
)

#--------- Self-report SSRI and SNRI response ---------------
response <- read.csv(PHENO_FILE)
response2 <- response %>%
mutate(SSRI_Responder = case_when(
  is.na(Sertraline) & is.na(Escitalopram) & is.na(Citalopram) & is.na(Fluoxetine) & is.na(Paroxetine) & is.na(Duloxetine) & is.na(Venlafaxine) & is.na(Desvenlafaxine) & is.na(Mirtazapine) & is.na(Amitriptyline) ~ NA_real_,
  Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 ~ 1,
  TRUE ~ 0
),
SSRI_SNRI_Responder = case_when(
  is.na(Sertraline) & is.na(Escitalopram) & is.na(Citalopram) & is.na(Fluoxetine) & is.na(Paroxetine) & is.na(Duloxetine) & is.na(Venlafaxine) & is.na(Desvenlafaxine) & is.na(Mirtazapine) & is.na(Amitriptyline) ~ NA_real_,
  Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Duloxetine == 1 | Venlafaxine == 1 | Desvenlafaxine == 1 ~ 1,
  TRUE ~ 0
)
) %>%
select(STUDYID, SSRI_Responder, SSRI_SNRI_Responder) %>%
distinct()


#---------- # Combine all datasets

#-- All data for those with a PC
dat <- pheno %>%
  full_join(groups, by = c("STUDYID" = "ParticipantID")) %>%
  full_join(pcs, by = c("IID" = "IID")) %>%
  full_join(response2, by = c("STUDYID")) %>%
  filter(!is.na(IID)) %>%
  mutate(
    AGE = as.numeric(AGE),
    SEX = ifelse(is.na(SEX), "NONE", SEX)
  )
  
#-- Filter for those of EUR ancestry
eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id")
dat_eur <- dat %>%
  filter(IID %in% eur$V2)

# Create and save covariate file
dat_eur %>%
  select(FID, IID, AGE, SEX, PC1, PC2, PC3) %>%
  write.table(
    file.path(OUTPUT_DIR, "AD_acceptability_covariates.txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )

# Create and save outcome file
dat_eur %>%
  select(FID, IID, DrugClass, SSRI_Acceptability, SSRI_SNRI_Acceptability, SSRI_Responder, SSRI_SNRI_Responder) %>%
  write.table(
    file.path(OUTPUT_DIR, "AD_acceptability_outcome.txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )


#-- Check case and control number for each phenotye
dat_eur %>%
select(STUDYID, SSRI_Acceptability, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_Acceptability)
  SSRI_Acceptability    n
1                  0 4722
2                  1 3566
3                 NA 1541

dat_eur %>%
select(STUDYID, SSRI_SNRI_Acceptability, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_SNRI_Acceptability)
  SSRI_SNRI_Acceptability    n
1                       0 2348
2                       1 5774
3                      NA 1707

dat_eur %>%
select(STUDYID, SSRI_Responder, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_Responder)

  SSRI_Responder    n
1              0 2744
2              1 6739
3             NA  346

dat_eur %>%
select(STUDYID, SSRI_SNRI_Responder, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_SNRI_Responder)

  SSRI_SNRI_Responder    n
1                   0  935
2                   1 8548
3                  NA  346

#--------------- Run SSRI acceptability GWAS using PLINK2 (common variants)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=SSRI_Acceptability_gwas
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/SSRI_Acceptability_gwas_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/SSRI_Acceptability_gwas_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

module load plink/2.00a3.6-gcc-11.3.0
wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/imputed_chr${c} \
--pheno ${wkdir}/AD_acceptability_outcome.txt \
--pheno-name SSRI_Acceptability \
--1 \
--covar ${wkdir}/AD_acceptability_covariates.txt \
--covar-variance-standardize \
--maf 0.01 \
--glm hide-covar cols=+a1freq \
--threads 8 \
--out SSRI_Acceptability_gwas_chr${c}_results

#---------------------------------------------

#--------------- Run SSRI+SNRI acceptability GWAS using PLINK2 (common variants)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=SSRI_SNRI_Acceptability_gwas
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/SSRI_SNRI_Acceptability_gwas_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/SSRI_SNRI_Acceptability_gwas_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

module load plink/2.00a3.6-gcc-11.3.0
wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/imputed_chr${c} \
--pheno ${wkdir}/AD_acceptability_outcome.txt \
--pheno-name SSRI_SNRI_Acceptability \
--1 \
--covar ${wkdir}/AD_acceptability_covariates.txt \
--covar-variance-standardize \
--maf 0.01 \
--glm hide-covar cols=+a1freq \
--threads 8 \
--out SSRI_SNRI_Acceptability_gwas_chr${c}_results

#--------------- Run SSRI self-report efficacy GWAS using PLINK2 (common variants)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=SSRI_Responder_gwas
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/SSRI_Responder_gwas_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/SSRI_Responder_gwas_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

module load plink/2.00a3.6-gcc-11.3.0
wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/imputed_chr${c} \
--pheno ${wkdir}/AD_acceptability_outcome.txt \
--pheno-name SSRI_Responder \
--1 \
--covar ${wkdir}/AD_acceptability_covariates.txt \
--covar-variance-standardize \
--maf 0.01 \
--glm hide-covar cols=+a1freq \
--threads 8 \
--out SSRI_Responder_gwas_chr${c}_results

#--------------- Run SSRI+SNRI self-report efficacy GWAS using PLINK2 (common variants)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=SSRI_SNRI_Responder_gwas
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/SSRI_SNRI_Responder_gwas_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/SSRI_SNRI_Responder_gwas_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

module load plink/2.00a3.6-gcc-11.3.0
wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/imputed_chr${c} \
--pheno ${wkdir}/AD_acceptability_outcome.txt \
--pheno-name SSRI_SNRI_Responder \
--1 \
--covar ${wkdir}/AD_acceptability_covariates.txt \
--covar-variance-standardize \
--maf 0.01 \
--glm hide-covar cols=+a1freq \
--threads 8 \
--out SSRI_SNRI_Responder_gwas_chr${c}_results


#----------- Concatenate GWAS results --------------------------------------

head -n1 SSRI_Acceptability_gwas_chr1_results.SSRI_Acceptability.glm.logistic.hybrid > SSRI_Acceptability_GWAS.txt

for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    awk 'NR>1 && $NF!="CONST_OMITTED_ALLELE"' SSRI_Acceptability_gwas_chr${chr}_results.SSRI_Acceptability.glm.logistic.hybrid >> SSRI_Acceptability_GWAS.txt
done
echo "Done! Results saved."

#----
head -n1 SSRI_SNRI_Acceptability_gwas_chr1_results.SSRI_SNRI_Acceptability.glm.logistic.hybrid > SSRI_SNRI_Acceptability_GWAS.txt

for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    awk 'NR>1 && $NF!="CONST_OMITTED_ALLELE"' SSRI_SNRI_Acceptability_gwas_chr${chr}_results.SSRI_SNRI_Acceptability.glm.logistic.hybrid >> SSRI_SNRI_Acceptability_GWAS.txt
done
echo "Done! Results saved."

#---
head -n1 SSRI_Responder_gwas_chr1_results.SSRI_Responder.glm.logistic.hybrid > SSRI_Responder_GWAS.txt

for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    awk 'NR>1 && $NF!="CONST_OMITTED_ALLELE"' SSRI_Responder_gwas_chr${chr}_results.SSRI_Responder.glm.logistic.hybrid >> SSRI_Responder_GWAS.txt
done
echo "Done! Results saved."

#---
head -n1 SSRI_SNRI_Responder_gwas_chr1_results.SSRI_SNRI_Responder.glm.logistic.hybrid > SSRI_SNRI_Responder_GWAS.txt

for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    awk 'NR>1 && $NF!="CONST_OMITTED_ALLELE"' SSRI_SNRI_Responder_gwas_chr${chr}_results.SSRI_SNRI_Responder.glm.logistic.hybrid >> SSRI_SNRI_Responder_GWAS.txt
done
echo "Done! Results saved."

#---------------------------- SSRI Acceptability ------------------------------------

# load libraries
library(data.table)
library(qqman)
library(dplyr)
library(openxlsx)

# genomic inflation factor (lambda)
calculate_lambda <- function(pvals) {
  chisq_stats <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq_stats) / 0.456
  return(lambda)
}

# Path to GWAS results
input_file <- "/scratch/user/uqawal15/SSRI_Acceptability_GWAS.txt"
output_dir <- "/scratch/user/uqawal15"
gwas <- fread(input_file, fill = TRUE)
              
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)

# Make the Manhattan plot
png(file.path(output_dir, "SSRI_Acceptability_GWAS_manhattan.png"), 
    bg = "white", 
    type = "cairo", 
    width = 2000,
    height = 1200,
    res = 300)

manhattan(gwas_cleaned, 
          chr="#CHROM", 
          bp="POS", 
          snp="ID", 
          p="P",
          cex = 0.6,
          cex.axis = 0.8,
          suggestiveline = -log10(5e-6),
          genomewideline = -log10(5e-8),
          main = "SSRI Acceptability GWAS")
dev.off()

# Make the qqplot
png(file.path(output_dir, "SSRI_Acceptability_GWAS_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI Acceptability")
dev.off()


#-- Convert to COJO format
cojo_temp <- gwas_cleaned %>%
  mutate(
    BETA = Z_STAT * `LOG(OR)_SE`,
    A2 = case_when(
      A1 == REF ~ ALT, 
      A1 == ALT ~ REF,
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    CHR = `#CHROM`,
    POS,
    SNP = ID, 
    REF,
    ALT,
    A1 = A1,      
    A2 = A2,      
    A1_FREQ,       
    b = BETA,       
    OR,
    SE = `LOG(OR)_SE`,  
    Z_STAT,
    P,                
    N = OBS_CT        
  )
  
#-- Save file for LD clumping
ld_format <- cojo_temp %>%
select(SNP, CHR, BP = POS, P)
fwrite(ld_format, file.path(output_dir, "SSRI_Acceptability_GWAS.temp"), sep = "\t")

#-- Save COJO format as tab-delimited
cojo_format <- cojo_temp %>%
select(-CHR, -POS,-REF, -ALT, -OR, -Z_STAT)
fwrite(cojo_format, file.path(output_dir, "SSRI_Acceptability_GWAS.ma"), sep = "\t")

#-- Format results for a Supp Table
sig <- cojo_temp %>%
    mutate_at(vars(SE, P, A1_FREQ), ~ signif(., 2)) %>%
    mutate_at(vars(OR, Z_STAT), ~ round(., 2)) %>%
    select(CHR, POS, SNP, REF, ALT, A1, A2, A1_FREQ, OR, SE, Z_STAT, P, N) %>%
    filter(P < 0.000005) %>%
    arrange(P)
    
wb <- createWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
#removeWorksheet(wb, "Table17")
addWorksheet(wb, "Table17")
writeData(wb, "Table17", sig)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)

# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 1.006742

#---------------------------- SSRI+SNRI Acceptability ------------------------------------

# load libraries
library(data.table)
library(qqman)
library(dplyr)
library(openxlsx)

# genomic inflation factor (lambda)
calculate_lambda <- function(pvals) {
  chisq_stats <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq_stats) / 0.456
  return(lambda)
}

# Path to GWAS results
input_file <- "/scratch/user/uqawal15/SSRI_SNRI_Acceptability_GWAS.txt"
output_dir <- "/scratch/user/uqawal15"
gwas <- fread(input_file, fill = TRUE)
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)
              
# Make the Manhattan plot
png(file.path(output_dir, "SSRI_SNRI_Acceptability_GWAS_manhattan.png"), 
    bg = "white", 
    type = "cairo", 
    width = 2000,
    height = 1200,
    res = 300)

manhattan(gwas_cleaned, 
          chr="#CHROM", 
          bp="POS", 
          snp="ID", 
          p="P",
          cex = 0.6,
          cex.axis = 0.8,
          suggestiveline = -log10(5e-6),
          genomewideline = -log10(5e-8),
          main = "SSRI+SNRI Acceptability GWAS")
dev.off()

# Make the qqplot
png(file.path(output_dir, "SSRI_SNRI_Acceptability_GWAS_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI+SNRI Acceptability")
dev.off()

#-- Convert to COJO format
cojo_temp <- gwas_cleaned %>%
  mutate(
    BETA = Z_STAT * `LOG(OR)_SE`,
    A2 = case_when(
      A1 == REF ~ ALT, 
      A1 == ALT ~ REF, 
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    CHR = `#CHROM`,
    POS,
    SNP = ID,      
    REF,
    ALT,
    A1 = A1,   
    A2 = A2,   
    A1_FREQ,    
    b = BETA,       
    OR,
    SE = `LOG(OR)_SE`,   
    Z_STAT,
    P,                 
    N = OBS_CT          
  )
  
#-- Save file for LD clumping
ld_format <- cojo_temp %>%
select(SNP, CHR, BP = POS, P)
fwrite(ld_format, file.path(output_dir, "SSRI_SNRI_Acceptability_GWAS.temp"), sep = "\t")

#-- Save COJO format as tab-delimited
cojo_format <- cojo_temp %>%
select(-CHR, -POS,-REF, -ALT, -OR, -Z_STAT)
fwrite(cojo_format, file.path(output_dir, "SSRI_SNRI_Acceptability_GWAS.ma"), sep = "\t")

#-- Format results for a Supp Table
sig <- cojo_temp %>%
    mutate_at(vars(SE, P, A1_FREQ), ~ signif(., 2)) %>%
    mutate_at(vars(OR, Z_STAT), ~ round(., 2)) %>%
    select(CHR, POS, SNP, REF, ALT, A1, A2, A1_FREQ, OR, SE, Z_STAT, P, N) %>%
    filter(P < 0.000005) %>%
    arrange(P)
    
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
#removeWorksheet(wb, "Table18")
addWorksheet(wb, "Table18")
writeData(wb, "Table18", sig)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)

# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 1.001765


#---------------------------- SSRI self-report efficacy ------------------------------------

# load libraries
library(data.table)
library(qqman)
library(dplyr)
library(openxlsx)

# genomic inflation factor (lambda)
calculate_lambda <- function(pvals) {
  chisq_stats <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq_stats) / 0.456
  return(lambda)
}

# Path to GWAS results
input_file <- "/scratch/user/uqawal15/SSRI_Responder_GWAS.txt"
output_dir <- "/scratch/user/uqawal15"
gwas <- fread(input_file, fill = TRUE)
              
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)

# Make the Manhattan plot
png(file.path(output_dir, "SSRI_Response_GWAS_manhattan.png"), 
    bg = "white", 
    type = "cairo", 
    width = 2000,
    height = 1200,
    res = 300)

manhattan(gwas_cleaned, 
          chr="#CHROM", 
          bp="POS", 
          snp="ID", 
          p="P",
          cex = 0.6,
          cex.axis = 0.8,
          suggestiveline = -log10(5e-6),
          genomewideline = -log10(5e-8),
          main = "SSRI Self-report Efficacy GWAS")
dev.off()

# Make the qqplot
png(file.path(output_dir, "SSRI_Response_GWAS_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI Self-report Efficacy")
dev.off()


#-- Convert to COJO format
cojo_temp <- gwas_cleaned %>%
  mutate(
    BETA = Z_STAT * `LOG(OR)_SE`,
    A2 = case_when(
      A1 == REF ~ ALT, 
      A1 == ALT ~ REF,
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    CHR = `#CHROM`,
    POS,
    SNP = ID,   
    REF,
    ALT,
    A1 = A1,   
    A2 = A2, 
    A1_FREQ,      
    b = BETA,   
    OR,
    SE = `LOG(OR)_SE`,
    Z_STAT,
    P,                   
    N = OBS_CT             
  )
  
#-- Save file for LD clumping
ld_format <- cojo_temp %>%
select(SNP, CHR, BP = POS, P)
fwrite(ld_format, file.path(output_dir, "SSRI_Responder_GWAS.temp"), sep = "\t")

#-- Save COJO format as tab-delimited
cojo_format <- cojo_temp %>%
select(-CHR, -POS,-REF, -ALT, -OR, -Z_STAT)
fwrite(cojo_format, file.path(output_dir, "SSRI_Responder_GWAS.ma"), sep = "\t")

#-- Format results for a Supp Table
sig <- cojo_temp %>%
    mutate_at(vars(SE, P, A1_FREQ), ~ signif(., 2)) %>%
    mutate_at(vars(OR, Z_STAT), ~ round(., 2)) %>%
    select(CHR, POS, SNP, REF, ALT, A1, A2, A1_FREQ, OR, SE, Z_STAT, P, N) %>%
    filter(P < 0.000005) %>%
    arrange(P)
    
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
#removeWorksheet(wb, "Table19")
addWorksheet(wb, "Table19")
writeData(wb, "Table19", sig)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)

# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 1.009734

#---------------------------- SSRI+SNRI self-report efficacy ------------------------------------

# load libraries
library(data.table)
library(qqman)
library(dplyr)
library(openxlsx)

# genomic inflation factor (lambda)
calculate_lambda <- function(pvals) {
  chisq_stats <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq_stats) / 0.456
  return(lambda)
}

# Path to GWAS results
input_file <- "/scratch/user/uqawal15/SSRI_SNRI_Responder_GWAS.txt"
output_dir <- "/scratch/user/uqawal15"
gwas <- fread(input_file, fill = TRUE)
              
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)

# Make the Manhattan plot
png(file.path(output_dir, "SSRI_SNRI_Response_GWAS_manhattan.png"), 
    bg = "white", 
    type = "cairo", 
    width = 2000,
    height = 1200,
    res = 300)

manhattan(gwas_cleaned, 
          chr="#CHROM", 
          bp="POS", 
          snp="ID", 
          p="P",
          cex = 0.6,
          cex.axis = 0.8,
          suggestiveline = -log10(5e-6),
          genomewideline = -log10(5e-8),
          main = "SSRI+SNRI Self-report Efficacy GWAS")
dev.off()

# Make the qqplot
png(file.path(output_dir, "SSRI_SNRI_Response_GWAS_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI+SNRI Self-report Efficacy")
dev.off()


#-- Convert to COJO format
cojo_temp <- gwas_cleaned %>%
  mutate(
    BETA = Z_STAT * `LOG(OR)_SE`,
    A2 = case_when(
      A1 == REF ~ ALT,
      A1 == ALT ~ REF, 
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    CHR = `#CHROM`,
    POS,
    SNP = ID,
    REF,
    ALT,
    A1 = A1,
    A2 = A2,
    A1_FREQ, 
    b = BETA,
    OR,
    SE = `LOG(OR)_SE`,
    Z_STAT,
    P,  
    N = OBS_CT
  )
  
#-- Save file for LD clumping
ld_format <- cojo_temp %>%
select(SNP, CHR, BP = POS, P)
fwrite(ld_format, file.path(output_dir, "SSRI_SNRI_Responder_GWAS.temp"), sep = "\t")

#-- Save COJO format as tab-delimited
cojo_format <- cojo_temp %>%
select(-CHR, -POS,-REF, -ALT, -OR, -Z_STAT)
fwrite(cojo_format, file.path(output_dir, "SSRI_SNRI_Responder_GWAS.ma"), sep = "\t")

#-- Format results for a Supp Table
sig <- cojo_temp %>%
    mutate_at(vars(SE, P, A1_FREQ), ~ signif(., 2)) %>%
    mutate_at(vars(OR, Z_STAT), ~ round(., 2)) %>%
    select(CHR, POS, SNP, REF, ALT, A1, A2, A1_FREQ, OR, SE, Z_STAT, P, N) %>%
    filter(P < 0.000005) %>%
    arrange(P)
    
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
#removeWorksheet(wb, "Table20")
addWorksheet(wb, "Table20")
writeData(wb, "Table20", sig)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)

# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 1.00191


############ PLINK LD CLUMPING ################################


#---- SSRI Acceptability

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=clumped_SSRI_Acceptability
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped_SSRI_Acceptability_chr%a.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped_SSRI_Acceptability_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

/home/uqawal15/bin/plink1.9/plink \
  --bfile ${bfiledir}/imputed_chr${c} \
  --clump ${wkdir}/SSRI_Acceptability_GWAS.temp \
  --clump-p1 5e-8 \
  --clump-r2 0.5 \
  --clump-kb 250 \
  --out clumped_snps_SSRI_Acceptability_chr${c}


#---- SSRI and SNRI Acceptability

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=clumped_SSRI_SNRI_Acceptability
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/clumped_SSRI_SNRI_Acceptability_chr%a.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/clumped_SSRI_SNRI_Acceptability_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

/home/uqawal15/bin/plink1.9/plink \
  --bfile ${bfiledir}/imputed_chr${c} \
  --clump ${wkdir}/SSRI_SNRI_Acceptability_GWAS.temp \
  --clump-p1 5e-8 \
  --clump-r2 0.5 \
  --clump-kb 250 \
  --out clumped_snps_SSRI_SNRI_Acceptability_chr${c}
  
#---- SSRI self-report efficacy

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=clumped_SSRI_Acceptability
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped_SSRI_Acceptability_chr%a.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped_SSRI_Acceptability_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

/home/uqawal15/bin/plink1.9/plink \
  --bfile ${bfiledir}/imputed_chr${c} \
  --clump ${wkdir}/SSRI_Acceptability_GWAS.temp \
  --clump-p1 5e-8 \
  --clump-r2 0.5 \
  --clump-kb 250 \
  --out clumped_snps_SSRI_Acceptability_chr${c}


##############################################################################


# Checking direction of effect for certain SNPs within rTMS study

# load libraries
library(data.table)
library(qqman)
library(dplyr)

# Path to your GWAS results
input_file <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Responder_gwas_concat.txt"
output_dir <- "/scratch/user/uqawal15/"

gwas <- fread(input_file, 
              fill = TRUE)
              
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P))
  
sum_stats <- gwas_cleaned %>%
select(SNP = ID, P)
fwrite(sum_stats, file.path(output_dir, "SSRI_Responder_gwas_sumstats.txt"), sep = "\t")

# Example data frame
rTMS_df <- data.frame(
  SNP_rTMS = c("rs960995", "rs17164813", "rs8035452", 
  "rs595562", "rs4648426", "rs12487861", 
  "rs198475", "rs560681", "rs11265485", 
  "rs1934115", "rs872994", "rs5995416", 
  "rs2189698", "rs4243296", "rs6899975",
  "rs1014129", "rs626904", "rs17673232",
  "rs11942069", "rs447347", "rs2271926",
  "rs4646515", "rs229526", "rs373521",
  "rs1283468", "rs11956034"),
  logOR_rTMS = c(log(0.18), log(0.173), log(7.25),
  log(0.1603), log(0.1489), log(16),
  log(0.1489), log(0.09502), log(0.08403),
  log(0.04528), log(0.1769), log(5.61),
  log(0.1839), log(0.1769), log(0.1729),
  log(0.1407), log(0.1674), log(0.07451),
  log(0.1489), log(0.1654), log(0.1674),
  log(0.2503), log(0.1324), log(0.1963),
  log(0.07051), log(0.09502)),
  A1_rTMS = c("G", "A", "G",
  "G", "G", "A",
  "A", "G", "G",
  "C", "A", "A",
  "C", "G", "A",
  "A", "A", "A",
  "A", "A", "G",
  "G", "G", "A",
  "A", "A")
)

mine_df <- gwas_cleaned %>%
  filter(ID %in% rTMS_df$SNP_rTMS) %>%
  mutate(logOR_mine = log(OR)) %>%
  select(
    SNP_mine = ID,
    logOR_mine,
    A1_mine = A1
  )

both_df <- left_join(rTMS_df, mine_df, by = c("SNP_rTMS" = "SNP_mine"))


# Align directions
both_df$aligned_logOR_rTMS <- ifelse(both_df$A1_mine != both_df$A1_rTMS,
                                 -both_df$logOR_rTMS,
                                 both_df$logOR_rTMS)

# Plot
outdir = "/scratch/user/uqawal15"
png(
  filename = file.path(outdir, "rTMS_effect_size_comparisons.png"),
  width = 250, height = 200, units = 'mm',
  bg = "white", res = 600, type = c("cairo")
)
plot(both_df$logOR_mine, both_df$aligned_logOR_rTMS,
     xlab = "log(OR) - SSRI Acceptability Study",
     ylab = "log(OR) - rTMS Study (Aligned)",
     main = "Effect Size Comparison",
     pch = 19)
abline(0, 1, col = "blue", lty = 2)
dev.off()


# alignd SNPs 
both_df$direction_match <- sign(both_df$logOR_mine) == sign(both_df$aligned_logOR_rTMS)

# View how many align
table(both_df$direction_match)

# Show which SNPs align
aligned_snps <- both_df[both_df$direction_match == TRUE, "SNP_rTMS"]
aligned_snps

[1] "rs960995"  "rs17164813" "rs11265485" "rs1934115"  "rs872994"
 [6] "rs626904"   "rs17673232" "rs11942069" "rs229526"   "rs11956034"


# ZNF471, THSD7A, LY9, LOC, LOC, LHFPL6, GTDC1, GRID2, EXOSC7, ADAMTS2

gwas_cleaned %>% filter(
  ID %in% aligned_snps
)

only one (ZNF471) reached suggestive significance (6.29500e-08)






