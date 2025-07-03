
#------------ Remove related individuals ---------------------------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --job-name=AGDS_unrelated
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /scratch/user/uqawal15/AGDS_unrelated.stdout
#SBATCH -e /scratch/user/uqawal15/AGDS_unrelated.stderr

module load plink/2.00a3.6-gcc-11.3.0
wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

plink2 --bfile ${bfiledir}/imputed_chr1 --king-cutoff 0.0884 --make-bed --out ${wkdir}/AGDS_chr22_unrelated

#------------ 

# Extract the list of individuals that were kept
cut -f1-2 AGDS_chr22_unrelated.fam > keep_individuals.txt

#----------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=AGDS_remove_related
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/AGDS_remove_related_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/AGDS_remove_related_chr%a.stderr

# Set the chromosome based on SLURM array index
chr=${SLURM_ARRAY_TASK_ID}

module load plink/2.00a3.6-gcc-11.3.0
wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

# Apply this filter to all chromosomes
for chr in {1..22}; do
    plink2 --bfile ${bfiledir}/imputed_chr${chr} --keep ${wkdir}/keep_individuals.txt --make-bed --out ${wkdir}/AGDS_unrelated_chr${chr}
done


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
  
#-- Filter for those of EUR ancestry 14603
eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id")
dat_eur <- dat %>%
  filter(IID %in% eur$V2)
  
#-- For those that are unrelated # 14026
unrelated <- read.table("/QRISdata/Q7280/pharmacogenomics/keep_individuals.txt", header = FALSE)
dat_eur_unrelated <- dat_eur %>%
filter(IID %in% unrelated$V2)

# Create and save covariate file
dat_eur_unrelated %>%
  select(FID, IID, AGE, SEX, PC1, PC2, PC3) %>%
  write.table(
    file.path(OUTPUT_DIR, "AD_acceptability_covariates.txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )

# Create and save outcome file
dat_eur_unrelated %>%
  select(FID, IID, DrugClass, SSRI_Acceptability, SSRI_SNRI_Acceptability, SSRI_Responder, SSRI_SNRI_Responder) %>%
  write.table(
    file.path(OUTPUT_DIR, "AD_acceptability_outcome.txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )


#-- Check case and control number for each phenotye
dat_eur_unrelated %>%
select(STUDYID, SSRI_Acceptability, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_Acceptability)
  SSRI_Acceptability    n
1                  0 4561
2                  1 3423
3                 NA 1479

dat_eur_unrelated %>%
select(STUDYID, SSRI_SNRI_Acceptability, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_SNRI_Acceptability)
  SSRI_SNRI_Acceptability    n
1                       0 2270
2                       1 5554
3                      NA 1639

dat_eur_unrelated %>%
select(STUDYID, SSRI_Responder, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_Responder)
  SSRI_Responder    n
1              0 2653
2              1 6474
3             NA  336

dat_eur_unrelated %>%
select(STUDYID, SSRI_SNRI_Responder, AGE, SEX, PC1, PC2, PC3) %>%
filter(!is.na(AGE) & SEX != "NONE" & !is.na(PC1) & !is.na(PC2) & !is.na(PC3)) %>%
count(SSRI_SNRI_Responder)

  SSRI_SNRI_Responder    n
1                   0  904
2                   1 8223
3                  NA  336



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
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/AGDS_unrelated_chr${c} \
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
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/AGDS_unrelated_chr${c} \
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
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/AGDS_unrelated_chr${c} \
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
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
cd ${wkdir}

plink2 \
--bfile ${bfiledir}/AGDS_unrelated_chr${c} \
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
input_file <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/SSRI_Acceptability_GWAS.txt"
output_dir <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability"
gwas <- fread(input_file, fill = TRUE)
              
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)
  
top_snps <- gwas_cleaned[gwas_cleaned$P < 5e-6, ]

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
          highlight = top_snps$ID,
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
    
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
#removeWorksheet(wb, "Table19")
addWorksheet(wb, "Table19")
writeData(wb, "Table19", sig)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)

# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 0.9992045

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
input_file <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/SSRI_SNRI_Acceptability_GWAS.txt"
output_dir <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability"
gwas <- fread(input_file, fill = TRUE)
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)

top_snps <- gwas_cleaned[gwas_cleaned$P < 5e-6, ]
            
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
          highlight = top_snps$ID,
          suggestiveline = -log10(5e-6),
          genomewideline = -log10(5e-8),
          main = "SSRI/SNRI Acceptability GWAS")
dev.off()

# Make the qqplot
png(file.path(output_dir, "SSRI_SNRI_Acceptability_GWAS_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI/SNRI Acceptability")
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
fwrite(ld_format, file.path(output_dir, "SSRI_SNRI_Acceptability.temp"), sep = "\t")

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
lambda <- calculate_lambda(gwas_cleaned$P) # 0.9974256


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
#removeWorksheet(wb, "Table21")
addWorksheet(wb, "Table21")
writeData(wb, "Table21", sig)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)

# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 1.004621

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
input_file <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Responder/SSRI_SNRI_Responder_GWAS.txt"
output_dir <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Responder"
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
          main = "SSRI/SNRI Self-report Efficacy GWAS")
dev.off()

# Make the qqplot
png(file.path(output_dir, "SSRI_SNRI_Response_GWAS_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI/SNRI Self-report Efficacy")
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
lambda <- calculate_lambda(gwas_cleaned$P) # 0.9992558


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
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped/clumped_SSRI_Acceptability_chr%a.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped/clumped_SSRI_Acceptability_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/"
outdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability/clumped"
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
cd ${wkdir}

/home/uqawal15/bin/plink1.9/plink \
  --bfile ${bfiledir}/AGDS_unrelated_chr${c} \
  --clump ${wkdir}/SSRI_Acceptability_GWAS.temp \
  --clump-p1 5e-8 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --out ${outdir}/clumped_snps_SSRI_Acceptability_chr${c}


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
r qu
wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/"
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
cd ${wkdir}

/home/uqawal15/bin/plink1.9/plink \
  --bfile ${bfiledir}/AGDS_unrelated_chr${c} \
  --clump ${wkdir}/SSRI_SNRI_Acceptability_GWAS.temp \
  --clump-p1 5e-8 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --out clumped_snps_SSRI_SNRI_Acceptability_chr${c}
  
#-------------- mBAT-combo on SSRI and SNRI Acceptability 

#!/bin/bash
#SBATCH --job-name=mBAT_all_genes_SSRI_SNRI_Acceptability
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-22
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/mBAT/mBAT_all_genes_SSRI_SNRI_Acceptability_chr%a.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/mBAT/mBAT_all_genes_SSRI_SNRI_Acceptability_chr%a.stderr

# Define chromosome from SLURM array task ID
c=${SLURM_ARRAY_TASK_ID}

# Define paths
mbat_version="/home/uqawal15/bin/gcta64_1.94.4/gcta64"
glist="/QRISdata/Q5059/mBAT/mBAT-main/hg19.pos.ensgid.all.gene.loc.extendedMHCexcluded"
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability"
# /scratch/project/genetic_data_analysis/uqali4/Repseudo_RNA.gene.loc.format.hT.autosome (brain expressed genes)

# Run mBAT for chromosome
$mbat_version \
--bfile ${bfiledir}/AGDS_unrelated_chr${c} \
--mBAT-combo ${wkdir}/SSRI_SNRI_Acceptability_GWAS.ma \
--mBAT-gene-list ${glist} \
--diff-freq 0.2 \
--mBAT-wind 50 \
--mBAT-svd-gamma 0.9 \
--mBAT-write-snpset \
--mBAT-print-all-p \
--out ${wkdir}/mBAT/mBAT_all_genes_SSRI_SNRI_Acceptability_chr${c}  \
--thread-num 6


#-------------- mBAT-combo on SSRI Acceptability 

#!/bin/bash
#SBATCH --job-name=mBAT_all_genes_SSRI_Acceptability
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-22
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/mBAT/mBAT_all_genes_SSRI_Acceptability_chr%a.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_SNRI_Acceptability/mBAT/mBAT_all_genes_SSRI_Acceptability_chr%a.stderr

# Define chromosome from SLURM array task ID
c=${SLURM_ARRAY_TASK_ID}

# Define paths
mbat_version="/home/uqawal15/bin/gcta64_1.94.4/gcta64"
glist="/QRISdata/Q5059/mBAT/mBAT-main/hg19.pos.ensgid.all.gene.loc.extendedMHCexcluded"
bfiledir="/QRISdata/Q7280/pharmacogenomics/genotypes/"
wkdir="/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Acceptability"
# /scratch/project/genetic_data_analysis/uqali4/Repseudo_RNA.gene.loc.format.hT.autosome (brain expressed genes)

# Run mBAT for chromosome
$mbat_version \
--bfile ${bfiledir}/AGDS_unrelated_chr${c} \
--mBAT-combo ${wkdir}/SSRI_Acceptability_GWAS.ma \
--mBAT-gene-list ${glist} \
--diff-freq 0.2 \
--mBAT-wind 50 \
--mBAT-svd-gamma 0.9 \
--mBAT-write-snpset \
--mBAT-print-all-p \
--out ${wkdir}/mBAT/mBAT_all_genes_SSRI_Acceptability_chr${c}  \
--thread-num 6
