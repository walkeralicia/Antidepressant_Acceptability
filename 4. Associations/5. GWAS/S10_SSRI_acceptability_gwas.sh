# First option (AD class diversity)
drug_mapping <- data.frame(
  ATCCode = c('N06AB06', 'N06AB10', 'N06AB04', 'N06AX16', 'N06AX21', 
              'N06AA09', 'N06AB05', 'N06AX11', 'N06AX23', 'N06AB03'),
  DrugClass = c('SSRI', 'SSRI', 'SSRI', 'SNRI', 'SNRI', 'TCA', 'SSRI', 'TeCA', 'SNRI', 'SSRI')
)
drugs_mapped <- left_join(drugs, drug_mapping, by = "ATCCode")

pharma <- drugs_mapped %>%
group_by(ParticipantID) %>%
mutate(
  num_class = length(unique((DrugClass))
)) %>%
select(ParticipantID, num_class) %>%
unique()


#---- Create a covariate and phenotype file
R
library(dplyr)
library(openxlsx)

# Load required libraries
fam <- read.table("/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink/imputed_chr1.fam", header=FALSE)

# Define file paths
PC_FILE <- "/QRISdata/Q5338/Ancestry_analysis/AGDS_R11_TOPMedr2_pruned.05.common_pca3.proj.eigenvec"
DEMO_FILE <- "/QRISdata/Q7280/pharmacogenomics/phenotypes/BMI_age_sex.csv"
TREATMENT_FILE <- "/QRISdata/Q7280/pharmacogenomics/treatment_groups/Final_Treatmentgroups_360days.csv"
PRESCRIPTION_FILE <- "/QRISdata/Q7280/pharmacogenomics/data/AGDSAcceptabilityTreatmentGroups_14122024.csv"
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

# Identify participants not in SSRI, BIP-L, or BIP+L groups
remaining <- groups %>% 
  filter(!(DrugClass %in% c("SSRI") | DrugName %in% c("BIP-L", "BIP+L")))

# Filter prescription data to include only those from the remaining participants
drugs_incl <- drugs %>% 
  filter(ParticipantID %in% remaining$ParticipantID)

# Identify participants with long-term SSRI use (>360 days)
long_term_users <- drugs_incl %>%
  filter(ATCCode %in% c("N06AB03", "N06AB04", "N06AB05", "N06AB06", "N06AB10") & PrescriptionDays > 360) %>%
  distinct(ParticipantID)

# Remove long-term users from the dataset
cleaned_data <- drugs_incl %>%
  filter(!ParticipantID %in% long_term_users$ParticipantID)

# Create SSRI responder phenotype
# 1 = SSRI responder, 0 = non-responder, NA = not applicable
groups$SSRI_Responder <- ifelse(
  groups$ParticipantID %in% cleaned_data$ParticipantID, 0,
  ifelse(groups$DrugClass == "SSRI", 1, NA)
)

#------ Second option (SSRI+SNRI acceptability) -------------------------

# Identify participants not in SSRI, SNRI, BIP-L, or BIP+L groups
remaining <- groups %>% 
  filter(!(DrugClass %in% c("SSRI", "SNRI") | DrugName %in% c("BIP-L", "BIP+L")))

# Filter prescriptions data to include only those from the remaining participants
drugs_incl <- drugs %>% 
  filter(ParticipantID %in% remaining$ParticipantID)

# Identify participants with long-term SSRI or SNRI use (>360 days)
long_term_users <- drugs_incl %>%
  filter(ATCCode %in% c("N06AB03", "N06AB04", "N06AB05", "N06AB06", "N06AB10", "N06AX16", "N06AX21", "N06AX23") & PrescriptionDays > 360) %>%
  distinct(ParticipantID)

# Remove long-term users from the dataset
cleaned_data <- drugs_incl %>%
  filter(!ParticipantID %in% long_term_users$ParticipantID)

# Create SSRI+SNRI responder phenotype
# 1 = SSRI+SNRI responder, 0 = non-responder, NA = not applicable
groups$SSRI_SNRI_Responder <- ifelse(
  groups$ParticipantID %in% cleaned_data$ParticipantID, 0,
  ifelse(groups$DrugClass %in% c("SSRI", "SNRI"), 1, NA)
)

#---------- # Combine all datasets

dat <- pheno %>%
  full_join(groups, by = c("STUDYID" = "ParticipantID")) %>%
  #full_join(pharma, by = c("STUDYID" = "ParticipantID")) %>%
  full_join(pcs, by = c("IID" = "IID")) %>%
  filter(!is.na(IID)) %>%
  mutate(
    # Ensure AGE is numeric
    AGE = as.numeric(AGE),
    # Replace NA in SEX with "NONE"
    SEX = ifelse(is.na(SEX), "NONE", SEX)
  )

# Create and save covariate file
dat %>%
  select(FID, IID, AGE, SEX, PC1, PC2, PC3) %>%
  write.table(
    file.path(OUTPUT_DIR, "AD_acceptability_covariates.txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )

# Create and save outcome file
dat %>%
  select(FID, IID, DrugClass, SSRI_Responder, num_class) %>%
  write.table(
    file.path(OUTPUT_DIR, "AD_acceptability_outcome.txt"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )


#--------------- Run SSRI responder vs SSRI non-responder GWAS using PLINK2 (common variants)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=SSRI_response_gwas
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/SSRI_response_gwas_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/SSRI_response_gwas_chr%a.stderr

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
--freq \
--maf 0.01 \
--glm hide-covar \
--threads 8 \
--out SSRI_response_gwas_chr${c}_results

#---------------------------------------------

# Concatenate gwas results

```{bash}
# First process chr1, keeping header

head -n1 SSRI_response_gwas_chr1_results.SSRI_Responder.glm.logistic.hybrid > SSRI_Responder_gwas_concat.txt
head -n1 SSRI_response_gwas_chr1_results.afreq > SSRI_Responder_afreq_concat.txt

# Process each chromosome
for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    # Filter out CONST_OMITTED_ALLELE and append (no header)
    awk 'NR>1 && $NF!="CONST_OMITTED_ALLELE"' SSRI_response_gwas_chr${chr}_results.SSRI_Responder.glm.logistic.hybrid >> SSRI_Responder_gwas_concat.txt
    # append (no header) AF
    awk 'NR>1' SSRI_response_gwas_chr${chr}_results.afreq >> SSRI_Responder_afreq_concat.txt
done
echo "Done! Results saved."

```

#################################################################

# load libraries
library(data.table)
library(qqman)
library(dplyr

# genomic inflation factor (lambda) function

calculate_lambda <- function(pvals) {
  # Convert p-values to chi-square statistics (1 df)
  chisq_stats <- qchisq(1 - pvals, df = 1)
  
  # Calculate the genomic inflation factor
  lambda <- median(chisq_stats) / 0.456
  
  return(lambda)
}

# Path to GWAS results
input_file <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Responder_gwas_concat.txt"
freq <- fread("/QRISdata/Q7280/pharmacogenomics/associations/GWAS/SSRI_Responder_afreq_concat.txt", fill = TRUE)
output_dir <- "/QRISdata/Q7280/pharmacogenomics/associations/GWAS"

gwas <- fread(input_file, fill = TRUE)
              
gwas_cleaned <- gwas %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P)) %>%
  arrange(P)

# Make the Manhattan plot
png(file.path(output_dir, "SSRI_Acceptability_gwas_manhattan.png"), 
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
png(file.path(output_dir, "SSRI_Acceptability_gwas_qqplot.png"), 
    bg = "white", 
    type = "cairo", 
    width = 1800,
    height = 1800,
    res = 300)

qq(gwas_cleaned$P, main = "QQ Plot: SSRI Acceptability")
dev.off()

#-- Save sig hits
sig <- gwas_cleaned %>%
filter(P < 0.000005) %>%
arrange(P)

#-- Join allele frequency information to GWAS sumstats
sig_freq <- sig %>%
  left_join(freq, by = c("#CHROM", "ID", "REF", "ALT"))
  
#-- Fix format of summary statistics
gwas_flipped <- sig_freq %>%
  mutate(
    # Identify rows to flip
    flip = REF != A1,
    
    # Flip A1 allele
    A1_new = case_when(
      flip & A1 == ALT ~ REF,
      flip & A1 == REF ~ ALT,
      TRUE ~ A1
    ),
    
    # Update statistics for flipped SNPs
    ALT_FREQS_new = ifelse(flip, 1 - ALT_FREQS, ALT_FREQS),
    OR_new = ifelse(flip, 1/OR, OR),
    Z_STAT_new = ifelse(flip, -Z_STAT, Z_STAT)
  ) %>%
  # Replace original columns
  mutate(
    A1 = A1_new,
    ALT_FREQS = ALT_FREQS_new,
    OR = OR_new,
    Z_STAT = Z_STAT_new
  ) %>%
  # Remove helper columns
  select(-flip, -A1_new, -ALT_FREQS_new, -OR_new, -Z_STAT_new)

gwas_final <- gwas_flipped %>%
    mutate_at(vars(`LOG(OR)_SE`, P, ALT_FREQS), ~ signif(., 2)) %>%
    mutate_at(vars(OR, Z_STAT), ~ round(., 2)) %>%
    select(`#CHROM`, POS, ID, REF, ALT, A1, OR, `LOG(OR)_SE`, Z_STAT, P, ALT_FREQS)
    
wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
removeWorksheet(wb, "Table17")
addWorksheet(wb, "Table17")
writeData(wb, "Table17", gwas_final)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)


fwrite(gwas_final, file.path(output_dir, "SSRI_Acceptability_gwas_top_hits.txt"))


# Lambda
lambda <- calculate_lambda(gwas_cleaned$P) # 1.006737



# Move results to QRISdata
/QRISdata/Q7280/pharmacogenomics/associations/GWAS




####################


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



############ PLINK LD CLUMPING ################################


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=clumped_SSRI_responder
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH --array=1-22
#SBATCH -o /scratch/user/uqawal15/clumped_SSRI_responder_chr%a.stdout
#SBATCH -e /scratch/user/uqawal15/clumped_SSRI_responder_chr%a.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

wkdir="/scratch/user/uqawal15"
bfiledir="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink"
cd ${wkdir}

/home/uqawal15/bin/plink1.9/plink \
  --bfile ${bfiledir}/imputed_chr${c} \
  --clump ${wkdir}/SSRI_Responder_gwas_sumstats.txt \
  --clump-p1 5e-8 \
  --clump-r2 0.5 \
  --clump-kb 250 \
  --out clumped_snps_SSRI_responder_chr${c}








