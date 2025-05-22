

#======= Extract path and file names for each PGS trait
library(data.table)
library(dplyr)

#-- Summary statistics collated by Tian
trait_list=read.csv("/QRISdata/Q6913/Predictor_Release_20250317.csv")
weights="/QRISdata/Q6913"

pred_list <- trait_list %>%
  select(Predictor, predictor_file)
pred_list <- pred_list[c(1:98, 149),]

write.csv(pred_list, "/QRISdata/Q7280/pharmacogenomics/pgs/weights/Predictor_Release_20250317.csv", quote = FALSE, row.names=FALSE)

#-- Summary statistics from AGDS-LOO and CNT_03
datdir="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
Trait <- c("ANO_LOO", "ANX_LOO", "BIP_LOO", "BMI_LOO","MDD_LOO","CNT_03", "OCD_2024")
predictor_files = c(file.path(datdir, Trait[1], "SBayesRC", "Anorexia_Watson2019_exclAus.txt_imp.ma.imputed_sbrc_gctb.snpRes"),
file.path(datdir,Trait[2], "SBayesRC", "daner_fullANX_v12_woQIMR_082224.txt_imp.ma.imputed_sbrc_gctb.snpRes"),
file.path(datdir, Trait[3], "SBayesRC", "PGC4_BIP_EUR_LOO_w23andMe_noAUS.trans.clean.sumstats_imp.ma.imputed_sbrc_gctb.snpRes"),
file.path(datdir, Trait[4], "SBayesRC", "BMI_yengo_2018_looQIMR.dat_imp.ma.imputed_sbrc_gctb.snpRes"),
file.path(datdir, Trait[5], "SBayesRC", "daner_pgc_mdd_w3_full_no.australia.txt_imp.ma.imputed_sbrc_gctb.snpRes"),
file.path("/QRISdata/Q7672/Restricted/Binary_Traits", Trait[6], "SBayesRC/GCTB_v2.5.2", "METAANALYSIS1.TBL_imp.ma.imputed_sbrc_gctb.snpRes"),
file.path(datdir, Trait[7], "SBayesRC", "ocs2024obsessive-compulsive_symptoms_daner_STR_NTR_SfS_TwinsUK_strometal.txt_imp.ma.imputed_sbrc_gctb.snpRes"))

pred_list_2 <- data.frame(
  Predictor = Trait,
  predictor_file = predictor_files)

#-- combine both sets of predictors for the AGDS  
agds_pred_list <- rbind(pred_list, pred_list_2)
write.csv(agds_pred_list, "/QRISdata/Q7280/pharmacogenomics/pgs/weights/AGDS_weights.csv", col.names = FALSE, row.names = FALSE, quote = FALSE)

#======= NEW: This bash job generates PGS for multiple traits =================

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=agds_pgs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/scores/jobs/agds_pgs.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/scores/jobs/agds_pgs.stderr

# Set the chromosome based on SLURM array index
c=${SLURM_ARRAY_TASK_ID}

# Load modules
module load plink/2.00a3.6-gcc-11.3.0

# Locations
outdir="/QRISdata/Q7280/pharmacogenomics/pgs/scores"
bfile="/QRISdata/Q5338/Genotypes/AGDS_R11_TOPMedr2/Plink_subset_for_PRS_profiling/AGDS_autosomes_7.3M_subset"
trait_list="/QRISdata/Q7280/pharmacogenomics/pgs/weights/AGDS_weights.csv"

# Process the trait list - use space as delimiter instead of comma
while IFS=' ' read -r trait_code gwas_path; do
  echo "Processing trait: $trait_code"
  
  # Create output directory if it doesn't exist
  mkdir -p "${outdir}/${trait_code}"
  
  echo "GWAS Filename: $gwas_path"
  
  # Define a clean output path
  out_path="${outdir}/${trait_code}/${trait_code}_agds_sbrc_gctb_plink2"
  
  echo "Output will be written to: $out_path"
  
  # Run plink2 with properly quoted paths
  plink2 --bfile "${bfile}" \
    --score "${gwas_path}" 2 5 8 header cols=+scoresums \
    --threads 1 \
    --out "${out_path}"
  
  # Check if plink2 was successful
  if [ $? -eq 0 ]; then
    echo "PGS generation successful for $trait_code"
  else
    echo "Error: PGS generation failed for $trait_code"
    echo "Command was: plink2 --bfile ${bfile} --score ${gwas_path} 2 5 8 header cols=+scoresums --threads 1 --out ${out_path}"
  fi
  
done < "$trait_list"

echo "All PGS calculations completed"
