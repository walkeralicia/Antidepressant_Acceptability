
#-- Set wkdir
workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

######################## ANOREXIA #################################

#-- COJO format

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --job-name=ANO_format
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANO_format.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANO_format.stderr

module load r/4.4.0-gfbf-2023a
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.4.0-gfbf-2023a

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

trait=ANO_LOO
gwas_file=Anorexia_Watson2019_exclAus.txt
ma_file=${trait}/${gwas_file}

Rscript  ${workingpath}/cojo_format_v8.R  \
  --file  ${workingpath}/${trait}/${gwas_file}  \
  --out  ${workingpath}/${trait}/${gwas_file}.ma   \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  FREQ   \
  --pvalue P  \
  --beta  OR  \
  --se  SE  \
  --samplesize  N
  

#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=ANO_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANO_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANO_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"

trait=ANO_LOO
gwas_file=Anorexia_Watson2019_exclAus.txt
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${ma_file}.ma --impute-summary --out ${ma_file}_imp.ma --thread 4


#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=ANO_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANO_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANO_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=ANO_LOO
gwas_file=Anorexia_Watson2019_exclAus.txt
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb


######################## BMI #################################

#-- COJO format

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --job-name=BMI_format
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BMI_format.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BMI_format.stderr

module load r/4.4.0-gfbf-2023a
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.4.0-gfbf-2023a

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

trait=BMI_LOO
gwas_file=BMI_yengo_2018_looQIMR.dat
ma_file=${trait}/${gwas_file}

Rscript  ${workingpath}/cojo_format_v8.R  \
  --file  ${workingpath}/${trait}/${gwas_file}  \
  --out  ${workingpath}/${trait}/${gwas_file}.ma   \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  Freq   \
  --pvalue P  \
  --beta  BETA  \
  --se  SE  \
  --samplesize  N

#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=BMI_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BMI_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BMI_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"

trait=BMI_LOO
gwas_file=BMI_yengo_2018_looQIMR.dat
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${workingpath}/${ma_file}.ma --impute-summary --out ${workingpath}/${ma_file}_imp.ma --thread 4


#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=BMI_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BMI_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BMI_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=BMI_LOO
gwas_file=BMI_yengo_2018_looQIMR.dat
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb

######################## BIP 2025 LOO #################################

#-- COJO format

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --job-name=BIP_format
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BIP_format.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BIP_format.stderr

module load r/4.4.0-gfbf-2023a
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.4.0-gfbf-2023a

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

trait=BIP_LOO
gwas_file=PGC4_BIP_EUR_LOO_w23andMe_noAUS.trans.clean.sumstats
ma_file=${trait}/${gwas_file}

Rscript  ${workingpath}/cojo_format_v8.R  \
  --file  ${workingpath}/${trait}/${gwas_file}  \
  --out  ${workingpath}/${trait}/${gwas_file}.ma   \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  FRQ_U_2309822   \
  --pvalue P  \
  --beta  OR  \
  --se  SE  \
  --samplesize  Nca,Nco
  
#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=BIP_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BIP_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BIP_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"

trait=BIP_LOO
gwas_file=PGC4_BIP_EUR_LOO_w23andMe_noAUS.trans.clean.sumstats
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${workingpath}/${ma_file}.ma --impute-summary --out ${workingpath}/${ma_file}_imp.ma --thread 4

#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=BIP_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BIP_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/BIP_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=BIP_LOO
gwas_file=PGC4_BIP_EUR_LOO_w23andMe_noAUS.trans.clean.sumstats
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb



######################## ANX 2025 LOO #################################

#-- COJO format

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --job-name=ANX_format
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANX_format.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANX_format.stderr

module load r/4.4.0-gfbf-2023a
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.4.0-gfbf-2023a

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

trait=ANX_LOO
gwas_file=daner_fullANX_v12_woQIMR_082224.txt
ma_file=${trait}/${gwas_file}

Rscript  ${workingpath}/cojo_format_v8.R  \
  --file  ${workingpath}/${trait}/${gwas_file}  \
  --out  ${workingpath}/${trait}/${gwas_file}.ma   \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  FRQ_U_709576   \
  --pvalue P  \
  --beta  OR  \
  --se  SE  \
  --samplesize  Nca,Nco

#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=ANX_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANX_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANX_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"

trait=ANX_LOO
gwas_file=daner_fullANX_v12_woQIMR_082224.txt
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${workingpath}/${ma_file}.ma --impute-summary --out ${workingpath}/${ma_file}_imp.ma --thread 4

#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=ANX_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANX_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/ANX_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=ANX_LOO
gwas_file=daner_fullANX_v12_woQIMR_082224.txt
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb


############## MDD 2024 LOO ##################################

#-- COJO format

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --job-name=MDD_format
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/MDD_format.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/MDD_format.stderr

module load r/4.4.0-gfbf-2023a
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.4.0-gfbf-2023a

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

trait=MDD_LOO
gwas_file=daner_pgc_mdd_w3_full_no.australia.txt
ma_file=${trait}/${gwas_file}

Rscript  ${workingpath}/cojo_format_v8.R  \
  --file  ${workingpath}/${trait}/${gwas_file}  \
  --out  ${workingpath}/${trait}/${gwas_file}.ma   \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  FRQ_U_3346858 \
  --pvalue P  \
  --beta  OR  \
  --se  SE  \
  --samplesize  Nca,Nco

#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=MDD_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/MDD_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/MDD_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/bin/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"

trait=MDD_LOO
gwas_file=daner_pgc_mdd_w3_full_no.australia.txt
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${workingpath}/${ma_file}.ma --impute-summary --out ${workingpath}/${ma_file}_imp.ma --thread 4

#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=MDD_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/MDD_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/MDD_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/bin/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=MDD_LOO
gwas_file=daner_pgc_mdd_w3_full_no.australia.txt
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb


############## Neuroticism 2018 ##################################

#-- COJO format
library(data.table)
library(dplyr)

dat <- fread("/QRISdata/Q7280/pharmacogenomics/pgs/sumstats/NEU_2018/sumstats_neuroticism_ctg_format.txt", header =TRUE)
dat <- dat %>%
mutate(BETA = Z / sqrt(2 * MAF_UKB * (1-MAF_UKB)* (N + Z^2)), 
SE = 1 / sqrt( 2 * MAF_UKB * ( 1 - MAF_UKB ) * ( N + Z^2 ) ))

cojo <- dat %>%
select(SNP = SNP,
A1 = A1,
A2 = A2,
freq = MAF_UKB,
beta = BETA,
se = SE,
pvalue = P,
samplesize = N)

fwrite(cojo, "/QRISdata/Q7280/pharmacogenomics/pgs/sumstats/NEU_2018/sumstats_neuroticism_ctg_format.txt.ma",
       sep = "\t", 
       quote = FALSE,
       row.names = FALSE,
       col.names = TRUE)
       
#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=NEU_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/NEU_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/NEU_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/bin/gctb_2.5.4_Linux"
ldm1="/QRISdata/Q6913/Pipeline/ukbEUR_Imputed/"

trait=NEU_2018
gwas_file=sumstats_neuroticism_ctg_format.txt
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${workingpath}/${ma_file}.ma --impute-summary --out ${workingpath}/${ma_file}_imp.ma --thread 4

#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=NEU_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/NEU_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/NEU_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/bin/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=NEU_2018
gwas_file=sumstats_neuroticism_ctg_format.txt
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb


############## OCD 2024 ##################################

#-- COJO format

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --job-name=OCD_format
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/OCD_format.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/OCD_format.stderr

module load r/4.4.0-gfbf-2023a
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.4.0-gfbf-2023a

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

trait=OCD_2024
gwas_file=ocs2024obsessive-compulsive_symptoms_daner_STR_NTR_SfS_TwinsUK_strometal.txt
ma_file=${trait}/${gwas_file}

Rscript  ${workingpath}/cojo_format_v8.R  \
  --file  ${workingpath}/${trait}/${gwas_file}  \
  --out  ${workingpath}/${trait}/${gwas_file}.ma   \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  FRQ_U_33943 \
  --pvalue P  \
  --beta  OR  \
  --se SE \
  --samplesize  Nca,Nco

#-- GCTB QC and Impute

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --job-name=OCD_impute
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/OCD_impute.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/OCD_impute.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"
cd $workingpath

path2gctb="/home/uqawal15/bin/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"

trait=OCD_2024
gwas_file=ocs2024obsessive-compulsive_symptoms_daner_STR_NTR_SfS_TwinsUK_strometal.txt
ma_file=${trait}/${gwas_file}

${path2gctb}/gctb --ldm-eigen $ldm1 --gwas-summary ${workingpath}/${ma_file}.ma --impute-summary --out ${workingpath}/${ma_file}_imp.ma --thread 4

#-- SBayesRC

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name=OCD_sbrc
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/OCD_sbrc.stdout
#SBATCH -e /QRISdata/Q7280/pharmacogenomics/pgs/sumstats/jobs/OCD_sbrc.stderr

workingpath="/QRISdata/Q7280/pharmacogenomics/pgs/sumstats"

path2gctb="/home/uqawal15/bin/gctb_2.5.2_Linux"
ldm1="/QRISdata/Q5514/For_Alicia./ukbEUR_Imputed"
annot="/QRISdata/Q6913/Pipeline/annot_baseline2.2.txt"

trait=OCD_2024
gwas_file=ocs2024obsessive-compulsive_symptoms_daner_STR_NTR_SfS_TwinsUK_strometal.txt
ma_file=${trait}/${gwas_file}

mkdir -p ${workingpath}/${trait}/SBayesRC

${path2gctb}/gctb \
--sbayes RC  \
--ldm-eigen   ${ldm1}   \
--gwas-summary   ${workingpath}/${ma_file}_imp.ma.imputed.ma   \
--annot  $annot  \
--thread 32 \
--out  ${workingpath}/${trait}/SBayesRC/${gwas_file}_imp.ma.imputed_sbrc_gctb
