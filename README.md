# Antidepressant_Acceptability

Scripts used to investigate the phenotypic and genetic heterogeneity of sustained antidepressant prescription patterns in individuals with Major Depressive Disorder (MDD). This repository accompanies the article *"Genetic and Phenotypic Associations with Sustained Antidepressant Use in Major Depressive Disorder"*.

## Overview

Using prescription dispensing records and survey, genotype, and polygenic score data from the Australian Genetics of Depression Study (AGDS), the pipeline:

1. Reconstructs antidepressant prescription episodes from raw dispensing records.
2. Assigns participants to sustained-use treatment groups (single drug, combination, lithium-treated bipolar, and mixed use) at defined duration thresholds.
3. Tests treatment groups and treatment characteristics for association with demographic, clinical, comorbidity, lifestyle, and symptom phenotypes.
4. Tests the same outcomes for association with a panel of polygenic scores.
5. Runs genome-wide association analyses of sustained SSRI and SSRI/SNRI use and of self-reported antidepressant response.

The ten most frequently dispensed N06A antidepressants are analysed, spanning SSRI, SNRI, TCA, and TeCA classes.

## Citation

If you use scripts from this repository, please cite the article: https://doi.org/10.1101/2025.08.19.25333784

## Repository structure

```
Antidepressant_Acceptability/
├── Drug_Reference_Table_AGDS.R     # ATC code to drug name and class reference table
├── scripts/
│   ├── S1_Extract_Prescriptions.R
│   ├── S2A_Duration_Groups.R
│   ├── S2B_Adherent_Duration_Groups.R
│   ├── S3_MDD_Sustained_AD_Use_Groups.R
│   ├── S4A_DTD.R
│   ├── S4A_DTD_SexStratified.R
│   ├── S4B_DTD_PGS.R
│   ├── S4B_DTD_PGS_SexStratified.R
│   ├── S5A_AD_Groups_Binary.R
│   ├── S5B_AD_Groups_Quantitative.R
│   ├── S6_AD_Groups_PGS.R
│   ├── S7_AD_Group_Multiple_Testing_Correction.R
│   ├── S8_BMI_adjusted.R
│   ├── S9_SSRI_acceptability_gwas.sh
│   └── S10_BMI_Mediation_Analyses.R
└── README.md
```

## Analysis pipeline

Scripts are numbered in the order they are run. Those sharing a number (e.g. S4A and S4B) are parallel analyses operating on common inputs.

### Reference table

**`Drug_Reference_Table_AGDS.R`** maps the ten antidepressant ATC codes of interest to drug names and drug classes (SSRI, SNRI, TCA, TeCA). It is sourced by the downstream scripts.

### Phenotype construction

**`S1_Extract_Prescriptions.R`** processes raw antidepressant dispensing records. It restricts to the ten antidepressants within the study window, estimates per-dispense supply duration, and aggregates consecutive dispenses into prescription episodes using a 90-day permissible gap.

**`S2A_Duration_Groups.R`** assigns participants to sustained-use treatment groups at duration thresholds of 360 and 600 days. It distinguishes single-antidepressant use, combination use (two or more drugs each dispensed for a sustained period, concatenated and resolved to a class where two or more drugs share a class), and a "Various" group for participants with no single sustained drug.

**`S2B_Adherent_Duration_Groups.R`** repeats the grouping in S2A under an adherent dispensing definition.

**`S3_MDD_Sustained_AD_Use_Groups.R`** restricts the sample to lifetime MDD cases and incorporates bipolar diagnosis (baseline and follow-up) together with lithium dispensing to define lithium-treated bipolar (BIP+L) and bipolar-without-lithium (BIP-L) groups. It links participants to QC-passed genotypes and European genetic ancestry, then produces the final treatment-group dataset mapped to drug names and classes across both duration and adherence definitions.

### Phenotypic and polygenic association analyses

**`S4A_DTD.R`** and **`S4A_DTD_SexStratified.R`** run phenotypic association analyses of per-drug and per-class treatment characteristics (including derived self-reported response variables for feeling well on, and stopping, each antidepressant) against demographic, clinical, comorbidity, lifestyle, and symptom variables, using generalised linear models and mixed models. The sex-stratified variant repeats the analyses within each sex.

**`S4B_DTD_PGS.R`** and **`S4B_DTD_PGS_SexStratified.R`** run the equivalent polygenic score association analyses, with principal-component adjustment, and a sex-stratified variant.

**`S5A_AD_Groups_Binary.R`** tests sustained-use treatment groups for association with binary phenotypes: clinical characteristics, risk factors, psychiatric and physical comorbidities, MDD symptoms, and concomitant medication dispenses.

**`S5B_AD_Groups_Quantitative.R`** tests treatment groups for association with quantitative phenotypes (e.g. age, BMI, education, prescription complexity measures).

**`S6_AD_Groups_PGS.R`** tests treatment groups for association with a panel of polygenic scores spanning psychiatric, anthropometric, metabolic, and related traits, restricted to European-ancestry participants with principal-component adjustment.

**`S7_AD_Group_Multiple_Testing_Correction.R`** applies multiple-testing correction jointly across the quantitative and binary association results.

### Sensitivity analyses

**`S8_BMI_adjusted.R`** repeats the treatment-group association analyses with additional adjustment for BMI and BMI polygenic score.

**`S10_BMI_Mediation_Analyses.R`** examines the relationship between BMI polygenic score, BMI, and treatment group, including the collider scenario in which conditioning on BMI polygenic score can induce spurious treatment-group associations.

### Genome-wide association analyses

**`S9_SSRI_acceptability_gwas.sh`** is an end-to-end workflow combining shell (SLURM) and R steps. It:

- Removes related individuals and prepares per-chromosome genotype files (PLINK2 KING cutoff).
- Defines four outcomes (sustained SSRI use, sustained SSRI/SNRI use, and self-reported SSRI and SSRI/SNRI response) and builds the corresponding covariate and phenotype files for European, unrelated participants.
- Runs association testing per outcome and per chromosome (PLINK2 logistic regression), then concatenates results.
- Computes genomic inflation, produces Manhattan and QQ plots, and exports formatted summary statistics (LD and COJO formats) and significant-hit tables.
- Performs LD clumping (PLINK 1.9) and gene-based testing (GCTA mBAT-combo).

## Software requirements

- **R** 
- **PLINK2** (association testing, relatedness filtering)
- **PLINK 1.9** (LD clumping)
- **GCTA** (mBAT-combo gene-based testing)
- An HPC environment with **SLURM** for the GWAS workflow

## Data availability

The scripts contain absolute file paths specific to the secure compute environment in which the analyses were run, and will need to be adapted to a local setup. AGDS prescription, phenotype, and genotype data are not publicly available and require approved access. `S1_Extract_Prescriptions.R` expects dispensing input in the column format of the referenced prescription data file (`ParticipantID`, `ATCCode`, `DateofSupply`).
