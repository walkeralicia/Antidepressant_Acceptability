

#-- Load R libraries
library(dplyr)
library(tidyr)

#-- set wkdir
wkdir="/QRISdata/Q7280/pharmacogenomics/"

#-- Set treatment duration
duration = 600

#-- Read in phenotypes
pheno <- read.csv(file.path(wkdir, "phenotypes/survey_phenotypes.csv"))

#-- Function to create new columns dynamically
create_treatment_columns <- function(data, prefix) {
  column_name <- paste0(prefix, "AD")
  
  # Dynamically create the new column using mutate and case_when
  data <- data %>%
    mutate(
      !!column_name := case_when(
        DrugName == "SSRI:Sertraline" & !!sym(paste0(prefix, "SERT")) == 1 ~ 1,
        DrugName == "SSRI:Sertraline" & !!sym(paste0(prefix, "SERT")) == 0 ~ 0,
        DrugName == "SSRI:Escitalopram" & !!sym(paste0(prefix, "ESCI")) == 1 ~ 1,
        DrugName == "SSRI:Escitalopram" & !!sym(paste0(prefix, "ESCI")) == 0 ~ 0,
        DrugName == "SSRI:Citalopram" & !!sym(paste0(prefix, "CITA")) == 1 ~ 1,
        DrugName == "SSRI:Citalopram" & !!sym(paste0(prefix, "CITA")) == 0 ~ 0,
        DrugName == "SSRI:Fluoxetine" & !!sym(paste0(prefix, "FLUO")) == 1 ~ 1,
        DrugName == "SSRI:Fluoxetine" & !!sym(paste0(prefix, "FLUO")) == 0 ~ 0,
        DrugName == "SSRI:Paroxetine" & !!sym(paste0(prefix, "PARO")) == 1 ~ 1,
        DrugName == "SSRI:Paroxetine" & !!sym(paste0(prefix, "PARO")) == 0 ~ 0,
        DrugName == "SNRI:Desvenlafaxine" & !!sym(paste0(prefix, "DESV")) == 1 ~ 1,
        DrugName == "SNRI:Desvenlafaxine" & !!sym(paste0(prefix, "DESV")) == 0 ~ 0,
        DrugName == "SNRI:Venlafaxine" & !!sym(paste0(prefix, "VENL")) == 1 ~ 1,
        DrugName == "SNRI:Venlafaxine" & !!sym(paste0(prefix, "VENL")) == 0 ~ 0,
        DrugName == "SNRI:Duloxetine" & !!sym(paste0(prefix, "DULO")) == 1 ~ 1,
        DrugName == "SNRI:Duloxetine" & !!sym(paste0(prefix, "DULO")) == 0 ~ 0,
        DrugName == "TCA:Amitriptyline" & !!sym(paste0(prefix, "AMIT")) == 1 ~ 1,
        DrugName == "TCA:Amitriptyline" & !!sym(paste0(prefix, "AMIT")) == 0 ~ 0,
        DrugName == "TeCA:Mirtazapine" & !!sym(paste0(prefix, "MIRT")) == 1 ~ 1,
        DrugName == "TeCA:Mirtazapine" & !!sym(paste0(prefix, "MIRT")) == 0 ~ 0,
        TRUE ~ NA_real_
      )
    )
  return(data)
}


#-- Load treatment data
ad <- read.csv(file.path(wkdir, paste0("treatment_groups/Final_Treatmentgroups_", duration, "days.csv")))

#=== Left merge treatment groups with phenotypes for summaries and analyses ===
inner <- left_join(ad, pheno, by = c("ParticipantID" = "STUDYID")) %>%
  mutate(
    WELLAD = case_when(
      DrugName == "SSRI:Sertraline" & Sertraline == 1 ~ 1,
      DrugName == "SSRI:Sertraline" & Sertraline == 0 ~ 0,
      DrugName == "SSRI:Escitalopram" & Escitalopram == 1 ~ 1,
      DrugName == "SSRI:Escitalopram" & Escitalopram == 0 ~ 0,
      DrugName == "SSRI:Citalopram" & Citalopram == 1 ~ 1,
      DrugName == "SSRI:Citalopram" & Citalopram == 0 ~ 0,
      DrugName == "SSRI:Fluoxetine" & Fluoxetine == 1 ~ 1,
      DrugName == "SSRI:Fluoxetine" & Fluoxetine == 0 ~ 0,
      DrugName == "SSRI:Paroxetine" & Paroxetine == 1 ~ 1,
      DrugName == "SSRI:Paroxetine" & Paroxetine == 0 ~ 0,
      DrugName == "SNRI:Desvenlafaxine" & Desvenlafaxine == 1 ~ 1,
      DrugName == "SNRI:Desvenlafaxine" & Desvenlafaxine == 0 ~ 0,
      DrugName == "SNRI:Venlafaxine" & Venlafaxine == 1 ~ 1,
      DrugName == "SNRI:Venlafaxine" & Venlafaxine == 0 ~ 0,
      DrugName == "SNRI:Duloxetine" & Duloxetine == 1 ~ 1,
      DrugName == "SNRI:Duloxetine" & Duloxetine == 0 ~ 0,
      DrugName == "TCA:Amitriptyline" & Amitriptyline == 1 ~ 1,
      DrugName == "TCA:Amitriptyline" & Amitriptyline == 0 ~ 0,
      DrugName == "TCA:Mirtazapine" & Mirtazapine == 1 ~ 1,
      DrugName == "TCA:Mirtazapine" & Mirtazapine == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WELLAD_any= case_when(
      DrugName == "SSRI:Sertraline" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Sertraline" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Escitalopram" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Escitalopram" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Citalopram" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Citalopram" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Fluoxetine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Fluoxetine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Paroxetine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Paroxetine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SNRI:Desvenlafaxine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SNRI:Desvenlafaxine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SNRI:Venlafaxine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SNRI:Venlafaxine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SNRI:Duloxetine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SNRI:Duloxetine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "TCA:Amitriptyline" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "TCA:Amitriptyline" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "TCA:Mirtazapine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "TCA:Mirtazapine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "BIP-L" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "BIP-L" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "BIP+L" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "BIP+L" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "Various" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "Various" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      TRUE ~ NA_real_
    ),
    STOPAD_any = case_when(
      DrugName == "SSRI:Sertraline" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Sertraline" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Escitalopram" & (STOPESCI == 1 | STOPSERT == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Escitalopram" & (STOPESCI != 1 & STOPSERT != 1 & STOPCITA != 1 & STOPFLUO!= 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Citalopram" & (STOPCITA == 1 |  STOPESCI == 1 | STOPSERT == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Citalopram" & (STOPCITA != 1 & STOPESCI != 1 & STOPSERT != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Fluoxetine" & (STOPFLUO == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPSERT == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Fluoxetine" & (STOPFLUO != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPSERT != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Paroxetine" & (STOPPARO == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPSERT == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Paroxetine" & (STOPPARO != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPSERT != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SNRI:Desvenlafaxine" & (STOPDESV == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPSERT == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SNRI:Desvenlafaxine" & (STOPDESV != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPSERT != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SNRI:Venlafaxine" & (STOPVENL == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPSERT == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SNRI:Venlafaxine" & (STOPVENL != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPSERT != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SNRI:Duloxetine" & (STOPDULO == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPSERT == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SNRI:Duloxetine" & (STOPDULO != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPSERT != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "TCA:Amitriptyline" & (STOPAMIT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPSERT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "TCA:Amitriptyline" & (STOPAMIT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPSERT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "TCA:Mirtazapine" & (STOPMIRT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPSERT == 1) ~ 1,
      DrugName == "TCA:Mirtazapine" & (STOPMIRT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPSERT != 1) ~ 0,
      DrugName == "BIP-L" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "BIP-L" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "BIP+L" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "BIP+L" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "Various" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "Various" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      TRUE ~ NA_real_
    )
  )

# Use the function to create multiple columns
inner <- inner %>%
  create_treatment_columns("DRYM") %>%
  create_treatment_columns("SWET") %>%
  create_treatment_columns("NAUS") %>%
  create_treatment_columns("VOM") %>%
  create_treatment_columns("DIAR") %>%
  create_treatment_columns("CONS") %>%
  create_treatment_columns("HEAD") %>%
  create_treatment_columns("DIZZ") %>%
  create_treatment_columns("SHAK") %>%
  create_treatment_columns("MUSC") %>%
  create_treatment_columns("DROW") %>%
  create_treatment_columns("WAKE") %>%
  create_treatment_columns("HANX") %>%
  create_treatment_columns("AGIT") %>%
  create_treatment_columns("FATI") %>%
  create_treatment_columns("WTGN") %>%
  create_treatment_columns("WTLS") %>%
  create_treatment_columns("RASH") %>%
  create_treatment_columns("RUN") %>%
  create_treatment_columns("RSEX") %>%
  create_treatment_columns("BLUR") %>%
  create_treatment_columns("SUIT") %>%
  create_treatment_columns("NOSE") %>%
  create_treatment_columns("STOP")

write.csv(inner, paste0("/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes/inner_data_", duration, "days.csv"), row.names = FALSE)


#=== Full merge treatment groups with phenotypes for some analyses ===

full <- full_join(ad, pheno, by = c("ParticipantID" = "STUDYID")) %>%
  mutate(
    WELLAD = case_when(
      DrugName == "SSRI:Sertraline" & Sertraline == 1 ~ 1,
      DrugName == "SSRI:Sertraline" & Sertraline == 0 ~ 0,
      DrugName == "SSRI:Escitalopram" & Escitalopram == 1 ~ 1,
      DrugName == "SSRI:Escitalopram" & Escitalopram == 0 ~ 0,
      DrugName == "SSRI:Citalopram" & Citalopram == 1 ~ 1,
      DrugName == "SSRI:Citalopram" & Citalopram == 0 ~ 0,
      DrugName == "SSRI:Fluoxetine" & Fluoxetine == 1 ~ 1,
      DrugName == "SSRI:Fluoxetine" & Fluoxetine == 0 ~ 0,
      DrugName == "SSRI:Paroxetine" & Paroxetine == 1 ~ 1,
      DrugName == "SSRI:Paroxetine" & Paroxetine == 0 ~ 0,
      DrugName == "SNRI:Desvenlafaxine" & Desvenlafaxine == 1 ~ 1,
      DrugName == "SNRI:Desvenlafaxine" & Desvenlafaxine == 0 ~ 0,
      DrugName == "SNRI:Venlafaxine" & Venlafaxine == 1 ~ 1,
      DrugName == "SNRI:Venlafaxine" & Venlafaxine == 0 ~ 0,
      DrugName == "SNRI:Duloxetine" & Duloxetine == 1 ~ 1,
      DrugName == "SNRI:Duloxetine" & Duloxetine == 0 ~ 0,
      DrugName == "TCA:Amitriptyline" & Amitriptyline == 1 ~ 1,
      DrugName == "TCA:Amitriptyline" & Amitriptyline == 0 ~ 0,
      DrugName == "TeCA:Mirtazapine" & Mirtazapine == 1 ~ 1,
      DrugName == "TeCA:Mirtazapine" & Mirtazapine == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WELLAD_any= case_when(
      DrugName == "SSRI:Sertraline" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Sertraline" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Escitalopram" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Escitalopram" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Citalopram" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Citalopram" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Fluoxetine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Fluoxetine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SSRI:Paroxetine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SSRI:Paroxetine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SNRI:Desvenlafaxine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SNRI:Desvenlafaxine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SNRI:Venlafaxine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SNRI:Venlafaxine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "SNRI:Duloxetine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "SNRI:Duloxetine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "TCA:Amitriptyline" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "TCA:Amitriptyline" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "TCA:Mirtazapine" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "TCA:Mirtazapine" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "BIP-L" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "BIP-L" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "BIP+L" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "BIP+L" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      DrugName == "Various" & (Sertraline == 1 | Escitalopram == 1 | Citalopram == 1 | Fluoxetine == 1 | Paroxetine == 1 | Desvenlafaxine == 1 | Venlafaxine == 1 | Duloxetine == 1 | Amitriptyline == 1 | Mirtazapine == 1) ~ 1,
      DrugName == "Various" & (Sertraline != 1 & Escitalopram != 1 & Citalopram != 1 & Fluoxetine != 1 & Paroxetine != 1 & Desvenlafaxine != 1 & Venlafaxine != 1 & Duloxetine != 1 & Amitriptyline != 1 & Mirtazapine != 1) ~ 0,
      TRUE ~ NA_real_
    ),
    STOPAD_any = case_when(
      DrugName == "SSRI:Sertraline" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Sertraline" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Escitalopram" & (STOPESCI == 1 | STOPSERT == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Escitalopram" & (STOPESCI != 1 & STOPSERT != 1 & STOPCITA != 1 & STOPFLUO!= 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Citalopram" & (STOPCITA == 1 |  STOPESCI == 1 | STOPSERT == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Citalopram" & (STOPCITA != 1 & STOPESCI != 1 & STOPSERT != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Fluoxetine" & (STOPFLUO == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPSERT == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Fluoxetine" & (STOPFLUO != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPSERT != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SSRI:Paroxetine" & (STOPPARO == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPSERT == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SSRI:Paroxetine" & (STOPPARO != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPSERT != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SNRI:Desvenlafaxine" & (STOPDESV == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPSERT == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SNRI:Desvenlafaxine" & (STOPDESV != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPSERT != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SNRI:Venlafaxine" & (STOPVENL == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPSERT == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SNRI:Venlafaxine" & (STOPVENL != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPSERT != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "SNRI:Duloxetine" & (STOPDULO == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPSERT == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "SNRI:Duloxetine" & (STOPDULO != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPSERT != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "TCA:Amitriptyline" & (STOPAMIT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPSERT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "TCA:Amitriptyline" & (STOPAMIT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPSERT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "TCA:Mirtazapine" & (STOPMIRT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPSERT == 1) ~ 1,
      DrugName == "TCA:Mirtazapine" & (STOPMIRT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPSERT != 1) ~ 0,
      DrugName == "BIP-L" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "BIP-L" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "BIP+L" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "BIP+L" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      DrugName == "Various" & (STOPSERT == 1 | STOPESCI == 1 | STOPCITA == 1 | STOPFLUO == 1 | STOPPARO == 1 | STOPDESV == 1 | STOPVENL == 1 | STOPDULO == 1 | STOPAMIT == 1 | STOPMIRT == 1) ~ 1,
      DrugName == "Various" & (STOPSERT != 1 & STOPESCI != 1 & STOPCITA != 1 & STOPFLUO != 1 & STOPPARO != 1 & STOPDESV != 1 & STOPVENL != 1 & STOPDULO != 1 & STOPAMIT != 1 & STOPMIRT != 1) ~ 0,
      TRUE ~ NA_real_
    )
  )

# Use the function to create multiple columns
full <- full %>%
  create_treatment_columns("DRYM") %>%
  create_treatment_columns("SWET") %>%
  create_treatment_columns("NAUS") %>%
  create_treatment_columns("VOM") %>%
  create_treatment_columns("DIAR") %>%
  create_treatment_columns("CONS") %>%
  create_treatment_columns("HEAD") %>%
  create_treatment_columns("DIZZ") %>%
  create_treatment_columns("SHAK") %>%
  create_treatment_columns("MUSC") %>%
  create_treatment_columns("DROW") %>%
  create_treatment_columns("WAKE") %>%
  create_treatment_columns("HANX") %>%
  create_treatment_columns("AGIT") %>%
  create_treatment_columns("FATI") %>%
  create_treatment_columns("WTGN") %>%
  create_treatment_columns("WTLS") %>%
  create_treatment_columns("RASH") %>%
  create_treatment_columns("RUN") %>%
  create_treatment_columns("RSEX") %>%
  create_treatment_columns("BLUR") %>%
  create_treatment_columns("SUIT") %>%
  create_treatment_columns("NOSE") %>%
  create_treatment_columns("STOP")

write.csv(full, paste0("/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes/full_data_", duration, "days.csv"), row.names = FALSE)



