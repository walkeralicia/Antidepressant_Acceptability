
#-- Load R libraries
library(dplyr)
library(data.table)

#-- Set path
wkdir <- "/QRISdata/Q7280/pharmacogenomics"

#-- Drug reference table
source("/QRISdata/Q7280/pharmacogenomics/Drug_Reference/Drug_Reference_Table.R")
atc_mapping <- setNames(as.list(drug_ref$DrugClass), drug_ref$ATCCode)


#-- Function to assign ATC class combination
get_atc_class_combo <- function(drugs) {
  if (drugs == "Miscellaneous" || !grepl("_", drugs)) return(NA)
  
  atc_codes <- strsplit(drugs, "_")[[1]]
  
  # Map each ATC code to its class
  classes <- sapply(atc_codes, function(code) {
    if (code %in% names(atc_mapping)) {
      return(atc_mapping[[code]])
    } else {
      return("Unknown")
    }
  })
  
  # Count occurrences of each class
  class_counts <- table(classes)
  
  # Check if any class has at least 2 occurrences
  for (class_name in drug_ref$DrugClass) {
    if (!is.na(class_counts[class_name]) && class_counts[class_name] >= 2) {
      return(class_name)
    }
  }
  
  # If no class has 2+ occurrences, return Miscellaneous
  return("Miscellaneous")
}

#-- Load treatment files and convert to Date format
groups_raw <- fread(file.path(wkdir, "data", "AGDSAcceptabilityTreatmentGroups_25082025.csv")) %>%
  mutate(
    EarliestPrescription = as.Date(EarliestPrescription, format = "%d/%m/%Y"),
    LatestPrescription = as.Date(LatestPrescription, format = "%d/%m/%Y")
  )

groups <- groups_raw %>%
  select(ParticipantID, ATCCode, PrescriptionDays, Adherent360)

groups_char <- groups_raw %>%
  select(ParticipantID, ATCCode, NumberOfPrescriptionEpisodes, AveragePrescriptionEpisodesDays, 
         PrescriptionDays, GA_score, AD_Window, All_AD_Window, 
         Augmentation_ADspecific, Augmentation, num_conditions, 
         ADHD, Analgesics, Anxiety, Asthma_and_COPD, Cancer, Cardiovascular, Dementia, Diabetes, 
         Dyslipidemia, Hepatic, Immunosuppressants, Siezures, Sleep, Thyroid)

#-- Split participants based on number of unique antidepressants dispensed
multiple_ads <- groups %>% 
  group_by(ParticipantID) %>%
  filter(n() > 1) %>% 
  ungroup() %>%
  as.data.frame()

#-- Participants with just one AD
one_ad <- groups[!groups$ParticipantID %in% multiple_ads$ParticipantID, ]


process_treatment_groups <- function(duration) {
  
  #-- Flag participants dispensing at least two antidepressants for a sustained period of time with 
  #-- both periods surpassing high adherence
  flagged <- multiple_ads %>%
    group_by(ParticipantID) %>%
    summarise(
      Flag = (sum(PrescriptionDays >= duration & Adherent360 == 1) >= 2) 
    ) %>%
    ungroup()
  
  #-- Annotate Participants as "Combination" or other based on flags
  multiple_ads_combo <- multiple_ads %>%
    left_join(flagged, by = "ParticipantID") %>%
    mutate(
      TreatmentGroup = ifelse(Flag, "Combination", NA)
    )
  
  #-- Process combination treatments (Concatenate drugs that were each dispensed for a sustained period of time
  #-- and were adherent)
  combination_participants <- multiple_ads_combo %>%
    filter(TreatmentGroup == "Combination") %>%
    group_by(ParticipantID) %>%
    mutate(
      ValidDrugs = PrescriptionDays >= duration & Adherent360 == 1,
      TreatmentDrugs = ifelse(
        sum(ValidDrugs) >= 2, 
        paste(ATCCode[ValidDrugs], collapse = "_"),
        NA
      )
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  #-- Summarise participants with sustained combination treatment
  sustained_combination <- combination_participants %>%
    filter(ValidDrugs) %>%
    group_by(ParticipantID) %>%
    summarise(
      PrescriptionDays = sum(PrescriptionDays),
      TreatmentGroup = first(TreatmentGroup),
      TreatmentDrugs = first(TreatmentDrugs)
    ) %>%
    as.data.frame()
  
  
  #-- Handle non-combination participants
  
  #-- Participants taking multiple antidepressants but only one for a sustained period of time with high adherence
  one_of_many_ad_sustained <- multiple_ads_combo %>%
    filter(is.na(TreatmentGroup) & PrescriptionDays >= duration & Adherent360 == 1) %>%
    mutate(TreatmentGroup = ATCCode,
           TreatmentDrugs = ATCCode) %>%
    select(-ATCCode, -Flag)
  
  #-- Participants taking multiple antidepressant with none for a sustained period of time
  unallocated_ids <- multiple_ads_combo %>%
    filter(is.na(TreatmentGroup) & (PrescriptionDays < duration | Adherent360 == 0)) %>%
    pull(ParticipantID) %>%
    unique()
  unallocated_ids <- setdiff(unallocated_ids, one_of_many_ad_sustained$ParticipantID)
  
  none_sustained <- multiple_ads_combo %>%
    filter(ParticipantID %in% unallocated_ids) %>%
    group_by(ParticipantID) %>%
    summarise(
      PrescriptionDays = sum(PrescriptionDays),
      TreatmentGroup = "Various",
      TreatmentDrugs = "Miscellaneous"
    ) %>%
    as.data.frame()
  
  #-- Process participants with only a exclusive single antidepressant dispensed
  one_ad_sustained <- one_ad %>%
    mutate(
      TreatmentGroup = ifelse(PrescriptionDays >= duration & Adherent360 == 1, ATCCode, "Various"),
      TreatmentDrugs = ifelse(PrescriptionDays >= duration & Adherent360 == 1, ATCCode, "Miscellaneous")
    ) %>%
    select(-ATCCode) %>%
    as.data.frame()
  
  #-- Combine all antidepressant treatment groups
  all_groups <- bind_rows(
    sustained_combination,
    one_of_many_ad_sustained,
    one_ad_sustained,
    none_sustained
  )
  
  #-- Add ATC class combinations
  all_groups <- all_groups %>%
    select(-Adherent360) %>%
    mutate(Combination_Class = sapply(TreatmentDrugs, get_atc_class_combo))
  
  #-- Add prescription characteristics for those within a sustained AD group
  sustained_char <- all_groups %>%
    select(-PrescriptionDays) %>%
    inner_join(groups_char, by = c("ParticipantID", "TreatmentGroup" = "ATCCode"))
  
  #-- Add characteristics for those within the combination or various AD group
  unsustained_char <- all_groups %>%
    filter(TreatmentGroup == "Combination" | TreatmentGroup == "Various")%>%
    inner_join(groups_char, by = c("ParticipantID")) %>%
    select(-PrescriptionDays.y, -AD_Window, -ATCCode, -Augmentation_ADspecific,
           -NumberOfPrescriptionEpisodes, -AveragePrescriptionEpisodesDays) %>%
    rename(PrescriptionDays = PrescriptionDays.x) %>%
    distinct()
  
  all_char <- bind_rows(sustained_char, unsustained_char)
  
  
  return(all_char)
}

# Run for both duration thresholds
results_360 <- process_treatment_groups(360)
results_600 <- process_treatment_groups(600)


#--Write results
write.csv(results_360, "/scratch/user/uqawal15/TEMP_TreatmentGroups_360days_Adherent.csv",  row.names = FALSE, quote = FALSE)
write.csv(results_600, "/scratch/user/uqawal15/TEMP_TreatmentGroups_600days_Adherent.csv",  row.names = FALSE, quote = FALSE)

