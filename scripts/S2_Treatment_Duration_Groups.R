
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

process_treatment_groups <- function(duration) {
  #-- Load treatment files and convert to Date format
  groups <- fread(file.path(wkdir, "data", "AGDSAcceptabilityTreatmentGroups_14072025.csv")) %>%
    mutate(
      EarliestPrescription = as.Date(EarliestPrescription, format = "%d/%m/%Y"),
      LatestPrescription = as.Date(LatestPrescription, format = "%d/%m/%Y")
    )
  
  #-- Split participants based on number of unique antidepressants dispensed
  multiple_ads <- groups %>% 
    group_by(ParticipantID) %>%
    filter(n() > 1) %>% 
    ungroup() %>%
    as.data.frame()
  
  one_ad <- groups[!groups$ParticipantID %in% multiple_ads$ParticipantID, ]
  
  #-- Flag participants dispensing at least two antidepressants for a sustained period of time
  flagged <- multiple_ads %>%
    group_by(ParticipantID) %>%
    summarise(
      Flag = sum(PrescriptionDays >= duration) >= 2
    ) %>%
    ungroup()
  
  #-- Annotate Participants as "Combination" or other based on flags
  multiple_ads <- multiple_ads %>%
    left_join(flagged, by = "ParticipantID") %>%
    mutate(
      TreatmentGroup = ifelse(Flag, "Combination", NA)
    )
  
  #-- Process combination treatments (Concatenate drugs that were each dispensed for a sustained period of time)
  combination_participants <- multiple_ads %>%
    filter(TreatmentGroup == "Combination") %>%
    group_by(ParticipantID) %>%
    mutate(
      ValidDrugs = PrescriptionDays >= duration,
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
    filter(TreatmentGroup == "Combination" & ValidDrugs) %>%
    group_by(ParticipantID) %>%
    summarise(
      NumberOfPrescriptionEpisodes = sum(NumberOfPrescriptionEpisodes),
      AveragePrescriptionEpisodeDays = mean(AveragePrescriptionEpisodeDays),
      EarliestPrescription = min(EarliestPrescription),
      LatestPrescription = max(LatestPrescription),
      PrescriptionDays = sum(PrescriptionDays),
      TreatmentGroup = first(TreatmentGroup),
      TreatmentDrugs = first(TreatmentDrugs)
    ) %>%
    as.data.frame()

  
  #-- Handle non-combination participants
  
  #-- Participants taking multiple antidepressants but only one for a sustained period of time
  one_of_many_ad_sustained <- multiple_ads %>%
    filter(is.na(TreatmentGroup) & PrescriptionDays >= duration) %>%
    mutate(TreatmentGroup = ATCCode,
           TreatmentDrugs = ATCCode) %>%
    select(-ATCCode, -Flag)
  
  #-- Participants taking multiple antidepressant with none for a sustained period of time
  unallocated_ids <- multiple_ads %>%
    filter(is.na(TreatmentGroup) & PrescriptionDays < duration) %>%
    pull(ParticipantID) %>%
    unique()
  unallocated_ids <- setdiff(unallocated_ids, one_of_many_ad_sustained$ParticipantID)
  
  none_sustained <- multiple_ads %>%
    filter(ParticipantID %in% unallocated_ids) %>%
    group_by(ParticipantID) %>%
    summarise(
      NumberOfPrescriptionEpisodes = sum(NumberOfPrescriptionEpisodes),
      AveragePrescriptionEpisodeDays = mean(AveragePrescriptionEpisodeDays),
      EarliestPrescription = min(EarliestPrescription),
      LatestPrescription = max(LatestPrescription),
      PrescriptionDays = sum(PrescriptionDays),
      TreatmentGroup = "Various",
      TreatmentDrugs = "Miscellaneous"
    ) %>%
    as.data.frame()
  
  #-- Process participants with only a single antidepressant dispensed
  one_ad_sustained <- one_ad %>%
    mutate(
      TreatmentGroup = ifelse(PrescriptionDays >= duration, ATCCode, "Various"),
      TreatmentDrugs = ifelse(PrescriptionDays >= duration, ATCCode, "Miscellaneous")
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
    mutate(Combination_Class = sapply(TreatmentDrugs, get_atc_class_combo))
  
  return(all_groups)
}

# Run for both duration thresholds
results_360 <- process_treatment_groups(360)
results_600 <- process_treatment_groups(600)


#--Write results
write.csv(results_360, "/scratch/user/uqawal15/TEMP_TreatmentGroups_360days.csv",  row.names = FALSE, quote = FALSE)
write.csv(results_600, "/scratch/user/uqawal15/TEMP_TreatmentGroups_600days.csv",  row.names = FALSE, quote = FALSE)

