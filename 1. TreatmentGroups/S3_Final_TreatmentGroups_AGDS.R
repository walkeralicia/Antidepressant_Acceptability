
#-- Load r libraries
library(dplyr)
library(data.table)

#-- Set path
wkdir <- "/QRISdata/Q7280/pharmacogenomics"

#-- Lithium prescription data
lithium <- read.table(file.path(wkdir, "data/AGDS_lithium_prescriptions.txt"), header=TRUE) %>%
  mutate(numPrescriptions = PrescriptionDays/30)
x3_lithium <- lithium %>% filter(numPrescriptions > 3)

#-- BIP diagnosis at baseline
baseline <- read.csv("/QRISdata/Q5338/Phenotypes/Combined_AGDS_2020_Freeze_20220620.csv")
baseline_bip <- baseline %>% 
  select(STUDYID, DXBPD) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(-STUDYID, ~ case_when(
    is.na(.) ~ 0,
    . == 1 ~ 1,
    . != 1 ~ 0
  )))

#-- Follow-up data to extract BIP diagnosis
followup_bip <- fread("/QRISdata/Q5338/Phenotypes/followUpData/dataforAlicia", fill = TRUE, header = TRUE) %>%
  select(Q3_2, StudyCode.0) %>%
  rename(BPD = Q3_2,
         STUDYID = StudyCode.0)

#-- Join BIP data frames
both_bip <- full_join(baseline_bip, followup_bip, by = "STUDYID") %>%
  mutate(DXBPD2 = ifelse(DXBPD == 1 | BPD == "Bipolar", 1, 0))
BIP <- both_bip %>% filter(DXBPD2 == 1)

#-- Linkage file for those with genetic data
link <- unique(read.table(file.path(wkdir,"data/StudyID_GenoID_link.txt"), header=T)) %>%
  rename("ParticipantID" = "StudyCodeID")

#-- AGDS participants with genotype data that passed QC
pgs <-read.table("/QRISdata/Q5338/Genotypes/AGDS_nonDup/AGDS_R11_TOPMedr2_imputed_nonDup_chr1.fam") %>%
  select(V2)

#-- AGDS participants of European genetic ancestry
eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id")

#--- AGDS participants with valid PBS data 
pbs <- read.csv(file.path(wkdir, "data/ALL_PBS.csv"))

#-- AGDS participants with at least one antidepressant dispense
filtered_pbs <- pbs %>% filter(grepl("^N06A", ATCCode))
ad_ids <- filtered_pbs$ParticipantID

##########################################################

#-- Vector of duration
durations <- c(360, 600)

for (duration in durations) {
  
  #-- Load treatment groups
  groups <- read.csv(paste0("/scratch/user/uqawal15/TEMP_TreatmentGroups_", duration, "days.csv"), header=TRUE)
  
  #-- Exclude those with self-report BIP from the AD acceptability treatment groups
  groups_noBIP <- groups %>%
    filter(!ParticipantID%in%BIP$STUDYID)
  
  #-- Those with BIP diagnosis and had more than three lithium prescriptions but at least one antidepressant dispense
  BIP_lithium <- full_join(both_bip, x3_lithium, by = c("STUDYID" = "ParticipantID"))
  BIP_lithium <- BIP_lithium %>%
    filter(DXBPD2 == 1 & numPrescriptions > 3) %>%
    select(STUDYID, PrescriptionDays) %>%
    rename("ParticipantID" = "STUDYID") %>%
    filter(ParticipantID %in% groups$ParticipantID)
  
  #-- Format the BIP_lithium dataframe into the same as groups_noBIP
  BIP_lithium_sorted <- BIP_lithium %>%
    mutate(
      NumberOfTreatmentPeriods = NA,
      AverageTreatmentPeriodDays = NA,
      MaxTreatmentPeriodDays = NA,
      MinTreatmentPeriodDays = NA,
      EarliestPrescription = NA,
      LatestPrescription = NA,
      PrescriptionDays = NA,
      TreatmentGroup = "BIP+L",
      TreatmentDrugs = "Lithium",
      Combination_Class = NA
    ) %>%
    select(names(groups_noBIP))
  
  
  #-- Those with a BIP diagnosis but not with at least 4 lithium prescriptions but at least one antidepressant dispense
  bip_ids <- BIP$STUDYID
  bip_wout_L_ids <- bip_ids[!bip_ids %in% BIP_lithium$ParticipantID] %>% 
    unique()
  bip_wout_L_ids <- bip_wout_L_ids[bip_wout_L_ids %in% groups$ParticipantID]
  bip_wout_L <- data.frame(ParticipantID = bip_wout_L_ids,
                           NumberOfTreatmentPeriods = NA,
                           AverageTreatmentPeriodDays = NA,
                           MaxTreatmentPeriodDays = NA,
                           MinTreatmentPeriodDays = NA,
                           EarliestPrescription = NA,
                           LatestPrescription = NA,
                           PrescriptionDays = NA,
                           TreatmentGroup = "BIP-L",
                           TreatmentDrugs = "Miscellaneous",
                           Combination_Class = NA)
  
  
  #-- Final AD treatment group dataframe
  final <- rbind(groups_noBIP, BIP_lithium_sorted, bip_wout_L)
  
  #-- Filter for participants with valid genetic data and of european ancestry
  final_link <- inner_join(link, final,by = "ParticipantID")
  final_link_pgs <- inner_join(final_link, pgs, by = c("IID" = "V2"))
  final_eur <- final_link_pgs %>%
    filter(IID %in% eur$V2)
  
  #-- Map ATC codes to drug names and ad classes
  ATCCodes <- c('N06AB06', 'N06AB10', 'N06AB04', 'N06AX16', 'N06AX21', 
                'N06AA09', 'N06AB05', 'N06AX11', 'N06AX23', 'N06AB03')
  DrugName <- c("SSRI:Sertraline", "SSRI:Escitalopram", "SSRI:Citalopram", "SNRI:Venlafaxine", "SNRI:Duloxetine", 
                "TCA:Amitriptyline", "SSRI:Paroxetine", "TeCA:Mirtazapine", "SNRI:Desvenlafaxine", "SSRI:Fluoxetine")
  DrugClass <- c('SSRI', 'SSRI', 'SSRI', 'SNRI', 'SNRI',
                 'TCA', 'SSRI', 'TeCA', 'SNRI', 'SSRI')
  atc_drug_ref = data.frame(DrugName = DrugName, ATCCodes = ATCCodes, DrugClass = DrugClass)
  
  final_mapped <- left_join(final_eur, atc_drug_ref, by = c("TreatmentGroup" = "ATCCodes"))
  
  
  #-- Final Drug and Drug Class groups
  final_mapped <- final_mapped %>%
    mutate(DrugName = ifelse(is.na(DrugName), TreatmentGroup, DrugName),
           DrugClass = ifelse(!is.na(Combination_Class) & Combination_Class != "Miscellaneous", Combination_Class,
                              ifelse(is.na(Combination_Class) & !is.na(DrugClass), DrugClass, 
                                     ifelse(TreatmentGroup == "BIP+L", "BIP+L", 
                                            ifelse(TreatmentGroup == "BIP-L", "BIP-L", "Various")))))
  
  #-- Save results
  write.csv(final_mapped, paste0("/scratch/user/uqawal15/Final_Treatmentgroups_", duration, "days.csv"), row.names=FALSE)

}

