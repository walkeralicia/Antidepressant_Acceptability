

# -- Load libraries
library(tidyverse)
library(lubridate)
library(tidyverse)
library(data.table)
library(stringr)

prescriptions <- read.csv("C://Users//walkera//OneDrive - Nexus365//Documents//PhD//AGDS//Antidepressant_Acceptability//Prescription_Data//Fake_Prescription_Data.csv")


#-- Prepare Data
#-- Make sure the prescriptions here are just for the top 10 N06A and within the 4.5 year window of interest
PrescriptionData <- prescriptions %>%
  filter(ATCCode %in% c("N06AA09", "N06AB04", "N06AB05", "N06AX16", "N06AX11",
                        "N06AB10", "N06AB06","N06AB03", "N06AX21", "N06AX23")) %>%
  arrange(ParticipantID, ATCCode, DateofSupply) %>%
  group_by(ParticipantID, ATCCode) %>%
  mutate(RowNum = row_number(),
         DateofSupply = as.Date(DateofSupply),
         PrescriptionDuration = 0,
         ExpectedEndDate = as.Date("2026-01-01")) %>%
  ungroup() %>%
  data.table()


# -- Calculate Prescription Duration and Expected End Date
PrescriptionData_Durations <- PrescriptionData %>%
  group_by(ParticipantID, ATCCode) %>%
  mutate(
    PrescriptionDuration = if_else(
      row_number() == n(),
      30,
      pmin(
        as.numeric(difftime(lead(DateofSupply, default = last(DateofSupply)), DateofSupply, units = "days")),
        30
      )
    ),
    PrescriptionDuration = replace_na(PrescriptionDuration, 0),
    ExpectedEndDate = DateofSupply + PrescriptionDuration
  ) %>%
  ungroup()

# -- Assign Prescription Episodes with 90-day gap condition
PrescriptionData_PE <- PrescriptionData_Durations %>%
  group_by(ParticipantID, ATCCode) %>%
  mutate(
    PrescriptionEpisode = cumsum(
      ifelse(
        as.numeric(difftime(DateofSupply, lag(ExpectedEndDate, default = first(ExpectedEndDate)), units = "days")) > 90,
        1, 
        0
      )
    ) + 1
  ) %>%
  ungroup()

#-- Calculate prescription episode duration
PrescriptionEpisodes <- PrescriptionData_PE %>%
  group_by(ParticipantID, ATCCode, PrescriptionEpisode) %>%
  summarize(
    FirstDispense = min(DateofSupply),
    ExpectedEndDate = max(ExpectedEndDate),
    EpisodeDuration = as.numeric(difftime(max(ExpectedEndDate), min(DateofSupply), units = "days")),
    .groups = "drop"
  )

#-- Classify a participant as adherent to their medication if at least one prescription episode is >= 360 days or >= 180 days long
Medication_Adherence <- PrescriptionEpisodes %>%
  group_by(ParticipantID, ATCCode) %>%
  mutate(
    EpisodeAdherence360 = case_when(
      EpisodeDuration >= 360 ~ 1,
      EpisodeDuration < 360 ~ 0,
      TRUE ~ NA
    ),
    EpisodeAdherence180 = case_when(
      EpisodeDuration >= 180 ~ 1,
      EpisodeDuration < 180 ~ 0,
      TRUE ~ NA
    )
  ) %>%
  summarise(
    SumAdherent360 = sum(EpisodeAdherence360),
    SumAdherent180 = sum(EpisodeAdherence180)
  ) %>%
  mutate(
    Adherent360 = case_when(
      SumAdherent360 >= 1 ~ 1,
      SumAdherent360 < 1 ~ 0,
      TRUE ~ NA
    ),
    Adherent180 = case_when(
      SumAdherent180 >= 1 ~ 1,
      SumAdherent180 < 1 ~ 0,
      TRUE ~ NA
    )
  ) %>%
  ungroup() %>%
  select(-SumAdherent360, -SumAdherent180)

#-- Calculate a general adherence score across all antidepressant medications
Adherent_Prescriptions <- PrescriptionData_Durations %>%
  arrange(ParticipantID, DateofSupply) %>%
  group_by(ParticipantID) %>%
  mutate(
    PrevExpectedEndDate = lag(ExpectedEndDate),
    WithinWindow = ifelse(
      !is.na(PrevExpectedEndDate) & 
        DateofSupply <= (PrevExpectedEndDate + 90), 
      1, 0
    )
  ) %>%
  select(-PrevExpectedEndDate) %>%
  ungroup()

General_Adherence <- Adherent_Prescriptions %>%
  filter(RowNum != 1) %>%
  group_by(ParticipantID) %>%
  mutate(
    GA_score = sum(WithinWindow)/n()
  ) %>%
  select(ParticipantID, GA_score) %>%
  distinct()


# -- Report number of prescription episodes, average prescription episode length, and max/min per individual-ATC combination
Summary_PE <- PrescriptionEpisodes %>%
  group_by(ParticipantID, ATCCode) %>%
  summarise(
    NumberOfPrescriptionEpisodes = max(PrescriptionEpisode),
    AveragePrescriptionEpisodesDays = mean(EpisodeDuration, na.rm = TRUE)) %>%
  arrange(desc(NumberOfPrescriptionEpisodes))

# -- Define a function to estimate individual-level treatment duration for a specific ATC code
process_duration <- function(data, atc_code) {
  filtered_data <- data %>%
    filter(ATCCode == atc_code) %>%
    group_by(ParticipantID) %>%
    summarise(
      EarliestPrescription = if (n() > 0) min(DateofSupply, na.rm = TRUE) else as.Date(NA),
      LatestPrescription = if (n() > 0) max(DateofSupply, na.rm = TRUE) else as.Date(NA),
      PrescriptionDays = sum(PrescriptionDuration)
    )
  return(filtered_data)
}

# -- Calculate overall treatment duration
prescription_tables <- list()
atc_codes <- unique(PrescriptionData$ATCCode)

for (atc_code in atc_codes) {
  prescription_tables[[atc_code]] <- process_duration(PrescriptionData_PE, atc_code)
}

prescription_tables <- Map(function(df, atc_code) {
  df$ATCCode <- atc_code
  return(df)
}, prescription_tables, names(prescription_tables))


all_prescriptions <- do.call(rbind, prescription_tables)

# -- Merge prescription episodes with all_prescriptions, and adherence measures
all_data <- full_join(Summary_PE, all_prescriptions, by = c("ParticipantID", "ATCCode")) %>%
  left_join(Medication_Adherence, by = c("ParticipantID", "ATCCode")) %>%
  left_join(General_Adherence, by = "ParticipantID") %>%
  mutate(
    AD_Window = as.numeric((LatestPrescription+30) - EarliestPrescription)
  ) 

#-- Define prescription dispensing windows
window <- all_data %>%
  group_by(ParticipantID) %>%
  summarise(
    All_AD_FirstDispense = min(EarliestPrescription),
    All_AD_LastDispense = max(LatestPrescription)
  ) %>%
  ungroup() %>%
  mutate(
    All_AD_Window = as.numeric((All_AD_LastDispense+30) - All_AD_FirstDispense)
  ) %>%
  select(ParticipantID, All_AD_Window, All_AD_FirstDispense, All_AD_LastDispense)

all_data <- all_data %>%
  left_join(window, by = c("ParticipantID"))


#-- Filter for Anti psychotics
AntipsychoticData <- prescriptions %>%
  filter(str_detect(ATCCode, "N05A") & ATCCode != "N05AN01") %>%
  arrange(ParticipantID, ATCCode, DateofSupply) %>%
  data.table()

#-- Check if any Anti psychotic Medication is dispensed within DateOfSupply and ExpectedEndDate of an Antidepressant prescription
antidep_antipsych <- PrescriptionData_Durations %>%
  left_join(AntipsychoticData, 
            by = c("ParticipantID"),
            relationship = "many-to-many",
            suffix = c("_antidep", "_antipsych")) %>%
  rename(ExpectedEndDate_antidep = ExpectedEndDate) %>%
  mutate(Overlap =
           case_when(
             DateofSupply_antipsych >= DateofSupply_antidep & DateofSupply_antipsych <= ExpectedEndDate_antidep ~ 1,
             TRUE ~ 0
           ),
         Overlap_buffer1 = 
           case_when(
             DateofSupply_antipsych >= DateofSupply_antidep - 7 & DateofSupply_antipsych <= ExpectedEndDate_antidep ~ 1,
             TRUE ~ 0
           ),
         Overlap_buffer2 = 
           case_when(
             DateofSupply_antipsych >= DateofSupply_antidep - 14 & DateofSupply_antipsych <= ExpectedEndDate_antidep ~ 1,
             TRUE ~ 0
           ))



antidep_antipsych_augmentation_AD <- antidep_antipsych %>%
  group_by(ParticipantID, ATCCode_antidep) %>%
  summarise(
    Augmentation_ADspecific = ifelse(sum(Overlap, na.rm = TRUE) >= 3, 1, 0),
    Augmentation_ADspecific_buffer1 = ifelse(sum(Overlap_buffer1, na.rm = TRUE) >= 3, 1, 0),
    Augmentation_ADspecific_buffer2 = ifelse(sum(Overlap_buffer2, na.rm = TRUE) >= 3, 1, 0)
  )

antidep_antipsych_augmentation_ALL <- antidep_antipsych %>%
  group_by(ParticipantID) %>%
  summarise(
    Augmentation = ifelse(sum(Overlap, na.rm = TRUE) >= 3, 1, 0),
    Augmentation_buffer1 = ifelse(sum(Overlap_buffer1, na.rm = TRUE) >= 3, 1, 0),
    Augmentation_buffer2 = ifelse(sum(Overlap_buffer2, na.rm = TRUE) >= 3, 1, 0)
  )

#-- Add in augmentation classifications
final <- all_data %>% 
  left_join(antidep_antipsych_augmentation_AD, by = c("ParticipantID", "ATCCode" = "ATCCode_antidep")) %>%
  left_join(antidep_antipsych_augmentation_ALL, by = c("ParticipantID"))


#--- Check for concomitant medication proxy for comorbidities
proxy_drugs <- prescriptions %>%
  filter(str_starts(ATCCode, paste(c("A10A", 
                                     "A10B",
                                     "L01",
                                     "N06D",
                                     "R03",
                                     "C10",
                                     "N05B",
                                     "N07BB",
                                     "J05AP",
                                     "L03AB",
                                     "A05",
                                     "B01A",
                                     "C01DA",
                                     "C02",
                                     "C03",
                                     "C07",
                                     "C08",
                                     "C09",
                                     "N06BA",
                                     "N03A",
                                     "N05C",
                                     "A08AA",
                                     "A08AB",
                                     "A10BJ",
                                     "N02",
                                     "L04A",
                                     "H03AA"
                                     ), collapse = "|")) 
         | ATCCode %in% c("C01AA05", "C02AC01", "C02AC02",
                          "N05CH01", "R06AD01", "R06AD02",
                          "A10BA02")) %>%
  arrange(ParticipantID, ATCCode, DateofSupply) %>%
  mutate(DateofSupply = as.Date(DateofSupply)) %>%
  data.table()

proxy_drugs_grouped <- proxy_drugs %>%
  mutate(
    Condition = case_when(
      str_detect(ATCCode, "^A10A") | str_detect(ATCCode, "^A10B") ~ "Diabetes",
      str_detect(ATCCode, "^L01") ~ "Cancer",
      str_detect(ATCCode, "^N06D") ~ "Dementia",
      str_detect(ATCCode, "^R03") ~ "Asthma_and_COPD",
      str_detect(ATCCode, "^C10") ~ "Dyslipidemia",
      str_detect(ATCCode, "^N05B") ~ "Anxiety",
      str_detect(ATCCode, "^N07BB") | str_detect(ATCCode, "^J05AP") |
        str_detect(ATCCode, "^L03AB") | str_detect(ATCCode, "^A05") ~ "Hepatic",
      str_detect(ATCCode, "^B01A") | str_detect(ATCCode, "^C01DA") | 
      ATCCode == "C01AA05" | str_detect(ATCCode, "^C02") | str_detect(ATCCode, "^C03") |
        str_detect(ATCCode, "^C07") | str_detect(ATCCode, "^C08") | str_detect(ATCCode, "^C09") ~ "Cardiovascular",
      str_detect(ATCCode, "^N06BA") | ATCCode == "C02AC01" | ATCCode == "C02AC02" ~ "ADHD",
      str_detect(ATCCode, "^N03A") ~ "Siezures",
      str_detect(ATCCode, "^N05C") | ATCCode == "N05CH01" | ATCCode == "R06AD01" |
        ATCCode == "R06AD02" ~ "Sleep",
      str_detect(ATCCode, "^A08AA") | str_detect(ATCCode, "^A08AB") | 
        str_detect(ATCCode, "^A10BJ") | ATCCode == "A10BA02" ~ "Weight_loss",
      str_detect(ATCCode, "^N02") ~ "Analgesics",
      str_detect(ATCCode, "^L04A") ~ "Immunosuppressants",
      str_detect(ATCCode, "^H03AA") ~ "Thyroid"
    )
  )

co_conditions <- final %>%
  select(ParticipantID, All_AD_FirstDispense, All_AD_LastDispense) %>%
  left_join(proxy_drugs_grouped, by = "ParticipantID", relationship = "many-to-many",
            suffix = c("_antidep", "_coconditions")) %>%
  filter((DateofSupply >= All_AD_FirstDispense - 90) & (DateofSupply < All_AD_LastDispense + 30)) %>%
  group_by(ParticipantID, Condition) %>%
  distinct() %>%
  summarise(
    num_prescriptions = n()
  ) %>%
  filter(
    num_prescriptions >= 3
  ) %>%
  ungroup()

num_conditions <- co_conditions %>%
  group_by(ParticipantID) %>%
  summarise(
    num_conditions = n_distinct(Condition)
  )
each_com <- co_conditions %>%
  setDT() %>%
  dcast(ParticipantID ~ Condition,
        value.var = "num_prescriptions",
        fill = 0,
        fun.aggregate = function(x) as.numeric(any(x > 0)))
  

original_cols <- names(final)
final_final <- final %>%
  left_join(num_conditions, by = "ParticipantID") %>%
  left_join(each_com, by = "ParticipantID")
new_cols <- setdiff(names(final_final), original_cols)
output <- final_final %>%
  mutate(across(all_of(new_cols), ~replace_na(.x, 0)))
  


# -- Save treatment groups as a .csv file
write.csv(output, "AGDSAcceptabilityTreatmentGroups_25082025.csv", row.names = FALSE, quote = FALSE)

