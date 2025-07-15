

# -- Load libraries
library(tidyverse)
library(lubridate)
library(tidyverse)
library(data.table)

# -- Sample prescription data frame
set.seed(42) 
num_individuals <- 40
num_prescriptions_per_individual <- 50

# -- Create a data frame to represent prescription data
prescriptions <- data.frame(
  ParticipantID = rep(1:num_individuals, each = num_prescriptions_per_individual),
  DateofSupply = rep(seq(as.Date('2020-01-01'), by = 'months', length.out = num_prescriptions_per_individual), num_individuals),
  ATCCode = c(
    rep(c('N06AB06', 'N06AB03'), each = 100, length.out = num_prescriptions_per_individual),
    rep(c('N06AB04', 'N06AB03'), each = 15, length.out = num_prescriptions_per_individual),
    rep(c('N06AB06', 'N06AB04'), each = 6, length.out = num_prescriptions_per_individual),
    rep(c('N06AX16'), length.out = num_prescriptions_per_individual), 
    rep(c('N06AX11'), length.out = num_prescriptions_per_individual),
    rep(c('N06AA09'), length.out = num_prescriptions_per_individual),
    rep(c('N06AB03'), length.out = num_prescriptions_per_individual), 
    rep(c('N06AB10'), length.out = num_prescriptions_per_individual), 
    rep(c('N06AX21'), length.out = num_prescriptions_per_individual), 
    rep(c('N06AB05'), length.out = num_prescriptions_per_individual) 
  )
)

## IMPORTANT ##
#-- Before proceeding:
#-- 1. Make sure to filter prescription data for only antidepressants and
#      participants with self-report/GP diagnosis of MDD and the 
#      4.5-year prescription window of interest
#-- 2. Select for columns 'ParticipantID', 'DateofSupply' and 'ATCCode'

# -- Prepare Data
PrescriptionData <- prescriptions %>%
  arrange(ParticipantID, ATCCode, DateofSupply) %>%
  group_by(ParticipantID, ATCCode) %>%
  mutate(RowNum = row_number()) %>%
  mutate(TreatmentPeriod = 1,
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

# -- Calculate treatment period duration
PrescriptionEpisodes <- PrescriptionData_PE %>%
  group_by(ParticipantID, ATCCode, PrescriptionEpisode) %>%
  summarize(
    FirstDispense = min(DateofSupply),
    ExpectedEndDate = max(ExpectedEndDate),
    EpisodeDuration = as.numeric(difftime(max(ExpectedEndDate), min(DateofSupply), units = "days")),
    .groups = "drop"
  )

# -- Report number of treatment periods, average treatment period length, and max/min per individual-ATC combination
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

# -- Merge treatment period with all_prescriptions
all_data <- full_join(Summary_PE, all_prescriptions, by = c("ParticipantID", "ATCCode"))

# -- Save treatment groups as a .csv file
write.csv(all_data, "AGDSAcceptabilityTreatmentGroups_14072025.csv", row.names = FALSE, quote = FALSE)

########################################################################################################
