

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

#====== Plot of antidepressants ranked by prescribing popularity each year ===

# Make sure to filter the 4.5-year period of prescription data for only antidepressants and
# participants with self-report/GP diagnosis of MDD

# -- Rank prescriptions by total dispenses by year
rank_prescriptions <- prescriptions %>%
  select(ParticipantID, ATCCode, DateofSupply) %>%
  group_by(year = format(as.Date(DateofSupply), "%Y"), ATCCode) %>%
  summarise(
    TotalPrescriptions = n()
  )  %>%
  group_by(year) %>%
  mutate(rank = as.integer(min_rank(-TotalPrescriptions))) %>%
  arrange(year, rank)

#-- Map ATC codes to drug names and classes
ATCCodes <- c('N06AB06', 'N06AB10', 'N06AB04', 'N06AX16', 'N06AX21', 
              'N06AA09', 'N06AB05', 'N06AX11', 'N06AX23', 'N06AB03')
DrugName <- c("Sertraline", "Escitalopram", "Citalopram", "Venlafaxine", "Duloxetine", 
              "Amitriptyline", "Paroxetine", "Mirtazapine", "Desvenlafaxine", "Fluoxetine")
atc_drug_ref = data.frame(DrugName = DrugName, ATCCodes = ATCCodes)
rank_prescriptions_mapped <- left_join(rank_prescriptions, atc_drug_ref, by = c("ATCCode" = "ATCCodes"))

#-- Plot code
rank_plot <- rank_prescriptions_mapped %>%
  ggplot(aes(x = year, y = rank, color = DrugName, group = DrugName)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(x = "Prescribing Year",
       y = "Rank",
       fill = "Antidepressant",
       title = "Antidepressant prescribing trends from 2020 to 2024") +
  scale_y_reverse(
    limits = c(max(rank_prescriptions_mapped$rank), 1),
    breaks = seq(1, max(rank_prescriptions_mapped$rank), 1)
  ) +
  scale_x_discrete(expand = c(0, 0.2)) + 
  theme_minimal(base_size = 15) +
  theme(axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.line = element_line(colour = 'black'),
        legend.position = "top",
        legend.title = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2))

#============= Begin Antidepressant Extraction =====================

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

# -- Assign Treatment Periods with 90-day gap condition
PrescriptionData_TP <- PrescriptionData_Durations %>%
  group_by(ParticipantID, ATCCode) %>%
  mutate(
    TreatmentPeriod = cumsum(
      ifelse(
        as.numeric(difftime(DateofSupply, lag(ExpectedEndDate, default = first(ExpectedEndDate)), units = "days")) > 90,
        1, 
        0
      )
    ) + 1
  ) %>%
  ungroup()

# -- Calculate treatment period duration
TreatmentPeriods <- PrescriptionData_TP %>%
  group_by(ParticipantID, ATCCode, TreatmentPeriod) %>%
  summarize(
    FirstDispense = min(DateofSupply),
    ExpectedEndDate = max(ExpectedEndDate),
    PeriodDuration = as.numeric(difftime(max(ExpectedEndDate), min(DateofSupply), units = "days")),
    .groups = "drop"
  )

# -- Report number of treatment periods, average treatment period length, and max/min per individual-ATC combination
Summary_TP <- TreatmentPeriods %>%
  group_by(ParticipantID, ATCCode) %>%
  summarise(
    NumberOfTreatmentPeriods = max(TreatmentPeriod),
    AverageTreatmentPeriodDays = mean(PeriodDuration, na.rm = TRUE),
    MaxTreatmentPeriodDays = max(PeriodDuration, na.rm = TRUE),
    MinTreatmentPeriodDays = min(PeriodDuration, na.rm = TRUE)) %>%
  arrange(desc(NumberOfTreatmentPeriods))

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
  prescription_tables[[atc_code]] <- process_duration(PrescriptionData_TP, atc_code)
}

prescription_tables <- Map(function(df, atc_code) {
  df$ATCCode <- atc_code
  return(df)
}, prescription_tables, names(prescription_tables))


all_prescriptions <- do.call(rbind, prescription_tables)

# -- Merge treatment period with all_prescriptions
all_data <- full_join(Summary_TP, all_prescriptions, by = c("ParticipantID", "ATCCode"))

# -- Save treatment groups as a .csv file
write.csv(all_data, "AcceptabilityTreatmentGroups.csv", row.names = FALSE, quote = FALSE)

########################################################################################################
