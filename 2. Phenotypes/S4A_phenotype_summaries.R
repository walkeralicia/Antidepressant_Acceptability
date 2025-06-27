
#-- Load R libraries
library(dplyr)
library(tidyr)
library(openxlsx)

#-- Vector of cumulative dispense thresholds
durations <- c(360, 600)

#-- Order of treatment groups
custom_order <- c('SSRI:Sertraline', 'SSRI:Escitalopram', 'SSRI:Citalopram', 'SSRI:Fluoxetine', 'SSRI:Paroxetine',
                  'SNRI:Desvenlafaxine', 'SNRI:Venlafaxine', 'SNRI:Duloxetine', 'TCA:Amitriptyline', 'TeCA:Mirtazapine',
                  'Combination', 'BIP+L', 'BIP-L', 'Various')


#============== Summaries ===============================

tables_list <- list()

for (duration in durations){
  
  #-- Load data
  dat <- read.csv(paste0("/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes/inner_data_", duration, "days.csv"))
  
  #=== Basic characteristics of MDD and prescriptions patterns ===
  # use non-missing response denominator 
  table1 <- dat %>%
    group_by(DrugName) %>%
    summarise(
      group_size = n(),
      perc_MDD = sum(MDD, na.rm = TRUE) / sum(!is.na(MDD)) * 100,
      median_MDD_score = median(CUM_SYM, na.rm = TRUE),
      sd_MDD_score = sd(CUM_SYM, na.rm = TRUE),
      mean_age = mean(AGE, na.rm = TRUE),
      sd_age = sd(AGE, na.rm = TRUE),
      perc_males = mean(SEX == 'Male', na.rm = TRUE) * 100,
      mean_age_of_MDDonset = mean(AGE2WKF, na.rm = TRUE), 
      sd_age_of_MDDonset = sd(AGE2WKF, na.rm = TRUE),
      median_MDD_eps = median(TIMES2WK, na.rm = TRUE),
      sd_MDD_eps = sd(TIMES2WK, na.rm = TRUE),
      perc_suicide = sum(SUICIDEA == 1, na.rm = TRUE) / sum(!is.na(SUICIDEA)) * 100,
      perc_selfharm = sum(SELFHARM == 1, na.rm = TRUE) / sum(!is.na(SELFHARM)) * 100,
      mean_bmi = mean(BMI, na.rm = TRUE),
      sd_bmi = sd(BMI, na.rm = TRUE),
      perc_T2D = sum(TYPE11B == 1, na.rm = TRUE) / sum(!is.na(TYPE11B)) * 100,
      perc_stomach_ulcers = sum(TYPE33C == 1, na.rm = TRUE) / sum(!is.na(TYPE33C)) * 100,
      medial_total_psych_comorbidities = median(TOTAL_COM, na.rm = TRUE),
      sd_total_psych_comorbidities = sd(TOTAL_COM, na.rm = TRUE),
      perc_backpain = sum(DXPHYS3, na.rm = TRUE) / sum(!is.na(DXPHYS3)) * 100,
      perc_chronicfatigue = sum(DXPHYS6, na.rm = TRUE) / sum(!is.na(DXPHYS6)) * 100,
      perc_epilepsy = sum(DXPHYS12, na.rm = TRUE) / sum(!is.na(DXPHYS12)) *100,
      perc_chronicpain = sum(DXPHYS34, na.rm = TRUE) / sum(!is.na(DXPHYS34)) * 100,
      perc_migraines = sum(MIGEVR, na.rm = TRUE) / sum(!is.na(MIGEVR)) * 100,
      perc_family_mentalhealthdisorder = sum(FAMMHD, na.rm = TRUE) / sum(!is.na(FAMMHD)) * 100,
      perc_family_MDD = sum(FAMMDD, na.rm = TRUE) / sum(!is.na(FAMMDD)) * 100,
      perc_family_BPD = sum(FAMBPD, na.rm = TRUE) / sum(!is.na(FAMBPD)) * 100,
      perc_FIBRO = sum(DXFIBRO, na.rm = TRUE)/sum(!is.na(DXFIBRO)) * 100,
      perc_ENDO = sum(DXENDO, na.rm = TRUE)/sum(!is.na(DXENDO)) * 100,
      perc_PCOS = sum(DXPCOS, na.rm = TRUE)/sum(!is.na(DXPCOS)) * 100,
      median_education_level = median(EDU, na.rm = TRUE),
      sd_education_level = sd(EDU, na.rm = TRUE),
      median_physical_health_level = median(PHYSHLTH, na.rm = TRUE),
      sd_physical_health_level = mean(PHYSHLTH, na.rm = TRUE),
      perc_regular_smoker = sum(REGSMK, na.rm = TRUE) / sum(!is.na(REGSMK)) * 100,
      mean_drink_3months = mean(DRK3FRQ, na.rm = TRUE),
      sd_drink_3months = sd(DRK3FRQ, na.rm = TRUE),
      perc_WELL = sum(WELLAD == 1, na.rm = TRUE) / sum(!is.na(WELLAD)) * 100,
      perc_WELL_any = sum(WELLAD_any == 1, na.rm = TRUE) / sum(!is.na(WELLAD_any)) * 100,
      perc_STOP = sum(STOPAD == 1, na.rm = TRUE) / sum(!is.na(STOPAD)) * 100,
      perc_STOP_any = sum(STOPAD_any == 1, na.rm = TRUE)/ sum(!is.na(STOPAD_any)) * 100
    )  %>%
    mutate(DrugName = factor(DrugName, levels = custom_order)) %>%
    arrange(DrugName) %>%
    as.data.frame()  %>%
    mutate_at(vars(-DrugName, -group_size), ~ round(., 1))
  
  transformed_table1 = setNames(data.frame(t(table1[,-1])), table1[,1])
  
  #=== Worst episode characteristics: Median time, cumulative_symptoms ===
  table2a <- dat %>%
    group_by(DrugName) %>%
    summarise(
      across(
        c(TIMEWKS),
        list(median = ~ median(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      )
    )  %>%
    mutate(DrugName = factor(DrugName, levels = custom_order)) %>%
    arrange(DrugName) %>%
    as.data.frame()  %>%
    mutate_at(vars(-DrugName), ~ round(., 1)) %>%
    select(DrugName, TIMEWKS_median)
  
  #=== Prevalence of MDD DSM-5 symptoms ===
  table2b <- dat %>%
    group_by(DrugName) %>%
    summarise(
      low_2weeks = sum(LOWINT2W == 1, na.rm = TRUE) / sum(!is.na(LOWINT2W)) * 100,
      depressed_2weeks = sum(DEP2WK == 1, na.rm = TRUE) / sum(!is.na(DEP2WK)) * 100,
      appetite_weight_changes = sum(APWTCHANGE == 1, na.rm = TRUE) / sum(!is.na(APWTCHANGE)) * 100,
      sleep_disturbances = sum(SLEEP == 1, na.rm = TRUE) / sum(!is.na(SLEEP)) * 100,
      movement_changes = sum(MOVEMENT == 1, na.rm = TRUE) / sum(!is.na(MOVEMENT)) * 100,
      fatigued = sum(FATIGUED.x == 1, na.rm = TRUE) / sum(!is.na(FATIGUED.x)) * 100,
      guilty = sum(GUILTY == 1, na.rm = TRUE) / sum(!is.na(GUILTY)) * 100,
      no_focus = sum(NOFOCUS == 1, na.rm = TRUE) / sum(!is.na(NOFOCUS)) * 100,
      death_thoughts = sum(DEATHTHK == 1, na.rm = TRUE) / sum(!is.na(DEATHTHK)) * 100
    )   %>%
    mutate(DrugName = factor(DrugName, levels = custom_order)) %>%
    arrange(DrugName) %>%
    as.data.frame() %>%
    mutate_at(vars(-DrugName), ~ round(., 1))
  
  #=== Prevalence of psychiatric comorbidities ===
  table2c <- dat %>%
    group_by(DrugName) %>%
    summarise(
      across(
        c(DXBPD2, DXPDMD, DXSCZ, DXANOR, DXADHD,
          DXSUD, DXPERSD, DXANX, DXSAD, DXOCD),
        ~ (sum(.x == 1, na.rm = TRUE) / sum(!is.na(.x)) * 100),
        .names = "{.col}"
      )
    ) %>%
    mutate(
      DrugName = factor(DrugName, levels = custom_order)
    ) %>%
    arrange(DrugName) %>%
    as.data.frame() %>%
    mutate(across(-DrugName, ~ round(., 1)))
  
  #=== Prevalence of antidepressant side effects ===
  table2d <- dat %>%
    group_by(DrugName) %>%
    summarise(
      across(
        c(DRYMAD, SWETAD, NAUSAD, VOMAD, DIARAD, CONSAD, HEADAD, DIZZAD, SHAKAD, MUSCAD, DROWAD,
          WAKEAD, HANXAD, AGITAD, FATIAD, WTGNAD, WTLSAD, RASHAD, RUNAD, RSEXAD, BLURAD, SUITAD, NOSEAD, STOPAD),
        ~ (sum(.x == 1, na.rm = TRUE) / sum(!is.na(.x))*100),
        .names = "{.col}"
      )
    ) %>%
    mutate(
      DrugName = factor(DrugName, levels = custom_order)
    ) %>%
    arrange(DrugName) %>%
    as.data.frame() %>%
    mutate(across(-DrugName, ~ round(., 1)))
  
  #=== Table 2 ===
  table2ab <- full_join(table2a, table2b, by = "DrugName")
  table2 <- full_join(table2ab, table2c, by = "DrugName")
  #table2 <- full_join(table2abc, table2d, by = "DrugName")
  transformed_table2 = setNames(data.frame(t(table2[,-1])), table2[,1])
  
  #== Combine Tables 1 and 2 ====
  all <- full_join(table1, table2, by = "DrugName")
  transformed_all = setNames(data.frame(t(all[,-1])), all[,1])
  
  #=== Write summaries ===
  if (duration == 360){
    wb <- createWorkbook()
    #removeWorksheet(wb, "Table5")
    addWorksheet(wb, "Table5")
    writeData(wb, "Table5", transformed_all, rowNames = TRUE)
    saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)
  } else {
    wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
    #removeWorksheet(wb, "Table7")
    addWorksheet(wb, "Table7")
    writeData(wb, "Table7", transformed_all, rowNames = TRUE)
    saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)
  }
  
  #=== Proportions that Self-reported that the antidepressant works well for them ===
  well <- dat %>% 
    select(ParticipantID, WELLAD, DrugName) %>%
    filter(DrugName != "BIP+L" & DrugName != "BIP-L" & DrugName != "Combination" & DrugName != "Various")
  
  well_modified <- well %>%
    mutate(across(-c(ParticipantID, DrugName), ~case_when(
      . == 0 ~ "Not at all well",
      . == 1 ~ "Moderately to very well",
      TRUE ~ NA_character_
    )))
  
  well_long <- well_modified %>%
    pivot_longer(-c(ParticipantID, DrugName), names_to = "Drug", values_to = "Response")
  
  response_proportions <- well_long %>%
    filter(!is.na(Response)) %>%
    group_by(DrugName, Response) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(DrugName) %>%
    mutate(Proportion = Count / sum(Count))
  
  write.csv(response_proportions, paste0("/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes/Drug_Response_Proportions_", duration, "days.csv", row.names = FALSE))
  
}

