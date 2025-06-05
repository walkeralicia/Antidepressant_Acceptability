# Complete Pharmacogenomics Analysis Script

# ========================= SETUP AND CONFIGURATION ===========================

#-- Load required libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggridges)
library(viridis)
library(patchwork)
library(forcats)
library(cowplot)

#-- Set working directories
wkdir <- "/QRISdata/Q7280/pharmacogenomics"
outdir <- "/scratch/user/uqawal15"

# ====================== FUNCTION TO SELECT GROUPING VARIABLE =================

# This function handles setup based on chosen grouping variable
setup_grouping <- function(group_var = "DrugName") {
  result <- list()
  
  if(group_var == "DrugClass") {
    # Define drug class settings
    result$ordered_values <- c("SSRI", "SNRI", "TCA", "TeCA")
    result$reversed_values <- c("SSRI", "SNRI", "TeCA", "TCA")
    result$color_map <- c(
      "SSRI" = "#2C7FB8",
      "SNRI" = "#BE3A8A",
      "TCA" = "#FE9929",
      "TeCA" = "#FFD92F"
    )
    result$group_label <- "Drug Class"
  } else {
    # Default to DrugName settings
    result$ordered_values <- c("SSRI:Sertraline", "SSRI:Escitalopram", "SSRI:Citalopram", "SSRI:Fluoxetine", 
                               "SSRI:Paroxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SNRI:Duloxetine", 
                               "TCA:Amitriptyline", "TeCA:Mirtazapine")
    result$reversed_values <- rev(c("SSRI:Paroxetine", "SSRI:Citalopram", "TCA:Amitriptyline", 
                                    "SNRI:Duloxetine", "SSRI:Fluoxetine", "TeCA:Mirtazapine",
                                    "SNRI:Venlafaxine", "SNRI:Desvenlafaxine", "SSRI:Sertraline", 
                                    "SSRI:Escitalopram"))
    result$color_map <- c(
      "SSRI:Sertraline" = "#2C7FB8", "SSRI:Escitalopram" = "#4DA6FF",
      "SSRI:Citalopram" = "#1A9850", "SSRI:Paroxetine" = "#7FCDBB",
      "SSRI:Fluoxetine" = "#225EA8", "SNRI:Desvenlafaxine" = "#BE3A8A",
      "SNRI:Venlafaxine" = "#984EA3", "SNRI:Duloxetine" = "#DC267F",
      "TCA:Amitriptyline" = "#FE9929", "TeCA:Mirtazapine" = "#FFD92F"
    )
    result$group_label <- "Antidepressant"
  }
  
  # Define thresholds consistently
  result$thresholds <- c(0, 90, 180, 360, 600, 800, 1000, 1200, 1500)
  result$threshold_labels <- c("1+ days", "90+ days", "180+ days", "360+ days", "600+ days", 
                               "800+ days", "1000+ days", "1200+ days", "1500+ days")
  
  # Pass through the grouping variable
  result$group_var <- group_var
  
  return(result)
}

# ====================== DATA LOADING AND PREPARATION ========================

#-- Function to load data
load_data <- function() {
  # Load treatment files
  ad <- fread(file.path(wkdir, "data/AGDSAcceptabilityTreatmentGroups_14122024.csv")) %>%
    mutate(
      EarliestPrescription = as.Date(EarliestPrescription, format = "%d/%m/%Y"),
      LatestPrescription = as.Date(LatestPrescription, format = "%d/%m/%Y")
    )
  
  #-- Load phenotype data
  pheno <- read.csv(file.path(wkdir, "phenotypes/survey_phenotypes.csv"))
  
  #-- Participants with phenotype and pharma data
  ad <- ad %>% filter(ParticipantID %in% pheno$STUDYID)
  
  # Drug reference data
  atc_drug_ref <- data.frame(
    DrugName = c("SSRI:Sertraline", "SSRI:Escitalopram", "SSRI:Citalopram", "SNRI:Venlafaxine", 
                 "SNRI:Duloxetine", "TCA:Amitriptyline", "SSRI:Paroxetine", "TeCA:Mirtazapine", 
                 "SNRI:Desvenlafaxine", "SSRI:Fluoxetine"),
    ATCCodes = c('N06AB06', 'N06AB10', 'N06AB04', 'N06AX16', 'N06AX21', 'N06AA09', 
                 'N06AB05', 'N06AX11', 'N06AX23', 'N06AB03'),
    DrugClass = c('SSRI', 'SSRI', 'SSRI', 'SNRI', 'SNRI', 'TCA', 'SSRI', 'TeCA', 'SNRI', 'SSRI')
  )
  
  # Join and prepare main dataset
  ad_mapped <- left_join(ad, atc_drug_ref, by = c("ATCCode" = "ATCCodes"))
  
  # Create combined dataset
  inner_combined <- inner_join(ad_mapped, pheno, by = c("ParticipantID" = "STUDYID"))
  
  return(list(
    ad_mapped = ad_mapped,
    pheno = pheno,
    inner_combined = inner_combined
  ))
}

# ====================== PLOT CREATION FUNCTIONS =============================

# Common theme elements
create_base_theme <- function() {
  theme_classic(base_size = 18) +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(color = "black", face = "bold", size = 14),
      axis.line = element_line(colour = 'black'),
      legend.position = "none",
      plot.title = element_text(size = 18)
    )
}

angle_x_theme <- function() {
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1)
  )
}

# ====================== INDIVIDUAL PLOT FUNCTIONS ==========================

#-- Function for Plot 1: Number of participants by grouping variable
create_participant_count_plot <- function(data, config) {
  group_var <- config$group_var
  
  # Count participants
  count_data <- data %>%
    group_by(!!sym(group_var)) %>%
    summarise(Total = n_distinct(ParticipantID)) %>%
    mutate(!!sym(group_var) := fct_reorder(!!sym(group_var), Total, .desc = TRUE))
  
  # Create plot
  plot1 <- count_data %>%
    ggplot(aes(x = Total, y = !!sym(group_var), 
               fill = !!sym(group_var), color = !!sym(group_var))) +
    geom_bar(stat = "identity") +
    labs(x = "Total AGDS Participants", y = config$group_label) +
    scale_y_discrete(position = "right") +
    scale_x_reverse() + 
    scale_color_manual(values = config$color_map) + 
    scale_fill_manual(values = config$color_map) +
    create_base_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16, color = "black"),
      plot.margin = unit(c(1, 1, 0, 0), "cm")
    )
  
  return(plot1)
}

#-- Function for Plot 2: Ridge plot of treatment duration
create_ridge_plot <- function(data, config) {
  group_var <- config$group_var
  
  # Prepare data with correct factor levels
  plot_data <- data %>%
    mutate(!!sym(group_var) := factor(!!sym(group_var), levels = config$reversed_values))
  
  # Create plot
  plot2 <- plot_data %>%
    ggplot(aes(x = PrescriptionDays, y = !!sym(group_var), 
               fill = !!sym(group_var), color = !!sym(group_var))) +
    geom_density_ridges(scale = 1) +
    geom_vline(xintercept = 360, linetype = "longdash", color = "black") +
    labs(x = "Cumulative Prescription\nDuration", y = config$group_label) +
    scale_x_continuous(limits = c(0, 1660)) +
    scale_color_manual(values = config$color_map) + 
    scale_fill_manual(values = config$color_map) +
    create_base_theme() +
    theme(
      plot.margin = unit(c(1, 1, 0, 0), "cm"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.title.y = element_blank())
  
  return(plot2)
}

#-- Function for Plot 3: Self-report efficacy by prescription duration
create_efficacy_plot <- function(data, config) {
  group_var <- config$group_var
  thresholds <- config$thresholds
  threshold_labels <- config$threshold_labels
  
  # Create WELLAD data
  wellad_data <- data %>%
    mutate(
      WELLAD = case_when(
        DrugName == "SSRI:Sertraline" & Sertraline == 1 ~ "Moderately to very well",
        DrugName == "SSRI:Sertraline" & Sertraline == 0 ~ "Not at all well",
        DrugName == "SSRI:Escitalopram" & Escitalopram == 1 ~ "Moderately to very well",
        DrugName == "SSRI:Escitalopram" & Escitalopram == 0 ~ "Not at all well",
        DrugName == "SSRI:Citalopram" & Citalopram == 1 ~ "Moderately to very well",
        DrugName == "SSRI:Citalopram" & Citalopram == 0 ~ "Not at all well",
        DrugName == "SSRI:Fluoxetine" & Fluoxetine == 1 ~ "Moderately to very well",
        DrugName == "SSRI:Fluoxetine" & Fluoxetine == 0 ~ "Not at all well",
        DrugName == "SSRI:Paroxetine" & Paroxetine == 1 ~ "Moderately to very well",
        DrugName == "SSRI:Paroxetine" & Paroxetine == 0 ~ "Not at all well",
        DrugName == "SNRI:Desvenlafaxine" & Desvenlafaxine == 1 ~ "Moderately to very well",
        DrugName == "SNRI:Desvenlafaxine" & Desvenlafaxine == 0 ~ "Not at all well",
        DrugName == "SNRI:Venlafaxine" & Venlafaxine == 1 ~ "Moderately to very well",
        DrugName == "SNRI:Venlafaxine" & Venlafaxine == 0 ~ "Not at all well",
        DrugName == "SNRI:Duloxetine" & Duloxetine == 1 ~ "Moderately to very well",
        DrugName == "SNRI:Duloxetine" & Duloxetine == 0 ~ "Not at all well",
        DrugName == "TCA:Amitriptyline" & Amitriptyline == 1 ~ "Moderately to very well",
        DrugName == "TCA:Amitriptyline" & Amitriptyline == 0 ~ "Not at all well",
        DrugName == "TeCA:Mirtazapine" & Mirtazapine == 1 ~ "Moderately to very well",
        DrugName == "TeCA:Mirtazapine" & Mirtazapine == 0 ~ "Not at all well",
        TRUE ~ NA_character_
      )
    ) %>%
    select(ParticipantID, !!sym(group_var), DrugName, PrescriptionDays, WELLAD)
  
  # Process each threshold manually
  all_responses_list <- list()
  
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    label <- threshold_labels[i]
    
    # Process this threshold
    threshold_data <- wellad_data %>%
      filter(PrescriptionDays >= threshold) %>%
      filter(!is.na(WELLAD)) %>%
      group_by(!!sym(group_var), WELLAD) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(!!sym(group_var)) %>%
      mutate(Proportion = Count / sum(Count), 
             Threshold = label)
    
    all_responses_list[[i]] <- threshold_data
  }
  
  all_responses <- bind_rows(all_responses_list)
  
  # Apply factor levels
  all_responses[[group_var]] <- factor(all_responses[[group_var]], levels = config$ordered_values)
  all_responses$Threshold <- factor(all_responses$Threshold, levels = threshold_labels)
  
  # Filter for well responses
  all_responses_well <- all_responses %>% 
    filter(WELLAD == "Moderately to very well")
  
  # Create plot
  plot3 <- ggplot(all_responses_well, 
                  aes(x = Threshold, y = Proportion, 
                      color = !!sym(group_var), group = !!sym(group_var))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(
      x = "Cumulative Prescription\nDuration Threshold", 
      y = "Proportion of participants",
      color = config$group_label,
      title = "Self-Report Responders"
    ) +
    scale_color_manual(values = config$color_map) +
    create_base_theme() + 
    angle_x_theme() +
    geom_vline(xintercept = "360+ days", linetype = "longdash", color = "black") +
    ylim(0, 1)
  
  return(plot3)
}

#-- Function for Plot 4: Self-report discontinuation due to side effects
create_discontinuation_plot <- function(data, config) {
  group_var <- config$group_var
  thresholds <- config$thresholds
  threshold_labels <- config$threshold_labels
  
  # Create STOPAD data
  stopad_data <- data %>%
    mutate(
      STOPAD = case_when(
        DrugName == "SSRI:Sertraline" & STOPSERT == 1 ~ "Yes",
        DrugName == "SSRI:Sertraline" & STOPSERT == 0 ~ "No",
        DrugName == "SSRI:Escitalopram" & STOPESCI == 1 ~ "Yes",
        DrugName == "SSRI:Escitalopram" & STOPESCI == 0 ~ "No",
        DrugName == "SSRI:Citalopram" & STOPCITA == 1 ~ "Yes",
        DrugName == "SSRI:Citalopram" & STOPCITA == 0 ~ "No",
        DrugName == "SSRI:Fluoxetine" & STOPFLUO == 1 ~ "Yes",
        DrugName == "SSRI:Fluoxetine" & STOPFLUO == 0 ~ "No",
        DrugName == "SSRI:Paroxetine" & STOPPARO == 1 ~ "Yes",
        DrugName == "SSRI:Paroxetine" & STOPPARO == 0 ~ "No",
        DrugName == "SNRI:Desvenlafaxine" & STOPDESV == 1 ~ "Yes",
        DrugName == "SNRI:Desvenlafaxine" & STOPDESV == 0 ~ "No",
        DrugName == "SNRI:Venlafaxine" & STOPVENL == 1 ~ "Yes",
        DrugName == "SNRI:Venlafaxine" & STOPVENL == 0 ~ "No",
        DrugName == "SNRI:Duloxetine" & STOPDULO == 1 ~ "Yes",
        DrugName == "SNRI:Duloxetine" & STOPDULO == 0 ~ "No",
        DrugName == "TCA:Amitriptyline" & STOPAMIT == 1 ~ "Yes",
        DrugName == "TCA:Amitriptyline" & STOPAMIT == 0 ~ "No",
        DrugName == "TeCA:Mirtazapine" & STOPMIRT == 1 ~ "Yes",
        DrugName == "TeCA:Mirtazapine" & STOPMIRT == 0 ~ "No",
        TRUE ~ NA_character_
      )
    ) %>%
    select(ParticipantID, !!sym(group_var), DrugName, PrescriptionDays, STOPAD)
  
  # Process each threshold manually
  all_stopad_list <- list()
  
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    label <- threshold_labels[i]
    
    # Process this threshold
    threshold_data <- stopad_data %>%
      filter(PrescriptionDays >= threshold) %>%
      filter(!is.na(STOPAD)) %>%
      group_by(!!sym(group_var), STOPAD) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(!!sym(group_var)) %>%
      mutate(Proportion = Count / sum(Count), 
             Threshold = label)
    
    all_stopad_list[[i]] <- threshold_data
  }
  
  all_stopad <- bind_rows(all_stopad_list)
  
  # Apply factor levels
  all_stopad[[group_var]] <- factor(all_stopad[[group_var]], levels = config$ordered_values)
  all_stopad$Threshold <- factor(all_stopad$Threshold, levels = threshold_labels)
  
  # Filter for "Yes" responses
  stopad_combined_yes <- all_stopad %>% filter(STOPAD == "Yes")
  
  # Create plot
  plot4 <- ggplot(stopad_combined_yes, 
                  aes(x = Threshold, y = Proportion, 
                      color = !!sym(group_var), group = !!sym(group_var))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(
      x = "Cumulative Prescription\nDuration Threshold", 
      y = "Proportion of participants",
      color = config$group_label,
      title = "Self-Report Discontinuation\ndue to Side Effects"
    ) +
    scale_color_manual(values = config$color_map) +
    create_base_theme() + 
    angle_x_theme() +
    geom_vline(xintercept = "360+ days", linetype = "longdash", color = "black") +
    ylim(0, 1)
  
  return(plot4)
}

#-- Function for Plot 5: Proportion of participants with >1 treatment period
create_repeat_prescription_plot <- function(data, config) {
  group_var <- config$group_var
  thresholds <- config$thresholds
  threshold_labels <- config$threshold_labels
  
  # Process each threshold manually
  recurrent_list <- list()
  
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    label <- threshold_labels[i]
    
    # Filter data by threshold
    filtered_data <- data %>% 
      filter(PrescriptionDays >= threshold)
    
    # Compute proportion with >1 treatment period
    summary_data <- filtered_data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        Total = n(),
        GreaterThanOne = sum(NumberOfTreatmentPeriods > 1),
        .groups = "drop"
      ) %>%
      mutate(Proportion = GreaterThanOne / Total, 
             Threshold = label)
    
    recurrent_list[[i]] <- summary_data
  }
  
  recurrent_combined <- bind_rows(recurrent_list)
  
  # Apply factor levels
  recurrent_combined[[group_var]] <- factor(recurrent_combined[[group_var]], levels = config$ordered_values)
  recurrent_combined$Threshold <- factor(recurrent_combined$Threshold, levels = threshold_labels)
  
  # Create plot
  plot5 <- ggplot(recurrent_combined, 
                  aes(x = Threshold, y = Proportion, 
                      color = !!sym(group_var), group = !!sym(group_var))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(
      x = "Cumulative Prescription\nDuration Threshold", 
      y = "Proportion of participants",
      color = config$group_label,
      title = "Repeat prescription episodes"
    ) +
    scale_color_manual(values = config$color_map) +
    create_base_theme() + 
    angle_x_theme() +
    geom_vline(xintercept = "360+ days", linetype = "longdash", color = "black") +
    ylim(0, 1)
  
  return(plot5)
}

#-- Function for Plot 6: Self-report number of depression episodes
create_depression_episodes_plot <- function(data, config) {
  group_var <- config$group_var
  thresholds <- config$thresholds
  threshold_labels <- config$threshold_labels
  
  # Process data for median episodes
  dep_data <- data %>% 
    select(ParticipantID, !!sym(group_var), DrugName, TIMES2WK, PrescriptionDays)
  
  # Process each threshold manually
  dep_list <- list()
  
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    label <- threshold_labels[i]
    
    # Filter and calculate median
    threshold_data <- dep_data %>%
      filter(PrescriptionDays >= threshold) %>%
      filter(!is.na(TIMES2WK)) %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        Median = median(TIMES2WK, na.rm = TRUE),
        IQR = IQR(TIMES2WK, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Threshold = label)
    
    dep_list[[i]] <- threshold_data
  }
  
  DEP_combined <- bind_rows(dep_list)
  
  # Apply factor levels
  DEP_combined[[group_var]] <- factor(DEP_combined[[group_var]], levels = config$ordered_values)
  DEP_combined$Threshold <- factor(DEP_combined$Threshold, levels = threshold_labels)
  
  # Create plot
  plot6 <- ggplot(DEP_combined, 
                  aes(x = Threshold, y = Median, 
                      color = !!sym(group_var), group = !!sym(group_var))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(
      x = "Cumulative Prescription\nDuration Threshold", 
      y = "Median",
      color = config$group_label,
      title = "Self-reported number of MDD episodes"
    ) +
    scale_color_manual(values = config$color_map) +
    create_base_theme() + 
    angle_x_theme() +
    geom_vline(xintercept = "360+ days", linetype = "longdash", color = "black") +
    ylim(3, 14) +
    theme(legend.position = "right")
  
  return(plot6)
}

#-- Function for Plot 7: Alternative for DrugClass and DrugName
create_plot7 <- function(data, config) {
  group_var <- config$group_var
  
  # Group participants by treatment periods
  ad_mapped_restricted <- data %>%
    filter(PrescriptionDays < 800) %>%
    mutate(NumberOfTreatmentPeriods_grouped = 
             ifelse(NumberOfTreatmentPeriods > 4, "5+", as.character(NumberOfTreatmentPeriods))) %>%
    mutate(NumberOfTreatmentPeriods_grouped = 
             factor(NumberOfTreatmentPeriods_grouped, 
                    levels = c("1", "2", "3", "4", "5+"), ordered = TRUE))
  
  if (group_var == "DrugName") {
    # For DrugName, use faceted plot with error bars
    summary_data <- ad_mapped_restricted %>%
      group_by(!!sym(group_var), NumberOfTreatmentPeriods_grouped) %>%
      summarise(
        Median_Duration = median(PrescriptionDays, na.rm = TRUE),
        IQR = IQR(PrescriptionDays, na.rm = TRUE)
      ) %>%
      arrange(!!sym(group_var), NumberOfTreatmentPeriods_grouped)
    
    plot7 <- ggplot(summary_data, 
                    aes(x = Median_Duration, y = NumberOfTreatmentPeriods_grouped, 
                        group = !!sym(group_var), color = !!sym(group_var))) +
      geom_line(orientation = "y", linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(xmin = Median_Duration - IQR/2, xmax = Median_Duration + IQR/2), 
                    width = 0.2, linewidth = 1.5) +
      facet_wrap(as.formula(paste("~", group_var)), nrow = 2) +
      scale_y_discrete() +
      theme_bw(base_size = 16) +
      scale_color_manual(values = config$color_map) +
      labs(
        x = "Median Cumulative Prescription Dispense (CPD) with IQR", 
        y = "Total Prescription Episodes (TPE)",
        title = "TPE across median TPD for those with <800 CPD", 
        color = config$group_label
      ) +
      angle_x_theme() +
      theme(
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14),
        legend.position = "none",
        plot.title = element_text(size = 18)
      )
  } else {
    # For DrugClass, create a stacked bar plot
    summary_data <- ad_mapped_restricted %>%
      group_by(DrugClass, NumberOfTreatmentPeriods_grouped) %>%
      summarise(Count = n(), .groups = "drop")
    
    plot7 <- ggplot(summary_data, 
                    aes(x = DrugClass, y = Count, fill = NumberOfTreatmentPeriods_grouped)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      labs(
        x = "Drug Class", 
        y = "Number of Participants",
        title = "Number of treatment episodes by drug class", 
        fill = "Treatment\nEpisodes"
      ) +
      theme_bw(base_size = 16) +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = 14),
        legend.position = "right",
        plot.title = element_text(size = 18)
      )
  }
  
  return(plot7)
}

#-- Function for Plot 8: Medication or Class specific effects
create_effects_plot <- function(config, effect_type = "medication") {
  # Determine which file and settings to use based on effect_type
  if(effect_type == "medication") {
    # Medication effects (using individual drugs)
    effects_file <- "LM_PrescriptionDays_MedicationEffects_Escitalopram_Reference.csv"
    reference_value <- "Reference:SSRI:Escitalopram"
    y_label <- "Medication"
    x_label <- "Cumulative Prescription\nDuration:\nDifference from\nreference SSRI:Escitalopram"
    filter_values <- c("Reference:SSRI:Escitalopram", "SSRI:Sertraline", "SNRI:Desvenlafaxine", 
                       "SNRI:Venlafaxine", "TeCA:Mirtazapine", "SSRI:Fluoxetine",
                       "SNRI:Duloxetine", "TCA:Amitriptyline", "SSRI:Citalopram", "SSRI:Paroxetine")
    factor_levels <- c("Reference:SSRI:Escitalopram", "SSRI:Sertraline", "SNRI:Desvenlafaxine", 
                       "SNRI:Venlafaxine", "TeCA:Mirtazapine", "SSRI:Fluoxetine",
                       "SNRI:Duloxetine", "TCA:Amitriptyline", "SSRI:Citalopram", "SSRI:Paroxetine")
    x_limits <- c(-400, 250)
    # For medication effects, hide y-axis text (will be shown in plot1)
    show_y_text <- FALSE
  } else {
    # Class effects
    effects_file <- "LM_PrescriptionDays_ClassEffects_SSRI_Reference.csv"
    reference_value <- "Reference:SSRI"
    y_label <- "Drug Class"
    x_label <- "Cumulative Prescription\nDuration:\nDifference from\nreference SSRI"
    filter_values <- c("Reference:SSRI", "SNRI", "TCA", "TeCA")
    factor_levels <- c("Reference:SSRI", "SNRI", "TeCA", "TCA")
    x_limits <- c(-400, 250)
    # For class effects, show y-axis text 
    show_y_text <- TRUE
  }
  
  # Read and process data
  effects_data <- read.csv(file.path("/QRISdata/Q7280/pharmacogenomics/pharma_summaries", effects_file)) %>%
    mutate(Dependent = "PrescriptionDays") %>%
    fill(Dependent, Total_N, .direction = "down")
  
  # Filter and calculate significance
  processed_data <- effects_data %>%
    filter(Term %in% filter_values) %>%
  mutate(
    Neglog10Pvalue = -log10(`Pr...t..`),
    BF_Sig = ifelse(`Pr...t..` < 0.05, "Significant", "Insignificant"),
    BF_Sig = ifelse(Term == reference_value, "Reference", BF_Sig),
    Term = factor(Term, levels = factor_levels)
  )

  # Create plot
  effects_plot <- ggplot(processed_data, aes(y = Term, x = Estimate, color = BF_Sig)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +  
  geom_errorbar(
    aes(xmin = Estimate - (1.96 * `Std..Error`), xmax = Estimate + (1.96 * `Std..Error`)), 
    position = position_dodge(width = 0.8), 
    linewidth = 1, width = 0.3
  ) +
  geom_vline(xintercept = 0) +
  ylab(y_label) +
  xlab(x_label) +
  theme_classic(base_size = 20) +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(color = "black", face = "bold", size = 16),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  xlim(x_limits) +
  scale_color_manual(
    values = c(
      "Significant" = "#0047AB",
      "Insignificant" = "#7F7F7F",
      "Reference" = "black"
    ),
    name = "Significance Level", na.translate = FALSE
  )

return(effects_plot)
}


# ====================== FIGURE COMBINATION FUNCTIONS ========================

#-- Function to combine and save plots
combine_and_save_plots <- function(plots_list, outdir, config) {
  group_var <- config$group_var
  
  # Determine appropriate filenames with group_var included
  if(group_var == "DrugClass") {
    main_figure_filename <- "MF_AGDS_by_DrugClass.png"
    efficacy_filename <- "SF_AGDS_efficacy_by_DrugClass.png"
    patterns_filename <- "SF_AGDS_patterns_by_DrugClass.png"
  } else {
    main_figure_filename <- "MF_AGDS_by_DrugName.png"
    efficacy_filename <- "SF_AGDS_efficacy_by_DrugName.png"
    patterns_filename <- "SF_AGDS_patterns_by_DrugName.png"
  }
  
  # Create main figure with plots 1, 2, and effects plots
  if(!is.null(plots_list$medication_effects)) {
    # Include medication effects for DrugName
    top_grid <- plot_grid(
      plots_list$plot1, plots_list$plot2, plots_list$medication_effects,
      nrow = 1,
      align = 'h',
      axis = 'tb',
      labels = c("A", "B", "C"),
      label_size = 18,
      hjust = -0.5,
      vjust = 1.5,
      rel_widths = c(0.6, 0.8, 0.6)
    ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  } else if(!is.null(plots_list$class_effects) && group_var == "DrugClass") {
    # Include class effects for DrugClass
    top_grid <- plot_grid(
      plots_list$plot1, plots_list$plot2, plots_list$class_effects,
      nrow = 1,
      align = 'h',
      axis = 'tb',
      labels = c("A", "B", "C"),
      label_size = 18,
      hjust = -0.5,
      vjust = 1.5,
      rel_widths = c(0.6, 0.8, 0.6)
    ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  } else {
    # Only include plots 1 and 2 if no effects plots
    top_grid <- plot_grid(
      plots_list$plot1, plots_list$plot2,
      nrow = 1,
      align = 'h',
      axis = 'tb',
      labels = c("A", "B"),
      label_size = 18,
      hjust = -0.5,
      vjust = 1.5,
      rel_widths = c(0.6, 0.8)
    ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  }
  
  # Save the top grid
  png(
    filename = file.path(outdir, main_figure_filename),
    width = 350, height = 200, units = 'mm',
    bg = "white", res = 600, type = c("cairo")
  )
  print(top_grid)
  dev.off()
  
  # Create efficacy figure with plots 3 and 4
  efficacy_grid <- plot_grid(
    plots_list$plot3, 
    plots_list$plot4 + theme(legend.position = "right"),
    nrow = 1,
    align = 'h',
    axis = 'tb',
    labels = c("D", "E"),
    label_size = 18,
    hjust = -2,
    vjust = 0,
    rel_widths = c(1, 1.5)
  ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # Save the efficacy grid
  png(
    filename = file.path(outdir, efficacy_filename),
    width = 350, height = 150, units = 'mm',
    bg = "white", res = 300, type = c("cairo")
  )
  print(efficacy_grid)
  dev.off()
  
  # Create patterns figure with plots 5, 6, and 7
  patterns_top_grid <- plot_grid(
    plots_list$plot5, plots_list$plot6,
    nrow = 1,
    align = 'h',
    axis = 'tb',
    labels = c("A", "B"),
    label_size = 18,
    hjust = -2,
    vjust = 0,
    rel_widths = c(1, 1.5)
  )
  
  # Combine with plot7
  patterns_grid <- plot_grid(
    patterns_top_grid, plots_list$plot7,
    nrow = 2,
    align = 'h',
    axis = 'tb',
    labels = c("", "C"),
    label_size = 18,
    hjust = -2,
    vjust = 0,
    rel_heights = c(1, 1),
    greedy = FALSE
  )
  
  # Add margins to ensure nothing gets cut off
  patterns_grid <- patterns_grid + theme(plot.margin = margin(10, 10, 10, 10, "pt"))
  
  # Save the patterns grid
  png(
    filename = file.path(outdir, patterns_filename),
    width = 250, height = 300, units = 'mm',
    bg = "white", res = 600, type = c("cairo")
  )
  print(patterns_grid)
  dev.off()
  
  # Combine top_grid and efficacy_grid vertically
  combined_grid <- plot_grid(
    top_grid, 
    efficacy_grid,
    nrow = 2,  # Stack vertically
    align = 'v',  # Vertical alignment
    axis = 'lr',  # Align left-right axes
    rel_heights = c(1.1, 1),  # Adjust relative heights as needed
    greedy = FALSE
  ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # Save the combined grid
  combined_filename <- if(group_var == "DrugClass") {
    "Combined_AGDS_by_DrugClass.png"
  } else {
    "Combined_AGDS_by_DrugName.png"
  }
  
  png(
    filename = file.path(outdir, combined_filename),
    width = 350, height = 350, units = 'mm',  # Increased height for combined plot
    bg = "white", res = 600, type = c("cairo")
  )
  print(combined_grid)
  dev.off()
  
  return(list(
    main_figure = top_grid, 
    efficacy_figure = efficacy_grid, 
    patterns_figure = patterns_grid,
    combined_figure = combined_grid
  ))
}

# ====================== MAIN ANALYSIS FUNCTION =============================

#-- Main function to run analysis with chosen grouping variable
run_analysis <- function(group_var = "DrugName", include_effects = TRUE) {
  # Get configuration for the chosen grouping variable
  config <- setup_grouping(group_var)
  
  # Load data
  data <- load_data()
  
  # Create all plots
  plots <- list()
  plots$plot1 <- create_participant_count_plot(data$ad_mapped, config)
  plots$plot2 <- create_ridge_plot(data$ad_mapped, config)
  plots$plot3 <- create_efficacy_plot(data$inner_combined, config)
  plots$plot4 <- create_discontinuation_plot(data$inner_combined, config)
  plots$plot5 <- create_repeat_prescription_plot(data$ad_mapped, config)
  plots$plot6 <- create_depression_episodes_plot(data$inner_combined, config)
  plots$plot7 <- create_plot7(data$ad_mapped, config)
  
  # Create effects plots based on group_var and options
  if(include_effects) {
    if(group_var == "DrugName") {
      # Create medication effects plot for DrugName
      plots$medication_effects <- create_effects_plot(config, effect_type = "medication")

    } else if(group_var == "DrugClass") {
      # Create class effects plot for DrugClass
      plots$class_effects <- create_effects_plot(config, effect_type = "class")
    }
  }
  
  # Combine and save plots
  figures <- combine_and_save_plots(plots, outdir, config)
  
  # Return both plots and combined figures
  return(list(individual_plots = plots, combined_figures = figures))
}

# ====================== USAGE EXAMPLES ====================================

# Run analysis with DrugName, including effects plots
drugname_results <- run_analysis("DrugName", include_effects = TRUE)

# Run analysis with DrugClass, including class effects plot
drugclass_results <- run_analysis("DrugClass", include_effects = TRUE)