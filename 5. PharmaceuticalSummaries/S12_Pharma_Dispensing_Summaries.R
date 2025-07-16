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

#-- Drug Reference Table
# load available antidepressants
ad <- fread(file.path(wkdir, "data/AGDSAcceptabilityTreatmentGroups_14122024.csv")) 
source("/QRISdata/Q7280/pharmacogenomics/Drug_Reference/Drug_Reference_Table.R")
drug_ref <- drug_ref %>%
  filter(ATCCode %in% ad$ATCCode)

# ====================== FUNCTION TO SELECT GROUPING VARIABLE =================

# This function handles setup based on chosen grouping variable
setup_grouping <- function(group_var = "DrugName") {
  result <- list()
  
  if(group_var == "DrugClass") {
    #result$ordered_values <- unique(drug_ref$DrugClass)
    #result$reversed_values <- unique(drug_ref$DrugClass)
    result$group_label <- "Drug Class"
  } else {
    #result$ordered_values <- drug_ref$DrugName
    #result$reversed_values <- rev(drug_ref$DrugName)
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
  
  # Join with drug reference table
  ad_mapped <- left_join(ad, drug_ref, by = c("ATCCode" = "ATCCode"))
  
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
  
  # Extract the ordered levels for consistency
  ordered_levels <- levels(count_data[[group_var]])
  
  # Create plot
  plot1 <- count_data %>%
    ggplot(aes(x = Total, y = !!sym(group_var), 
               fill = !!sym(group_var), color = !!sym(group_var))) +
    geom_bar(stat = "identity") +
    labs(x = "Total Participants", y = config$group_label) +
    scale_y_discrete(position = "right") +
    scale_x_reverse() + 
    create_base_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16, color = "black"),
      plot.margin = unit(c(1, 1, 0, 0), "cm")
    )
  
  return(list(plot = plot1, levels = ordered_levels))
}

#-- Function for Plot 2: Ridge plot of treatment duration
create_ridge_plot <- function(data, config, ordered_levels = NULL) {
  group_var <- config$group_var
  
  # Prepare data with correct factor levels
  plot_data <- data %>%
    mutate(!!sym(group_var) := factor(!!sym(group_var), levels = ordered_levels))
  
  # Create plot
  plot2 <- plot_data %>%
    ggplot(aes(x = PrescriptionDays, y = !!sym(group_var), 
               fill = !!sym(group_var), color = !!sym(group_var))) +
    geom_density_ridges(scale = 1) +
    geom_vline(xintercept = 360, linetype = "longdash", color = "black") +
    labs(x = "Cumulative Prescription\nDuration", y = config$group_label) +
    scale_x_continuous(limits = c(0, 1660)) +
    create_base_theme() +
    theme(
      plot.margin = unit(c(1, 1, 0, 0), "cm"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.title.y = element_blank())
  
  return(plot2)
}

#-- Function for Plot 3: Proportion of participants with >1 prescription episode
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
    
    # Compute proportion with >1 prescription episode
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
  recurrent_combined[[group_var]] <- factor(recurrent_combined[[group_var]])
  recurrent_combined$Threshold <- factor(recurrent_combined$Threshold, levels = threshold_labels)
  
  # Create plot
  plot3 <- ggplot(recurrent_combined, 
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
    create_base_theme() + 
    angle_x_theme() +
    geom_vline(xintercept = "360+ days", linetype = "longdash", color = "black") +
    ylim(0, 1)
  
  return(plot3)
}

#-- Function for Plot 4: Self-report number of depression episodes
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
  DEP_combined[[group_var]] <- factor(DEP_combined[[group_var]])
  DEP_combined$Threshold <- factor(DEP_combined$Threshold, levels = threshold_labels)
  
  # Create plot
  plot4 <- ggplot(DEP_combined, 
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
    create_base_theme() + 
    angle_x_theme() +
    geom_vline(xintercept = "360+ days", linetype = "longdash", color = "black") +
    ylim(3, 14) +
    theme(legend.position = "right")
  
  return(plot4)
}


# ====================== FIGURE COMBINATION FUNCTIONS ========================

#-- Function to combine and save plots
combine_and_save_plots <- function(plots_list, outdir, config) {
  group_var <- config$group_var
  
  # Determine appropriate file names with group_var included
  if(group_var == "DrugClass") {
    top_figure_filename <- "Top_grid_patterns_by_DrugClass.png"
    bottom_figure_filename <- "Bottom_grid_patterns_by_DrugClass.png"
  } else {
    top_figure_filename <- "Top_grid_patterns_by_DrugName.png"
    bottom_figure_filename <- "Bottom_grid_patterns_by_DrugName.png"
  }
    
  # Top grid: Plots 1 and 2
  top_grid <- plot_grid(
    plots_list$plot1, plots_list$plot2,
    nrow = 1,
    align = 'h',
    axis = 'tb',
    labels = c("A", "B"),
    label_size = 18,
    hjust = -0.5,
    vjust = 1.5,
    rel_widths = c(1, 1.2)
  ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # Save the top grid
  png(
    filename = file.path(outdir, top_figure_filename),
    width = 260, height = 170, units = 'mm',
    bg = "white", res = 600, type = c("cairo")
  )
  print(top_grid)
  dev.off()
  
  # Bottom grid: plots 3 and 4
  bottom_grid <- plot_grid(
    plots_list$plot3, plots_list$plot4,
    nrow = 1,
    align = 'h',
    axis = 'tb',
    labels = c("C", "D"),
    label_size = 18,
    hjust = -2,
    vjust = 0,
    rel_widths = c(1, 1.5)
  )
  
  # Add margins to ensure nothing gets cut off
  bottom_grid <- bottom_grid + theme(plot.margin = margin(10, 10, 10, 10, "pt"))
  
  # Save the patterns grid
  png(
    filename = file.path(outdir, bottom_figure_filename),
    width = 270, height = 180, units = 'mm',
    bg = "white", res = 600, type = c("cairo")
  )
  print(bottom_grid )
  dev.off()
  
  # Combine top_grid and bottom_grid vertically
  combined_grid <- plot_grid(
    top_grid, 
    bottom_grid,
    nrow = 2,  # Stack vertically
    align = 'v',  # Vertical alignment
    axis = 'lr',  # Align left-right axes
    rel_heights = c(1.1, 1),  # Adjust relative heights as needed
    greedy = FALSE
  ) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # Save the combined grid
  combined_filename <- if(group_var == "DrugClass") {
    "All_patterns_by_DrugClass.png"
  } else {
    "All_patterns_by_DrugName.png"
  }
  
  png(
    filename = file.path(outdir, combined_filename),
    width = 350, height = 350, units = 'mm',
    bg = "white", res = 600, type = c("cairo")
  )
  print(combined_grid)
  dev.off()
  
  return(list(
    top_figure = top_grid, 
    bottom_figure = bottom_grid,
    combined_figure = combined_grid
  ))
}
  

# ====================== MAIN ANALYSIS FUNCTION =============================

#-- Main function to run analysis with chosen grouping variable
run_analysis <- function(group_var = "DrugName") {
  # Get configuration for the chosen grouping variable
  config <- setup_grouping(group_var)
  
  # Load data
  data <- load_data()
  
  # Create all plots
  plots <- list()
  plot1_result <- create_participant_count_plot(data$ad_mapped, config)
  plots$plot1  <- plot1_result$plot
  plots$plot2 <- create_ridge_plot(data$ad_mapped, config, plot1_result$levels)
  plots$plot3 <- create_repeat_prescription_plot(data$ad_mapped, config)
  plots$plot4 <- create_depression_episodes_plot(data$inner_combined, config)
  
  # Combine and save plots
  figures <- combine_and_save_plots(plots, outdir, config)
  
  # Return both plots and combined figures
  return(list(individual_plots = plots, combined_figures = figures))
}

# ====================== FINAL USAGE ====================================

# Run analysis with DrugName
drugname_results <- run_analysis("DrugName")

# Run analysis with DrugClass
drugclass_results <- run_analysis("DrugClass")