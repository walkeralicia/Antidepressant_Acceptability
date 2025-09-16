# -- Load R libraries
library(tidyverse)  # Includes dplyr, ggplot2, tidyr
library(cowplot)
library(readxl)

# -- Read in association results and process data in a single pipeline
processed_data <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet = "Table4") %>%
  # Fill in NA values
  fill(Analysis, Dependent, PGS, .direction = "down") %>%
  # Filter for full samples
  filter(Analysis == "Full_Sample") %>%
  # Filter for std_pgs terms only
  filter(Term == "std_pgs") %>%
  select(Dependent, PGS, estimate, std.error, P.value, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf) %>%
  # Add significance labels
  mutate(
    P.value = as.numeric(P.value),
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Sig_Nom = case_when(
      P.value < 0.05 ~ "*",
      P.value >= 0.05 ~ NA_character_
    )
  )

# -- Get ordered PGS based on Class Diversity estimates
ordered_pgs <- processed_data %>%
  filter(Dependent == "Class Diversity") %>%
  arrange(estimate) %>%
  pull(PGS)

# -- Create dataset by outcome, with consistent PGS ordering
plot_data <- list(
  duration = processed_data %>% 
    filter(Dependent == "Cumulative Prescription Dispense (days)") %>%
    mutate(PGS = factor(PGS, levels = ordered_pgs)),
  
  diversity = processed_data %>%
    filter(Dependent == "Medication Diversity") %>%
    mutate(PGS = factor(PGS, levels = ordered_pgs)),
  
  class_diversity = processed_data %>%
    filter(Dependent == "Class Diversity") %>%
    mutate(PGS = factor(PGS, levels = ordered_pgs))
)

# -- Create a function for plot generation
create_plot <- function(data, title, x_offset = 0.01, show_y_axis = TRUE) {
  dodge_width <- 0.8
  
  data$point_color <- ifelse(data$Sig_Nom == "*", "#2C7FB8", "white")
  
  p <- ggplot(data, aes(y = PGS, x = estimate)) +
    geom_point(position = position_dodge(width = dodge_width), 
               size = 3, fill = data$point_color, shape = 21, color = "#2C7FB8") +  
    geom_errorbar(aes(xmin = estimate - (1.96 * std.error), 
                      xmax = estimate + (1.96 * std.error)), 
                  position = position_dodge(width = dodge_width), 
                  linewidth = 0.5, width = 0.3, color = "#2C7FB8") +
    geom_vline(xintercept = 0, linetype = "longdash") +
    xlab(title) +
    theme_classic(base_size = 16) +
    geom_text(aes(label = Sig_label, 
                  y = PGS, 
                  x = (estimate + (1.96 * std.error) + x_offset)),
              color = "black",
              size = 6,
              position = position_dodge(width = dodge_width)) +  
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing = unit(1, "lines"),
      panel.heights = unit(c(5, 8, 12, 10), "null"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(color = "black", face = "bold", size = 12),
      legend.position = "none"
    )
  
  # Configure y-axis display based on parameter
  if (show_y_axis) {
    p <- p + 
      ylab("Polygenic Risk Scores") + 
      theme(
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold")
      )
  } else {
    p <- p + 
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  return(p)
}

# -- Create plots with different configurations
# Calculate appropriate x-offsets based on the range of data
x_offset1 <- max(plot_data$duration$estimate + 1.96 * plot_data$duration$std.error, na.rm = TRUE) * 0.1
x_offset2 <- max(plot_data$diversity$estimate + 1.96 * plot_data$diversity$std.error, na.rm = TRUE) * 0.15
x_offset3 <- max(plot_data$class_diversity$estimate + 1.96 * plot_data$class_diversity$std.error, na.rm = TRUE) * 0.15

p1 <- create_plot(plot_data$duration, "Cumulative AD Dispense", x_offset1, TRUE)
p2 <- create_plot(plot_data$diversity, "AD Diversity", x_offset2, FALSE)
p3 <- create_plot(plot_data$class_diversity, "Class Diversity", x_offset3, FALSE)

# -- Combine plots
plots <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  align = 'h',
  axis = 'tb',
  rel_widths = c(1.8, 1, 1)
)

# -- Save combined plot
ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\Figure_2.png", 
       plot = plots, device = "png", width = 230, height = 140, units = "mm", dpi = 600)

#====== Plot for sensitivity analysis
# -- Read in association results and process data in a single pipeline
processed_data <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet = "Table4") %>%
  # Fill in NA values
  fill(Analysis, Dependent, PGS, .direction = "down") %>%
  # Filter for BIP Excluded
  filter(Analysis == "BIP_Excluded") %>%
  # Filter for std_pgs terms only
  filter(Term == "std_pgs") %>%
  select(Dependent, PGS, estimate, std.error, P.value, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf) %>%
  # Add significance labels
  mutate(
    P.value = as.numeric(P.value),
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Sig_Nom = case_when(
      P.value < 0.05 ~ "*",
      P.value >= 0.05 ~ NA_character_
    )
  )

# -- Get ordered PGS based on Class Diversity estimates
ordered_pgs <- processed_data %>%
  filter(Dependent == "Class Diversity") %>%
  arrange(estimate) %>%
  pull(PGS)

# -- Create dataset by outcome, with consistent PGS ordering
plot_data <- list(
  duration = processed_data %>% 
    filter(Dependent == "Cumulative Prescription Dispense (days)") %>%
    mutate(PGS = factor(PGS, levels = ordered_pgs)),
  
  diversity = processed_data %>%
    filter(Dependent == "Medication Diversity") %>%
    mutate(PGS = factor(PGS, levels = ordered_pgs)),
  
  class_diversity = processed_data %>%
    filter(Dependent == "Class Diversity") %>%
    mutate(PGS = factor(PGS, levels = ordered_pgs))
)

# -- Create a function for plot generation
create_plot <- function(data, title, x_offset = 0.01, show_y_axis = TRUE) {
  dodge_width <- 0.8
  
  data$point_color <- ifelse(data$Sig_Nom == "*", "purple", "white")
  
  p <- ggplot(data, aes(y = PGS, x = estimate)) +
    geom_point(position = position_dodge(width = dodge_width), 
               size = 3, fill = data$point_color, shape = 21, color = "purple") +  
    geom_errorbar(aes(xmin = estimate - (1.96 * std.error), 
                      xmax = estimate + (1.96 * std.error)), 
                  position = position_dodge(width = dodge_width), 
                  linewidth = 0.5, width = 0.3, color = "purple") +
    geom_vline(xintercept = 0, linetype = "longdash") +
    xlab(title) +
    theme_classic(base_size = 16) +
    geom_text(aes(label = Sig_label, 
                  y = PGS, 
                  x = (estimate + (1.96 * std.error) + x_offset)),
              color = "black",
              size = 6,
              position = position_dodge(width = dodge_width)) +  
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing = unit(1, "lines"),
      panel.heights = unit(c(5, 8, 12, 10), "null"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(color = "black", face = "bold", size = 12),
      legend.position = "none"
    )
  
  # Configure y-axis display based on parameter
  if (show_y_axis) {
    p <- p + 
      ylab("Polygenic Risk Scores") + 
      theme(
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold")
      )
  } else {
    p <- p + 
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  return(p)
}

# -- Create plots with different configurations
# Calculate appropriate x-offsets based on the range of data
x_offset1 <- max(plot_data$duration$estimate + 1.96 * plot_data$duration$std.error, na.rm = TRUE) * 0.1
x_offset2 <- max(plot_data$diversity$estimate + 1.96 * plot_data$diversity$std.error, na.rm = TRUE) * 0.15
x_offset3 <- max(plot_data$class_diversity$estimate + 1.96 * plot_data$class_diversity$std.error, na.rm = TRUE) * 0.15

p1 <- create_plot(plot_data$duration, "Cumulative AD Dispense", x_offset1, TRUE)
p2 <- create_plot(plot_data$diversity, "AD Diversity", x_offset2, FALSE)
p3 <- create_plot(plot_data$class_diversity, "Class Diversity", x_offset3, FALSE)

# -- Combine plots
plots_bip_excluded <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  align = 'h',
  axis = 'tb',
  rel_widths = c(1.8, 1, 1)
)

# -- Save combined plot
ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\Figure_2_Sensitivity.png", 
       plot = plots_bip_excluded, device = "png", width = 230, height = 140, units = "mm", dpi = 600)
