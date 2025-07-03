# -- Load R libraries
library(tidyverse)
library(cowplot)
library(readxl)
library(viridis)

# -- Read in association results
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet = "Table3")

# -- Pre-processing: fill missing values, filter, and categorize in a single pipeline
processed_data <- dat %>%
  # Fill in NA values
  fill(Dependent, Total_N,Independent, .direction = "down") %>%
  # Filter relevant terms
  filter(!Term %in% c("Intercept")) %>%
  filter(Term != "Age" & Term != "Sex" | 
           (Term == "Age" & Independent == "Age") | 
           (Term == "Sex" & Independent == "Sex")) %>%
  # Select needed columns  
  select(Dependent, Independent, Total_N, Estimate, `Std. Error`, `Pr(>|t|)`, FDR_P, Bonf_P, Sig_FDR, Sig_Bonf, Total_Cases, Total_Controls)

cases <- dat %>%
  filter(!is.na(Total_N) & !is.na(Total_Cases) & Dependent == "Cumulative Prescription Dispense (days)") %>%
  select(Independent,Total_Cases, Total_Controls) %>%
  mutate(
    Cases_perc = round(Total_Cases/(Total_Cases+Total_Controls)*100, 0)
  ) %>%
  select(-Total_Controls, -Total_Cases)
  

processed_data <- processed_data %>%
  left_join(cases, by = c("Independent")) %>%
  # Create label with participant count
  mutate(
    Label = paste0(Independent, " (", Total_N, ";", Cases_perc, "%)"),
    # Significance labeling
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA_character_
    ),
    # Categorize variables
    Category = case_when(
      Independent %in% c("Age", "Sex", "Age of MDD Onset", "Recurrent 2-Weeks of MDD", 
                         "Body Mass Index", "Suicidal Ideation", "Self harm", "Education level", 
                         "Physical health", "Regular Smoker", "Drinks over 3 Months") ~ "Risk Factors",
      Independent %in% c("Type 2 Diabetes", "Stomach Ulcers", "Back pain", "Chronic pain", 
                         "Migraines or Headaches", "Endometriosis", "Fibroids (uterus)", 
                         "PCOS", "Chronic Fatigue Syndrome", "Epilepsy") ~ "Physical",
      Independent %in% c("Anxiety Disorder", "Bipolar Disorder", "Schizophrenia", "Anorexia Nervosa", "Personality Disorder", 
                         "Substance Use Disorder", "ADHD", "Obsessive-compulsive Disorder",
                         "Seasonal Affective Disorder", "Chronic pain", 
                         "Premenstrual Dysphoric Disorder") ~ "Neuropsychiatric",
      TRUE ~ "Worst MDD Episode"
    )
  ) %>%
  # Convert category to factor with specified levels
  mutate(Category = factor(Category, levels = c("Risk Factors", "Physical", "Neuropsychiatric", "Worst MDD Episode")))

# -- Function to order labels by estimate within each category
order_labels <- function(data) {
  # Get ordered labels for each category
  ordered_labels <- data %>%
    filter(Dependent == "Cumulative Prescription Dispense (days)") %>%
    group_by(Category) %>%
    arrange(Estimate) %>%
    pull(Label)
  
  # Return unique labels in the order they appeared
  return(unique(ordered_labels))
}

# Get ordered labels
ordered_labels <- order_labels(processed_data)

# -- Split data by dependent variable for plotting
plot_data <- list(
  TotalPrescriptionDays = processed_data %>% 
    filter(Dependent == "Cumulative Prescription Dispense (days)") %>%
    mutate(Label = factor(Label, levels = ordered_labels)),
  NumberOfATC = processed_data %>% 
    filter(Dependent == "Medication Diversity") %>%
    mutate(Label = factor(Label, levels = ordered_labels)),
  NumberADClass = processed_data %>% 
    filter(Dependent == "Class Diversity") %>%
    mutate(Label = factor(Label, levels = ordered_labels))
)

# -- Function to create a plot
create_plot <- function(data, title, x_limits, show_y_axis = TRUE, show_facet_labels = TRUE) {
  BASE_SIZE <- 10
  dodge_width <- 0.8
  
  # Base plot
  p <- ggplot(data, aes(y = Label, x = Estimate, color = Category)) +
    geom_point(position = position_dodge(width = dodge_width), size = 8, shape = 18) +  
    geom_errorbar(aes(xmin = Estimate - (1.96 * `Std. Error`), 
                      xmax = Estimate + (1.96 * `Std. Error`)), 
                  position = position_dodge(width = dodge_width), 
                  linewidth = 1.5, width = 0.3) +
    geom_vline(xintercept = 0, linetype = "longdash") +
    xlab(title) +
    theme_classic(base_size = BASE_SIZE) +
    scale_color_viridis(option = "D", discrete = TRUE) +
    geom_text(aes(label = Sig_label, 
                  y = Label, 
                  x = (Estimate + (1.96 * `Std. Error`) + x_limits[2]/10)),
              color = "black",
              size = 12,
              position = position_dodge(width = dodge_width)) +
    xlim(x_limits) +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.heights = unit(c(5, 8, 12, 10), "null"),
      axis.text.x = element_text(color = "black", size = 18),
      axis.title.x = element_text(color = "black", size = 17, face = "bold"),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  # Add facet_grid with or without switch
  if (show_facet_labels) {
    p <- p + facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme(
        strip.text.y = element_text(size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        strip.background = element_rect(fill = "white")
      )
  } else {
    p <- p + facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
      theme(
        strip.text = element_blank(),
        strip.background = element_blank()
      )
  }
  
  # Modify y-axis display
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  return(p)
}

# -- Create plots using the function
p1 <- create_plot(plot_data$TotalPrescriptionDays, "Cumulative AD Dispense", c(-130, 270), TRUE, TRUE)
p2 <- create_plot(plot_data$NumberOfATC, "AD Diversity", c(-0.25, 0.8), FALSE, FALSE)
p3 <- create_plot(plot_data$NumberADClass, "Class Diversity", c(-0.25, 0.4), FALSE, FALSE)

# -- Combine plots
plots <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  align = 'h',
  axis = 'tb',
  rel_widths = c(2.5, 1, 1)
)

# -- Save combined plot
ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\1. TreatmentDuration\\Results\\PrescriptionDuration_Associations.png", 
       plot = plots, device = "png", width = 400, height = 500, units = "mm", dpi = 600)