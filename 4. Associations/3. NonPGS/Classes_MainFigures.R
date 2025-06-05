
# Load required packages
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(patchwork)
library(readxl)
library(tidyr)

#=================== Scaled Quantitative Features =========================================

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table10")

#-- Drug order
drug_order <- c("SNRI", "TeCA", "TCA", "BIP+L", "BIP-L", "Various")

#-- Fill in NA
dat_filled <- dat %>%
  fill(threshold, reference_term, outcome, total_n, .direction = "down")

#-- Filter for test statistics for medications groups in which the reference drug is SSRI
dat_ssri <- dat_filled %>% 
  filter(term != "AGE" & 
           term != "SEX" & 
           term != "SSRI" & 
           reference_term == "SSRI", 
         threshold == "360 days",
         (outcome != "Age" & outcome != "BMI"))

#-- Incorporate Number of Participants for each Outcome Variable in the association
dat_labeled <- dat_ssri %>%
  mutate(Label = paste0(outcome, " (", total_n, ")"))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_labeled %>%
  mutate(
    Sig_label = case_when(
      sig_bonf == "*" ~ "**",
      sig_fdr == "*" & is.na(sig_bonf) ~ "*",
      TRUE ~ NA
    )
  )

#-- Order Dependent Variables
dat_reordered <- dat_labeled %>%
  mutate(Outcome_label = factor(Label, levels = c("Total Prescription Dispense (6102)", 
                                            "Average Prescription Episode Length (6102)",
                                            "Number of Prescription Episodes (6102)", 
                                            "Number of Unique Antidepressants (6102)",
                                            "Number of Unique AD Class (6102)",
                                            "Age of MDD Onset (8976)",
                                            "Times 2-Weeks of MDD (8948)",
                                            "Education Level (9817)",
                                            "Physical health (8746)")))


#-- Put the Dependent Variables into Two Categories
dat_category <- dat_reordered %>%
  mutate(
    Category = ifelse(Outcome_label %in% 
                        c("Age of MDD Onset (8976)", "Times 2-Weeks of MDD (8948)", "Education Level (9817)", "Physical health (8746)"), "Clinical","Pharmaceutical Metrics")
  )


#-- Process the data for the heatmap
heatmap_data_scaled <- dat_category %>%
  mutate(
    # Create a label with estimate and significance
    label = ifelse(!is.na(Sig_label), paste0(estimate, Sig_label), NA),
    Term = factor(term, levels = drug_order)
  )

create_label_breaks <- function(x) {
  # Create temporary copies without the (N=XXXX) part
  temp_x <- x
  # Extract (N=XXXX) parts to add back later
  n_counts <- regmatches(x, regexpr("\\(N=[0-9]+\\)", x))
  
  # Remove counts temporarily for easier text manipulation
  temp_x <- gsub(" \\(N=[0-9]+\\)", "", temp_x)
  
  # Apply specific breaks based on content
  temp_x <- gsub("Age of MDD Onset", "Age of\nMDD Onset", temp_x)
  temp_x <- gsub("Times 2-Weeks of MDD", "Times 2-Weeks\nof MDD", temp_x)
  temp_x <- gsub("Education Level", "Education\nLevel", temp_x)
  temp_x <- gsub("Physical health", "Physical\nhealth", temp_x)
  
  temp_x <- gsub("Total Prescription Dispense", "Total Prescription\nDispense", temp_x)
  temp_x <- gsub("Average Prescription Episode Length", "Average Prescription\nEpisode Length", temp_x)
  temp_x <- gsub("Number of Prescription Episodes", "Number of\nPrescription Episodes", temp_x)
  temp_x <- gsub("Number of Unique Antidepressants", "Number of Unique\nAntidepressants", temp_x)
  temp_x <- gsub("Number of Unique AD Class", "Number of\nUnique AD Class", temp_x)
  
  # If we extracted N counts, add them back at the end
  if (length(n_counts) > 0) {
    for (i in 1:length(temp_x)) {
      if (i <= length(n_counts)) {
        temp_x[i] <- paste0(temp_x[i], "\n", n_counts[i])
      }
    }
  }
  
  return(temp_x)
}

#-- Create plot
scaled_quantitative <- ggplot(heatmap_data_scaled, aes(y = Label, x = term, fill = estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#0047AB", mid = "white", high = "#DC143C", 
    midpoint = 0
  ) +
  geom_text(aes(label = label),
            color = "black", size = 3) +
  scale_x_discrete(limits = drug_order) +
  scale_y_discrete(labels = create_label_breaks) +
  labs(y = NULL, x = NULL, fill = "Estimate(in SD units)\n(SSRI as reference)") +
  theme_classic() +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(
    axis.text.x = element_text(color="black", face = "bold"),
    axis.text.y = element_text(color = "black", size = 10),
    plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "black"),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.title = element_text(face = "bold")
  )

scaled_quantitative

#=== Unscaled Quantitative Features

#-- Filter for test statistics for medications groups in which the reference drug is SSRI
dat_ssri <- dat_filled %>% 
  filter(term != "AGE" & term != "SEX",
         term != "SSRI" & 
           reference_term == "SSRI" &
           (outcome == "Age" |
              outcome == "BMI") & threshold == "360 days")


#-- Incorporate Number of Participants for each Outcome Variable in the association
dat_labeled <- dat_ssri %>%
  mutate(Label = paste0(outcome, " (", total_n, ")"))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_labeled %>%
  mutate(
    Sig_label = case_when(
      sig_bonf == "*" ~ "**",
      sig_fdr == "*" & is.na(sig_bonf) ~ "*",
      TRUE ~ NA
    )
  )

#-- Order Outcome Variables
dat_reordered <- dat_labeled %>%
  mutate(Label = factor(Label, levels = c("Age (9829)", 
                                          "BMI (9316)")))

#-- Process the data for the heatmap
heatmap_data_unscaled <- dat_reordered %>%
  mutate(
    # Create a label with estimate and significance
    label = ifelse(!is.na(Sig_label), paste0(estimate, Sig_label), NA),
    # Create a new grouping variable for combined plots
    plot_group = case_when(
      outcome %in% c("Age") ~ "Age",
      outcome %in% c("BMI") ~ "BMI",
      TRUE ~ as.character(outcome)
    ),
    term = factor(term, levels = drug_order)
  )

#-- Function to create a plot for a group of outcomes
create_group_plot <- function(data, group) {
  filtered_data <- data %>% filter(plot_group == group)
  
  p <- ggplot(filtered_data, aes(x = term, y = Label, fill = estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#0047AB", mid = "white", high = "#DC143C", 
      midpoint = 0
    ) +
    geom_text(aes(label = label),
              color = "black", size = 3) +
    scale_x_discrete(limits = drug_order) +
    scale_y_discrete(labels = create_label_breaks) +
    labs(title = NULL, y = NULL, x = NULL) +
    facet_grid(plot_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_classic() +
    theme(
      axis.text.x = element_text(color="black", face = "bold"),
      axis.text.y = element_text(color="black",size = 10),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
      legend.position = "none",
      strip.background = element_rect(fill = "white")
    )
  
  return(p)
}

# Get unique plot groups
plot_groups <- unique(heatmap_data_unscaled$plot_group)
num_groups <- length(plot_groups)

# Create plots for each group
plot_list <- list()
for (i in 1:num_groups) {
  plot_list[[i]] <- create_group_plot(heatmap_data_unscaled, plot_groups[i])
}

# Combine with patchwork
unscaled_quantitative<- (plot_list[[1]] + theme(axis.text.x = element_blank())) / (plot_list[[2]] + labs(x="Estimate\n(SSRI as reference)") + theme(axis.title.x = element_text(face = "bold")))
  plot_layout(heights = c(1, 1))


unscaled_quantitative_expanded <- unscaled_quantitative + 
  # Increase right margin - adjust the value (currently 3cm) as needed
  theme(plot.margin = margin(5, 5, 50, 5, "pt"))  # top, right, bottom, left

scaled_quantitative_expanded <- scaled_quantitative + 
  # Increase right margin - adjust the value (currently 3cm) as needed
  theme(plot.margin = margin(5, 5, 130, 5, "pt"))  # top, right, bottom, left

combined_plot <- plot_grid(
  unscaled_quantitative_expanded, scaled_quantitative_expanded,
  ncol = 1,
  rel_heights = c(2, 5),
  labels = c("A", "C"),
  axis = "lr",
  align = "v"
)

combined_plot
#====================== Unscaled Binary Features

library(dplyr)
library(ggplot2)
library(tidyr)

#-- Load in data
dat_bin <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table8")

#-- Fill in NA
dat_filled <- dat_bin %>%
  fill(Threshold, Reference, Dependent, Total_N, .direction = "down")

#-- Drug order
drug_order <- c("SNRI", "TeCA", "TCA", "BIP+L", "BIP-L", "Various")

#-- Filter for test statistics for medications groups in which the reference drug is SSRI
dat_ssri <- dat_filled %>% 
  filter(Term != "AGE" & 
           Term != "SEX1" & 
           Term != "SSRI" & 
           Reference == "SSRI", 
         Threshold == 360)
  
cases <- dat_ssri %>%
  group_by(Dependent) %>%
  summarise(
    Total_Controls = sum(Class_N_0),
    Total_Cases = sum(Class_N_1)
  ) %>%
  mutate(
    Case_perc = round((Total_Cases/(Total_Cases + Total_Controls)*100),0)
  )

dat_ssri <- dat_ssri %>%
  left_join(cases, by = "Dependent")

#-- Incorporate Number of Participants for each Outcome Variable in the association
dat_labeled <- dat_ssri %>%
  mutate(
    Dependent = ifelse(Dependent == "Migraines or Headaches", "Migraines", 
                       ifelse(Dependent == "Chronic Fatigue Syndrome", "CFS", 
                              ifelse(Dependent == "Premenstrual Dysphoric Disorder", "PMDD",
                                     ifelse(Dependent == "Obsessive-compulsive Disorder", "OCD",
                                            ifelse(Dependent == "Seasonal Affective Disorder", "SAD", Dependent)))))
  ) %>% 
  mutate(Label = paste0(Dependent, " (", Total_N, ";", Case_perc,  "%)"))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_labeled %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  )


#-- Put the Outcome Variables into Two Categories
dat_category <- dat_labeled %>%
  mutate(
    Category = ifelse(Dependent %in% 
                        c("Sex",
                          "Self-report Responders",
                          "Self-report Discontinuation", 
                          "Suicidal Ideation",
                          "Regular Smoker"), "Clinical", 
                      ifelse(Dependent %in% 
                               c("Type 2 Diabetes", 
                                 "Stomach Ulcers",
                                 "Migraines",
                                 "Back pain",
                                 "Chronic pain",
                                 "CFS",
                                 "Epilepsy",
                                 "Endometriosis",
                                 "Fibroids (uterus)",
                                 "PCOS"), "Physical", 
                             ifelse(Dependent %in% 
                                      c("PMDD",
                                        "Anxiety Disorder",
                                        "Personality Disorder",
                                        "Substance Use Disorder",
                                        "ADHD",
                                        "OCD",
                                        "SAD"), "Psychiatric", "Worst Episode Symptom")))
  )

#-- Process the data for the heatmap
heatmap_data <- dat_category %>%
  mutate(
    # Create a label with estimate and significance
    label = ifelse(!is.na(Sig_label), paste0(OR, Sig_label), NA),
    Term = factor(Term, levels = drug_order)
  )


# Create the main drugs plot with a reasonable scale
main_plot <- ggplot(heatmap_data, aes(x = Term, y = Label, fill = OR)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#0047AB", mid = "white", high = "#DC143C", 
    midpoint = 1,
    limits = c(0, 5)  # Limit the scale to make differences visible
  ) +
  geom_text(aes(label = label),
            color = "black", size = 3) +
  labs(title = NULL, 
       y = NULL, 
       x = NULL, 
       fill = "OR\n(SSRI as reference)") +
  theme_classic() +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(
    axis.text.x = element_text(face = "bold", color="black"),
    axis.text.y = element_text(color = "black",size = 10),
    plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "black"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "white")
  )


final_plot <- plot_grid(
  combined_plot, main_plot,
  ncol = 2,
  rel_widths = c(1, 1.25),
  labels = c("", "B")
)

final_plot


ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\3. NonPGS\\Results\\Figures\\Class_Heatmap.png", 
       plot = final_plot, device = "png", width = 250, height = 250, units = "mm", dpi=1000)

