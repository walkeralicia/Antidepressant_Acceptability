
#=================== Only Self-Report Characteristics as Predictors =========================

#-- Load R libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(viridis)
library(readxl)

#-- Load data
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet = "Table18") %>%
  fill(Model, y.level, .direction = "down") %>%
  filter(Model == "Self Report Predictors") %>%
  filter(term != "PDMD")

# First, transform your data to better interpret continuous predictors
plot_data <- dat %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Identify variable type
    var_type = case_when(
      grepl("PGS", term) ~ "Continuous",
      term %in% c("Age", "BMI", "Number of MDD Episodes", "Physical Health", 
                  "Education Level") ~ "Continuous",
      TRUE ~ "Binary"
    )
  ) 

# Improved forest plot
p1 <- ggplot(plot_data, aes(x = estimate, y = reorder(term, estimate), color = p.value)) +
  geom_vline(xintercept = 1, linetype = "longdash", color = "black") +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0) +
  facet_wrap(~ y.level, scales = "free_x") +
  labs(
    title = "Multinomial Regression: Self-Report Predictors of Treatment Class",
    subtitle = "Reference: SSRI",
    x = "Odds Ratio (with 95% CI)", 
    y = NULL
  ) +
  theme_classic() +
  scale_color_viridis(option = "D") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black")
  )


#Note: For continuous predictors (_scaled), OR represents\nthe change in odds for a 1 SD increase in the predictor
#==================== Include PGS as Predictors =====================
#-- Load data
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet ="Table18") %>%
  fill(Model, y.level, .direction = "down") %>%
  filter(Model == "All Predictors")

# First, transform your data to better interpret continuous predictors
plot_data <- dat %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Identify variable type
    var_type = case_when(
      grepl("PGS", term) ~ "Continuous",
      term %in% c("Age", "BMI", "Number of MDD Episodes", "Physical Health", 
                  "Education Level", "PUD PGS") ~ "Continuous",
      TRUE ~ "Binary"
    )
  ) 


# Improved forest plot
p3 <- ggplot(plot_data, aes(x = estimate, y = reorder(term, estimate), color = p.value)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0) +
  facet_wrap(~ y.level, scales = "free_x") +
  labs(
    title = "Multinomial Regression: Self-Report and PGS Predictors of Treatment Class",
    subtitle = "Reference: SSRI",
    x = "Odds Ratio (with 95% CI)", 
    y = NULL
  ) +
  theme_classic() +
  scale_color_viridis(option = "D") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black")
  ) 



#============= Only PGS as Predictors ====================

#-- Load data
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet ="Table18") %>%
  fill(Model, y.level, .direction = "down") %>%
  filter(Model == "PGS Predictors")

# First, transform your data to better interpret continuous predictors
plot_data <- dat %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Identify variable type
    var_type = case_when(
      grepl("PGS", term) ~ "Continuous",
      TRUE ~ "Binary"
    )
  ) 



# Improved forest plot
p2 <- ggplot(plot_data, aes(x = estimate, y = reorder(term, estimate), color = p.value)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0) +
  facet_wrap(~ y.level, scales = "free_x") +
  labs(
    title = "Multinomial Regression: PGS Predictors of Treatment Class",
    subtitle = "Reference: SSRI",
    x = "Odds Ratio (with 95% CI)", 
    y = NULL
  ) +
  theme_classic() +
  scale_color_viridis(option = "D") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black")
  )

# Create a combined figure with one column
combined_plot <- plot_grid(
  p1, p2, p3,
  labels = c("A", "B", "C"),
  ncol = 1,
  align = "v",
  axis = "lr", 
  rel_heights = c(1.7, 1, 1.7)
)



# Save the combined figure
ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\4. Multinomial_Regression\\combined_multinomial_regression_models.png", 
       combined_plot, width = 10, height = 10, dpi = 300)

