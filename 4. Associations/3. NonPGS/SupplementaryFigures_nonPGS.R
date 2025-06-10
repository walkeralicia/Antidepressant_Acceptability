
#-- Load R libraries
library(dplyr)
library(ggplot2)
library(patchwork)

custom_colors <- c("360 days" = "#4169E1",
                   "600 days" = "#004400")

#====================== Pharmaceutical Plot =========================================

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table11") %>%
  fill(threshold, reference_term, outcome, total_n, .direction = "down")

dat <- dat %>% filter(term != "AGE" & term != "SEX" & reference_term == "SSRI:Sertraline" & term != "SSRI:Sertraline",
                        outcome %in% c("Average Prescription Episode Length", "Number of Prescription Episodes",
                                       "Total Prescription Dispense", "Number of Unique Antidepressants", "Number of Unique AD Class"))

#-- Renaming treatment and PGS labels
dat$term <- factor(dat$term, levels = rev(c("TeCA:Mirtazapine", "TCA:Amitriptyline", "SNRI:Duloxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SSRI:Paroxetine", "SSRI:Fluoxetine", "SSRI:Escitalopram", "SSRI:Citalopram")))
dat$threshold <- factor(dat$threshold, levels = c("360 days", "600 days"))

variable_order <- c("Total Prescription Dispense",
                    "Number of Prescription Episodes",
                    "Average Prescription Episode Length",
                    "Number of Unique Antidepressants",
                    "Number of Unique AD Class")
dat_factored <- dat %>%
  mutate(outcome = factor(outcome, levels = variable_order))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_factored %>%
  mutate(
    Sig_label = case_when(
      sig_bonf == "*" ~ "**",
      sig_fdr == "*" & is.na(sig_bonf) ~ "*",
      TRUE ~ NA
    )
  )


# Plotting code
dodge_width <- 0.6
p <- ggplot(dat_labeled, aes(x = term, y = estimate, group = threshold, colour = threshold)) +
  geom_point(position = position_dodge(width = dodge_width), size = 5) +  
  geom_errorbar(aes(ymin = estimate - (1.96 *  std.error), ymax = estimate + (1.96 * std.error)), 
                position = position_dodge(width = dodge_width), 
                linewidth = 1.5, width = 0) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = Sig_label, 
                x = term, 
                y = (estimate + (1.96 * std.error) + 0.01)), 
            size = 10, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~outcome, ncol = 3, scales = "free_y") +
  xlab("Treatment group") +
  ylab("Estimate difference from Sertraline") +
  theme_bw(base_size = 27) +
  scale_color_manual(values = custom_colors) +
  theme(
    strip.text = element_text(color = "black"), 
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.7),
    axis.text.y = element_text(color = "black"),
    legend.position = "top"
  )

p

ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\3. NonPGS\\Results\\Figures\\PharmaceuticalMetrics_QuantitativeAssociations.png", 
       plot = p, device = "png", width = 500, height = 400, units = "mm")



#======================== Quantitative survey phenotypes =================================

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table11") %>%
  fill(threshold, reference_term, outcome, total_n, .direction = "down")

dat <- dat %>% filter(term != "AGE" & term != "SEX" & reference_term == "SSRI:Sertraline" & term != "SSRI:Sertraline",
                      outcome %in% c("Age", "Age of MDD Onset", "Times 2-Weeks of MDD", "BMI", "Education Level",
                                     "Physical health", "Drinks over 3 Months"))

#-- Renaming treatment and PGS labels
dat$term <- factor(dat$term, levels = rev(c( "Various", "Combination", "BIP-L", "BIP+L", "TeCA:Mirtazapine", "TCA:Amitriptyline", "SNRI:Duloxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SSRI:Paroxetine", "SSRI:Fluoxetine", "SSRI:Escitalopram", "SSRI:Citalopram")))
dat$threshold <- factor(dat$threshold, levels = c("360 days", "600 days"))

variable_order <- c("Age",
                    "Age of MDD Onset",
                    "Times 2-Weeks of MDD",
                    "BMI",
                    "Education Level",
                    "Physical health",
                    "Drinks over 3 Months")
dat_factored <- dat %>%
  mutate(outcome = factor(outcome, levels = variable_order))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_factored %>%
  mutate(
    Sig_label = case_when(
      sig_bonf == "*" ~ "**",
      sig_fdr == "*" & is.na(sig_bonf) ~ "*",
      TRUE ~ NA
    )
  )


# Plotting code
dodge_width <- 0.5
p <- ggplot(dat_labeled, aes(x = term, y = estimate, group = threshold, colour = threshold)) +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +  
  geom_errorbar(aes(ymin = estimate - (1.96 *  std.error), ymax = estimate + (1.96 * std.error)), 
                position = position_dodge(width = dodge_width), 
                linewidth = 1, width = 0) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = Sig_label, 
                x = term, 
                y = (estimate + (1.96 * std.error) + 0.01)), 
            size = 10, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~outcome, nrow = 2, scales = "free_y") +
  xlab("Treatment group") +
  ylab("Estimate difference from Sertraline") +
  theme_bw(base_size = 27) +
  scale_color_manual(values = custom_colors) +
  theme(
    strip.text = element_text(color = "black"), 
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.7),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

p
ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\3. NonPGS\\Results\\Figures\\ClinicalFeatures_QuantitativeAssociations.png", 
       plot = p, device = "png", width = 500, height = 400, units = "mm")


#=================== Binary Worst MDD Episode Symptoms ===============================

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table9") %>%
  fill(Threshold, Reference, Dependent, Total_N, .direction = "down")

variable_order <- c("Low Interest 2-Weeks",
                    "Depressed 2-Weeks",
                    "Fatigued",
                    "Guilty Feelings",
                    "No Focus",
                    "Death Thoughts",
                    "Appetite/Weight Change",
                    "Sleep Disturbances",
                    "Movement Changes")

dat <- dat %>% filter(Term != "AGE" & Term != "SEX1" & Reference == "SSRI:Sertraline" & Term != "SSRI:Sertraline",
                      Dependent %in% variable_order)

#-- Renaming treatment and PGS labels
dat$Term <- factor(dat$Term, levels = rev(c( "Various", "Combination", "BIP-L", "BIP+L", "TeCA:Mirtazapine", "TCA:Amitriptyline", "SNRI:Duloxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SSRI:Paroxetine", "SSRI:Fluoxetine", "SSRI:Escitalopram", "SSRI:Citalopram")))
dat <- dat %>%
  mutate(Threshold = ifelse(Threshold == 360, "360 days", "600 days")) %>%
  mutate(Threshold = factor(Threshold, levels = c("360 days", "600 days")))


dat_factored <- dat %>%
  mutate(Dependent = factor(Dependent, levels = variable_order))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_factored %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  )


# Plotting code
dodge_width <- 0.7
p <- ggplot(dat_labeled, aes(x = Term, y = OR, group = Threshold, colour = Threshold)) +
  geom_point(position = position_dodge(width = dodge_width), size = 1.5) +  
  geom_errorbar(aes(ymin = LCI, ymax = UCI), 
                position = position_dodge(width = dodge_width), 
                linewidth = 0.3, width = 0) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  geom_text(aes(label = Sig_label, 
                x = Term, 
                y = UCI + 0.01), 
            size = 7, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~Dependent, nrow = 3, scales = "free_y") +
  xlab("Treatment group") +
  ylab("Odds Ratio (95% CI) with SSRI:Sertraline as Reference") +
  scale_color_manual(values = custom_colors) +
  theme(
    strip.text = element_text(color = "black"), 
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.7),
    axis.text.y = element_text(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", face = "bold")
  )

p
ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\3. NonPGS\\Results\\Figures\\MDDSymptoms_BinaryAssociations.png", 
       plot = p, device = "png", width = 200, height = 200, units = "mm", dpi = 600)


############################### Clinical Binary Survey Phenotypes ############################################

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table9") %>%
  fill(Threshold, Reference, Dependent, Total_N, .direction = "down")


variable_order <- c("Sex",
                    "Self-report Responders",
                    "Self-report Discontinuation", 
                    "Suicidal Ideation",
                    "Regular Smoker")

dat <- dat %>% filter(Term != "AGE" & Term != "SEX1" & Reference == "SSRI:Sertraline" & Term != "SSRI:Sertraline",
                      Dependent %in% variable_order)

#-- Renaming treatment and PGS labels
dat$Term <- factor(dat$Term, levels = rev(c( "Various", "Combination", "BIP-L", "BIP+L", "TeCA:Mirtazapine", "TCA:Amitriptyline", "SNRI:Duloxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SSRI:Paroxetine", "SSRI:Fluoxetine", "SSRI:Escitalopram", "SSRI:Citalopram")))
dat <- dat %>%
  mutate(Threshold = ifelse(Threshold == 360, "360 days", "600 days")) %>%
  mutate(Threshold = factor(Threshold, levels = c("360 days", "600 days")))


dat_factored <- dat %>%
  mutate(Dependent = factor(Dependent, levels = variable_order))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_factored %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  )


# Plotting code
dodge_width <- 0.7
p <- ggplot(dat_labeled, aes(x = Term, y = OR, group = Threshold, colour = Threshold)) +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +  
  geom_errorbar(aes(ymin = LCI, ymax = UCI), 
                position = position_dodge(width = dodge_width), 
                linewidth = 1, width = 0) +
  geom_hline(yintercept = 1) +
  geom_text(aes(label = Sig_label, 
                x = Term, 
                y = UCI + 0.01), 
            size = 10, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~Dependent, nrow = 2, scales = "free_y") +
  xlab("Treatment group") +
  ylab("Odds Ratio (95% CI) with Sertraline as Reference") +
  theme_bw(base_size = 27) +
  scale_color_manual(values = custom_colors) +
  theme(
    strip.text = element_text(color = "black"), 
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.7),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )


ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\3. NonPGS\\Results\\Figures\\ClinicalFeatures_BinaryAssociations.png", 
       plot = p, device = "png", width = 500, height = 400, units = "mm")

#================= Binary Comorbidities ================================================

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table9") %>%
  fill(Threshold, Reference, Dependent, Total_N, .direction = "down")


variable_order <- c("Type 2 Diabetes", 
                    "Stomach Ulcers",
                    "PDMD",
                    "Anxiety Disorder",
                    "Personality Disorder",
                    "Schizophrenia",
                    "Anorexia Nervosa",
                    "Substance Use Disorder",
                    "ADHD",
                    "Obsessive-compulsive Disorder",
                    "Seasonal Affective Disorder",
                    "Premenstrual Dysphoric Disorder",
                    "Back pain",
                    "Chronic Fatigue Syndrome",
                    "Epilepsy",
                    "Chronic pain",
                    "Migraines or Headaches",
                    "Endometriosis",
                    "PCOS",
                    "Fibroids (uterus)")

dat <- dat %>% filter(Term != "AGE" & Term != "SEX1" & Reference == "SSRI:Sertraline" & Term != "SSRI:Sertraline",
                      Dependent %in% variable_order)

#-- Renaming treatment and PGS labels
dat$Term <- factor(dat$Term, levels = rev(c( "Various", "Combination", "BIP-L", "BIP+L", "TeCA:Mirtazapine", "TCA:Amitriptyline", "SNRI:Duloxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SSRI:Paroxetine", "SSRI:Fluoxetine", "SSRI:Escitalopram", "SSRI:Citalopram")))
dat <- dat %>%
  mutate(Threshold = ifelse(Threshold == 360, "360 days", "600 days")) %>%
  mutate(Threshold = factor(Threshold, levels = c("360 days", "600 days")))


dat_factored <- dat %>%
  mutate(Dependent = factor(Dependent, levels = variable_order))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat_factored %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  )


# Plotting code
dodge_width <- 0.7
p <- ggplot(dat_labeled, aes(x = Term, y = OR, group = Threshold, colour = Threshold)) +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +  
  geom_errorbar(aes(ymin = LCI, ymax = UCI), 
                position = position_dodge(width = dodge_width), 
                linewidth = 1, width = 0) +
  geom_hline(yintercept = 1) +
  geom_text(aes(label = Sig_label, 
                x = Term, 
                y = UCI + 0.01), 
            size = 10, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~Dependent, nrow = 7, scales = "free_y") +
  xlab("Treatment group") +
  ylab("Odds Ratio (95% CI) with Sertraline as Reference") +
  theme_bw(base_size = 27) +
  scale_color_manual(values = custom_colors) +
  theme(
    strip.text = element_text(color = "black"), 
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.7),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )


ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\3. NonPGS\\Results\\Figures\\Comorbidities_BinaryAssociations.png", 
       plot = p, device = "png", width = 500, height = 650, units = "mm")


