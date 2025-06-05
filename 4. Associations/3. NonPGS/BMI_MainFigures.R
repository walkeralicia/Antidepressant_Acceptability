
library(ggplot2)
library(dplyr)

custom_colors <- c("No" = "black",
                   "Yes" = "#2C7FB8")

#-- Antidepressant Treatment Class Level analysis

#-- BMI unadjusted
bmi_unadj <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet = "Table10")
bmi_unadj <- bmi_unadj %>%
  rename(
    Threshold = threshold, 
    ReferenceTerm = reference_term, 
    Outcome = outcome, 
    Term = term, 
    Estimate = estimate, 
    `Std..Error` = std.error,
    `t.value` = statistic,
    `Pr...t..` = p.value,
    Sig_FDR = sig_fdr,
    Sig_Bonf = sig_bonf
  ) %>%
  fill(Outcome, Threshold, ReferenceTerm, .direction = "down") %>% 
  filter(Term != "AGE" & Term != "SEX") %>%
  filter(Threshold == "360 days") %>%
  filter(Outcome == "BMI") %>%
  filter(ReferenceTerm == "SSRI") %>%
  select(Outcome, ReferenceTerm, Term, Estimate, `Std..Error`, `t.value`, `Pr...t..`, Sig_FDR, Sig_Bonf) %>%
  mutate(Adjusted = "No") %>%
  filter(Term != "SSRI")

#-- BMI adjusted by BMI PGS
bmi_adj <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_Results.xlsx", sheet = "Table14")
bmi_adj <- bmi_adj %>%
  fill(Outcome, Threshold, ReferenceTerm, .direction = "down") %>% 
  filter(Term != "sd_pgs" & Term != "(Intercept)") %>%
  filter(Term != "AGE" & Term != "SEX") %>%
  filter(Threshold == "360 days") %>%
  filter(Outcome == "BMI") %>%
  filter(ReferenceTerm == "SSRI") %>%
  filter(Term != "SSRI") %>%
  select(Outcome, ReferenceTerm, Term, Estimate, `Std..Error`, `t.value`, `Pr...t..`, Sig_FDR, Sig_Bonf) %>%
  mutate(Adjusted = "Yes")

#-- Combine BMI adjusted and unadjusted analysis into one dataframe
dat <- rbind(bmi_adj, bmi_unadj)

#-- Prepare data for plotting
dat$Term <- factor(dat$Term, levels = rev(c( "Various","BIP-L" ,"BIP+L", "TeCA", "TCA", "SNRI")))
dodge_width <- 0.5

#-- Main Plot 
main_p <- ggplot(dat, aes(x = Term, y = Estimate, group = Adjusted, color = Adjusted)) +
  geom_point(position = position_dodge(width = dodge_width), size = 3, shape = 18) +  
  geom_errorbar(aes(ymin = Estimate - (1.96 * `Std..Error`), ymax = Estimate + (1.96 * `Std..Error`)), 
                position = position_dodge(width = dodge_width), 
                linewidth = 0.5, width = 0) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.5, color = "#333333") +
  xlab("Antidepressant Class") +
  ylab("BMI calculated from self-reported\nheight and weight\n(SSRI as reference)") +
  labs(color = "Covariate\nBMI PGS") +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black", face = "bold"),
    legend.position = "right", 
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

main_p

###################################################


# Combine the quantitative plots

main_p_expanded <- main_p + 
  # Increase right margin - adjust the value (currently 3cm) as needed
  theme(plot.margin = margin(5, 150, 5, 40, "pt"))  # top, right, bottom, left

combined_plot <- plot_grid(
  p, main_p_expanded,
  ncol = 1,
  labels = c("A", "B"),
  rel_heights = c(2, 1.2)
)

combined_plot

ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\2. PGS\\Results\\AGDS_MF_PGS_plots.png", 
       plot = combined_plot, device = "png", width = 230, height = 200, units = "mm", dpi=600)

#####################################################

#-- Antidepressant Treatment Level analysis

custom_colors <- c("No" = "black",
                   "Yes" = "#2C7FB8")

#-- BMI unadjusted
bmi_unadj <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Pharmacogenomics\\All_Results.xlsx", sheet = 9)
bmi_unadj <- bmi_unadj %>%
  fill(Outcome, Threshold, ReferenceTerm, .direction = "down") %>% 
  filter(Term != "AGE" & Term != "SEX") %>%
  filter(Outcome == "BMI") %>%
  filter(ReferenceTerm == "SSRI:Sertraline") %>%
  select(Threshold, Outcome, ReferenceTerm, Term, Estimate, `Std..Error`, `t.value`, `Pr...t..`, Sig_FDR, Sig_Bonf) %>%
  mutate(Adjusted = "No") %>%
  filter(Term != "SSRI:Sertraline")

#-- BMI adjusted by BMI PGS
bmi_adj <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Pharmacogenomics\\All_Results.xlsx", sheet = 13)
bmi_adj <- bmi_adj %>%
  fill(Outcome, Threshold, ReferenceTerm, .direction = "down") %>% 
  filter(Term != "sd_pgs" & Term != "(Intercept)") %>%
  filter(Term != "AGE" & Term != "SEX") %>%
  filter(Outcome == "BMI") %>%
  filter(ReferenceTerm == "SSRI:Sertraline") %>%
  select(Threshold, Outcome, ReferenceTerm, Term, Estimate, `Std..Error`, `t.value`, `Pr...t..`, Sig_FDR, Sig_Bonf) %>%
  mutate(Adjusted = "Yes")

#-- Combine BMI adjusted and unadjusted analysis into one dataframe
dat <- rbind(bmi_adj, bmi_unadj)

#-- Prepare data for plotting
dat$Term <- factor(dat$Term, levels = rev(c( "Various", "Combination", "BIP-L", "BIP+L", "TeCA:Mirtazapine", "TCA:Amitriptyline", "SNRI:Duloxetine", "SNRI:Desvenlafaxine", "SNRI:Venlafaxine", "SSRI:Paroxetine", "SSRI:Fluoxetine", "SSRI:Escitalopram", "SSRI:Citalopram")))
dat$Threshold <- factor(dat$Threshold, levels = c("360 days", "600 days"))

#-- Label for FDR or Bonferroni significance
dat_labeled <- dat %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  )

dodge_width <- 0.5

#-- Main Plot 
main_p <- ggplot(dat_labeled, aes(x = Term, y = Estimate, group = Adjusted, color = Adjusted)) +
  geom_point(position = position_dodge(width = dodge_width), size = 3, shape = 18) +  
  geom_errorbar(aes(ymin = Estimate - (1.96 * `Std..Error`), ymax = Estimate + (1.96 * `Std..Error`)), 
                position = position_dodge(width = dodge_width), 
                linewidth = 0.5, width = 0) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.5, color = "#333333") +
  geom_text(aes(label = Sig_label, 
                x = Term, 
                y = (Estimate + (1.96 * `Std..Error`) + 0.01)), 
            size = 6, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(. ~Threshold, nrow = 2) +
  xlab("Treatment Group") +
  ylab(expression(atop("BMI (SSRI:Sertraline as reference)", paste("(", pm, " 95% CI)")))) +
  labs(color = "Covariate\nBMI PGS") +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black", face = "bold", angle = 60, hjust = 1),
    legend.position = "right", 
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

main_p

ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Pharmacogenomics\\4. Associations\\2. PGS\\Results\\AGDS_SF_PGS_plots.png", 
       plot = main_p, device = "png", width = 150, height = 200, units = "mm", dpi=600)

