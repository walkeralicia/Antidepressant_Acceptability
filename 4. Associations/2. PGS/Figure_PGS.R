
#----------------- PGS ~ ATC Class -------------------------------------------------------------------------
options(scipen = 999)

#-- Load R libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(readxl)
library(viridis)
dodge_width <- 0.8

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table13")

#-- Drug order
source("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\Drug_Reference_Table.R")
drug_order <- intersect(c(unique(drug_ref$DrugClass), "BIP+L", "BIP-L", "Various"), dat$Term)

#-- Fill in NA
dat_filled <- dat %>%
  fill(Threshold, Reference, PGS, Total_N, .direction = "down")

#-- Filter for test statistics for medications groups in which the reference drug is SSRI
dat_ssri <- dat_filled %>% 
  filter(Term != "AIIa" & 
        Reference == "AIIa", 
         Threshold == "360 days") %>%
  mutate(PGS = ifelse(PGS == "MDD", "MD", PGS))

dat_ssri$PGS <- factor(dat_ssri$PGS, levels = rev(c("MD", "BIP", "SCZ", "ADHD", "ANO", "ANX", "OCD", "Neuroticism",
                                                "BMI", "T2D", "SBP", "CNT", "LRA", "Migraine", "PUD", "CRP", "LDL-c")))
dat_ssri$Term <- factor(dat_ssri$Term, levels = drug_order)
                       

dat_ssri <- dat_ssri %>%
  mutate(PGS_label = paste0(PGS, " (N=", Total_N, ")")) %>%
  mutate(Med_label = paste0(Term, "\n(N=", Group_N, ")"))


# Create a named vector 
med_labels <- setNames(
  dat_ssri$Med_label,
  dat_ssri$Term
)
# Remove duplicates to get a clean mapping
med_labels <- med_labels[!duplicated(names(med_labels))]

#-- Label for FDR or Bonferroni significance
dat_ssri <- dat_ssri %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  ) %>%
  mutate(
    P.value = as.numeric(P.value)
  )


p <- ggplot(dat_ssri, aes(y = PGS, x = estimate, color = P.value)) +
  geom_point(size = 3, shape = 18) +  
  geom_errorbar(aes(xmin = estimate - (1.96 * `Std..Error`), xmax = estimate + (1.96 * `Std..Error`)), 
                linewidth = 0.5, width = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", linewidth = 0.5, color = "#333333") +
  geom_text(aes(label = Sig_label, 
                y = PGS, 
                x = (estimate + (1.96 * `Std..Error`))),
            color = "black",
            size = 5,
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~Term, nrow = 1, labeller = labeller(Term = med_labels)) +
  ylab("PGS Trait") +
  xlab("Estimate (SD units; AIIa as reference [N=3,573])") +
  theme_classic(base_size = 12) +
  scale_color_viridis(option = "D") +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black", face="bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p

ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\2. PGS\\Results\\AGDS_360days_PGS_Associations_Reference_SSRI.png", 
       plot = p, device = "png", width = 250, height = 120, units = "mm", dpi = 1000)



#----------------- PGS ~ ATC Drug -------------------------------------------------------------------------

custom_colors <- c("360 days" = "#4169E1",
                   "600 days" = "#004400")

#-- Load R libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(readxl)
library(viridis)
dodge_width <- 0.8

#-- Read in results in which the dependent variable is quantitative
dat <- read_excel("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\All_results.xlsx", sheet = "Table14")

#-- Fill in NA
dat_filled <- dat %>%
  fill(Threshold, Reference, PGS, Total_N, .direction = "down")

#-- Filter for test statistics for medications groups in which the reference drug is SSRI
dat_ssri <- dat_filled %>% 
  filter(Term != "AIIa:Sertraline" & 
           Reference== "AIIa:Sertraline")

dat_ssri$PGS <- factor(dat_ssri$PGS, levels = rev(c("MDD", "BIP", "SCZ", "ADHD", "ANO", "ANX", "OCD", "Neuroticism",
                                                    "BMI", "T2D", "SBP", "CNT", "LRA", "Migraine", "PUD", "CRP", "LDL-c")))
dat_ssri$Term <- factor(dat_ssri$Term, levels = rev(c( "Various", "Combination", "BIP-L", "BIP+L", drug_ref$DrugName)))
dat_ssri$Threshold <- factor(dat_ssri$Threshold, levels = c("360 days", "600 days"))

#-- Label for FDR or Bonferroni significance
dat_ssri <- dat_ssri %>%
  mutate(
    Sig_label = case_when(
      Sig_Bonf == "*" ~ "**",
      Sig_FDR == "*" & is.na(Sig_Bonf) ~ "*",
      TRUE ~ NA
    )
  ) %>%
  mutate(
    P.value = as.numeric(P.value)
  )



# Plotting code
dodge_width <- 0.5
p <- ggplot(dat_ssri, aes(x = Term, y = estimate, group = Threshold, colour = Threshold)) +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +  
  geom_errorbar(aes(ymin = estimate - (1.96 *  `Std..Error`), ymax = estimate + (1.96 * `Std..Error`)), 
                position = position_dodge(width = dodge_width), 
                linewidth = 1, width = 0) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = Sig_label, 
                x = Term, 
                y = (estimate + (1.96 * `Std..Error`) + 0.01)), 
            size = 10, color = "black", 
            position = position_dodge(width = dodge_width)) +  
  facet_wrap(~PGS, nrow = 5, scales = "free_y") +
  xlab("Treatment group") +
  ylab("Estimate difference from AIIa:Sertraline") +
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

ggsave("C:\\Users\\walkera\\OneDrive - Nexus365\\Documents\\PhD\\AGDS\\Antidepressant_Acceptability\\4. Associations\\2. PGS\\Results\\Supp_PGS_groups.png", 
       plot = p, device = "png", width = 500, height = 650, units = "mm")


