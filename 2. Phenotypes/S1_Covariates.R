
#========== Load required R libraries =========
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(viridis)
library(ggpubr)
library(gridExtra)
library(forcats)

#=== Read in Phenotypes ========================
pheno <- read.csv("/QRISdata/Q5338/Phenotypes/Combined_AGDS_2020_Freeze_20220620.csv") 
wkdir="/QRISdata/Q7280/pharmacogenomics"

#=== BMI Quality Control ===================================
phys <- pheno %>% 
  select("STUDYID","SEX", "AGE", "HTFTTXT", "HTINTXT", "HTCMTXT", "WTKG") %>%
  mutate(SEX = as.factor(SEX)) %>%
  mutate(AGE = as.numeric(AGE)) %>%
  filter(SEX != "3") %>%
  mutate(SEX =  fct_recode(SEX, "Male" = "1", "Female" = "2")) %>%
  mutate(HTFTCM = ifelse(!is.na(HTFTTXT), HTFTTXT*30.48, NA)) %>%
  mutate(HTINCM = ifelse(!is.na(HTINTXT), HTINTXT*2.54, NA)) %>%
  mutate(HTCM = ifelse(!is.na(HTFTCM), HTFTCM+HTINCM, 
                       ifelse(!is.na(HTCMTXT), HTCMTXT, NA))) %>%
  mutate(BMI_raw = WTKG/(HTCM/100)^2) %>%
  mutate(HTIN = ifelse(SEX == "Female" & HTCM>=140 & HTCM <=190, HTCM, 
                       ifelse(SEX == "Male" & HTCM >= 140 & HTCM <= 210, HTCM, NA))) %>%
  mutate(WTIN = ifelse(SEX == "Female" & WTKG>= 38 & WTKG <=200, WTKG, 
                       ifelse(SEX == "Male" & WTKG >= 40 & WTKG <= 200, WTKG, NA))) %>%
  mutate(BMIIN = WTIN/(HTIN/100)^2) %>%
  mutate(BMI = as.numeric(ifelse(BMIIN >= 14 & BMIIN <= 60, BMIIN, NA)))

write.csv(phys, file.path(wkdir, "phenotypes/BMI_age_sex.csv"), row.names = FALSE)

#== Split data by gender
fphys <- phys %>% filter(SEX == "Female")
mphys <- phys %>% filter(SEX == "Male")

#== Distribution of raw height and weight across gender
p_height_f <- ggplot(fphys, aes(x = HTCM)) +
  geom_density(size = 2, color = "deeppink", fill = "deeppink", alpha = 0.2) + 
  geom_vline(xintercept=140, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  geom_vline(xintercept=190, color='deeppink', linetype='dashed', alpha=2, linewidth=1.5) +
  ggtitle("Female Height") +
  labs(x = "Raw Height (cm)") +
  theme_light(base_size = 50)
p_height_m <- ggplot(mphys, aes(x = HTCM)) +
  geom_density(size = 2, color = "blue", fill = "blue", alpha = 0.2) + 
  geom_vline(xintercept=140, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  geom_vline(xintercept=210, color='blue', linetype='dashed', alpha=2, linewidth=1.5) +
  ggtitle("Male Height") +
  labs(x = "Raw Height (cm)") +
  theme_light(base_size = 50)
p_weight_f <- ggplot(fphys, aes(x = WTKG)) +
  geom_density(size = 2, color = "deeppink", fill = "deeppink", alpha = 0.2) + 
  geom_vline(xintercept=200, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  geom_vline(xintercept=38, color='deeppink', linetype='dashed', alpha=2, linewidth=1.5) +
  ggtitle("Female Weight") +
  labs(x = "Raw Weight (kg)") +
  theme_light(base_size = 50)
p_weight_m <- ggplot(mphys, aes(x = WTKG)) +
  geom_density(size = 2, color = "blue", fill = "blue", alpha = 0.2) + 
  geom_vline(xintercept=200, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  geom_vline(xintercept=40, color='blue', linetype='dashed', alpha=2, linewidth=1.5) +
  ggtitle("Male Weight") +
  labs(x = "Raw Weight (kg)") +
  theme_light(base_size = 50)

phys_gender_plots <- ggarrange(p_height_f, p_height_m,
                          p_weight_f, p_weight_m,
                       ncol = 2, nrow = 2, widths=c(1,1,1,1))

#== Height and weight values for those with extreme BMI values
out <- phys %>% filter(BMIIN >= 60) # 44 individuals
out$BMIIN <- as.numeric(as.character(out$BMIIN))

p_out <- ggplot(out, aes(x= WTIN, y = HTIN, color = BMIIN)) +
  geom_point(size = 5, aes(shape = SEX)) +
  ggtitle(expression("BMI of at least 60 kg/m"^2*"")) + 
  labs(x = "Weight (kg)", y = "Height (cm)") +
  scale_color_gradient(low = "green", high = "blue") +
  theme_light(base_size = 35) +
  guides(color=guide_legend(title="BMI"))
png(file.path(wkdir, "phenotypes/BMI_QC/Qced_biplot_extremeBMI.png"), bg="white", type=c("cairo"),width=15, height=10, units='in', res = 500)
p_out
dev.off()

#--- Define ancestries
ancestries <- c("ANCEST1", "ANCEST2", "ANCEST3", "ANCEST4", "ANCEST5", "ANCEST6", "ANCEST7", "ANCEST8", "ANCEST9", "ANCEST10",
                "ANCEST11", "ANCEST12", "ANCEST13", "ANCEST14", "ANCEST15", "ANCEST16", "ANCEST17", "ANCEST18", "ANCEST20")

#--- Select ID, AGE, SEX and ancestries
pheno <- read.csv(file.path(wkdir, "/data/Combined_AGDS_2020_Freeze_20220620.csv"))
pheno_select <- pheno %>%
  select(STUDYID, AGE, SEX, all_of(ancestries)) %>%
  mutate(SEX = case_when(
    SEX == 2 ~ "Female",
    SEX == 1 ~ "Male"
  ))

#--- Select the first non-NA ancestry value for each individual
agds_anc <- c()
for (r in 1:nrow(pheno_select)){
  print(r/nrow(pheno_select))
  anc <- pheno_select[r,ancestries]
  if (length(which(!is.na(anc))) == 0){
    anc_first <- NA
  } else {
    anc_idx <- which(!is.na(anc))[1]
    anc_first <- ancestries[anc_idx]
  }
  agds_anc <- c(agds_anc, anc_first)
}

#--- Map ancestry codes to names
pheno_ancestry <- pheno_select %>%
  mutate(CODE = agds_anc) %>%
  mutate(ANCESTRY = case_when(
    CODE == "ANCEST1" ~ "England, Ireland, Scotland or Wales",
    CODE == "ANCEST2" ~ "Australia - not of Aboriginal or Torres Strait Islander descent",
    CODE == "ANCEST3" ~ "Australia - of Aboriginal or Torres Strait Islander descent",
    CODE == "ANCEST4" ~ "New Zealand - not of Maori descent",
    CODE == "ANCEST5" ~ "New Zealand - of Maori descent",
    CODE == "ANCEST6" ~ "Northern Europe including Sweden, Norway, Finland and surrounding countries",
    CODE == "ANCEST7" ~ "Western Europe including France, Germany, the Netherlands and surrounding countries",
    CODE == "ANCEST8" ~ "Southern Europe including Italy, Greece, Spain, Portugal and surrounding countries",
    CODE == "ANCEST9" ~ "Eastern Europe including Russia, Poland, Hungary and surrounding countries",
    CODE == "ANCEST10" ~ "Middle East including Lebanon, Turkey and surrounding countries",
    CODE == "ANCEST11" ~ "Eastern Asia including China, Japan, South Korea, North Korea, Taiwan and Hong Kong",
    CODE == "ANCEST12" ~ "South-East Asia including Thailand, Malaysia, Indonesia, Singapore and surrounding countries",
    CODE == "ANCEST13" ~ "South Asia including India, Pakistan, Sri Lanka and surrounding countries",
    CODE == "ANCEST14" ~ "Polynesia, Micronesia or Melanesia including Tonga, Fiji, Papua New Guinea and surrounding countries",
    CODE == "ANCEST15" ~ "Africa",
    CODE == "ANCEST16" ~ "North America - not of First Nations, Native American, Inuit or Métis descent",
    CODE == "ANCEST17" ~ "North America - of First Nations, Native American, Inuit or Métis descent",
    CODE == "ANCEST18" ~ "Caribbean, Central or South America",
    CODE =="ANCEST20" ~ "Other"
  )) %>%
  select(STUDYID, AGE, SEX, ANCESTRY)



