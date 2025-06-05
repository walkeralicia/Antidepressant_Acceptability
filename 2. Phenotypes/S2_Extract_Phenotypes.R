

#========== Load required R libraries =========
library(dplyr)
library(tibble)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(RColorBrewer)

#== phenotypes
wkdir="/QRISdata/Q7280/pharmacogenomics/"
outdir="/scratch/user/uqawal15"
pheno_raw <- read.csv("/QRISdata/Q5338/Phenotypes/Combined_AGDS_2020_Freeze_20220620.csv") 
covariates <- read.csv(file.path(wkdir, "phenotypes/BMI_age_sex.csv"))

invalid_entries <- pheno_raw[!grepl("^MDD", pheno_raw$STUDYID), ]
pheno <- pheno_raw %>% filter(!STUDYID %in% invalid_entries$STUDYID)

#== baseline and follow-up bip diagnosis
followup_bip <- fread("/QRISdata/Q5338/Phenotypes/followUpData/dataforAlicia", fill = TRUE, header = TRUE) %>%
  select(Q3_2, StudyCode.0) %>%
  rename(BPD = Q3_2,
         STUDYID = StudyCode.0)

#-- Join BIP data frames
pheno <- left_join(pheno, followup_bip, by = "STUDYID") %>%
  mutate(DXBPD2 = ifelse(DXBPD == 1 | BPD == "Bipolar", 1, 0))

#=== Depression ===
dep <- pheno %>% select(STUDYID, AGE2WKF, TIMEWKS, TIMES2WK) %>%
  mutate(across(-STUDYID, as.numeric))

#=== Proportion of patients with DSM-5 symptoms during worst MDD episode ===
sym <- pheno %>% select(STUDYID, LOWINT2W, DEP2WK, APCHANGE, WTCHANGE, DIFFFALL, 
                        SLEPMORE, FIDGETY, SLOWMOVE, FATIGUED, GUILTY, NOFOCUS,
                        DEATHTHK)  %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(APWTCHANGE = ifelse(APCHANGE == 2 | 
                               WTCHANGE%in%c(1:3), 2, 1),
         SLEEP = ifelse(DIFFFALL == 2 | SLEPMORE == 2, 2, 1),
         MOVEMENT = ifelse(FIDGETY == 2 | SLOWMOVE == 2, 2, 1)) %>%
  select(-APCHANGE, -WTCHANGE, -DIFFFALL, -SLEPMORE, -FIDGETY, -SLOWMOVE) %>%
  mutate(across(-STUDYID, ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1,
    TRUE ~ .
  ))) %>% 
  mutate(CUM_SYM = rowSums(select(., -STUDYID), na.rm=FALSE)) %>%
  mutate(MDD = case_when(
    (is.na(LOWINT2W) & is.na(DEP2WK)) | is.na(CUM_SYM) ~ NA_real_,
    (LOWINT2W == 1 | DEP2WK == 1) & CUM_SYM >= 5 ~ 1,
    TRUE ~ 0)
  )

#-- Percentage of the cohort with Lifetime MDD
perc_MDD <- sum(sym$MDD, na.rm = TRUE) / sum(!is.na(sym$MDD)) * 100

#=== Weight change ===
wtchange <- pheno %>% select(STUDYID, WTCHANGE, WTKILO) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(WTLOSS = ifelse(WTCHANGE == 2, 1, 0),
         WTGAIN = ifelse(WTCHANGE == 1, 1, 0)) %>%
  mutate(WTLOSS_KG = ifelse(WTLOSS == 1, WTKILO, NA),
         WTGAIN_KG = ifelse(WTGAIN == 1, WTKILO, NA)) %>%
  select(-WTCHANGE, -WTKILO)

#=== Number of co-existing co-morbidity ===
com <- pheno %>% select(STUDYID, DXBPD2, DXPDMD, DXSCZ, DXANOR, DXBUL, DXADHD,
                        DXASD, DXSUD, DXPERSD, DXAGORA, DXANX, DXSAD, DXPHOB,
                        DXPTSD, DXHOARD, DXOCD, DXPANIC, DXTOUR) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(-STUDYID, ~ case_when(
    . == 1 ~ 1,
    . != 1 ~ 0,
    is.na(.) ~ 0
  ))) %>%
  mutate(TOTAL_COM = rowSums(select(., -STUDYID), na.rm=TRUE))


#=== Comorbidity (T2D, stomach ulcers) prevalence ===
pre <- pheno %>% select(STUDYID, TYPE11B, TYPE33C)  %>%
  mutate(across(-STUDYID, as.numeric)) %>% 
  mutate(across(TYPE11B, ~ case_when(
    . == 1 ~ 1,
    is.na(.) ~ 0
  ))) %>% 
  mutate(across(TYPE33C, ~ case_when(
    . == 1 ~ 1,
    is.na(.) ~ 0
  )))

#=== Suicidal Ideation and self-harm ===
sui <- pheno %>% select(STUDYID, SUICIDEA, SELFHARM) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(SUICIDEA, ~ case_when(
    . == 2 ~ 0,
    . == 1 ~ 1
  ))) %>%
  mutate(across(SELFHARM, ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1
  ))) 

#=== Familial mental health history

fam <- pheno %>% select(STUDYID, FAMMHD, FAMMDD, FAMBPD) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(FAMMHD, ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1
  )))

#=== Overall education and physical health level
edu <- pheno %>% select(STUDYID, EDU, PHYSHLTH) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(EDU, ~ case_when(
    . == 8 ~ NA,
    TRUE ~ .
  )))

#=== Drug use
drugs <- pheno %>% select(STUDYID, REGSMK, REGTOB, COC10, AMP10, KET10, ECS10, DRK3FRQ, DRK5FRQ) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(-c(STUDYID, DRK3FRQ, DRK5FRQ), ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1
  )))

#=== Medical conditions
# back pain, chronic fatigue syndrome, epilepsy, chronic pain 
med <- pheno %>% select(STUDYID, DXPHYS3, DXPHYS6, DXPHYS12, DXPHYS34) %>%
  mutate(across(-STUDYID, as.numeric)) %>% 
  mutate(across(-STUDYID, ~ case_when(
    . == 1 ~ 1,
    is.na(.) ~ 0
  )))

#=== Associated with headaches
head <- pheno %>% select(STUDYID, MIGEVR, MIGSTOM, MIGNVD, MIGVIS) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(-STUDYID, ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1
  )))

#=== Female health
fem <- pheno %>% select(STUDYID, DXFIBRO, DXPCOS, DXENDO) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(-STUDYID, ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1
  )))

#=== self-reported that the antidepressant works well for them

# Define the columns to rename
columns_to_rename <- c("WELLSERT" = "Sertraline",
                       "WELLESCI" = "Escitalopram",
                       "WELLCITA" = "Citalopram",
                       "WELLFLUO" = "Fluoxetine",
                       "WELLPARO" = "Paroxetine",
                       "WELLDESV" = "Desvenlafaxine",
                       "WELLVENL" = "Venlafaxine",
                       "WELLDULO" = "Duloxetine",
                       "WELLAMIT" = "Amitriptyline",
                       "WELLMIRT" = "Mirtazapine")

well <- pheno %>% 
  select(STUDYID, all_of(names(columns_to_rename))) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  rename(!!!setNames(names(columns_to_rename),columns_to_rename)) %>%
  mutate(across(-STUDYID, ~case_when(
    . == 1 ~ 0,
    . == 2 ~ 1,
    . == 3 ~ 1,
    TRUE ~ NA
  )))

#=== self reported side effects =========

sideeffs <- pheno %>%
  select(STUDYID, 
         c(starts_with("DRYM"), starts_with("SWET"), starts_with("NAUS"), 
           starts_with("VOM"), starts_with("DIAR"), starts_with("CONS"),
           starts_with("HEAD"), starts_with("DIZZ"), starts_with("SHAK"), 
           starts_with("MUSC"), starts_with("DROW"), starts_with("WAKE"), 
           starts_with("HANX"), starts_with("AGIT"), starts_with("FATI"), 
           starts_with("WTGN"), starts_with("WTLS"), starts_with("RASH"), 
           starts_with("RUN"), starts_with("RSEX"), starts_with("BLUR"), 
           starts_with("SUIT"), starts_with("NOSE"), starts_with("STOP"))) %>%
  mutate(across(-STUDYID, as.numeric)) %>%
  mutate(across(-c(STUDYID, starts_with("STOP")), ~ case_when(
    . == 1 ~ 1,
    TRUE ~ 0
  ))) %>%
  mutate(across(starts_with("STOP"), ~ case_when(
    . == 1 ~ 0,
    . == 2 ~ 1,
    TRUE ~ NA_real_
  )))


#=== Covariates ===
cov <- covariates %>% select(STUDYID, AGE, SEX, BMI)

#=== Combine all AGDS phenotypes ===
datasets_to_join <- list(cov, dep, sym, com, wtchange, pre, sui, fam , edu, drugs, med, head, fem, well, sideeffs)
tab <- Reduce(function(x, y) full_join(x, y, by = "STUDYID"), datasets_to_join)
write.csv(tab, file.path(wkdir, "phenotypes/survey_phenotypes.csv"), row.names = FALSE)


######## Proportion of self-responders plot ################

well_modified <- well %>%
  mutate(across(-STUDYID, ~case_when(
    . == 0 ~ "Not at all well",
    . == 1 ~ "Moderately to very well",
    TRUE ~ NA_character_
  )))

# Reshape the data to long format
well_long <- well_modified %>%
  pivot_longer(-STUDYID, names_to = "Drug", values_to = "Response")

response_proportions <- well_long %>%
  filter(!is.na(Response)) %>%  # Exclude rows where Response is NA
  group_by(Drug, Response) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Drug) %>%
  mutate(Proportion = Count / sum(Count))

write.csv(response_proportions, file.path(wkdir, "pharma_summaries/response_proportions_usage_basis.csv"), row.names = FALSE)

#-- Group by antidepressant for total N
anti_N <- well_long %>%
  filter(!is.na(Response)) %>%
  group_by(Drug) %>%
  summarise(TotalParticipants = n_distinct(STUDYID), .groups = "drop")

#-- Merge group sizes with self-responders proportions and create drug label with N
response_N <- left_join(response_proportions, anti_N, by = "Drug") %>%
  mutate(
    DrugLabel = paste0(Drug, "\n(N = ", TotalParticipants, ")")
  )

ordered_drugs <- c("Sertraline\n(N = 9152)", "Escitalopram\n(N = 6996)", "Citalopram\n(N = 3933)", 
                   "Paroxetine\n(N = 2404)", "Fluoxetine\n(N = 5817)", "Desvenlafaxine\n(N = 4013)",
                   "Venlafaxine\n(N = 6290)", "Duloxetine\n(N = 3145)", 
                   "Amitriptyline\n(N = 2498)", "Mirtazapine\n(N = 3058)") 

response_N$DrugLabel <- factor(response_N$DrugLabel, 
                                    levels = ordered_drugs)

# Create a bar plot
base_theme <- theme_minimal(base_size = 20) +
  theme(text = element_text(family = "Helvetica"))

p1 <- ggplot(response_N, aes(x = DrugLabel, y = Proportion, fill = Response)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  #geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
           # position = position_dodge(), size = 5, angle = 90) +
  base_theme +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Antidepressants", y = "Percentage of Participants (%)", 
       title = "Self-responders", fill = "Self-response") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set1")

png(file.path(outdir, "SelfResponse_by_antidepressant.png"), bg="white", type=c("cairo"),
    width=250, height=180, units='mm', res = 500)
p1
dev.off()

###### Proportion of discontinuation due to any side-effect #################

stopad <- pheno %>%
  select(STUDYID, starts_with("STOP"))

stopad_modified <- stopad %>%
  mutate(across(-STUDYID, ~case_when(
    . == 1 ~ "Yes",
    TRUE ~ "No"
  )))

# Define the columns to rename
columns_to_rename <- c("STOPSERT" = "Sertraline",
                       "STOPESCI" = "Escitalopram",
                       "STOPCITA" = "Citalopram",
                       "STOPFLUO" = "Fluoxetine",
                       "STOPPARO" = "Paroxetine",
                       "STOPDESV" = "Desvenlafaxine",
                       "STOPVENL" = "Venlafaxine",
                       "STOPDULO" = "Duloxetine",
                       "STOPAMIT" = "Amitriptyline",
                       "STOPMIRT" = "Mirtazapine")


# Reshape the data to long format
stopad_long <- stopad_modified %>%
  select(STUDYID, all_of(names(columns_to_rename))) %>%
  rename(!!!setNames(names(columns_to_rename),columns_to_rename)) %>%
  pivot_longer(-STUDYID, names_to = "Drug", values_to = "Stopped")

stopped_proportions <- stopad_long %>%
  group_by(Drug, Stopped) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Drug) %>%
  mutate(Proportion = Count / sum(Count))

write.csv(stopped_proportions, file.path(wkdir, "pharma_summaries/stopped_proportions_usage_basis.csv"), row.names = FALSE)

ordered_drugs <- c("Sertraline", "Escitalopram", "Citalopram", 
                   "Paroxetine", "Fluoxetine","Desvenlafaxine",
                   "Venlafaxine", "Duloxetine", 
                   "Amitriptyline", "Mirtazapine") 

stopped_proportions$Drug <- factor(stopped_proportions$Drug, 
                                    levels = ordered_drugs)

colors <- brewer.pal(n = 9, name = "Set1")
third_color <- colors[3]
fourth_color <- colors[4]

# Create a bar plot
base_theme <- theme_minimal(base_size = 20) +
  theme(text = element_text(family = "Helvetica"))

p2 <- ggplot(stopped_proportions, aes(x = Drug, y = Proportion, fill = Stopped)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  #geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
            #position = position_stack(vjust = 0.5), size = 5, angle = 90) +
  base_theme +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Antidepressants", y = "Percentage of Participants (N = 23211)", 
       title = "Reported discontinuation due to side-effects", fill = "Discontinuation") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
        legend.title = element_blank(),
        legend.position = "top") +
  scale_fill_manual(values = c(third_color, fourth_color))

png(file.path(outdir, "Discontinuation_by_antidepressant.png"), bg="white", type=c("cairo"),
    width=250, height=180, units='mm', res = 500)
p2
dev.off()









