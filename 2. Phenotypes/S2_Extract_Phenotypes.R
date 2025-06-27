

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

##########
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
  mutate(
    # First calculate if all symptoms are missing
    CORE_MISSING = is.na(LOWINT2W) & is.na(DEP2WK),
    ALL_MISSING = rowSums(is.na(select(., -STUDYID))) == 9,
    # Calculate symptom sum
    CUM_SYM = case_when(
      CORE_MISSING ~ NA_real_,
      TRUE ~ rowSums(select(., -STUDYID), na.rm=TRUE)
    ),
    # Then calculate MDD
    MDD = case_when(
      CORE_MISSING ~ NA_real_,
      (LOWINT2W == 1 | DEP2WK == 1) & CUM_SYM >= 5 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  select(-ALL_MISSING, -CORE_MISSING)  # Remove helper variable if desired

#-- Percentage of the cohort with Lifetime MDD (excluding NAs)
perc_MDD <- sum(sym$MDD, na.rm = TRUE) / sum(!is.na(sym$MDD)) * 100 # 88.
# including NAs
perc_MDD <- sum(sym$MDD, na.rm = TRUE) / length(sym$MDD) * 100

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

