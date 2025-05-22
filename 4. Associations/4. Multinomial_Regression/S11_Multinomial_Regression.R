

#-- Load R libraries
library(dplyr)
library(tibble)

#-- set wkdir
wkdir="/QRISdata/Q7280/pharmacogenomics/"

#-- Load data
duration = 360

dat <- read.csv("/QRISdata/Q7280/pharmacogenomics/phenotypes/treatment_phenotypes/inner_data_360days.csv") %>%
  filter(DrugClass != "BIP-L" & DrugClass != "BIP+L" & DrugClass != "Various")

#-- Set SSRI as the reference group
dat$DrugClass <- factor(dat$DrugClass)
dat$DrugClass <- relevel(factor(dat$DrugClass), ref = "SSRI")

#-- Scale quantitative variables
dat <- dat %>%
  mutate(AGE_scaled = scale(AGE),
         BMI_scaled = scale(BMI),
         AGE2WKF_scaled = scale(AGE2WKF),
         TIMES2WK_scaled = scale(TIMES2WK),
         EDU_scaled = scale(EDU),
         PHYSHLTH_scaled = scale(PHYSHLTH)
         )


rename_mapping <- c("AGE_scaled" = "Age",
                    "BMI_scaled" = "BMI",
                    "TIMES2WK_scaled" = "Number of MDD Episodes",
                    "PHYSHLTH_scaled" = "Physical Health",
                    "SUICIDEA" = "Suicidal Ideation",
                    "DEATHTHK" = "Death Thoughts",
                    "TYPE11B" = "T2D",
                    "DXANX" = "Anxiety Disorder",
                    "DXPERSD" = "Personality Disorder",
                    "DXOCD" = "OCD",
                    "DXADHD" = "ADHD",
                    "DXPHYS34" = "Chronic Pain",
                    "DXPCOS" = "PCOS",
                    "DXPDMD" = "PDMD",
                    "MDD_LOO_pgs_std" = "MDD PGS",
                    "PUD_01_pgs_std" = "PUD PGS",
                    "T2D_03_pgs_std" = "T2D PGS",
                    "EDU_scaled" = "Education Level",
                    "DXPHYS6" = "CFS",
                    "ADHD_01_pgs_std" = "ADHD PGS",
                    "BMI_LOO_pgs_std" = "BMI PGS",
                    "Migraine_01_pgs_std" = "Migraine PGS", 
                    "BIP_LOO_pgs_std" = "BIP PGS"
)


#=================== Only Self-Report Characteristics ===============================

#-- Complete dataset
#-- na_count <- colSums(is.na(dat))
#-- do not include REGSMK due to high level of missingness (N=2617)
#-- do not use sex because of its dominance in SNRIs
dat_full <- dat %>%
  select("ParticipantID", "IID", "DrugClass", "AGE_scaled", "BMI_scaled", "AGE2WKF_scaled", "TIMES2WK_scaled", "EDU_scaled", "PHYSHLTH_scaled","SUICIDEA",
         "TYPE11B", "TYPE33C", "DXPDMD", "DXANX", "DXPERSD", "DXSUD", "DXADHD", "DXOCD", "DXSAD","DXPHYS3", "DXPHYS6", "DXPHYS12", "DXPHYS34", "MIGEVR", "DXFIBRO", "DXPCOS", "DXENDO",
         "FATIGUED.x", "GUILTY", "NOFOCUS", "DEATHTHK", "APWTCHANGE", "SLEEP", "MOVEMENT") %>%
  na.omit()

table(dat_full$DrugClass) # 360 days
#SSRI SNRI  TCA TeCA
#1752 1061   62   61
counts <- dat_full %>%
  group_by(DrugClass) %>%
  summarise(
    Group_N = n()
  )


# Multinomial logistic regression
# Use first level as reference category
multinomial_model <- nnet::multinom(
  DrugClass ~ AGE_scaled + BMI_scaled + AGE2WKF_scaled + TIMES2WK_scaled + EDU_scaled + PHYSHLTH_scaled + SUICIDEA +
    TYPE11B + TYPE33C + DXPDMD + DXANX + DXPERSD + DXSUD + DXADHD + DXOCD + DXSAD + 
    DXPHYS3 + DXPHYS6 + DXPHYS12 + DXPHYS34 + 
    MIGEVR + DXFIBRO + DXPCOS + DXENDO +
    FATIGUED.x + GUILTY + NOFOCUS + 
    DEATHTHK + APWTCHANGE + SLEEP + MOVEMENT,
  data = dat_full
)

coef_df <- broom::tidy(multinomial_model, exponentiate = TRUE, conf.int = TRUE)

#-- Backward Step wise Variable Selection
reduced_model <- step(multinomial_model, direction = "backward", trace = TRUE)

AIC(multinomial_model, reduced_model)
anova(multinomial_model, reduced_model, test = "Chisq")

simple_coef_df <- broom::tidy(reduced_model, exponentiate = TRUE, conf.int = TRUE)

df <- left_join(simple_coef_df, counts, by = c("y.level" = "DrugClass"))

#-- Rename Dependent Variables
df_renamed <- df %>%
  mutate(term = recode(term, !!!rename_mapping)) %>%
  mutate_at(vars(p.value, std.error), ~ signif(., 2)) %>%
  mutate_at(vars(estimate, statistic, conf.low, conf.high), ~ round(., 2)) %>%
  add_column(Model = rep("Self Report Predictors", nrow(df)), .before = 1)

#=============== Include PGS into Multinomial regression model =================================

dat_copy <- dat

#-- Genetic data
pgslist <- system('find /QRISdata/Q7280/pharmacogenomics/pgs -name "*agds_sbrc_gctb_plink2.sscore" -type f -exec ls {} \\;', intern = TRUE)
pgs_codes <- c("PUD_01", "T2D_03", "CNT_03", "ADHD_01", "BIP_LOO", "BMI_LOO", "MDD_LOO", "MDD_07_hsa04726", "SCZ_02", "Migraine_01", "SBP_01", "ANX_LOO", "ANO_LOO", "UKB_35BM_2021_LDL_direct_adjstatins", "CRP_01", "OCD_2024", "Neuroticism_01", "LRA_01")
num <- length(pgs_codes)
selected_pgs <- pgslist[grep(paste(pgs_codes, collapse = "|"), pgslist)]

#-- List of Europeans
eur <- read.table("/QRISdata/Q5338/Ancestry_analysis/PCA/AGDS_EUR.id") # 14,603 individuals
ids_to_keep <- as.character(eur$V2)

#-- For each PGS trait fit a linear model with treatment group as the independent variable
for (j in 1:length(selected_pgs)){
  
  if (j != 16){
    #-- Read in PGS file for the selected trait
    pgs_file_path <- selected_pgs[j]
    pgs <- read.table(selected_pgs[j], header=FALSE)
    
    #-- Extract PGS trait name
    pgs_filename <- basename(pgs_file_path)
    trait_parts <- unlist(strsplit(unlist(strsplit(pgs_filename, "\\."))[1], "_"))
    trait <- paste0(trait_parts[1], "_", trait_parts[2])
    cat(trait, j, '\n')
    
    #-- Filter PGS for Europeans (loss of 806 individuals)
    ids_to_keep <- as.character(eur$V2)
    pgs_e <- pgs[pgs$V2 %in% ids_to_keep, ]
    
    #-- Standardize PGS
    pgs_e$std_pgs <- scale(pgs_e$V6) 
    
    #-- Join PGS file with ad treatment group data
    pgs_e <- pgs_e %>% select(V2, std_pgs)
    colnames(pgs_e) <- c("V2", paste0(trait, "_pgs_std"))
    dat_copy <- inner_join(dat_copy, pgs_e, by = c("IID" = "V2"))
  } else {
    
    #-- Read in PGS file for the selected trait
    pgs_file_path <- selected_pgs[j]
    pgs <- read.table(selected_pgs[j], header=TRUE)
    
    #-- Extract PGS trait name
    pgs_filename <- basename(pgs_file_path)
    trait_parts <- unlist(strsplit(unlist(strsplit(pgs_filename, "\\."))[1], "_"))
    trait <- paste0(trait_parts[1], "_", trait_parts[2], "_", trait_parts[3])
    cat(trait, j, '\n')
    
    #-- Filter PGS for Europeans (loss of 806 individuals)
    ids_to_keep <- as.character(eur$V2)
    pgs_e <- pgs[pgs$IID %in% ids_to_keep, ]
    
    #-- Standardize PGS
    pgs_e$std_pgs <- scale(pgs_e$SCORESUM) 
    
    #-- Join PGS file with ad treatment group data
    pgs_e <- pgs_e %>% select(IID, std_pgs)
    colnames(pgs_e) <- c("IID", paste0(trait, "_pgs_std"))
    dat_copy <- inner_join(dat_copy, pgs_e, by = c("IID" = "IID"))
  }
}


#-- Complete dataset
pgs_atc_full <- dat_copy %>%
  select("ParticipantID", "IID", "DrugClass", "AGE_scaled", "BMI_scaled", "AGE2WKF_scaled", "TIMES2WK_scaled", "EDU_scaled", "PHYSHLTH_scaled","SUICIDEA",
         "TYPE11B", "TYPE33C", "DXPDMD", "DXANX", "DXPERSD", "DXSUD", "DXADHD", "DXOCD", "DXSAD","DXPHYS3", "DXPHYS6", "DXPHYS12", "DXPHYS34", "MIGEVR", "DXFIBRO", "DXPCOS", "DXENDO",
         "FATIGUED.x", "GUILTY", "NOFOCUS", "DEATHTHK", "APWTCHANGE", "SLEEP", "MOVEMENT", tail(colnames(dat_copy), 18)) %>%
  na.omit()

table(pgs_atc_full$DrugClass)
#SSRI SNRI  TCA TeCA
#1030  675   31   41

counts <- pgs_atc_full %>%
  group_by(DrugClass) %>%
  summarise(
    Group_N = n()
  )


# Multinomial logistic regression
# Use first level as reference category
multinomial_model <- nnet::multinom(
  DrugClass ~ AGE_scaled + BMI_scaled + AGE2WKF_scaled + TIMES2WK_scaled + EDU_scaled + PHYSHLTH_scaled + SUICIDEA +
    TYPE11B + TYPE33C + DXPDMD + DXANX + DXPERSD + DXSUD + DXADHD + DXOCD + DXSAD + 
    DXPHYS3 + DXPHYS6 + DXPHYS12 + DXPHYS34 + 
    MIGEVR + DXFIBRO + DXPCOS + DXENDO +
    FATIGUED.x + GUILTY + NOFOCUS + 
    DEATHTHK + APWTCHANGE + SLEEP + MOVEMENT +
    PUD_01_pgs_std + T2D_03_pgs_std + CNT_03_pgs_std + ADHD_01_pgs_std + BIP_LOO_pgs_std + BMI_LOO_pgs_std + Neuroticism_01_pgs_std +
    MDD_LOO_pgs_std + SCZ_02_pgs_std + Migraine_01_pgs_std + SBP_01_pgs_std + ANX_LOO_pgs_std + ANO_LOO_pgs_std + LRA_01_pgs_std +
    OCD_2024_pgs_std + CRP_01_pgs_std + UKB_35BM_pgs_std,
  data = pgs_atc_full
)

coef_df <- broom::tidy(multinomial_model, exponentiate = TRUE, conf.int = TRUE)

#-- Backward Step wise Variable Selection
reduced_model <- step(multinomial_model, direction = "backward", trace = TRUE)

AIC(multinomial_model, reduced_model)
anova(multinomial_model, reduced_model, test = "Chisq")

simple_coef_df <- broom::tidy(reduced_model, exponentiate = TRUE, conf.int = TRUE)

df <- left_join(simple_coef_df, counts, by = c("y.level" = "DrugClass"))

#-- Rename Dependent Variables
df_renamed_2 <- df %>%
  mutate(term = recode(term, !!!rename_mapping)) %>%
  mutate_at(vars(p.value, std.error), ~ signif(., 2)) %>%
  mutate_at(vars(estimate, statistic, conf.low, conf.high), ~ round(., 2)) %>%
  add_column(Model = rep("All Predictors", nrow(simple_coef_df)), .before = 1)


#===================== Only PGS as Predictors ==================================

#-- Complete dataset
pgs_full <- dat_copy %>%
  select("ParticipantID", "IID", "DrugClass", tail(colnames(dat_copy), 18)) %>%
  na.omit()

table(pgs_full$DrugClass)
#SSRI SNRI  TCA TeCA
#1593 1119   56   82

counts <- pgs_full %>%
  group_by(DrugClass) %>%
  summarise(
    Group_N = n()
  )

# Multinomial logistic regression
# Use first level as reference category
multinomial_model <- nnet::multinom(
  DrugClass ~ PUD_01_pgs_std + T2D_03_pgs_std + CNT_03_pgs_std + ADHD_01_pgs_std + BIP_LOO_pgs_std + BMI_LOO_pgs_std + Neuroticism_01_pgs_std +
    MDD_LOO_pgs_std + SCZ_02_pgs_std + Migraine_01_pgs_std + SBP_01_pgs_std + ANX_LOO_pgs_std + ANO_LOO_pgs_std + LRA_01_pgs_std +
    OCD_2024_pgs_std + CRP_01_pgs_std + UKB_35BM_pgs_std,
  data = pgs_full
)

coef_df <- broom::tidy(multinomial_model, exponentiate = TRUE, conf.int = TRUE)

#-- Backward Step wise Variable Selection
reduced_model <- step(multinomial_model, direction = "backward", trace = TRUE)

AIC(multinomial_model, reduced_model)
anova(multinomial_model, reduced_model, test = "Chisq")

simple_coef_df <- broom::tidy(reduced_model, exponentiate = TRUE, conf.int = TRUE)

df <- left_join(simple_coef_df, counts, by = c("y.level" = "DrugClass"))

#-- Rename Dependent Variables
df_renamed_3 <- df %>%
  mutate(term = recode(term, !!!rename_mapping)) %>%
  mutate_at(vars(p.value, std.error), ~ signif(., 2)) %>%
  mutate_at(vars(estimate, statistic, conf.low, conf.high), ~ round(., 2)) %>%
  add_column(Model = rep("PGS Predictors", nrow(simple_coef_df)), .before = 1)


df_all <- rbind(df_renamed, df_renamed_2, df_renamed_3)
df_all <- df_all %>%
  group_by(Model, y.level) %>%
  mutate(
    Model = if_else(Model != lag(Model, default = ""), Model, NA_character_),
    y.level = if_else(y.level != lag(y.level, default = ""), y.level, NA_character_)
  )


wb <- loadWorkbook("/scratch/user/uqawal15/All_Results.xlsx")
addWorksheet(wb, "Table17")
writeData(wb, "Table17", df_all)
saveWorkbook(wb, file.path("/scratch/user/uqawal15", "All_Results.xlsx"), overwrite = TRUE)




