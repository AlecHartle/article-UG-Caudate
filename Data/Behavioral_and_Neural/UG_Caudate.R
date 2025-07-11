library(car)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyverse)
library(sjPlot)
library(knitr) 
library(kableExtra)
library(magick)
library(tidyr) 
library(dplyr) 
library(ggplotify)  
library(ggplot2)


# Set the correct working directory (replace with the correct path)
setwd('')
file_name = "UG_Caudate_Table_R.csv" ###set .csv file

dat <- read.csv(file_name, na = "NaN", header = TRUE)

dat$subjectNum = as.factor(dat$subjectNum)
dat$PEtype <- as.factor(dat$PEtype)
dat$PatientType <- as.factor(dat$PatientType)
dat$emoRating <- as.numeric(dat$emoRating)
dat$offerDecision = as.factor(dat$offerDecision)

setwd('')
file_name = 'UG_Caudate_Clinical_Measures_R.csv' ###set .csv file for clinical measures

dat_coord <- read.csv(file_name, na = "NaN", header = TRUE)

dat_coord$subjectNum = as.factor(dat_coord$subjectNum)
dat_coord$UPDP_score = as.numeric(dat_coord$UPDP_score)
dat_coord$non_motor_score = as.numeric(dat_coord$non_motor_score)
dat_coord$PD_duration = as.numeric(dat_coord$PD_duration)


dat <- dat %>%
  left_join(dat_coord[, c("subjectNum", "anterior", "UPDP_score", "non_motor_score", "PD_duration")], by = "subjectNum")

dat_ET = subset(dat, PatientType ==1)
dat_PD = subset(dat, PatientType ==2)

lower_cutoff <- quantile(dat_coord$UPDP_score, 0.25, na.rm = TRUE)
upper_cutoff <- quantile(dat_coord$UPDP_score, 0.75, na.rm = TRUE)

# Create new column for tail category
dat_coord$UPDP_tail <- NA  # initialize with NA
dat_coord$UPDP_tail[dat_coord$UPDP_score <= lower_cutoff] <- "low"
dat_coord$UPDP_tail[dat_coord$UPDP_score > lower_cutoff & dat_coord$UPDP_score < upper_cutoff] <- "mid"
dat_coord$UPDP_tail[dat_coord$UPDP_score >= upper_cutoff] <- "high"


# Create new data frames for high and low tails
PD_UPDP_low <- dat_coord[dat_coord$UPDP_tail == "low", ]
PD_UPDP_mid <- dat_coord[dat_coord$UPDP_tail == "mid", ]
PD_UPDP_high <- dat_coord[dat_coord$UPDP_tail == "high", ]
# Get subjectNum from low tail (excluding NAs)
PD_UPDP_low_subjects <- na.omit(PD_UPDP_low$subjectNum)
PD_UPDP_mid_subjects <- na.omit(PD_UPDP_mid$subjectNum)
PD_UPDP_high_subjects <- na.omit(PD_UPDP_high$subjectNum)

PD_UPDP_low_subjects=as.numeric(PD_UPDP_low_subjects)
PD_UPDP_mid_subjects=as.numeric(PD_UPDP_mid_subjects)
PD_UPDP_high_subjects = as.numeric(PD_UPDP_high_subjects)

# Define the cutoffs
lower_cutoff <- quantile(dat_coord$non_motor_score, 0.25, na.rm = TRUE)
upper_cutoff <- quantile(dat_coord$non_motor_score, 0.75, na.rm = TRUE)

# Create new column for tail category
dat_coord$non_motor_tail <- NA  # initialize with NA
dat_coord$non_motor_tail[dat_coord$non_motor_score <= lower_cutoff] <- "low"
dat_coord$non_motor_tail[dat_coord$non_motor_score > lower_cutoff & dat_coord$non_motor_score < upper_cutoff] <- "mid"
dat_coord$non_motor_tail[dat_coord$non_motor_score >= upper_cutoff] <- "high"


# Create new data frames for high and low tails
PD_non_motor_low <- dat_coord[dat_coord$non_motor_tail == "low", ]
PD_non_motor_mid <- dat_coord[dat_coord$non_motor_tail == "mid", ]
PD_non_motor_high <- dat_coord[dat_coord$non_motor_tail == "high", ]
# Get subjectNum from low tail (excluding NAs)
PD_non_motor_low_subjects <- na.omit(PD_non_motor_low$subjectNum)
PD_non_motor_mid_subjects <- na.omit(PD_non_motor_mid$subjectNum)
PD_non_motor_high_subjects <- na.omit(PD_non_motor_high$subjectNum)


PD_non_motor_low_subjects=as.numeric(PD_non_motor_low_subjects)
PD_non_motor_mid_subjects=as.numeric(PD_non_motor_mid_subjects)
PD_non_motor_high_subjects = as.numeric(PD_non_motor_high_subjects)



####create the labels
dat_PD <- dat_PD %>%
  mutate(UPDP_tail = case_when(
    subjectNum %in% PD_UPDP_low_subjects    ~ 1L,
    subjectNum %in% PD_UPDP_mid_subjects    ~ 2L,
    subjectNum %in% PD_UPDP_high_subjects ~ 3L,
    TRUE                       ~ NA_integer_
  ))
dat_PD$UPDP_tail <- as.factor(dat_PD$UPDP_tail)

dat_PD <- dat_PD %>%
  mutate(non_motor_tail = case_when(
    subjectNum %in% PD_non_motor_low_subjects    ~ 1L,
    subjectNum %in% PD_non_motor_mid_subjects    ~ 2L,
    subjectNum %in% PD_non_motor_high_subjects ~ 3L,
    TRUE                       ~ NA_integer_
  ))
dat_PD$non_motor_tail <- as.factor(dat_PD$non_motor_tail)

# Define the cutoffs
lower_cutoff <- quantile(dat_coord$PD_duration, 0.25, na.rm = TRUE)
upper_cutoff <- quantile(dat_coord$PD_duration, 0.75, na.rm = TRUE)

# Create new column for tail category
dat_coord$PD_duration_tail <- NA  # initialize with NA
dat_coord$PD_duration_tail[dat_coord$PD_duration <= lower_cutoff] <- "low"
dat_coord$PD_duration_tail[dat_coord$PD_duration > lower_cutoff & dat_coord$PD_duration < upper_cutoff] <- "mid"
dat_coord$PD_duration_tail[dat_coord$PD_duration >= upper_cutoff] <- "high"

# View the new column
table(dat_coord$PD_duration_tail)

# Create new data frames for high and low tails
PD_duration_low <- dat_coord[dat_coord$PD_duration_tail == "low", ]
PD_duration_mid <- dat_coord[dat_coord$PD_duration_tail == "mid", ]
PD_duration_high <- dat_coord[dat_coord$PD_duration_tail == "high", ]
# Get subjectNum from low tail (excluding NAs)
PD_duration_low_subjects <- na.omit(PD_duration_low$subjectNum)
PD_duration_mid_subjects <- na.omit(PD_duration_mid$subjectNum)
PD_duration_high_subjects <- na.omit(PD_duration_high$subjectNum)


PD_duration_low_subjects=as.numeric(PD_duration_low_subjects)
PD_duration_mid_subjects=as.numeric(PD_duration_mid_subjects)
PD_duration_high_subjects = as.numeric(PD_duration_high_subjects)

dat_PD <- dat_PD %>%
  mutate(PD_duration_tail = case_when(
    subjectNum %in% PD_duration_low_subjects    ~ 1L,
    subjectNum %in% PD_duration_mid_subjects    ~ 2L,
    subjectNum %in% PD_duration_high_subjects ~ 3L,
    TRUE                       ~ NA_integer_
  ))
dat_PD$PD_duration_tail <- as.factor(dat_PD$PD_duration_tail)

dat_long_PD_score <- dat_PD %>%
  # Pivot longer for AUC columns
  pivot_longer(cols = c(AUCDA, AUCNE, AUC5HT),
               names_to = "Neurotransmitter",
               values_to = "AUC_Measurement") %>%
  # Add NTtype classification for the AUC neurotransmitters
  mutate(NTtype = case_when(
    Neurotransmitter == "AUCDA" ~ 1,
    Neurotransmitter == "AUCNE" ~ 2,
    Neurotransmitter == "AUC5HT" ~ 3
  )) %>%
  # Remove the temporary Neurotransmitter column
  select(-Neurotransmitter)
dat_long_PD_score$NTtype = as.factor(dat_long_PD_score$NTtype)


dat_long_AUC <- dat %>%
  # Pivot longer for AUC columns
  pivot_longer(cols = c(AUCDA, AUCNE, AUC5HT),
               names_to = "Neurotransmitter",
               values_to = "AUC_Measurement") %>%
  # Add NTtype classification for the AUC neurotransmitters
  mutate(NTtype = case_when(
    Neurotransmitter == "AUCDA" ~ 1,
    Neurotransmitter == "AUCNE" ~ 2,
    Neurotransmitter == "AUC5HT" ~ 3
  )) %>%
  # Remove the temporary Neurotransmitter column
  select(-Neurotransmitter)



PD_no_SSRI = c(6,8,10,11,13)+6
PD_SSRI = c(2,3,4,5,9,12)+6
PD_all_Rx = c(1,2,3,4,5,7,9,12)+6

ET_SSRI = c(1,2,4)
ET_no_SSRI = c(3,5,6)

####create the labels
dat_PD <- dat_PD %>%
  mutate(SSRI = case_when(
    subjectNum %in% PD_SSRI    ~ 1L,
    subjectNum %in% PD_no_SSRI ~ 2L,
    TRUE                       ~ NA_integer_
  ))
dat_PD$SSRI <- as.factor(dat_PD$SSRI)

dat_PD <- dat_PD %>%
  mutate(all_Rx = case_when(
    subjectNum %in% PD_all_Rx    ~ 1L,
    subjectNum %in% PD_no_SSRI ~ 2L,
    TRUE                       ~ NA_integer_
  ))
dat_PD$all_Rx <- as.factor(dat_PD$all_Rx)

dat_ET <- dat_ET %>%
  mutate(SSRI = case_when(
    subjectNum %in% ET_SSRI    ~ 1L,
    subjectNum %in% ET_no_SSRI ~ 2L,
    TRUE                       ~ NA_integer_
  ))
dat_ET$SSRI <- as.factor(dat_ET$SSRI)


dat_long_PD <- dat_PD %>%
  # Pivot longer for AUC columns
  pivot_longer(cols = c(AUCDA, AUCNE, AUC5HT),
               names_to = "Neurotransmitter",
               values_to = "AUC_Measurement") %>%
  # Add NTtype classification for the AUC neurotransmitters
  mutate(NTtype = case_when(
    Neurotransmitter == "AUCDA" ~ 1,
    Neurotransmitter == "AUCNE" ~ 2,
    Neurotransmitter == "AUC5HT" ~ 3
  )) %>%
  # Remove the temporary Neurotransmitter column
  select(-Neurotransmitter)


dat_long_ET <- dat_ET %>%
  # Pivot longer for AUC columns
  pivot_longer(cols = c(AUCDA, AUCNE, AUC5HT),
               names_to = "Neurotransmitter",
               values_to = "AUC_Measurement") %>%
  # Add NTtype classification for the AUC neurotransmitters
  mutate(NTtype = case_when(
    Neurotransmitter == "AUCDA" ~ 1,
    Neurotransmitter == "AUCNE" ~ 2,
    Neurotransmitter == "AUC5HT" ~ 3
  )) %>%
  # Remove the temporary Neurotransmitter column
  select(-Neurotransmitter)

#
rejection_by_PEtype_PD <- dat_PD %>%
  group_by(PEtype, subjectNum) %>%
  summarize(
    total_decisions = n(),
    rejections = sum(offerDecision == 0, na.rm = TRUE),
    rejection_rate = (rejections / total_decisions) * 100  # Rejection rate as a percentage
  )

# Add PatientType column for PD
rejection_by_PEtype_PD <- rejection_by_PEtype_PD %>%
  mutate(PatientType = "PD")

# Calculate rejection rate per subjectNum within each PEtype for ET
rejection_by_PEtype_ET <- dat_ET %>%
  group_by(PEtype, subjectNum) %>%
  summarize(
    total_decisions = n(),
    rejections = sum(offerDecision == 0, na.rm = TRUE),
    rejection_rate = (rejections / total_decisions) * 100  # Rejection rate as a percentage
  )

# Add PatientType column for ET
rejection_by_PEtype_ET <- rejection_by_PEtype_ET %>%
  mutate(PatientType = "ET")

combined_rejection_by_PEtype <- bind_rows(rejection_by_PEtype_PD, rejection_by_PEtype_ET)



PD_0_rejects_pos <- combined_rejection_by_PEtype[
  combined_rejection_by_PEtype$PatientType   == "PD" &
    combined_rejection_by_PEtype$PEtype        == 1   &
    combined_rejection_by_PEtype$rejection_rate == 0, 
]


# subset the rows you want
PD_non_0_rejects_pos <- combined_rejection_by_PEtype[
  combined_rejection_by_PEtype$PatientType   == "PD" &
    combined_rejection_by_PEtype$PEtype        == 1   &
    combined_rejection_by_PEtype$rejection_rate > 0, 
]


# subset the rows you want
PD_100_rejects_neg <- combined_rejection_by_PEtype[
  combined_rejection_by_PEtype$PatientType   == "PD" &
    combined_rejection_by_PEtype$PEtype        == 2   &
    combined_rejection_by_PEtype$rejection_rate == 100, 
]

# subset the rows you want
PD_non_100_rejects_neg <- combined_rejection_by_PEtype[
  combined_rejection_by_PEtype$PatientType   == "PD" &
    combined_rejection_by_PEtype$PEtype        == 2   &
    combined_rejection_by_PEtype$rejection_rate < 100, 
]

PD_0_reject_pos_subjects = PD_0_rejects_pos$subjectNum
PD_non_0_reject_pos_subjects = PD_non_0_rejects_pos$subjectNum
PD_100_rejects_neg_subjects = PD_100_rejects_neg$subjectNum
PD_non_100_rejects_neg_subjects = PD_non_100_rejects_neg$subjectNum


##creating new column to run stats
dat_PD$reject_rate_0 <- ifelse(
  dat_PD$subjectNum %in% PD_0_reject_pos_subjects ,
  1,
  2
)



##creating new column to run stats. label subjects that 100% rejected negative PEs with a 1, all others with a 2
dat_PD$reject_rate_100 <- ifelse(
  dat_PD$subjectNum %in% PD_100_rejects_neg_subjects ,
  1,
  2
)


dat_ET$reject_rate_0 = 1
dat_ET$reject_rate_100 = 2

# 1. Find the extra columns in dat_PD
extra_cols <- setdiff(names(dat_PD), names(dat_ET))

# 2. Add those columns to dat_ET with NA values
for (col in extra_cols) {
  dat_ET[[col]] <- NA
}

# 3. (Optional) Reorder dat_ET columns to match dat_PD
dat_ET <- dat_ET[, names(dat_PD)]


dat_rejection <- rbind(dat_PD, dat_ET)

dat_rejection_long_AUC <- dat_rejection %>%
  # Pivot longer for AUC columns
  pivot_longer(cols = c(AUCDA, AUCNE, AUC5HT),
               names_to = "Neurotransmitter",
               values_to = "AUC_Measurement") %>%
  # Add NTtype classification for the AUC neurotransmitters
  mutate(NTtype = case_when(
    Neurotransmitter == "AUCDA" ~ 1,
    Neurotransmitter == "AUCNE" ~ 2,
    Neurotransmitter == "AUC5HT" ~ 3
  )) %>%
  # Remove the temporary Neurotransmitter column
  select(-Neurotransmitter)


######----------emotional rating by patienttype and PEtype

# Calculate rejection rate per subjectNum within each PEtype for PD
emorating_by_subject_PD <- dat_PD %>%
  group_by(PEtype, subjectNum) %>%
  summarize(
    total_decisions = n(),
    ratings = mean(emoRating, na.rm = TRUE),
  )

# Add PatientType column for PD
emorating_by_subject_PD <- emorating_by_subject_PD %>%
  mutate(PatientType = "PD")

# Calculate rejection rate per subjectNum within each PEtype for ET
emorating_by_subject_ET <- dat_ET %>%
  group_by(PEtype, subjectNum) %>%
  summarize(
    total_decisions = n(),
    ratings = mean(emoRating, na.rm = TRUE),
  )

# Add PatientType column for ET
emorating_by_subject_ET <- emorating_by_subject_ET %>%
  mutate(PatientType = "ET")

# Combine both datasets
combined_emorating_by_subject <- bind_rows(emorating_by_subject_PD, emorating_by_subject_ET)

combined_emorating_by_subject <- combined_emorating_by_subject %>%
  mutate(PEtype = ifelse(PEtype == 1, "PosPE", "NegPE"))

# Calculate mean and SEM for each PEtype and PatientType
summary_stats_emo <- dat %>%
  group_by(PEtype, PatientType) %>%
  summarize(
    mean_emo_rate = mean(emoRating, na.rm = TRUE),
    sem_emo_rate = sd(emoRating, na.rm = TRUE) / sqrt(n())  # SEM calculation
  )
summary_stats_emo <- summary_stats_emo %>%
  mutate(PEtype = ifelse(PEtype == 1, "PosPE", "NegPE"))


##########--------------------------------------------------------------------Behavior Model ID: B-M1 Reaction

model <- lmer(reactionTime ~ PatientType*PEtype+ (1 | subjectNum), data = dat)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                    Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# PatientType        54.394  54.394     1  17.40  1.2748 0.2742
# PEtype             69.100  69.100     1 555.12  1.6194 0.2037
# PatientType:PEtype 76.194  76.194     1 555.12  1.7857 0.1820

##########--------------------------------------------------------------------Behavior Model ID: B-M2 Rejection Rate 

model = lmer(rejection_rate ~ 1 + PatientType*PEtype + (1 | subjectNum),  data = combined_rejection_by_PEtype)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                     Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# PatientType          174.6   174.6     1    17  0.4507    0.5110    
# PEtype             15930.1 15930.1     1    17 41.1128 6.441e-06 ***
# PatientType:PEtype   190.5   190.5     1    17  0.4917    0.4926 

##########--------------------------------------------------------------------Behavior Model ID: B-M3 Emotional Rating 

model_PE <- lmer(emoRating ~ PatientType*PEtype + (1 | subjectNum), data = dat)
as.data.frame(anova(model_PE, ddf = "Kenward-Roger"))
# Sum Sq      Mean Sq NumDF     DenDF      F value       Pr(>F)
# PatientType          0.05897642   0.05897642     1  17.29108   0.02357177 8.797571e-01
# PEtype             361.03548804 361.03548804     1 323.69246 144.29911768 9.696604e-28
# PatientType:PEtype   1.15120786   1.15120786     1 323.69246   0.46011620 4.980543e-01

##########--------------------------------------------------------------------Neural Model ID: N-M1 
model_full <- lmer(AUC_Measurement ~ factor(NTtype) *factor(PEtype)*factor(PatientType) + (1 | subjectNum) , data = dat_long_AUC)
anova(model_full, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                                                    Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# factor(NTtype)                                     114.95   57.47     2 1669.02  0.2089 0.81150  
# factor(PEtype)                                      32.42   32.42     1 1213.90  0.1178 0.73146  
# factor(PatientType)                                299.57  299.57     1   19.38  1.0889 0.30957  
# factor(NTtype):factor(PEtype)                     2522.20 1261.10     2 1669.02  4.5838 0.01035 *
# factor(NTtype):factor(PatientType)                1208.20  604.10     2 1669.02  2.1958 0.11160  
# factor(PEtype):factor(PatientType)                 283.76  283.76     1 1213.90  1.0314 0.31003  
# factor(NTtype):factor(PEtype):factor(PatientType) 2122.21 1061.10     2 1669.02  3.8568 0.02132 *

##########--------------------------------------------------------------------Neural Model ID: N-M2 Essential Tremor
model_ET <- lmer(AUC_Measurement ~ factor(NTtype) *PEtype + (1 | subjectNum ),  data = subset(dat_long_AUC, PatientType == 1))
anova(model_ET, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                       Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# factor(NTtype)         723.9  361.96     2 526.01  1.4019 0.247036   
# PEtype                 172.3  172.34     1 353.06  0.6675 0.414476   
# factor(NTtype):PEtype 3242.5 1621.23     2 526.01  6.2793 0.002018 **

###post hocs
pairs_AUC_DA_PEtype_ET <- as.data.frame(pairs(emmeans(model_ET, ~ factor(NTtype)*PEtype, at = list(NTtype = 1), lmer.df = "satterthwaite", adjust = "bonferroni")))
pairs_AUC_NE_PEtype_ET = as.data.frame(pairs(emmeans(model_ET, ~  factor(NTtype)*PEtype, at = list(NTtype = 2), lmer.df = "satterthwaite", adjust = "bonferroni")))
pairs_AUC_5HT_PEtype_ET = as.data.frame(pairs(emmeans(model_ET, ~  factor(NTtype)*PEtype, at = list(NTtype = 3), lmer.df = "satterthwaite", adjust = "bonferroni")))

combined_pairs_AUC_ET <- bind_rows(
  pairs_AUC_DA_PEtype_ET %>% mutate(Comparison = "DA"),
  pairs_AUC_NE_PEtype_ET %>% mutate(Comparison = "NE"),
  pairs_AUC_5HT_PEtype_ET %>% mutate(Comparison = "5HT")
) 
# 1 NTtype1 PEtype1 - NTtype1 PEtype2  5.695859 2.54415 531  2.238807 0.02558195         DA
# 2 NTtype2 PEtype1 - NTtype2 PEtype2 -2.551761 2.54415 531 -1.002992 0.31632163         NE
# 3 NTtype3 PEtype1 - NTtype3 PEtype2 -6.849079 2.54415 531 -2.692090 0.00732470        5HT

##########--------------------------------------------------------------------Neural Model ID: N-M3 Parkinson's disease
model_PD <- lmer(AUC_Measurement ~ factor(NTtype) *PEtype + (1 | subjectNum ),  data = subset(dat_long_AUC, PatientType == 2))
anova(model_PD, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                       Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# factor(NTtype)        511.09 255.545     2 1143.01  0.9033 0.4055
# PEtype                105.82 105.816     1  933.09  0.3740 0.5410
# factor(NTtype):PEtype 100.96  50.482     2 1143.01  0.1784 0.8366

##########--------------------------------------------------------------------Neural Model ID: N-M4 Dopamine
model_DA <- lmer(AUCDA ~ PEtype *PatientType + (1 | subjectNum) ,  data = dat)
anova(model_DA, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                     Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# PEtype             1073.10 1073.10     1 558.31  4.0082 0.04576 *
# PatientType         202.56  202.56     1  18.52  0.7566 0.39553  
# PEtype:PatientType  681.50  681.50     1 558.31  2.5455 0.11117 

##########--------------------------------------------------------------------Neural Model ID: N-M5 Serotonin
model <- lmer(AUC5HT ~ factor(PEtype) *PatientType + (1 | subjectNum) , data = dat)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                     Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# PEtype             1362.84 1362.84     1 545.01  4.9560 0.02641 *
# PatientType         604.94  604.94     1  19.03  2.1999 0.15440  
# PEtype:PatientType 1226.86 1226.86     1 545.01  4.4615 0.03512 *

###post hocs
pairs(emmeans(model, ~  factor(PEtype)*PatientType,at = list(PEtype = 1), lmer.df = "satterthwaite",adjust = "bonferroni"))
# contrast                                    estimate   SE   df t.ratio p.value
# PEtype1 PatientType1 - PEtype1 PatientType2    -5.87 2.58 81.3  -2.281  0.0252

pairs(emmeans(model, ~  PEtype*PatientType,at = list(PEtype = 2), lmer.df = "satterthwaite",adjust = "bonferroni"))
# contrast                                    estimate   SE   df t.ratio p.value
# PEtype2 PatientType1 - PEtype2 PatientType2    0.809 2.03 36.7   0.398  0.6929

##########--------------------------------------------------------------------Neural Model ID: N-M6 Noradrenaline
model <- lmer(AUCNE ~ PEtype *PatientType + (1 | subjectNum) , data = dat)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                    Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# PEtype              42.57   42.57     1 529.91  0.1546 0.6944
# PatientType        454.58  454.58     1  19.38  1.6504 0.2140
# PEtype:PatientType 414.47  414.47     1 529.91  1.5048 0.2205



##########--------------------------------------------------------------------Supplementary Table 3 | Medication status

###ETs on SSRIs (n=3) or not (n=3)
model = lmer(AUC_Measurement ~ 1 + SSRI*NTtype*PEtype + (1 | subjectNum),  data = dat_long_ET)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                     Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
# SSRI                 20.88   20.88     1 169.52  0.0806 0.776858    
# NTtype              577.29  577.29     1 525.01  2.2276 0.136165    
# PEtype             2165.00 2165.00     1 527.95  8.3540 0.004007 ** 
# SSRI:NTtype           0.09    0.09     1 525.01  0.0004 0.984887    
# SSRI:PEtype           0.91    0.91     1 527.95  0.0035 0.952855    
# NTtype:PEtype      3138.55 3138.55     1 525.01 12.1107 0.000543 ***
# SSRI:NTtype:PEtype    0.26    0.26     1 525.01  0.0010 0.974969  

###PDs on SSRIs (n=5) or not (n=5)
model = lmer(AUC_Measurement ~ 1 + SSRI*NTtype*PEtype + (1 | subjectNum),  data = dat_long_PD)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                    Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# SSRI                34.03   34.03     1 331.10  0.1192 0.7301
# NTtype             137.04  137.04     1 964.01  0.4802 0.4885
# PEtype             265.65  265.65     1 973.00  0.9309 0.3349
# SSRI:NTtype         80.09   80.09     1 964.01  0.2806 0.5964
# SSRI:PEtype         76.09   76.09     1 973.00  0.2666 0.6057
# NTtype:PEtype      221.13  221.13     1 964.01  0.7749 0.3789
# SSRI:NTtype:PEtype 389.47  389.47     1 964.01  1.3647 0.2430

###PDs on SSRIs + Buproprion (n=7) or not (n=5)
model = lmer(AUC_Measurement ~ 1 + all_Rx*NTtype*PEtype + (1 | subjectNum),  data = dat_long_PD)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                      Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# all_Rx                82.28   82.28     1  403.72  0.2914 0.5896
# NTtype               244.52  244.52     1 1142.01  0.8659 0.3523
# PEtype               161.18  161.18     1 1153.00  0.5708 0.4501
# all_Rx:NTtype        161.49  161.49     1 1142.01  0.5719 0.4497
# all_Rx:PEtype        193.12  193.12     1 1153.00  0.6839 0.4084
# NTtype:PEtype        147.00  147.00     1 1142.01  0.5206 0.4707
# all_Rx:NTtype:PEtype 607.15  607.15     1 1142.01  2.1501 0.1428

##########--------------------------------------------------------------------Supplementary Table 5 | PD clinical measures.
##UPDRS-III score
model = lmer(AUC_Measurement ~ 1 + UPDP_tail*NTtype*PEtype + (1 | subjectNum),  data = dat_long_PD_score)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# UPDP_tail               322.39 161.195     2   8.02  0.5594 0.5923
# NTtype                  352.43 176.215     2 955.01  0.6116 0.5427
# PEtype                   13.44  13.439     1 884.08  0.0466 0.8291
# UPDP_tail:NTtype        652.28 163.069     4 955.01  0.5659 0.6874
# UPDP_tail:PEtype        268.21 134.107     2 880.42  0.4654 0.6280
# NTtype:PEtype           214.50 107.252     2 955.01  0.3722 0.6893
# UPDP_tail:NTtype:PEtype 875.93 218.983     4 955.01  0.7600 0.5514

##Non-motor scores (NMS-Quest)
model = lmer(AUC_Measurement ~ 1 + non_motor_tail*NTtype*PEtype + (1 | subjectNum),  data = dat_long_PD_score)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                               Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# non_motor_tail                271.58  135.79     2   6.86  0.4970 0.62871  
# NTtype                        130.46   65.23     2 875.01  0.2387 0.78767  
# PEtype                         62.93   62.93     1 542.55  0.2303 0.63148  
# non_motor_tail:NTtype        1226.51  306.63     4 875.01  1.1222 0.34462  
# non_motor_tail:PEtype        1297.33  648.67     2 539.38  2.3739 0.09409 .
# NTtype:PEtype                 281.44  140.72     2 875.01  0.5150 0.59766  
# non_motor_tail:NTtype:PEtype 1000.16  250.04     4 875.01  0.9151 0.45440 

##Disease duration
model = lmer(AUC_Measurement ~ 1 + PD_duration_tail*NTtype*PEtype + (1 | subjectNum),  data = dat_long_PD_score)
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                                Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# PD_duration_tail               610.80 305.401     2   10.39  1.0725 0.3772
# NTtype                         447.22 223.608     2 1133.01  0.7853 0.4562
# PEtype                          77.70  77.705     1 1033.24  0.2729 0.6015
# PD_duration_tail:NTtype        470.47 117.618     4 1133.01  0.4131 0.7993
# PD_duration_tail:PEtype        114.00  57.002     2 1012.22  0.2002 0.8186
# NTtype:PEtype                   96.09  48.045     2 1133.01  0.1687 0.8448
# PD_duration_tail:NTtype:PEtype  66.75  16.688     4 1133.01  0.0586 0.9936
##########--------------------------------------------------------------------Supplementary Figure 4 | Electrode positioning
model_full_anterior <- lmer(AUC_Measurement ~ factor(NTtype) *factor(PEtype)*anterior + (1 | subjectNum) , data = dat_long_AUC)
anova(model_full_anterior, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# factor(NTtype)                          941.72  470.86     2 1669.02  1.7134 0.18057  
# factor(PEtype)                          403.78  403.78     1 1161.98  1.4693 0.22570  
# anterior                                900.74  900.74     1   16.77  3.2777 0.08818 .
# factor(NTtype):factor(PEtype)          1800.85  900.43     2 1669.02  3.2765 0.03800 *
# factor(NTtype):anterior                1341.53  670.76     2 1669.02  2.4408 0.08740 .
# factor(PEtype):anterior                 441.36  441.36     1 1443.39  1.6060 0.20525  
# factor(NTtype):factor(PEtype):anterior 1024.65  512.32     2 1669.02  1.8643 0.15533 

##########--------------------------------------------------------------------Supplementary Figure 5 | PD Rejection rate sub groups

#####PD patients that had 0% rejection rate (n =5) compared to PD patients that had >0% rejection rate (n=7) during positive NPEs 
model <- lmer( AUC_Measurement ~ factor(NTtype) * factor(reject_rate_0)*PEtype + (1 | subjectNum),  data = subset(dat_rejection_long_AUC,  PatientType == 2))
anova(model, ddf = "Kenward-Roger")
#                                               Sum Sq   Mean Sq NumDF      DenDF    F value    Pr(>F)
# factor(NTtype)                              804.98959 402.49480     2 1138.01268 1.42066624 0.2419811
# factor(reject_rate_0)                       122.44415 122.44415     1   10.79944 0.43218515 0.5246825
# PEtype                                      118.31870 118.31870     1  866.64958 0.41762373 0.5182962
# factor(NTtype):factor(reject_rate_0)        820.25724 410.12862     2 1138.01268 1.44761098 0.2355640
# factor(NTtype):PEtype                        55.35323  27.67662     2 1138.01268 0.09768880 0.9069387
# factor(reject_rate_0):PEtype                 12.21233  12.21233     1  866.64958 0.04310527 0.8355760
# factor(NTtype):factor(reject_rate_0):PEtype 354.70480 177.35240     2 1138.01268 0.62599211 0.5349147



############ PD patients with a 0% rejection rate during positive prediction errors, and the ET patient group.

model <- lmer( AUC_Measurement ~ factor(NTtype) * factor(PatientType) + (1 | subjectNum),
               data = subset(dat_rejection_long_AUC, PEtype  == 1 & reject_rate_0  == 1))
anova(model, ddf = "Kenward-Roger")
# Type III Analysis of Variance Table with Kenward-Roger's method
#                                    Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# factor(NTtype)                      245.9  122.93     2 373.03  0.4828 0.617463   
# factor(PatientType)                 536.5  536.51     1   7.62  2.1069 0.186534   
# factor(NTtype):factor(PatientType) 3506.1 1753.05     2 373.03  6.8843 0.001159 **

####post hocs
pairs_AUC_DA_reject=as.data.frame(pairs(emmeans(model, ~ factor(NTtype)*factor(PatientType), at = list( NTtype =1), lmer.df = "satterthwaite", adjust = "bonferroni")))
pairs_AUC_NE_reject=as.data.frame(pairs(emmeans(model, ~ factor(NTtype)*factor(PatientType), at = list( NTtype =2), lmer.df = "satterthwaite", adjust = "bonferroni")))
pairs_AUC_5HT_reject=as.data.frame(pairs(emmeans(model, ~ factor(NTtype)*factor(PatientType), at = list(NTtype =3), lmer.df = "satterthwaite", adjust = "bonferroni")))

combined_pairs_ET_PD_0_reject <- bind_rows(
  pairs_AUC_DA_reject %>% mutate(Comparison = "DA"),
  pairs_AUC_NE_reject %>% mutate(Comparison = "NE"),
  pairs_AUC_5HT_reject %>% mutate(Comparison = "5HT")
) 
#                                      contrast  estimate       SE  df   t.ratio     p.value Comparison
# 1 NTtype1 PatientType1 - NTtype1 PatientType2  5.739647 2.816828 381  2.037627 0.042277624         DA
# 2 NTtype2 PatientType1 - NTtype2 PatientType2 -4.284593 2.816828 381 -1.521070 0.129071747         NE
# 3 NTtype3 PatientType1 - NTtype3 PatientType2 -8.680337 2.816828 381 -3.081599 0.002208588        5HT



##########--------------------------------------------------------------------Plotting Fig. 1f


emo_plot <- ggplot(summary_stats_emo, aes(x = PEtype, y = mean_emo_rate, fill = PatientType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_emo_rate - sem_emo_rate, ymax = mean_emo_rate + sem_emo_rate),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(data = combined_emorating_by_subject, 
             aes(x = PEtype, y = ratings, color = PatientType),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
             size = 2, 
             shape = 21, 
             fill = "black") +  # Points for each individual emoRating
  geom_hline(yintercept = 0, color = "black", size = 0.8) +
  labs(
    title = "emotional rating by PEtype and Patient Type",
    x = "",
    y = "Emotional Rating"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("1" = rgb(80/255, 233/255, 145/255), "2" = rgb(155/255, 25/255, 245/255))) +  
  scale_color_manual(values = c("1" = rgb(80/255, 233/255, 145/255), "2" = rgb(155/255, 25/255, 245/255))) +  
  theme(
    axis.text.x = element_text(size = 25), 
    axis.text.y = element_text(size = 25), 
    axis.title = element_text(size = 14),
    legend.position = "none", 
    panel.border = element_blank(),
    axis.line.y = element_line(color = "black", size = 1)
  ) +
  scale_x_discrete(labels = c("NegPE" = "Neg", "PosPE" = "Pos")) +
  scale_y_continuous(
    breaks = c(0,3, 6, 9), 
    limits = c(0, 9)   
  )
print(emo_plot)

##########--------------------------------------------------------------------Plotting Fig. 2c

####updated AUC averages line plots rather than bar graphs.

# Function to calculate SEM
sem <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

# Filter data for PatientType == 1 (ET)
dat_ET <- subset(dat, PatientType == 1)

# Compute group-level means and SEMs
mean_AUCDA <- tapply(dat_ET$AUCDA, dat_ET$PEtype, mean, na.rm = TRUE)
sem_AUCDA <- tapply(dat_ET$AUCDA, dat_ET$PEtype, sem)

mean_AUCNE <- tapply(dat_ET$AUCNE, dat_ET$PEtype, mean, na.rm = TRUE)
sem_AUCNE <- tapply(dat_ET$AUCNE, dat_ET$PEtype, sem)

mean_AUC5HT <- tapply(dat_ET$AUC5HT, dat_ET$PEtype, mean, na.rm = TRUE)
sem_AUC5HT <- tapply(dat_ET$AUC5HT, dat_ET$PEtype, sem)

# Create plotting data
data_plot <- data.frame(
  Neurotransmitter = rep(c("AUCDA", "AUCNE", "AUC5HT"), each = 2),
  PEtype = rep(c("posPE", "negPE"), times = 3),
  Mean = c(mean_AUCDA, mean_AUCNE, mean_AUC5HT),
  SEM = c(sem_AUCDA, sem_AUCNE, sem_AUC5HT)
)
# ─────────────────────────────────────────────────────────────────────────────
# (1) Make sure PEtype is ordered so that negPE comes first, then posPE:
data_plot$PEtype <- factor(
  data_plot$PEtype,
  levels = c("posPE", "negPE")
)

ggplot() +
  
  # ── (a) Group‐level lines + error bars ───────────────────────────────────────
  geom_line(
    data = data_plot,
    aes(
      x = PEtype,
      y = Mean,
      group = Neurotransmitter,
      color = Neurotransmitter
    ),
    size = 1
  ) +
  geom_point(
    data = data_plot,
    aes(
      x = PEtype,
      y = Mean,
      color = Neurotransmitter
    ),
    size = 3
  ) +
  geom_errorbar(
    data = data_plot,
    aes(
      x = PEtype,
      ymin = Mean - SEM,
      ymax = Mean + SEM,
      color = Neurotransmitter
    ),
    width = 0.1,
    size = 0.8
  ) +
  
  # ── (c) Color mapping: black = DA, cyan = NE, magenta = 5HT ────────────────
  scale_color_manual(
    values = c(
      "AUCDA"  = "black",
      "AUCNE"  = "cyan",
      "AUC5HT" = "magenta"
    )
  ) +
  scale_fill_manual(
    values = c(
      "AUCDA"  = "black",
      "AUCNE"  = "cyan",
      "AUC5HT" = "magenta"
    )
  ) +
  
  # ── (d) Axis labels / theme tweaks ─────────────────────────────────────────
  scale_x_discrete(
    labels = c("negPE" = "Neg", "posPE" = "Pos")
  ) +
  labs(
    x = "NPE Type",
    y = "AUC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    axis.title       = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.15, "cm"), 
    aspect.ratio      = 1
    
  )+
  scale_x_discrete(expand = c(0.12, 0.12)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  coord_cartesian(ylim = c(-6.8, 6.5))  # Preferred over `ylim()` to avoid dropping data


##########--------------------------------------------------------------------Plotting Fig. 2f

# Function to calculate SEM
sem <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

# Filter data for PatientType == 1 (ET)
dat_PD <- subset(dat, PatientType == 2)

# Compute group-level means and SEMs
mean_AUCDA <- tapply(dat_PD$AUCDA, dat_PD$PEtype, mean, na.rm = TRUE)
sem_AUCDA <- tapply(dat_PD$AUCDA, dat_PD$PEtype, sem)

mean_AUCNE <- tapply(dat_PD$AUCNE, dat_PD$PEtype, mean, na.rm = TRUE)
sem_AUCNE <- tapply(dat_PD$AUCNE, dat_PD$PEtype, sem)

mean_AUC5HT <- tapply(dat_PD$AUC5HT, dat_PD$PEtype, mean, na.rm = TRUE)
sem_AUC5HT <- tapply(dat_PD$AUC5HT, dat_PD$PEtype, sem)

# Create plotting data
data_plot <- data.frame(
  Neurotransmitter = rep(c("AUCDA", "AUCNE", "AUC5HT"), each = 2),
  PEtype = rep(c("posPE", "negPE"), times = 3),
  Mean = c(mean_AUCDA, mean_AUCNE, mean_AUC5HT),
  SEM = c(sem_AUCDA, sem_AUCNE, sem_AUC5HT)
)
# ─────────────────────────────────────────────────────────────────────────────
# (1) Make sure PEtype is ordered so that negPE comes first, then posPE:
data_plot$PEtype <- factor(
  data_plot$PEtype,
  levels = c("posPE", "negPE")
)

ggplot() +
  
  # ── (a) Group‐level lines + error bars ───────────────────────────────────────
  geom_line(
    data = data_plot,
    aes(
      x = PEtype,
      y = Mean,
      group = Neurotransmitter,
      color = Neurotransmitter
    ),
    size = 1
  ) +
  geom_point(
    data = data_plot,
    aes(
      x = PEtype,
      y = Mean,
      color = Neurotransmitter
    ),
    size = 3
  ) +
  geom_errorbar(
    data = data_plot,
    aes(
      x = PEtype,
      ymin = Mean - SEM,
      ymax = Mean + SEM,
      color = Neurotransmitter
    ),
    width = 0.1,
    size = 0.8
  ) +
  
  # ── (c) Color mapping: black = DA, cyan = NE, magenta = 5HT ────────────────
  scale_color_manual(
    values = c(
      "AUCDA"  = "black",
      "AUCNE"  = "cyan",
      "AUC5HT" = "magenta"
    )
  ) +
  scale_fill_manual(
    values = c(
      "AUCDA"  = "black",
      "AUCNE"  = "cyan",
      "AUC5HT" = "magenta"
    )
  ) +
  
  # ── (d) Axis labels / theme tweaks ─────────────────────────────────────────
  scale_x_discrete(
    labels = c("negPE" = "Neg", "posPE" = "Pos")
  ) +
  labs(
    x = "NPE Type",
    y = "AUC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    axis.title       = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    aspect.ratio      = 1
  )+
  scale_x_discrete(expand = c(0.12, 0.12)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  coord_cartesian(ylim = c(-6.8, 6.5))  # Preferred over `ylim()` to avoid dropping data


##########--------------------------------------------------------------------Plotting Supplementary Figure 2c
######----------reaction time  by patienttype and PEtype

# Calculate rejection rate per subjectNum within each PEtype for PD
reaction_by_subject_PD <- dat_PD %>%
  group_by(PEtype, subjectNum) %>%
  summarize(
    total_decisions = n(),
    reaction = mean(reactionTime, na.rm = TRUE),
  )

# Add PatientType column for PD
reaction_by_subject_PD <- reaction_by_subject_PD %>%
  mutate(PatientType = "PD") %>%
  mutate(PEtype = ifelse(PEtype == 1, "Pos", "Neg"))

# Calculate rejection rate per subjectNum within each PEtype for ET
reaction_by_subject_ET <- dat_ET %>%
  group_by(PEtype, subjectNum) %>%
  summarize(
    total_decisions = n(),
    reaction = mean(reactionTime, na.rm = TRUE),
  )

# Add PatientType column for ET
reaction_by_subject_ET <- reaction_by_subject_ET %>%
  mutate(PatientType = "ET") %>%
  mutate(PEtype = ifelse(PEtype == 1, "Pos", "Neg"))

# Combine both datasets
combined_reactionTime_by_subject <- bind_rows(reaction_by_subject_PD, reaction_by_subject_ET)

# Calculate mean and SEM for each PEtype and PatientType
summary_stats_reaction <- dat %>%
  group_by(PEtype, PatientType) %>%
  summarize(
    mean_reaction = mean(reactionTime, na.rm = TRUE),
    sem_reaction = sd(reactionTime, na.rm = TRUE) / sqrt(n())  # SEM calculation
  )

summary_stats_reaction <- summary_stats_reaction %>%
  mutate(PEtype = ifelse(PEtype == 1, "Pos", "Neg"))
# Create bar plot with custom colors, no legend, and scatter dots for each subject
reaction_plot <- ggplot(summary_stats_reaction, aes(x = PEtype, y = mean_reaction, fill = PatientType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_reaction - sem_reaction, ymax = mean_reaction + sem_reaction),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(data = combined_reactionTime_by_subject, 
             aes(x = PEtype, y = reaction, color = PatientType),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
             size = 2, 
             shape = 21, 
             fill = "black") +  # Points for each individual rejection rate
  labs(
    title = "reaction time by PEtype and Patient Type",
    x = "Prediction Error",
    y = "Reaction time (seconds)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("1" = rgb(80/255, 233/255, 145/255), "2" = rgb(155/255, 25/255, 245/255))) +  # Custom colors
  scale_color_manual(values = c("1" = rgb(80/255, 233/255, 145/255), "2" = rgb(155/255, 25/255, 245/255))) +  # Point colors
  theme(
    axis.text.x = element_text(size = 25), 
    axis.text.y = element_text(size = 25), 
    axis.title = element_text(size = 14),
    legend.position = "none",   # Remove the legend
    plot.margin = unit(c(0, 0, 0, 0), "lines"),  # Remove the margin around the plot
    panel.background = element_rect(color = "black", size = 1),  # Border around the panel
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

print(reaction_plot)


##########--------------------------------------------------------------------Plotting Supplementary Figure 2d

combined_rejection_by_PEtype <- combined_rejection_by_PEtype %>%
  mutate(PEtype = recode(PEtype, `1` = "Pos", `2` = "Neg"))

# Calculate mean and SEM for each PEtype and PatientType
summary_stats_rejection <- combined_rejection_by_PEtype %>%
  group_by(PEtype, PatientType) %>%
  summarize(
    mean_rejection_rate = mean(rejection_rate, na.rm = TRUE),
    sem_rejection_rate = sd(rejection_rate, na.rm = TRUE) / sqrt(n())  # SEM calculation
  )
# Ensure PEtype is a factor with the desired order
summary_stats_rejection$PEtype <- factor(summary_stats_rejection$PEtype, levels = c("Neg", "Pos" ))
combined_rejection_by_PEtype$PEtype <- factor(combined_rejection_by_PEtype$PEtype, levels = c("Neg", "Pos" ))

# Create bar plot with custom colors, no legend, and scatter dots for each subject
rejection_rate_plot <- ggplot(summary_stats_rejection, aes(x = PEtype, y = mean_rejection_rate, fill = PatientType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rejection_rate - sem_rejection_rate, ymax = mean_rejection_rate + sem_rejection_rate),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(data = combined_rejection_by_PEtype, 
             aes(x = PEtype, y = rejection_rate, color = PatientType),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
             size = 2, 
             shape = 21, 
             fill = "black") +  # Points for each individual rejection rate
  labs(
    title = "Rejection Rate by PEtype and Patient Type",
    x = "Prediction Error",
    y = "Rejection Rate (%)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("ET" = rgb(80/255, 233/255, 145/255), "PD" = rgb(155/255, 25/255, 245/255))) +  # Custom colors
  scale_color_manual(values = c("ET" = rgb(80/255, 233/255, 145/255), "PD" = rgb(155/255, 25/255, 245/255))) +  # Point colors
  theme(
    axis.text.x = element_text(size = 25), 
    axis.text.y = element_text(size = 25), 
    axis.title = element_text(size = 14),
    legend.position = "none",  # Remove the legend
    plot.margin = unit(c(0, 0, 0, 0), "lines"),  # Remove the margin around the plot
    panel.background = element_rect(color = "black", size = 1),  # Border around the panel
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
print(rejection_rate_plot)
##########--------------------------------------------------------------------Plotting Supplementary Figure 4


ggplot(dat_coord, aes(x = anterior, fill = PatientType)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("ET" = "#50E991", "PD" = "#9B19F5")) +
  labs(title = "Distribution of Anterior Values by Patient Type", x = "Anterior relative to Midcommissural Point (mm)", y = "Density") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),  # Removes the box around the panel
    axis.line = element_line(color = "black", size = 0.4),  # Adds lines on the x and y axes
    axis.ticks.x = element_line(size = .5),                            # Adds x-axis ticks
    axis.ticks.y = element_line(size = .5),                            # Adds x-axis ticks
    
    axis.ticks.length = unit(0.1, "cm")                                # Sets length of tick marks
  ) 


##########--------------------------------------------------------------------Plotting Supplementary Figure 5
######subset by PD comparing 0% rejection rates

dat_PD_subset <- dat_PD[
  dat_PD$subjectNum %in% PD_0_reject_pos_subjects,
]

# Compute group-level means and SEMs
mean_AUCDA <- tapply(dat_PD_subset$AUCDA, dat_PD_subset$PEtype, mean, na.rm = TRUE)
sem_AUCDA <- tapply(dat_PD_subset$AUCDA, dat_PD_subset$PEtype, sem)

mean_AUCNE <- tapply(dat_PD_subset$AUCNE, dat_PD_subset$PEtype, mean, na.rm = TRUE)
sem_AUCNE <- tapply(dat_PD_subset$AUCNE, dat_PD_subset$PEtype, sem)

mean_AUC5HT <- tapply(dat_PD_subset$AUC5HT, dat_PD_subset$PEtype, mean, na.rm = TRUE)
sem_AUC5HT <- tapply(dat_PD_subset$AUC5HT, dat_PD_subset$PEtype, sem)

# Create plotting data
data_plot <- data.frame(
  Neurotransmitter = rep(c("AUCDA", "AUCNE", "AUC5HT"), each = 2),
  PEtype = rep(c("posPE", "negPE"), times = 3),
  Mean = c(mean_AUCDA, mean_AUCNE, mean_AUC5HT),
  SEM = c(sem_AUCDA, sem_AUCNE, sem_AUC5HT)
)
# ─────────────────────────────────────────────────────────────────────────────
# (1) Make sure PEtype is ordered so that negPE comes first, then posPE:
data_plot$PEtype <- factor(
  data_plot$PEtype,
  levels = c("posPE", "negPE")
)

ggplot() +
  
  # ── (a) Group‐level lines + error bars ───────────────────────────────────────
  geom_line(
    data = data_plot,
    aes(
      x = PEtype,
      y = Mean,
      group = Neurotransmitter,
      color = Neurotransmitter
    ),
    size = 1
  ) +
  geom_point(
    data = data_plot,
    aes(
      x = PEtype,
      y = Mean,
      color = Neurotransmitter
    ),
    size = 3
  ) +
  geom_errorbar(
    data = data_plot,
    aes(
      x = PEtype,
      ymin = Mean - SEM,
      ymax = Mean + SEM,
      color = Neurotransmitter
    ),
    width = 0.1,
    size = 0.8
  ) +
    scale_color_manual(
    values = c(
      "AUCDA"  = "black",
      "AUCNE"  = "cyan",
      "AUC5HT" = "magenta"
    )
  ) +
  scale_fill_manual(
    values = c(
      "AUCDA"  = "black",
      "AUCNE"  = "cyan",
      "AUC5HT" = "magenta"
    )
  ) +
  
  # ── (d) Axis labels / theme tweaks ─────────────────────────────────────────
  scale_x_discrete(
    labels = c("posPE" = "Pos", "negPE" = "Neg")
  ) +
  labs(
    x = "NPE Type",
    y = "AUC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    axis.title       = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    aspect.ratio      = 1
  )+
  scale_x_discrete(expand = c(0.12, 0.12)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  coord_cartesian(ylim = c(-6.8, 6.5))  




