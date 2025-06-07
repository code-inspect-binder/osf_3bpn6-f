#################### Description ####################

# Author: Andrius Pašukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-05-08

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Stats for poison frog space use


#################### Load packages ####################

library(ggplot2) #plots
library(dplyr) #deal with data easy
library(lme4) #GLMMs
library(emmeans) #pairwise comparisons
library(multcomp)
library(lmerTest) #pvalues for GLMMs
library(sjPlot) #for model plots
library(car) #for qqp
library(coefplot) #for model coef plots
library(MuMIn) #for model selection using AICc

#################### Import data ####################

# Set your working directory 
setwd("~/")

frog_daily <- read.csv(file="all_spaceuse_daily.csv", header = T, sep =",") # cumulative daily distancess with covariates 
frog_hr <- read.csv(file="all_spaceuse_hr.csv", header = T, sep =",") # home range with covariates
injuries <- read.csv(file="all_spaceuse_injuries.csv", header = T, sep =",") # tag-related and natural injuries

frog_hr <- left_join(frog_hr, injuries)
frog_daily <- left_join(frog_daily, injuries)

# subset by species
af_hr <- frog_hr[frog_hr$species=="afemo",]
af_daily <- frog_daily[frog_daily$species=="afemo",]

# Because injury status has effect on space use parameters in A. femoralis
# we remove two individuals that had skin lesions (injury_cat_2 = 2)
af_hr <- filter(af_hr, injury_cat_2 < 2 | is.na(injury_cat_2)) #filter out the frogs with deep wounds
af_daily <- filter(af_daily, injury_cat_2 < 2 | is.na(injury_cat_2)) #filter out the frogs with deep wounds

dt_hr <- frog_hr[frog_hr$species=="dtinc",]
dt_daily <- frog_daily[frog_daily$species=="dtinc",]

os_hr <- frog_hr[frog_hr$species=="osylv",]
os_daily <- frog_daily[frog_daily$species=="osylv",]

#################### Models: Home Range (UD95_HPI) ####################
#################### All species ####################

## Data exploration
ggplot(frog_hr, aes(x=log(ud95_hpi), fill=species)) + geom_density(alpha=.3) # log transformed distributions
cor.test(frog_hr$ud95_hpi, frog_hr$relocs) # no correlation between total relocs and ud95
cor.test(frog_hr$ud95_hpi, frog_hr$days_tracked) # no correlation between number of days tracked and ud95
cor.test(frog_hr$ud95_hpi, frog_hr$duration) # positive corelation between tracking span (including days not tracked) and ud95
# use duration as model weight where appropriate

## Full LM for log transformed data
m11_allud <- lm(log(ud95_hpi) ~ species*sex + duration,
               data = frog_hr)
# Diagnostics: based on Cook's d Point 66 is an outlier

# Excluding the outlier pt 66 (D. tinctorius t08) because it showed long unidirectional movement inflating HR values (no site fidelity)
# Check if the results change when removing these
frog_hr_outliers <- frog_hr[-c(66),]

## LM for log transformed data and outliers removed
m11_allud_1 <- lm(log(ud95_hpi) ~ species*sex + duration,
            data = frog_hr_outliers)
# Diagnostics look good

## Compare model results with and without outliers
Anova(m11_allud)
Anova(m11_allud_1) #removing "outlier" does not change the results
frog_hr[66,] #This is the crazy runner D. tinctorius t08 Jane (no site fidelity). Proceed with Jane excluded.

# Summary of the final best model
Anova(m11_allud_1) #species and sex interaction very significant
summary(m11_allud_1)
#################### A. femoralis ####################

cor.test(af_hr$ud95_hpi, af_hr$relocs) # no correlation between total relocs and ud95
cor.test(af_hr$ud95_hpi, af_hr$days_tracked) # no correlation between number of days tracked and ud95
cor.test(af_hr$ud95_hpi, af_hr$duration) # no correlation between tracking span (including days not tracked) and ud95

## Full LM for log transformed data
m12_afud <- lm(log(ud95_hpi) ~ sex + duration,
            data = af_hr) 
# Diagnostics: look good

# Summary of the final best model
summary(m12_afud) # sex is significant
# males have larger HR than females
#################### D. tinctorius ####################

cor.test(dt_hr$ud95_hpi, dt_hr$relocs) # no correlation between total relocs and ud95
cor.test(dt_hr$ud95_hpi, dt_hr$days_tracked) # no correlation between number of days tracked and ud95
cor.test(dt_hr$ud95_hpi, dt_hr$duration) # no correlation between tracking span (including days not tracked) and ud95

## Full LM for log transformed data
m13_dtud <- lm(log(ud95_hpi) ~ sex + duration,
           data = dt_hr) 
# Diagnostics: unequal variance. based ok Cook's d pt 8 is a clear pt of influence
# Excluding the outlier pt 66 (D. tinctorius t08) because it showed long unidirectional movement inflating HR values (no site fidelity).
dt_hr_outliers <- filter(dt_hr, !id == "t08")

## LM for log transformed data with outliers removed
m13_dtud_1 <- lm(log(ud95_hpi) ~ sex + duration,
           data = dt_hr_outliers)
# Diagnostics look good
# Compare model results with and without outliers
summary(m13_dtud)
summary(m13_dtud_1) #removing "outlier" does not change the results

# Summary of the final best model
summary(m13_dtud_1) 
# sex has no influence, R2 is close to 0

#################### O. sylvatica #################### 

cor.test(os_hr$ud95_hpi, os_hr$relocs) # no correlation between total relocs and ud95
cor.test(os_hr$ud95_hpi, os_hr$days_tracked) # no correlation between number of days tracked and ud95
cor.test(os_hr$ud95_hpi, os_hr$duration) # no correlation between tracking span (including days not tracked) and ud95

# Full LM for log transformed data
m14_osud <- lm(log(ud95_hpi) ~ sex + duration,
            data = os_hr) 
# Diagnostics look good
# Summary of the final best model
summary(m14_osud) 
# sex has a strong influence on HR

#### Home range model tables ####
tab_model(m12_afud, m13_dtud_1, m14_osud, collapse.ci = TRUE,
          dv.labels = c("A. femoralis home rage",
                            "D. tinctorus home rage",
                            "O. sylvatica home range"))
#################### Models: movement extent (MCP100) #################### 
#################### All species ####################

cor.test(frog_hr$mcp100, frog_hr$relocs) # no correlation between total relocs and mcp100
cor.test(frog_hr$mcp100, frog_hr$days_tracked) # no correlation between number of days tracked and mcp100
cor.test(frog_hr$mcp100, frog_hr$duration) # positive corelation between tracking span (including days not tracked) and mcp100

# Excluding the outlier pt 66 (D. tinctorius t08) because it showed long unidirectional movement inflating HR values
frog_hr_outliers <- filter(frog_hr, !id == "t08")
## Full LM for log transformed data and pt66 removed
m15_allmcp <- lm(log(mcp100) ~ species*sex + duration,
           data = frog_hr_outliers)
# Diagnostics look good

# Summary of the final model
summary(m15_allmcp)
# model anova with car package
Anova(m15_allmcp) 
# species and sex interaction very significant
# proceed for separate analysis per species

#################### A. femoralis ####################

cor.test(af_hr$mcp100, af_hr$relocs) # no correlation between total relocs and mcp100
cor.test(af_hr$mcp100, af_hr$days_tracked) # correlation between number of days tracked and mcp100
cor.test(af_hr$mcp100, af_hr$duration) # correlation between tracking span (including days not tracked) and mcp100

## Full LM for log transformed data
m16_afmcp <- lm(log(mcp100) ~ sex + duration,
            data = af_hr) 
# Diagnostics look good
# Summary
summary(m16_afmcp) 
# sex has significant influence on movement extent

#################### D. tinctorius ####################

cor.test(dt_hr$mcp100, dt_hr$relocs) # no correlation between total relocs and mcp100
cor.test(dt_hr$mcp100, dt_hr$days_tracked) # no correlation between number of days tracked and mcp100
cor.test(dt_hr$mcp100, dt_hr$duration) # no correlation between tracking span (including days not tracked) and mcp100

## Full LM for log transformed data
m17_dtmcp <- lm(log(mcp100) ~ sex + duration,
          data = dt_hr) 
# Diagnostics: based ok Cook's d pt 8 is a clear pt of influence
# Excluding the outlier pt 66 (D. tinctorius t08) because it showed long unidirectional movement inflating HR values (no site fidelity).
# Check if the results change when removing these
dt_hr_outliers <- filter(dt_hr, !id == "t08")

## LM for log transformed data and outliers removed
m17_dtmcp_1 <- lm(log(mcp100) ~ sex + duration,
             data = dt_hr_outliers)

## Compare model results with and without outliers
summary(m17_dtmcp)
summary(m17_dtmcp_1)
# removing "outlier" does not change the results

## Summary of the final best model
summary(m17_dtmcp_1) 
# sex has no influence, R2 is close to 0

#################### O. sylvatica ####################

cor.test(os_hr$mcp100, os_hr$relocs) # no correlation between total relocs and mcp100
cor.test(os_hr$mcp100, os_hr$days_tracked) # no correlation between number of days tracked and mcp100
cor.test(os_hr$mcp100, os_hr$duration) # no correlation between tracking span (including days not tracked) and mcp100

## Full LM for log transformed data
m18_osmcp <- lm(log(mcp100) ~ sex + duration,
            data = os_hr) 
# Diagnostics look good
# Summary 
summary(m18_osmcp) 
# sex has very a strong influence on mcp

#### Movement extent model tables ####
tab_model(m16_afmcp, m17_dtmcp_1, m18_osmcp, collapse.ci = TRUE,
          dv.labels = c("A. femoralis Log movement extent",
                        "D. tinctorus Log movement extent",
                        "O. sylvatica Log movement extent"))

#################### Models for daily cumulative distance #################### 

#################### All species #################### 

## Data exploration
ggplot(frog_daily, aes(x=daily_dist, fill=species)) + geom_density(alpha=.3) # data distribution by species
ggplot(frog_daily, aes(x=log(daily_dist + 1), fill=species)) + geom_density(alpha=.3) # log transform (add 1 to avoid negative log) 
nrow(frog_daily[frog_daily$daily_dist==0,])  # n rows with zeros
cor.test(frog_daily$daily_dist, frog_daily$relocs) # weak positive correlation between number of relocs per day and daily distance

## Full LMM for log transformed daily distance
m1_alldl <- lmer(log(daily_dist + 1) ~ sex*species +  (1|id) + (1|date) + relocs,
                 REML = TRUE, data = frog_daily) 
# Model summaries
summary(m1_alldl)
# model Anova with car package
Anova(m1_alldl) # sex and species interaction significant
# The effect of sex varies strongly between species (i.e. significant interaction sex * species)
# We analyse each species separately to understand the sex difference within each species

#################### A. femoralis ####################

cor.test(af_daily$daily_dist, af_daily$relocs) 
# no correlation between relocs and daily dist
# we will not use relocs as covariate for A. femoralis

## Full LMM for log transformed daily distance
m2_afdl_1 <- lmer(log(daily_dist + 1) ~ sex + behavior + temp + (1|id) + (1|date),
                  REML = TRUE, data = af_daily) 
## AICc model selection
# No temperature
m2_afdl_2 <- lmer(log(daily_dist + 1) ~ sex + behavior +  (1|id) + (1|date),
                  REML = TRUE, data = af_daily)
# No sex
m2_afdl_3 <- lmer(log(daily_dist + 1) ~ behavior + temp + (1|id) + (1|date),
                  REML = TRUE, data = af_daily)
# No behavior
m2_afdl_4 <- lmer(log(daily_dist + 1) ~ sex + temp + (1|id) + (1|date),
                  REML = TRUE, data = af_daily)
# No sex and behavior
m2_afdl_5 <- lmer(log(daily_dist + 1) ~ temp + (1|id) + (1|date),
                  REML = TRUE, data = af_daily)
# No sex and temperature
m2_afdl_6  <- lmer(log(daily_dist + 1) ~ behavior +  (1|id) + (1|date),
                   REML = TRUE, data = af_daily)
# No behavior and temperature
m2_afdl_7  <- lmer(log(daily_dist + 1) ~ sex + (1|id) + (1|date),
                   REML = TRUE, data = af_daily)
# Null model: no sex, behavior and temperature
m2_afdl_8  <- lmer(log(daily_dist + 1) ~ (1|id) + (1|date),
                   REML = TRUE, data = af_daily)

# model selection 
model.sel(m2_afdl_1, m2_afdl_2, m2_afdl_3, m2_afdl_4,
          m2_afdl_5,  m2_afdl_6,  m2_afdl_7, m2_afdl_8)
# Best model (by far) is m2_afdl_1 (weight 0.9)
# Sex, behavior and temperature all have an effect on daily distance
# Males move more than females. 
# Frogs movements differ by behavior. 
# Frogs move less with higher temp.

# Split by sex for analyzing behavior influence on movement
# because tadpoles transport only occured once in female A. femoralis
afm <- subset(af_daily, af_daily$sex == "m")
aff <- subset(af_daily, af_daily$sex == "f")

## LMM for log transformed daily speed: males by behavior 
m3_afmdl <- lmer(log(daily_dist + 1) ~ behavior + temp + (1|id) + (1|date),
                 REML = TRUE, data = afm) 

# Summary for males
summary(m3_afmdl)
pairs(emmeans(m3_afmdl, ~behavior)) 
# pool_move > others, base = w_op_sex  

## LMM for log transformed daily distance: females by behavior 
m4_affdl<- lmer(log(daily_dist + 1) ~ behavior + temp + (1|id) + (1|date),
                REML = TRUE, data = aff)

# Summary for females
summary(m4_affdl)
pairs(emmeans(m4_affdl, ~behavior)) 
# Mating > others, pool_move only occured once in females

#################### D. tinctorius ####################

cor.test(dt_daily$daily_dist, dt_daily$relocs)
# positive correlation between relocs and daily dist
# we will use relocs as a covariate in all models

## Full LMM for log transformed daily speed
m5_dtdl_1 <- lmer(log(daily_dist + 1) ~ behavior + sex +  temp + (1|id) + (1|date) + relocs,
                   REML = TRUE, data = dt_daily) 
# Diagnostics: random effect residuals look skewed but no better alternative. 
Anova(m5_dtdl_1)
summary(m5_dtdl_1)

## AICc model selection
# No temperature
m5_dtdl_2 <- lmer(log(daily_dist + 1) ~ sex + behavior +  (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)
# No sex
m5_dtdl_3 <- lmer(log(daily_dist + 1) ~ behavior + temp + (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)
# No behavior
m5_dtdl_4 <- lmer(log(daily_dist + 1) ~ sex + temp + (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)
# No sex and behavior
m5_dtdl_5 <- lmer(log(daily_dist + 1) ~ temp + (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)
# No sex and temperature
m5_dtdl_6 <- lmer(log(daily_dist + 1) ~ behavior +  (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)
# No behavior and temperature
m5_dtdl_7 <- lmer(log(daily_dist + 1) ~ sex + (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)
# Null model: no sex, behavior and temperature
m5_dtdl_8 <- lmer(log(daily_dist + 1) ~ sex + (1|id) + (1|date), weight = relocs, REML = TRUE, data = dt_daily)

model.sel(m5_dtdl_1, m5_dtdl_2, m5_dtdl_3, m5_dtdl_4,
          m5_dtdl_5, m5_dtdl_6, m5_dtdl_7, m5_dtdl_8)
# m5_dtdl_6 is the best (weight 0.8)
# Frogs movements differ by behavior. 
# Temperature and sex have no influence

# Split by sex for analyzing behavior influence on movement
# no tadpole transport in females
dtm <- subset(dt_daily, dt_daily$sex == "m")
dtf <- subset(dt_daily, dt_daily$sex == "f")

## LMM for log transformed daily distance: males by behavior 
m6_dtmdl <- lmer(log(daily_dist + 1) ~ behavior + (1|id) + (1|date) + relocs,
                  REML = TRUE, data = dtm) 
# Summary for males
summary(m6_dtmdl)
pairs(emmeans(m6_dtmdl, ~behavior)) 
# Pool move significantly increases the daily distance

## LMM for log transformed daily distance: females by behavior 
m7_dtfdl <- lmer(log(daily_dist + 1) ~ behavior + (1|id) + relocs,
                REML = TRUE, data = dtf) 
# Diagnostics look good

# Summary for females
summary(m7_dtfdl)
pairs(emmeans(m7_dtfdl , ~behavior)) 
# No effect of mating behavior 

#################### O. sylvatica #################### 

cor.test(os_daily$daily_dist, os_daily$relocs)
# weak positive correlation between relocs and daily dist
# weigh the model by relocs

## Full LMM for log transformed daily speed
m8_osdl  <- lmer(log(daily_dist + 1) ~ sex*behavior +  temp + (1|id) + (1|date) + relocs,
                  REML = TRUE, data = os_daily) 
# Diagnostics: random effecct residuals skewed but no better alternative available 
# model Anova with car package
Anova(m8_osdl) # sex*behavior interaction is not significant
# continue without interaction

## LMM for log transformed daily speed without interaction
m8_osdl_1 <-  lmer(log(daily_dist + 1) ~ sex + behavior + temp + (1|id) + (1|date) + relocs,
                   REML = TRUE, data = os_daily) # weights controls for n relocs

## AICc model selection
# No temperature
m8_osdl_2 <- lmer(log(daily_dist + 1) ~ sex + behavior +  (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)
# No sex
m8_osdl_3 <- lmer(log(daily_dist + 1) ~ behavior + temp + (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)
# No behavior
m8_osdl_4 <- lmer(log(daily_dist + 1) ~ sex + temp + (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)
# No sex and behavior
m8_osdl_5 <- lmer(log(daily_dist + 1) ~ temp + (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)
# No sex and temperature
m8_osdl_6 <- lmer(log(daily_dist + 1) ~ behavior +  (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)
# No behavior and temperature
m8_osdl_7 <- lmer(log(daily_dist + 1) ~ sex + (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)
# Null model: no sex, behavior and temperature
m8_osdl_8 <- lmer(log(daily_dist + 1) ~ sex + (1|id) + (1|date), weight = relocs, REML = TRUE, data = os_daily)


model.sel(m8_osdl_1, m8_osdl_2, m8_osdl_3, m8_osdl_4,
          m8_osdl_5, m8_osdl_6, m8_osdl_7, m8_osdl_8 )
# Best model (by far) is m8_osdl_6 (weight 0.85)
# Temprature and sex has no effect on daily movement
# Behavior has an effect on daily movement

## Split by sex for behavior
## Model for behavior by sex
osm <- subset(os_daily, os_daily$sex == "m")
osf <- subset(os_daily, os_daily$sex == "f")

## LMM for log transformed daily speed: males by behavior (temp has no influence and so dropped)
m9_osmdl <- lmer(log(daily_dist + 1) ~ behavior + (1|id) + (1|date),
                 weight = relocs, REML = TRUE, data = osm) # weights controls for n relocs
# Singular fit. Random effect structure probably too complex. Drop "date" as random factor.
m9_osmdl <- lmer(log(daily_dist + 1) ~ behavior + (1|id) + relocs,
                 REML = TRUE, data = osm) # weights controls for n relocs
# Diagnostics: all look good

# Summary for males
summary(m9_osmdl)
pairs(emmeans(m9_osmdl, ~behavior)) 
# No significant differences in males

## LMM for log transformed daily distance: females by behavior (temp has no influence and so dropped)
m10_osfdl<- lmer(log(daily_dist + 1) ~ behavior + (1|id) + (1|date) + relocs,
                REML = TRUE, data = osf) # weights controls for n relocs
# model did not converge. Random effect structure probably too complex. Drop "date" as random factor. 
m10_osfdl <- lmer(log(daily_dist + 1) ~ behavior + (1|id) + relocs,
                  REML = TRUE, data = osf) # weights controls for n relocs
## Diagnostics: some outliers and random effects residuals skewed to the left. 

# Summary for females
summary(m10_osfdl)
pairs(emmeans(m10_osfdl, ~behavior)) 
# No effect of mating behavior but effect of "pool_move"
# pool_move > others

