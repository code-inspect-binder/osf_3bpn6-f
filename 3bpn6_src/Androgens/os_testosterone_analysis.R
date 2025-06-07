#################### Description ####################

# Author: Andrius Pašukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-05-11

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Plots and models for testosterone data

#################### Load packages ####################

library(dplyr) # for data managments
library(tidyr) # for pivoting data
library(MuMIn) # for model selection
library(betareg) # regression models for values between 0 - 1
library(lme4) #for running mixed models
library(lmerTest) #add p-values to mixed models
library(sjPlot) #for model plots
library(car) #for leveneTest, model Anova, and distribution plots
library(ggplot2) #plots

#################### Import and clean up data ####################

# set your working directory
setwd("")

# read in, display, and summarize 'csv' data
homing <- read.table(file="os_homing_export.csv", header = T, sep =",")
t_conc_corr <- read.table(file="os_tconc.csv", header = T, sep = ",")
size <- read.table(file="os_baseline_info.csv", header = T, sep = ",")

#Join navigation data 
hormones <- t_conc_corr %>%
  left_join(size, by = "id") %>%
  left_join(homing, by = "id") %>%
  dplyr::select(id, species = species.x, sex = sex.x, trans_group, condition, size, weight, homing_bin,
                track.time, track.daytime, area, sc, mean_temp, cumul_rain, t_conc_corr, cv)

# remove not translocated frogs
hormones <- hormones %>%
  filter(!is.na(track.time))%>%
  filter(cv < 15.5) # filter by CV 
hormones$homing_bin <- as.factor(hormones$homing_bin)

# Subset the data set to use base2 instead of base1 for the individuals
# that have the 2nd baseline
# Extract values of frogs that have two baselines
base2_frogs <-filter(hormones, condition == "base_02")
# Select base1 rows of frogs with two baselines
remove <- hormones %>%
  filter(id %in% base2_frogs$id)%>%
  filter(condition == "base_01")

# Subset hormone dataset by removing base1 of frogs that have base2
hormones <- anti_join(hormones, remove)
# Create grouped "baseline" variable and rename variables
hormones$condition <- recode_factor(hormones$condition, 
                                    base_01 = "baseline", base_02 = "baseline", base_03 = "baseline")
# Select only baseline data 
base <- hormones %>%
  filter(condition == "baseline")
# Select only homing frogs
homing_base <- filter(base, homing_bin == "1")

# Select delta change in testosterone
# pivot datasets for delta values
# calculate delta_t
hormones_wide <-hormones %>% 
  pivot_wider(id_cols = c("id", "sex"), names_from = condition, values_from = t_conc_corr)%>% 
  mutate(delta_t = back_home - baseline)

delta_t <- left_join(hormones, hormones_wide)
delta_t <- delta_t %>%
  filter(condition == "baseline")%>%
  filter(!is.na(delta_t))


#################### Plot sex and time point ####################
## Boxplot for sex differences
t_condition <- hormones %>%
  ggplot(aes(x = condition, y = log(t_conc_corr), fill = sex)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  #geom_jitter(aes(position=position_jitterdodge(0.1), size = 2, stroke = 1, shape = 21, alpha = 0.6)  +
  #geom_dotplot(position = position_dodge(0.4), binwidth = 0.06, binaxis='y', stackdir='center', dotsize = 2, stroke = 1, alpha = 0.6) +
  geom_point(position=position_jitterdodge(0.2), size = 2, stroke = 1, shape = 21, alpha = 0.6) +
  labs(x="Time point", y= expression(log~Androgen~(pg~mL^-1))) +
  scale_x_discrete(breaks=c("baseline","back_home"), labels=c("baseline", "back home"), expand = expansion(add = 0.5)) +
  theme_bw(24) +
  theme(aspect.ratio = 1.5, legend.position = "none") +
  ylim(1, NA)
t_condition
#ggsave(file = "os_t_sexcondition_boxplot_v3.pdf", width = 12, height = 12, units = "cm")

#################### Models: sex and time point ####

# data distribution
ggplot(hormones, aes(x=t_conc_corr, fill=sex)) + geom_density(alpha=.3) # data distribution by sex
ggplot(hormones, aes(x=log(t_conc_corr), fill=sex)) + geom_density(alpha=.3) # log transformed data distribution by sex

## LMER testosterone by sex and condition with interaction
m1_osT <- lmer(log(t_conc_corr) ~ sex*condition + (1|id),
               data = hormones)
# No interaction
m1_osT_1 <- lmer(log(t_conc_corr) ~ sex + condition + (1|id),
                 data = hormones)
# Diagnostics ok
anova(m1_osT, m1_osT_1, test = "LRT") # no interaction

summary(m1_osT_1) 
# Sex difference in T-level

# Model plot wit sjPlot
plot_model(m1_osT_1, title = "",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Time point\n (back home)", "Sex \n(male)"),
           dot.size = 5,
           line.size = 2,
           transform = NULL,
           sort.est = TRUE, colors = "gs", vline.color = "darkgrey")+
  ylim(-0.5, 1) +
  theme_sjplot2(24) 
#ggsave(file = "os_m1_osT_coefplot_v2.eps", width = 12, height = 12, units = "cm")

#################### Models: baseline T and movement ####
# Calculate standardized (centered and scaled) predictors
base$z.t_conc_corr <- stdize(log(base$t_conc_corr))
base$z.track.time <- stdize(base$track.time)
base$z.weight <- stdize(base$weight)
base$z.mean_temp <- stdize(base$mean_temp)

### Homing success ###
m2_osT <- glm(as.factor(homing_bin) ~ z.t_conc_corr + sex + trans_group +
                z.weight + z.mean_temp, family = binomial(link = "logit"), data = base)
# Drop weight?
drop1(m2_osT, test = "Chisq") 
# drop mean_temp but not weight

m2_osT <- glm(as.factor(homing_bin) ~ z.t_conc_corr + sex + trans_group +
                z.weight, family = binomial(link = "logit"), data = base)
summary(m2_osT)
# no effect of T

### Explored area ###

## LM for area explored with standardized continuous predictors 
m3_osT <- lm(log(area) ~ z.t_conc_corr + sex + trans_group + as.factor(homing_bin) + 
               z.weight + z.mean_temp,
             na.action = na.fail, data = base) 
# Diagnostics ok 

# Drop weight and mean_temp?
drop1(m3_osT, test = "F") 
# yes, both
## reduced model
m3_osT <- lm(log(area) ~ z.t_conc_corr + sex + trans_group + as.factor(homing_bin),
             na.action = na.fail, data = base) 
# Diagnostics ok
summary(m3_osT)
# no effect of T-conc  

# Model plot wit sjPlot
plot_model(m3_osT , title = "",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Homing\n(yes)", "Distance\n(50 m)",
                           "Sex \n(male)", "Andr\nlevel"),
           dot.size = 5,
           line.size = 2,
           transform = NULL,
           show.p = TRUE, colors = "bw", vline.color = "darkgrey", grid = FALSE) +
  theme_sjplot2(24) 
#ggsave(file = "os_m3_osT_coefplot_v2.eps", width = 12, height = 12, units = "cm")

### Homing duration ###

# Calculate standardized (centered and scaled) predictors
homing_base$z.t_conc_corr <- stdize(log(homing_base$t_conc_corr))
homing_base$z.track.time <- stdize(homing_base$track.time)
homing_base$z.weight <- stdize(homing_base$weight)
homing_base$z.mean_temp <- stdize(homing_base$mean_temp)

## LM for homing duration with standardized continuous predictors 
m4_osT<- lm(track.daytime ~ z.t_conc_corr + sex + z.weight +
               z.mean_temp, na.action = na.fail, data = homing_base)

# Drop weight and mean_temp?
drop1(m4_osT, test = "F") 
# drop both
## reduced model
m4_osT<- glm(track.daytime ~ z.t_conc_corr + sex,
             na.action = na.fail, family = Gamma(link = "inverse"), data = homing_base)
# Diagnostics ok
summary(m4_osT) 
# t_conc has no effect

### Trajectory straightness ###
## GLM for trajectory straigghtness (avPlots not available for betaregression)
m5_osT <- glm(sc ~ sex +  z.t_conc_corr +
                z.mean_temp + z.weight, family = Gamma(link = "identity"), data = homing_base)
# Drop weight and mean_temp?
drop1(m5_osT, test = "F") 
# drop both
## reduced model
m5_osT <- glm(sc ~ z.t_conc_corr + sex, family = Gamma(link = "identity"), data = homing_base)
# Diagnostics okish
summary(m5_osT) 
# t_conc has strong positive effect 

plot_model(m5_osT, title = "",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Sex\n(male) ", "Andr\nlevel"),
           dot.size = 5,
           line.size = 2,
           transform = NULL,
           colors = "bw", vline.color = "darkgrey") +
  ylim(-0.5, 0.5) +
  theme_sjplot2(24) 
#ggsave(file = "os_m5_osT_coefplot_v2.eps", width = 12, height = 12, units = "cm")


#pdf("os_m5_osT_regplot.pdf")
avPlots(m5_osT, ~z.t_conc_corr, id = FALSE, col = "black",
        col.lines = "red", lwd = 6, pch = 1, grid = FALSE,
        xlab = "Resid. Testosterone | others",
        ylab = "Resid. SC | others",
        cex=4, cex.axis = 2.5, cex.lab = 2.5)
#dev.off()

### model table for all dependent variables ####
tab_model(m3_osT, m4_osT, m5_osT,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("O. sylvatica Log explored area",
                        "O. sylvatica  homing duration",
                        "O. sylvatica  trajectory straightness"))

#################### Models: delta_T and movement ####

# Calculate standardized (centered and scaled) predictors
delta_t$z.area <- stdize(delta_t$area)
delta_t$z.track.time <- stdize(delta_t$track.time)
delta_t$z.sex <- stdize(delta_t$sex)
delta_t$z.homing_bin <- stdize(delta_t$homing_bin)
delta_t$z.weight <- stdize(delta_t$weight)
delta_t$z.mean_temp <- stdize(delta_t$mean_temp)

## LM for area explored with standardized continuous predictors 
m6_osT <- lm(delta_t ~  z.area + z.sex + z.track.time + z.homing_bin +
               z.weight + z.mean_temp,
             na.action = na.fail, data = delta_t)
# Drop weight and mean_temp?
drop1(m6_osT, test = "F") 
# yes, both
## reduced model
m6_osT <- lm(delta_t ~  z.area + sex + z.track.time + homing_bin,
             na.action = na.fail, data = delta_t)
# Diagnostics ok
# summary
summary(m6_osT)
# no significant effect

plot_model(m6_osT, title = "Delta androgen concentration",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Sex (male)", "Area explored", "Tracking duration", "Homing (yes)"),
           dot.size = 3, 
           show.values = TRUE, 
           sort.est = TRUE, colors = "bw", vline.color = "darkgrey") +
  #font_size(axis_title.x = 18, title = 18, axis_title.y = 18, labels.x = 18, labels.y = 18) +
  theme_sjplot2(16) 
#save_plot("os_m6_osT_coefplot.svg")

#### summary ####

## Summary notes for effects of testosterone
# Sex and Time point: strong sex difference and no time point effect
# Baseline: T has effect on SC
# Delta_t: no significant effects on delta_T