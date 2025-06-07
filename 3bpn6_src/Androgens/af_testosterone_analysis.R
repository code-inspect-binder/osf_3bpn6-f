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
homing <- read.table(file="af_homing_export.csv", header = T, sep =",")
t_conc_corr <- read.table(file="af_tconc.csv", header = T, sep = ",")
size <- read.table(file="af_baseline_info.csv", header = T, sep = ",")

# Replace one missing weight value with approximate value from size ~ weight regression line
plot(size$size, size$weight)
abline(lm(weight ~ size, data = size))
lm(weight ~ size, data = size )
# y = -1.83 + 0.14 * x
# for size af20_28 size = 26mm, estimated weight = -1.83 + 0.14 * 26 = 1.81 g
size$weight[67] <- 1.81

# Join data with covariates
hormones <- t_conc_corr %>%
  left_join(size, by = "id") %>%
  left_join(homing, by = "id") %>%
  dplyr::select(id, species = species.x, sex = sex.x, trans_group, condition, size, weight, homing_bin,
                track.time, track.daytime, area, sc, mean_temp, cumul_rain, t_conc_corr, cv)

# remove not translocated frogs
hormones <- hormones %>%
  filter(!is.na(track.time)) 
hormones$homing_bin <- as.factor(hormones$homing_bin)

# filter by CV 
hormones <- filter(hormones, cv < 15.5)

## Subset the data set to use base2 instead of base1 for the individuals
## that have the 2nd baseline
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

## Calculate delta change in testosterone
# pivot data sets for delta values
hormones_wide <-hormones %>% 
  pivot_wider(id_cols = c("id", "sex"), names_from = condition, values_from = t_conc_corr) %>% 
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
#ggsave(file = "af_t_sexcondition_boxplot_v3.pdf", width = 12, height = 12, units = "cm")

#################### Models: sex and time point ####

# data distribution
ggplot(hormones, aes(x=t_conc_corr, fill=sex)) + geom_density(alpha=.3) # data distribution by sex
ggplot(hormones, aes(x=log(t_conc_corr), fill=sex)) + geom_density(alpha=.3) # log transformed data distribution by sex

## GLMM for testosterone by sex and condition with interaction 
m1_afT <- lmer(log(t_conc_corr) ~ sex*condition + (1|id),
             data = hormones)
# Diagnostics ok
# No interaction
m1_afT_1 <- lmer(log(t_conc_corr) ~ sex + condition + (1|id),
             data = hormones)
# Diagnostics ok

anova(m1_afT, m1_afT_1, test = "LRT") # no interaction
summary(m1_afT_1) 
# Male have higher Testosterone
# No difference in time point

# Model plot wit sjPlot
plot_model(m1_afT_1, title = "",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Time point\n (back home)", "Sex \n(male)"),
           dot.size = 5,
           line.size = 2,
           show.p = TRUE, transform = NULL,
           sort.est = TRUE, colors = "gs", vline.color = "darkgrey") +
  theme_sjplot2(24) 
#ggsave(file = "af_m1_afT_1_coefplot_v2.eps", width = 12, height = 12, units = "cm")


#################### Models: baseline T and movement ####

# Calculate standardized (centered and scaled) predictors
base$z.t_conc_corr <- stdize(log(base$t_conc_corr))
base$z.track.time <- stdize(base$track.time)
base$z.weight <- stdize(base$weight)
base$z.mean_temp <- stdize(base$mean_temp)


### Homing success ###
m2_afT <- glm(homing_bin ~ z.t_conc_corr + sex + trans_group 
             + z.mean_temp + z.weight, family = binomial(link = "logit"), data = base)
# Drop weight and mean_temp?
drop1(m2_afT, test = "Chisq") 
# no, keep both 
summary(m2_afT)
#no effect of T

### Explored area ###

## LM for area explored with standardized continuous predictors 
m3_afT <- lm(log(area) ~ z.t_conc_corr + sex + trans_group + as.factor(homing_bin)
           + z.mean_temp + z.weight, data = base)

# Drop weight and mean_temp?
drop1(m3_afT, test = "F") 
# yes, better model without weight and mean_temp

## LM for area explored with standardized continuous predictors 
m3_afT <- lm(log(area) ~ z.t_conc_corr + sex + trans_group + as.factor(homing_bin),
             na.action = na.fail, data = base)
# Diagnostics ok
summary(m3_afT) 
# No influence of T-levels on exploration
plot_model(m3_afT)

# Model plot wit sjPlot
plot_model(m3_afT, title = "",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Homing\n(yes)", "Distance\n(50 m)",
                           "Sex \n(male)", "Andr\nlevel"),
           dot.size = 5,
           line.size = 2,
           transform = NULL,
           show.p = TRUE, colors = "bw", vline.color = "darkgrey", grid = FALSE) +
  theme_sjplot2(24) 
#ggsave(file = "af_m3_afT_coefplot_v2.eps", width = 12, height = 12, units = "cm")

### Homing duration ###

# Calculate standardized (centered and scaled) predictors
homing_base$z.t_conc_corr <- stdize(log(homing_base$t_conc_corr))
homing_base$z.track.time <- stdize(homing_base$track.time)
homing_base$z.weight <- stdize(homing_base$weight)
homing_base$z.mean_temp <- stdize(homing_base$mean_temp)

## LM for area explored with standardized continuous predictors 
m4_afT <- glm(track.daytime ~ z.t_conc_corr + trans_group + z.mean_temp + z.weight, family = Gamma(link = "inverse"), data = homing_base)

# Drop weight and mean_temp?
drop1(m4_afT, test = "F") 
# drop weight but not mean temp
## LM for area explored with standardized continuous predictors 
m4_afT <- glm(track.daytime ~ z.t_conc_corr + trans_group + z.mean_temp,
              na.action = na.fail, family = Gamma(link = "inverse"), data = homing_base)
summary(m4_afT) 
# No influence of T-cont

### Trajectory straightness ###

## GLM for trajectory straigghtness (avPlots not available for betaregression)
m5_afT <- glm(sc ~ z.t_conc_corr + trans_group + z.mean_temp + z.weight,
              family = Gamma(link = "identity"), data = homing_base)
# No sex as fixed factor, cause only 1 female with T samples homing
# Drop weight and mean_temp?
drop1(m5_afT, test = "F") 
# drop weight but not temp
# Final reduced model
m5_afT <- glm(sc ~ z.t_conc_corr + trans_group + z.mean_temp,
              na.action = na.fail, family = Gamma(link = "identity"), data = homing_base)
# Diangostics indicate outliers, remove pt 11 and 9
homing_base_outliers <-  homing_base[-c(9, 11), ]
# Model without outliers
m5_afT_out <- glm(sc ~ z.t_conc_corr + trans_group + z.mean_temp,
              na.action = na.fail, family = Gamma(link = "identity"), data = homing_base_outliers)
# Diagnostics not so much better

# Compare with and without outlier
summary(m5_afT)
summary(m5_afT_out) 
# No effect of T-conc on trajectory straightness either way

plot_model(m5_afT, title = "",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Temp.",  "Distance\n(50 m) ", "Andr\nlevel"),
           dot.size = 5,
           line.size = 2,
           transform = NULL,
           colors = "bw", vline.color = "darkgrey") +
  ylim(-0.25, 0.25) +
  theme_sjplot2(24) 
#ggsave(file = "af_m5_afT_coefplot_v2.eps", width = 12, height = 12, units = "cm")

#pdf("af_sct5_regplot.pdf")
avPlots(m4_afT, ~z.t_conc_corr, id = FALSE, col = "black",
        col.lines = "red", lwd = 6, pch = 1, grid = FALSE,
        xlab = "Resid. Testosterone | others",
        ylab = "Resid. SC | others",
        cex=4, cex.axis = 2.5, cex.lab = 2.5)
#dev.off()

### model table for all dependent variables ####
tab_model(m3_afT, m4_afT, m5_afT,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("A. femoralis Log explored area",
                        "A. femoralis homing duration",
                        "A. femoralis trajectory straightness"))


#################### Models: delta_T and movement ####
# Calculate standardized (centered and scaled) predictors
delta_t$z.area <- stdize(delta_t$area)
delta_t$z.track.time <- stdize(delta_t$track.time)
delta_t$z.weight <- stdize(delta_t$weight)
delta_t$z.mean_temp <- stdize(delta_t$mean_temp)

## LM for delta_T explored with standardized continuous predictors
m6_afT <- lm(delta_t ~  z.area + z.weight + sex + z.track.time+ z.mean_temp + homing_bin, data = delta_t)
# Diagnostics okish 
# Drop weight and mean_temp?
drop1(m6_afT, test = "F") 
# yes, better model without weight and mean_temp

## LM for delta_T explored reduced
m6_afT <- lm(delta_t ~  z.area +  sex + z.track.time+  homing_bin,
             na.action = na.fail, data = delta_t)

# Final model summary
summary(m6_afT) 
# Expored area has significant positive effect on delta_t
# homing success has significant negative effect on delta_t

plot_model(m6_afT, title = "Delta androgen concentration",
           axis.title = "Fixed effect estimates",
           axis.labels = c("Homing (yes)", "Tracking duration", "Sex (male)", "Area explored"),
           dot.size = 3, 
           show.values = TRUE, 
           ort.est = TRUE, colors = "bw", vline.color = "darkgrey") +
  #font_size(axis_title.x = 18, title = 18, axis_title.y = 18, labels.x = 18, labels.y = 18) +
  theme_sjplot2(16) 
#save_plot("af_m6_afT_coefplot.svg")

#### summary ####

## Summary notes for effects of testosterone
# Sex and Time point: sex difference in testosterone and no time point effect
# Baseline: no effect on T
# Delta_t: strong effect of exploration on delta_T, negative effect of homing success on delta_T