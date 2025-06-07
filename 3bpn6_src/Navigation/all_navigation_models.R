#################### Description ####################

# Author: Andrius Pasšukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-05-08

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Stats and plots for poison frog navigation

#################### Load packages ####################

library(ggplot2) #plots
library(dplyr) #deal with data easy
library(circular) #circular data and stats
library(CircStats) #circular data and stats
library(lme4) #GLMMs
library(emmeans) #pairwise comparisons
library(sjPlot) #for model plots
library(car) #for qqp and model Anova
library(coefplot) #for model coef plots
library(MuMIn) #for model selection using AICc
library(betareg) #models for 0-1 data 

#################### Import and clean data ####################

# Set your working directory
setwd("~/")

# import data
af_sum <- read.csv(file="af_homing_export.csv", header = T, sep =",") # navigation variables and covariates
af_size <- read.csv(file="af_nav_info.csv", header = T, sep =",") # size covariates
af_ang <- read.csv(file="af_angular.csv", header = T, sep =",") # orientation angles at difference thresholds

dt_sum <- read.csv(file="dt_homing_export.csv", header = T, sep =",")
dt_size <- read.csv(file="dt_nav_info.csv", header = T, sep =",")
dt_ang <- read.csv(file="dt_angular.csv", header = T, sep =",")

os_sum <- read.csv(file="os_homing_export.csv", header = T, sep =",")
os_size <- read.csv(file="os_nav_info.csv", header = T, sep =",")
os_ang <- read.csv(file="os_angular.csv", header = T, sep =",")

injuries <- read.csv(file="all_navigation_injuries.csv", header = T, sep =",") 
injuries$injuries_cat_2 <- as.character(injuries$injuries_cat_2)
injuries$injuries_cat_2 <- as.numeric(injuries$injuries_cat_2)

## Clean, split and merge data 

# Afemoralis
af_sum  <- af_sum %>%
  left_join(af_size, by = c("id", "sex")) %>%
  left_join(injuries)%>%
  dplyr::select(id, species, sex, size, weight, injuries_cat_2, trans_group, homing_bin, track.time,
                track.daytime, area, sc, mean_temp, cumul_rain)
  #filter(injuries_cat_2 < 2 | is.na(injuries_cat_2)) #filter out the frogs with injuries


# subset by translocation distance
af50 <- filter(af_sum, trans_group == "50m")
af200 <- filter(af_sum, trans_group == "200m")

# Reorder factors for "trans_group" for more intuitive plotting
af_sum$trans_group <- factor(af_sum$trans_group, levels= c("50m", "200m"))

# Dtinctorius
dt_sum  <- dt_sum %>%
  mutate(species = "dt")%>%
  left_join(dt_size)%>%
  left_join(injuries)%>%
  dplyr::select(id, species, sex, size, weight, injuries_cat_2, trans_group, homing_bin, track.time,
                track.daytime, area, sc, mean_temp, cumul_rain)
  #filter(injuries_cat_2  | is.na(injuries_cat_2)) #filter out the frogs with injuries

# subset by translocation distance
dt50 <- filter(dt_sum, trans_group == "50m")
dt200 <- filter(dt_sum, trans_group == "200m")

# Reorder factors for "trans_group" for more intuitive plotting
dt_sum$trans_group <- factor(dt_sum$trans_group, levels= c("50m", "200m"))

# Osylvatica
os_sum  <- os_sum %>%
  mutate(species = "os")%>%
  left_join(os_size)%>%
  left_join(injuries)%>%
  dplyr::select(id, species, sex, size, weight, injuries_cat_2, trans_group, homing_bin, track.time,
                track.daytime, area, sc, mean_temp, cumul_rain)%>%
  filter(injuries_cat_2 < 2 | is.na(injuries_cat_2)) #filter out the frogs with injuries

# subset for translocation distance
os50 <- filter(os_sum, trans_group == "50m")
os200 <- filter(os_sum, trans_group == "200m")

# Reorder factors for "trans_group" for more intuitive plotting
os_sum$trans_group <- factor(os_sum$trans_group, levels= c("50m", "200m"))

# merge all species
allfrogs <- bind_rows(af_sum, dt_sum, os_sum)

#subset only for homing frogs
af50_home <- filter(af50, homing_bin == "1")
af200_home <- filter(af200, homing_bin == "1")

dt50_home <- filter(dt50, homing_bin == "1")
dt200_home <- filter(dt200, homing_bin == "1")

os50_home <- filter(os50, homing_bin == "1")

afdt_200home <- bind_rows(af_sum, dt_sum)
afdt_200home <- filter(afdt_200home, homing_bin == "1" & trans_group == 
                         "200m")

#################### Circlular data subsets #####

## Afemo 50m
af50_10 <- af_ang%>%
  filter(trans_group == "50m")%>%
  filter(!(is.na(orient_angle1)))%>%
  dplyr::select(id, sex, homing_bin, lat = latency1, dist = orient_dist1,
                angle = orient_angle1)
af50_10$angle <- circular(af50_10$angle, units = "degrees")
af50_10$ang_rad <- rad(af50_10$angle)
af50_10$ang_rad <- circular(af50_10$ang_rad , units = "radians")

#subset for each sex
af50_males <- subset(af50_10, sex == "male")
af50_females <- subset(af50_10, sex == "female")

## Afemo 200m 
af200_40 <- af_ang%>%
  filter(trans_group == "200m")%>%
  filter(!(is.na(orient_angle2)))%>%
  dplyr::select(id, sex, homing_bin, lat = latency2, dist = orient_dist2,
                angle = orient_angle2)
af200_40$angle <- circular(af200_40$angle, units = "degrees")
af200_40$ang_rad <- rad(af200_40$angle)
af200_40$ang_rad <- circular(af200_40$ang_rad , units = "radians")

# subset for each sex
af200_males <- subset(af200_40, sex == "male")
af200_females <- subset(af200_40, sex == "female")

## Dtinc 50m
dt50_10 <- dt_ang%>%
  filter(trans_group == "50m")%>%
  filter(!(is.na(orient_angle1)))%>%
  dplyr::select(id, sex, homing_bin, lat = latency1, dist = orient_dist1,
                angle = orient_angle1)
dt50_10$angle <- circular(dt50_10$angle, units = "degrees")
dt50_10$ang_rad <- rad(dt50_10$angle)
dt50_10$ang_rad <- circular(dt50_10$ang_rad , units = "radians")

# subset for each sex
dt50_males <- subset(dt50_10, sex == "male")
dt50_females <- subset(dt50_10, sex == "female")

## Dtinc 200m
dt200_40 <- dt_ang%>%
  filter(trans_group == "200m")%>%
  filter(!(is.na(orient_angle2)))%>%
  dplyr::select(id, sex, homing_bin, lat = latency2, dist = orient_dist2,
                angle = orient_angle2)
dt200_40$angle <- circular(dt200_40$angle, units = "degrees")
dt200_40$ang_rad <- rad(dt200_40$angle)
dt200_40$ang_rad <- circular(dt200_40$ang_rad , units = "radians")

#subset for each sex
dt200_males <- subset(dt200_40, sex == "male")
dt200_females <- subset(dt200_40, sex == "female")

## Osylv 50m
os50_10 <- os_ang%>%
  filter(trans_group == "50m")%>%
  filter(!(is.na(orient_angle1)))%>%
  dplyr::select(id, sex, homing_bin, lat = latency1, dist = orient_dist1,
                angle = orient_angle1)
os50_10$angle <- circular(os50_10$angle, units = "degrees")
os50_10$ang_rad <- rad(os50_10$angle)
os50_10$ang_rad <- circular(os50_10$ang_rad , units = "radians")

#subset for each sex
os50_males <- subset(os50_10, sex == "male")
os50_females <- subset(os50_10, sex == "female")

## osinc 200m
os200_40 <- os_ang%>%
  filter(trans_group == "200m")%>%
  filter(!(is.na(orient_angle2)))%>%
  dplyr::select(id, sex, homing_bin, lat = latency2, dist = orient_dist2,
                angle = orient_angle2)
os200_40$angle <- circular(os200_40$angle, units = "degrees")
os200_40$ang_rad <- rad(os200_40$angle)
os200_40$ang_rad <- circular(os200_40$ang_rad , units = "radians")

#subset for each sex
os200_males <- subset(os200_40, sex == "male")
os200_females <- subset(os200_40, sex == "female")

#################### Navigation boxplots #################### 

## Explored area
# 50m trans
area50_log <-allfrogs%>%
  filter(trans_group == "50m")%>%
  ggplot(aes(x= sex, y=log(area), fill = species)) + 
  theme_bw(20) +
  geom_boxplot(aes(group = sex), width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex, fill = as.factor(homing_bin)),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, color = "black", size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("NA", "black", "#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = expression(Explored~area~(Log(m^2)))) +
  facet_wrap(~species, labeller = labeller(species=c("af" = "A. femoralis",
                                                     "dt" = "D. tinctorius",
                                                     "os" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(4,10)
area50_log
#ggsave(area50_log, file = "expl_50m_log_boxplot.pdf", device = pdf)

# 200m trans
area200_log <-allfrogs%>%
  filter(trans_group == "200m")%>%
  ggplot(aes(x= sex, y=log(area), fill = species)) + 
  theme_bw(20) +
  geom_boxplot(aes(color = sex), width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex, fill = as.factor(homing_bin)),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, color = "black", size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("NA", "black", "#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = expression(Explored~area~(Log(m^2)))) +
  facet_wrap(~species, labeller = labeller(species=c("af" = "A. femoralis",
                                                     "dt" = "D. tinctorius",
                                                     "os" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(4,10)
area200_log
#ggsave(area200_log, file = "expl_200m_log_boxplot.pdf", device = pdf)

## Trajectory straightness
# 50m trans
sc50_3 <-allfrogs%>%
  filter(trans_group == "50m" & homing_bin == 1)%>%
  ggplot(aes(x= sex, y=sc, fill = species)) + 
  theme_bw(20) +
  geom_boxplot(width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = "Trajectory straightness") +
  facet_wrap(~species, labeller = labeller(species=c("af" = "A. femoralis",
                                                     "dt" = "D. tinctorius",
                                                     "os" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(0, 1)
sc50_3
#ggsave(sc50_3, file = "sc_50m_boxplot.pdf", device = pdf)

# 200m trans
sc200_3 <-allfrogs%>%
  filter(trans_group == "200m" & homing_bin == 1)%>%
  ggplot(aes(x= sex, y=sc, fill = species)) + 
  theme_bw(20) +
  geom_boxplot(width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = "Trajectory straightness") +
  facet_wrap(~species, labeller = labeller(species=c("af" = "A. femoralis",
                                                     "dt" = "D. tinctorius",
                                                     "os" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(0, 1)
sc200_3
#ggsave(sc200_3, file = "sc_200m_boxplot.pdf", device = pdf)

## Homing duration
# 50m trans
duration50_3 <-allfrogs%>%
  filter(trans_group == "50m" & homing_bin == 1)%>%
  ggplot(aes(x= sex, y=track.daytime, fill = species)) + 
  theme_bw(20) +
  geom_boxplot(width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = "Homing duration (h)") +
  facet_wrap(~species, labeller = labeller(species=c("af" = "A. femoralis",
                                                     "dt" = "D. tinctorius",
                                                     "os" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(0, 75)
duration50_3
#ggsave(duration50_3, file = "duration_50m_boxplot.pdf", device = pdf)

# 200m trans
duration200_3 <-allfrogs%>%
  filter(trans_group == "200m" & homing_bin == 1)%>%
  ggplot(aes(x= sex, y=track.daytime, fill = species)) + 
  theme_bw(20) +
  geom_boxplot(width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = "Homing duration (h)") +
  facet_wrap(~species, labeller = labeller(species=c("af" = "A. femoralis",
                                                     "dt" = "D. tinctorius",
                                                     "os" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(0, 75)
duration200_3
#ggsave(duration200_3, file = "duration_200m_boxplot.pdf", device = pdf)

#################### Weather colinearity ####################

# Colinearity  of rain and temp
cor.test(allfrogs$cumul_rain, allfrogs$mean_temp, na.rm = TRUE) 
# r = -0.5. Temp and rain are collinear -> remove rain from models and use temp
# alsobecause rain has gaps
# Thomas (2015): Correlations between independent variables are generally regarded as excessive collinearity if r is substantially greater than 0.3

#################### Models: Homing success ####################
#################### A. femoralis ####################
## 50m trans
m1_afprob <- glm(homing_bin ~ sex +  mean_temp + weight, family = binomial(link = "logit"), data = af50)
summary(m1_afprob)
# Males more likely to home from 50m
# Negative influence of temperature
tab_model(m1_afprob, transform = NULL)

## 200m trans no females return, clearly a difference 

#################### D. tinctorius ####################
## 50m trans clearly no difference (1 male and 1 female did not return)

## 200m trans 
m2_dtprob <- glm(homing_bin ~ sex + mean_temp + weight , family = binomial(link = "logit"), data = dt200)
# model summary
summary(m2_dtprob)
# no difference
tab_model(m2_dtprob, transform = NULL)
#################### O. sylvatica ####################

## 50m trans
m3_osprob <- glm(homing_bin ~ sex + mean_temp, family = binomial(link = "logit"), data = os50)
# model does not converge when weight is included
# weight is evenly distributed among groups and clearly has no effect. 
# model summary
summary(m3_osprob)
# no difference

## 200m trans no frogs returned

#################### homing success model tables ####
tab_model(m1_afprob, m2_dtprob , m3_osprob,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("A. femoralis 50m homing success",
                        "D. tinctorius 200m homing success",
                        "O. sylvatica 50m homing success"))
#################### Models: Explored area ####################
# Because homing behavior is qualitatively very different between species and
# translocation distance, we will proceed separate models per species per translocation distance

#################### A. femoralis ####################

## 50m translocation
## LM for area explored with log transformation
m4_af50ex <- lm(log(area) ~ sex + mean_temp + weight, data = af50)
# Diagnostics look ok.
# Point 32 is an outlier. Check if the results change when removing 
af50_outliers <- af50[-c(32), ]
## LM for area explored with outliers removed
m4_af50ex_1 <- lm(log(area) ~ sex + mean_temp  + weight, data = af50_outliers)
summary(m4_af50ex)
summary(m4_af50ex_1) 
# removing "outliers" does not change the results
# males explore more than females
# higher temp reduced exploration
## Model plot
sjPlot::plot_model(m4_af50ex)

## 200m trans
## Lm with log transformation of dependent variable
m5_af200ex <-lm(log(area) ~ sex  + mean_temp + weight, data = af200)
# Diagnostics look good
# Model summary
summary(m5_af200ex)
# males explore more than females
# higher temp reduced exploration
## Model plot
sjPlot::plot_model(m5_af200ex)

#################### D. tinctorius ####################

## 50m trans
## Lm with log transformation of dependent variable
m6_dt50ex <-lm(log(area) ~ sex + mean_temp + weight, data = dt50)
# Diagnostics ok
# Model summary
summary(m6_dt50ex)
# no differences
# Model plot
sjPlot::plot_model(m6_dt50ex)

## 200m
## Lm with log transformation of dependent variable
m7_dt200ex <-lm(log(area) ~ sex  + mean_temp + weight, data = dt200)
# Diagnostics ok
# Model summary
summary(m7_dt200ex)
# Males explored more than females
# Frogs size had positive influence on explortaion
# no differences
# Model plot
sjPlot::plot_model(m7_dt200ex)

#################### O. sylvatica ####################

## 50m trans
## Lm with log transformation of dependent variable
m8_os50ex <-lm(log(area) ~ sex + mean_temp + weight, data = os50)
# Diagnostics ok
# Model summary
summary(m8_os50ex)
# no differences
# Model plot
sjPlot::plot_model(m8_os50ex)

## 200m trans
## Lm with log transformation of dependent variable
m9_os200ex <-lm(log(area) ~ sex + mean_temp + weight, data = os200)
# Diagnostics ok, but Point 16 close to outlier
# Try excluding the outlier Points 16. Check if the results change when removing these
os200_outliers <- os200[-c(13), ]
## Lm with log transformation of dependent variable and outliers removed
m9_os200ex_1 <-lm(log(area) ~ sex + mean_temp + weight, data = os200_outliers)
# Compare model results with and without outliers
summary(m9_os200ex)
summary(m9_os200ex_1) 
# removing "outlier" has strong effect on sex effect from P = 0.05 to P = 0.004
# This is the crazy female runner "Sylvia"
# males tend to me more explorative, no other influence
# Model plot
sjPlot::plot_model(m9_os200ex_1)

#################### exploration model tables ####
tab_model(m4_af50ex, m6_dt50ex, 
          m8_os50ex,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("Af 50m Log explored area",
                        "Dt 50m Log explored area",
                        "Os 50m Log explored area"))

tab_model(m5_af200ex, m7_dt200ex,
          m9_os200ex,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("Af 200m Log explored area",
                        "Dt 200m Log explored area",
                        "Os 200m Log explored area"))

#################### Models: homing trajectory straightness ####

#################### A. femoralis #####################

## 50 m trans
## beta regression for sc 
m10_af50sc <- betareg(sc ~ sex  + mean_temp,  data = af50_home)
# Diagnostics ok
# model summary
summary(m10_af50sc) 
# Males have more straight trajectories
# Temp has negative influence on trajectory straightness

## 200m trans: only males return

#################### D. tinctorius #####################

## 50m trans
## beta regression for sc 
m11_dt50sc <- betareg(sc ~ sex  + mean_temp,  data = dt50_home)
# Diagnostics ok
# model summary
summary(m11_dt50sc) 
# no differences whatsoever

## 200m  trans
## beta regression for sc 
m12_dt200sc <- betareg(sc ~ sex  + mean_temp,  data = dt200_home)
# Diagnostics: skewed residuals, One strong outlier point based on Cook's d.
# Idenfity and remove influencial points based on Cooks'D. Not necessary here, based on diagnostic plots
cooksD <- cooks.distance(m12_dt200sc)
n <- nrow(dt200_home)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 4/n, lty = 2, col = "steelblue")
# identify influenial points (e.g., cooksD > 6/n)
influential_obs <- as.numeric(names(cooksD)[(cooksD > (6/n))])
# Remove outliers and run the model again
dt200_home_outliers <- dt200_home[-influential_obs, ]

## beta regression for sc without one outlier point
m12_dt200sc_1 <- betareg(sc ~ sex  + mean_temp,  data = dt200_home_outliers)
# Diagnostics: residual plots look much better.
# model summaries with and without outlier
summary(m12_dt200sc)
summary(m12_dt200sc_1)
# In both cases males have less direct homing trajectories (explore more)

#################### O. sylvatica #####################

## 50m trans
## beta regression for sc 
m13_os50sc <- betareg(sc ~ sex  + mean_temp,  data = os50_home)
# Diagnostics ok
# model summary
summary(m13_os50sc)
# no sex difference
# temperature has negative influence on straightness

## 200m trans: no homing

#################### trajectory straightness model tables ####
tab_model(m10_af50sc, m11_dt50sc, m12_dt200sc, m13_os50sc,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("Af 50-m homing staightness",
                        "Dt 50-m homing staightness",
                        "Dt 200-m homing staightness",
                        "Os 50-m homing staightness"))

#################### Models: homing duration ####
#################### A. femoralis #####################

##  50m trans LM
## LM with log-tranfromed homing duration
m14_af50dr <- lm(log(track.daytime) ~ sex +  mean_temp + weight,
                 data =  af50_home)
# Diagnostics ok
# Model summary
summary(m14_af50dr)
sjPlot::plot_model(m14_af50dr)
# males home faster
# temp increases homing speed
# Back transform the coefficient from log scale into "percent" change
(exp(coef(m14_af50dr)[2]) - 1) * 100

## 200 m trans: only males returned

#################### D. tinctorius #####################

## 50 m trans LM
## LM with log-tranfromed homing duration
m15_dt50dr <- lm(log(track.daytime) ~ sex + mean_temp + weight,
                 data =  dt50_home)
# Diasgnostics marginal but better than LM 
# model summary
summary(m15_dt50dr)
sjPlot::plot_model(m15_dt50dr)
# No difference whatsovever 

## 200 m trans LM
m16_dt200dr <- lm(log(track.daytime) ~ sex + mean_temp + weight,
                  data =  dt200_home)
# Diagnostics not great but no better alternative
# Model summary
summary(m16_dt200dr)
sjPlot::plot_model(m16_dt200dr)
# No difference whatsoever

#################### O. sylvatica #####################

## 50 m trans LM 
m17_os50dr <- lm(log(track.daytime) ~ sex  + mean_temp + weight,
                 data = os50_home)
# Diagnostics ok
# Model summary
summary(m17_os50dr)
sjPlot::plot_model(m17_os50dr)
# No sex difference
# Higher temp reducedes homing speed

#################### homing duration model tables ####
tab_model(m14_af50dr, m15_dt50dr, m16_dt200dr, m17_os50dr,
          collapse.ci = TRUE, 
          transform = NULL,
          dv.labels = c("Af 50-m homing duration",
                        "Dt 50-m homing duration",
                        "Dt 200-m homing duration",
                        "Os 50-m homing duration"))

#################### Circular statistics for angular deviations ####

#################### A. femoralis ####

## 50 m translocation 
# angles measured at ~ 10m from start
# histogram of orientation angles
ggplot(af50_10, aes(x = angle)) +
  geom_histogram(aes(y = ..density..,), binwidth = 22.5, col = "black", fill = "grey") +
  geom_density(aes(col = sex), size = 1) +
  labs(x = "delta angle") +
  theme_bw(20)

#Significance of orientation towards home
circ.summary(af50_males$ang_rad)
circ.summary(af50_females$ang_rad)
home<- circular(0, units = "radians")
rayleigh.test(af50_males$ang_rad, home) # significant homeward orientation
rayleigh.test(af50_females$ang_rad, home) # no significant homeward orientation 

# Manova approach for comparing two circular distributions based on Landler et al. 2021
m18_af50cr<-manova(cbind(cos(af50_10$ang_rad), sin(af50_10$ang_rad)) ~ af50_10$sex)
summary(m18_af50cr)
# males are better orientated than females

## 200 m translocation 
# angles measured at ~ 40 m from start
# histogram of orientation angles
ggplot(af200_40, aes(x = angle)) +
  geom_histogram(aes(y = ..density..,), binwidth = 22.5, col = "black", fill = "grey") +
  geom_density(aes(col = sex), size = 1) +
  labs(x = "delta angle") +
  theme_bw(20)

# Significance of orientation towards home
# too few females to measure
circ.summary(af200_males$ang_rad)
home<- circular(0, units = "radians")
rayleigh.test(af200_males$ang_rad, home) #Signficicantly home oriented
# Too few females to do a group comparison

#################### D. tinctorius ####

## 50 m translocation 
# angles measured at ~ 10m from start
# histogram of orientation angles
ggplot(dt50_10, aes(x = angle)) +
  geom_histogram(aes(y = ..density..,), binwidth = 22.5, col = "black", fill = "grey") +
  geom_density(aes(col = sex), size = 1) +
  labs(x = "delta angle") +
  theme_bw(20)

# Significance of orientation towards home
circ.summary(dt50_males$ang_rad)
circ.summary(dt50_females$ang_rad)
home<- circular(0, units = "radians")
rayleigh.test(dt50_males$ang_rad, home) # significant homeward orientation
rayleigh.test(dt50_females$ang_rad, home) # significant homeward orientation

# Manova approach for comparing two circular distributions based on Landler et al. 2021
m19_dt50cr <-manova(cbind(cos(dt50_10$ang_rad),
                  sin(dt50_10$ang_rad)) ~ dt50_10$sex)
summary(m19_dt50cr)
# No sex differences

## 200 m translocation 
# angles measured at ~ 40 m from start
# histogram of orientation angles
ggplot(dt200_40, aes(x = angle)) +
  geom_histogram(aes(y = ..density..,), binwidth = 22.5, col = "black", fill = "grey") +
  geom_density(aes(col = sex), size = 1) +
  labs(x = "delta angle") +
  theme_bw(20)

#Significance of orientation towards home
circ.summary(dt200_males$ang_rad)
circ.summary(dt200_females$ang_rad)
home<- circular(0, units = "radians")
rayleigh.test(dt200_males$ang_rad, home) # no significant orientation 
rayleigh.test(dt200_males$ang_rad, home) # no significant orientation 

# Manova approach for comparing two circular distributions based on Landler et al. 2021
m20_dt200cr <- manova(cbind(cos(dt200_40$ang_rad),
                  sin(dt200_40$ang_rad)) ~ dt200_40$sex)
summary(m20_dt200cr)
# No sex difference

#################### O. sylvatica ######

## 50 m translocation 
# angles measured at ~ 10m from start
# histogram of orientation angles
ggplot(os50_10, aes(x = angle)) +
  geom_histogram(aes(y = ..density..,), binwidth = 22.5, col = "black", fill = "grey") +
  geom_density(aes(col = sex), size = 1) +
  labs(x = "delta angle") +
  theme_bw(20)

# Significance of orientation towards home
circ.summary(os50_males$ang_rad)
circ.summary(os50_females$ang_rad)
home<- circular(0, units = "radians")
rayleigh.test(os50_males$ang_rad, home) # significant howeward orientation 
rayleigh.test(os50_females$ang_rad, home) # no significant homeward orientation 

# Manova approach for comparing two circular distributions based on Landler et al. 2021
m21_os50cr <-manova(cbind(cos(os50_10$ang_rad),
                  sin(os50_10$ang_rad)) ~ os50_10$sex)
summary(m21_os50cr)
# No sex difference

## Osyl 200 m: no homing 
#Significance of orientation towards home
circ.summary(os200_males$ang_rad)
circ.summary(os200_females$ang_rad)
home<- circular(0, units = "radians")
rayleigh.test(os200_males$ang_rad, home) # no significant orientation 
rayleigh.test(os200_males$ang_rad, home) # no significant orientation 
