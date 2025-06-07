#################### Description ####################

# Author: Andrius Pašukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-03-24

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Space use comparison for O. sylvatica tracked in enclosures and at a natural site

###########
# libraries
library(sp)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)

# set working folder
setwd("~/")

# read in, display, and summarize 'csv' data
data_raw <- read.table("osylv_spaceuse_raw.csv", header = T, sep =",")
head(data_raw)

# Create timestamp 
timestamp<-paste(data_raw$date, data_raw$time)
tracks <- as.data.frame(data_raw)

# create a POSIXct object from timestamp and attach it do the dataframe
tracks$dt <- as.POSIXct(timestamp,  format="%d.%m.%y %H:%M")
head(tracks)

# Round up the coordinates to 0.1m 
tracks$x_utm <- round(tracks$x_utm, digits = 0.1)
tracks$y_utm <- round(tracks$y_utm, digits = 0.1)

# Add collumns with delta distance and delta time 
tracks <- tracks %>%
  group_by(id) %>%
  mutate(time_diff = difftime(dt, lag(dt, n = 1L), units = "min"), delta_x = x_utm - lag(x_utm, n = 1L), delta_y = y_utm - lag(y_utm, n = 1L), dist = sqrt(delta_x^2+delta_y^2))

## Tracking duration
# first pt of trajectory
tracks %>%
  group_by(id) %>%
  slice_head() -> start 
# last pt of trajectory
tracks %>%
  group_by(id) %>%
  slice_tail() -> end 
duration <- left_join(start, end, by = c("id", "sex", "species")) #join first and last pt of trajectory

# calculate durations and clean up dataframe
duration <- duration %>%
  mutate(duration = difftime(dt.y,dt.x, units = "days")) %>%
  dplyr::select(id, sex, species, start_time = dt.x, end_time = dt.y, duration)

## Daily movement 
# Daily cumulative distances
extras <-  dplyr::distinct(tracks, date, id, .keep_all=TRUE)
extras_perid <-  dplyr::distinct(tracks, id, .keep_all=TRUE)

# Tally the number relocations for each individual per day
tracks %>%
  group_by(date) %>%
  dplyr::count(id) -> relocs.perday

daily <- tracks %>%
  group_by(date,id) %>%
  filter(!dist == "NA") %>%
  dplyr::summarize(daily_dist = sum(dist), na.rm=TRUE)

daily <- full_join(daily, extras)
daily <- full_join(daily, relocs.perday)

# Maximum individual distances moved
splitted_id <- split(tracks[,c("x_utm","y_utm")], tracks$id)
distances <- lapply(splitted_id, dist)
max_dist_perid <-  lapply(distances, max)
daily$max_dist_perid <- unsplit(max_dist_perid, daily$id)

# Remove first row (tagging date) and last row (untagging day) per individual 
# Sort the daily, group by id and add index of id+date as a new variable
daily %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(id_day = paste(id, date)) -> daily_grouped

# select last row of each id group and make new variable
daily_grouped %>%  
  top_n(1, date) -> last_rows

# select first row of each id group and make new variable
daily_grouped %>%  
  top_n(-1, date) -> first_rows  

# filter daily_grouped by removing rows that match "id_day" of "last_rows" and "first_rows"
daily_select <- filter(daily_grouped, !(id_day %in% last_rows$id_day), !(id_day %in% first_rows$id_day))
daily_select <- dplyr::select(daily_select, id, date, sex, site, relocs = n, daily_dist, max_dist_perid)

# summarize by individual
daily_perid <- daily_select %>%
  group_by(id) %>%
  dplyr::summarize(mean_cumul = mean(daily_dist), max_dist_perid = first(max_dist_perid),
                   sex = first(sex), site = first(site)) %>%
  left_join(duration)
  
## Plots 
# boxplot of mean cumulative distance by sex
daily_plot <-
  ggplot(daily_perid, aes(x= sex, y=mean_cumul, fill = sex)) + 
  theme_bw(20) +
  geom_boxplot(aes(fill = sex), outlier.shape = NA) +
  labs(y = "Daily travel (m)") +
  geom_jitter(aes(fill=sex), position=position_jitterdodge(0.2), shape = 21, size = 4, alpha=0.4) +
  facet_wrap(~site, labeller = labeller(
    site=c(can ="Natural site", oto = "Enclosures"))) +
  scale_x_discrete(labels= c("F", "M")) +
  theme(legend.position= "none", axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        aspect.ratio = 3) +
  ylim(0, 25)
daily_plot

#boxplot of maximum dist per id
maxdist_plot <-
  ggplot(daily_perid, aes(x= sex, y=max_dist_perid, fill = sex)) + 
  theme_bw(20) +
  geom_boxplot(aes(fill = sex), outlier.shape = NA) +
  labs(y = "Longest movement (m)") +
  geom_jitter(aes(fill=sex), position=position_jitterdodge(0.2), shape = 21, size = 4, alpha=0.4) +
  facet_wrap(~site, labeller = labeller(
    site=c(can ="Natural site", oto = "Enclosures"))) +
  scale_x_discrete(labels= c("F", "M")) +
  theme(legend.position= "none", axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        aspect.ratio = 3) +
  ylim(0, 25) # excludes one male point at >40m from Canande
maxdist_plot

## Models
# Max dist
can <- filter(daily_perid, site == "can")
oto <- filter(daily_perid, site == "oto")

m23_oscan_max <- lm(log(max_dist_perid) ~ sex + duration, data = can)
plot(m23_oscan_max)
summary(m23_oscan_max)

m24_osoto_max <- lm(log(max_dist_perid) ~ sex + duration, data = oto)
plot(m24_osoto_max)
summary(m24_osoto_max)

# Daily dist
can_daily <- filter(daily_select, site == "can")
oto_daily <- filter(daily_select, site == "oto")

m25_oscan_dl <- lmer(log(daily_dist + 1) ~ sex + (1|id), data = can_daily)
plot(m25_oscan_dl)
summary(m25_oscan_dl)

m26_osoto_dl <- lmer(log(daily_dist + 1) ~ sex + (1|id), data = oto_daily)
plot(m26_oscan_dl)
summary(m26_osoto_dl)