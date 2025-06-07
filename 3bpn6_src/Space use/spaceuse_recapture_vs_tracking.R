#################### Description ####################

# Author: Andrius Pašukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-03-24

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Sex difference in movement extent comparison for tracking and capture-recapture data 

###########

# libaries
library(ggplot2) 
library(dplyr)
library(lubridate)

# Import data
# Set your working directory 
setwd("~/")

frog_tracks <- read.table("all_spaceuse_raw.csv", header = T, sep =",")
dt_recap <- read.table("dtinc_recap_2009-2011.csv", header = T, sep =",")
af_recap <- read.table("afemo_recap_2014-2019.csv", header = T, sep =",")

# Clean up data
af_recap$date <- ymd(af_recap$date)
dt_recap$date <- ymd(dt_recap$date)

af_tracks <- frog_tracks %>%
  filter(species == "afemo") %>%
  mutate(dt = paste(date, time))
af_tracks$sex <- dplyr::recode(af_tracks$sex, "f" = "female", "m" = "male" )
af_tracks$dt <- dmy_hm(af_tracks$dt)

dt_tracks <- frog_tracks %>%
  filter(species == "dtinc") %>%
  mutate(dt = paste(date, time))
dt_tracks$sex <- dplyr::recode(dt_tracks$sex, "f" = "female", "m" = "male" )
dt_tracks$dt <- dmy_hm(dt_tracks$dt)


# Max linear distance per ID 

# A. femoralis recap
# summarize per ID
af_recap_sum <- af_recap %>%
  group_by(id)%>%
  add_count() %>%
  arrange(date)%>%
  summarise(sex = first(sex), period = last(date) - first(date), count = first(n)) 

# Maximum individual distances moved
splitted_id <- split(af_recap[,c("x","y")], af_recap$id)
distances <- lapply(splitted_id, dist)
max_dist_perid <-  lapply(distances, max)
af_recap_sum$max_dist_perid <- unsplit(max_dist_perid, af_recap_sum$id)

#Filter all ids with observation time span of less than 30 days
af_recap_sum <- filter(af_recap_sum, period > 30)

# D. tinctorius recap
# summarize per ID
dt_recap_sum <- dt_recap %>%
  group_by(id)%>%
  add_count() %>%
  arrange(date)%>%
  summarise(sex = first(sex), period = last(date) - first(date), count = first(n)) 

# Maximum individual distances moved
splitted_id <- split(dt_recap[,c("x","y")], dt_recap$id)
distances <- lapply(splitted_id, dist)
max_dist_perid <-  lapply(distances, max)
dt_recap_sum$max_dist_perid <- unsplit(max_dist_perid, dt_recap_sum$id)

#Filter all ids with observation time span of less than 30 days
dt_recap_sum <- filter(dt_recap_sum, period > 30)

# A. femoralis track
# summarize per ID
af_track_sum <- af_tracks %>%
  group_by(id)%>%
  add_count() %>%
  arrange(dt)%>%
  summarise(sex = first(sex), period = last(dt) - first(dt), count = first(n)) 

# Maximum individual distances moved
splitted_id <- split(af_tracks[,c("x_utm","y_utm")], af_tracks$id)
distances <- lapply(splitted_id, dist)
max_dist_perid <-  lapply(distances, max)
af_track_sum $max_dist_perid <- unsplit(max_dist_perid, af_track_sum$id)

# D. tinctorius track
# summarize per ID
dt_track_sum <- dt_tracks %>%
  group_by(id)%>%
  add_count() %>%
  arrange(dt)%>%
  summarise(sex = first(sex), period = last(dt) - first(dt), count = first(n)) 

# Maximum individual distances moved
splitted_id <- split(dt_tracks[,c("x_utm","y_utm")], dt_tracks$id)
distances <- lapply(splitted_id, dist)
max_dist_perid <-  lapply(distances, max)
dt_track_sum$max_dist_perid <- unsplit(max_dist_perid, dt_track_sum$id)

# Join data 
af_recap_sum$data <- "recapture"
af_track_sum$data <- "tracking"
dt_recap_sum$data <- "recapture"
dt_track_sum$data <- "tracking"

af_maxdist <- bind_rows(af_recap_sum, af_track_sum)
dt_maxdist <- bind_rows(dt_recap_sum, dt_track_sum)

## Boxplots grouped 

# boxplot of maximum dist per id
af_maxdist_plot <-
  ggplot(af_maxdist, aes(x= sex, y=max_dist_perid, fill = sex)) + 
  theme_bw(20) +
  geom_boxplot(aes(fill = sex), outlier.shape = NA) +
  labs(y = "Longest movement (m)") +
  geom_jitter(aes(fill= sex), position=position_jitterdodge(0.2), shape = 21, size = 4, alpha=0.4) +
  facet_wrap(~data, labeller = labeller(
    data=c(recapture ="Recapture", tracking = "Tracking"))) +
  scale_x_discrete(labels= c("F", "M")) +
  theme(legend.position= "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 18), aspect.ratio = 3) +
  ylim(0, 200)
af_maxdist_plot

dt_maxdist_plot <-
  ggplot(dt_maxdist, aes(x= sex, y=max_dist_perid, fill = sex)) + 
  theme_bw(20) +
  geom_boxplot(aes(fill = sex), outlier.shape = NA) +
  labs(y = "Longest movement (m)") +
  geom_jitter(aes(fill= sex), position=position_jitterdodge(0.2), shape = 21, size = 4, alpha=0.4) +
  facet_wrap(~data, labeller = labeller(
    data=c(recapture ="Recapture", tracking = "Tracking"))) +
  scale_x_discrete(labels= c("F", "M")) +
  theme(legend.position= "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 18), aspect.ratio = 3) #+
  #ylim(0, 200)
dt_maxdist_plot

## Statistics
# A. femo track
m19_aftr <- lm(log(max_dist_perid + 1) ~ sex + period, data = af_track_sum)
plot(m19_aftr) # look okish
summary(m19_aftr) # sex and period have effect on distance 

# A. femo recap
m20_afcr <- lm(log(max_dist_perid + 1) ~ sex + period, data = af_recap_sum)
plot(m20_afcr) # look good
summary(m20_afcr) # sex and period have effect on distance

# D. tinc track
m21_dttr <- lm(log(max_dist_perid + 1) ~ sex + period, data = dt_track_sum)
plot(m21_dttr) # look okish
summary(m21_dttr) # No effect of sex nor period on distance

# D. tinc recap
m22_dtcr <- lm(log(max_dist_perid + 1) ~ sex + period, data = dt_recap_sum)
plot(m22_dtcr) # look good
summary(m22_dtcr) # sex and period have effect on distance

