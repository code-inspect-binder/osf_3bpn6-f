#################### Description ####################

# Author: Andrius Pašukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-03-25

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Homing data processing and export for O. sylvatica

#### IMPORT DATA ####
library(ggplot2) #plots
library(ggforce) #plot circles for ggplot
library(dplyr) #deal with data easy
library(stringr) #extracting an replacing strings 
library(adehabitatLT) #split trajectory into bursts
library(lubridate) #deal with dates and times easy
library(circular) #circular data and stats
library(sf) #plot spatial 
library(viridis) #color palette 

#choose your working directory
setwd("~/")

trajectory.df <- read.csv(file="os_homing_traj.csv", header = T, sep =",")
territory.df <- read.csv(file="os_territorypts.csv", header = T, sep =",")

head(trajectory.df)
head(territory.df)

#changing the format of the date time to POSIXct format
trajectory.df$dt  <- ymd_hms(trajectory.df$dt)
territory.df$dt  <- ymd_hms(territory.df$dt)

#Round up the coordinates to 0.1m to include exclude the measurment error
trajectory.df$x_utm <- round(trajectory.df$x_utm, digits = 1)
trajectory.df$y_utm <- round(trajectory.df$y_utm, digits = 1)

#### FROGS TO EXCLUDE OR TRAJECTORIES TO TRIM ####

#os19_23f 
#Exluded in the final output but included for homing success

#Remove "os19_44f" because translocation stopped early due to injury
trajectory.df <- filter(trajectory.df, id != "os19_44f")
territory.df <- filter(territory.df, id != "os19_44f")

#check and remove duplicates
summary(duplicated(trajectory.df[1:2]))
which(duplicated(trajectory.df[1:2]))
trajectory.df <- distinct(trajectory.df, id, dt, .keep_all = TRUE)

#Check sample size
sample_size <- trajectory.df %>%
  group_by(id)%>%
  summarise(id = first(id))

#### DATA DOWNSAMPLING ####
#add collumns with delta distance and delta time for downsampling
trajectory.df <- trajectory.df %>%
  group_by(id) %>%
  mutate(time_diff = difftime(dt, lag(dt, n = 1L), units = "min"), delta_x = x_utm - lag(x_utm, n = 1L), delta_y = y_utm - lag(y_utm, n = 1L), dist = sqrt(delta_x^2+delta_y^2))

#Store last pt of trajectory to add back later
trajectory.df %>%
  group_by(id) %>%
  slice_tail() -> end

#Loop to downsample data based on minimum time between points
#Time parameter is set in "while(diff <= time in minutes)"
trajectory_dsmpl <- data.frame()
ids <- unique(trajectory.df$id)
for(i in ids) {
  traj = subset(trajectory.df, trajectory.df$id == i)
  for(i in 1:nrow(traj)) {
    diff <- difftime(traj$dt[i+1], traj$dt[i], units = "mins")
    if(is.na(diff)) {break}
    while(diff < 15) {
      traj <- traj[-(i+1),]
      diff <- difftime(traj$dt[i+1], traj$dt[i], units = "mins")
      if(is.na(diff)) {break}
    }
  }
  trajectory_dsmpl <- bind_rows(trajectory_dsmpl, traj)
}

trajectory_dsmpl <- trajectory_dsmpl %>%
  group_by(id) %>%
  mutate(time_diff = difftime(dt, lag(dt, n = 1L), units = "min"), delta_x = x_utm - lag(x_utm, n = 1L), delta_y = y_utm - lag(y_utm, n = 1L), dist = sqrt(delta_x^2+delta_y^2))

# #check which exceptions were made based on distance
# table(trajectory_dsmpl$id)
# check <- subset(trajectory_dsmpl, trajectory_dsmpl$time_diff <= 120)
trajectory.df <- trajectory_dsmpl %>%
  dplyr::select(id, dt, x_utm, y_utm, sex, time_diff, dist)


#add back last points to downsampled dataset in case they were removed because those are the "retun home" points
trajectory.df <- bind_rows(trajectory.df, end)
trajectory.df <- arrange(trajectory.df, id, dt)

#check and remove duplicates
summary(duplicated(trajectory.df[1:2]))
which(duplicated(trajectory.df[1:2]))
trajectory.df <- distinct(trajectory.df, id, dt, .keep_all = TRUE)

#Clean up df
trajectory.df <- trajectory.df%>%
  dplyr::select(id, dt, x_utm, y_utm, sex)

#### ADDING HYPOTHETICAL MORNING POINTS FOR EACH FROG EACH DAY ####

#Extract all evening points except for the last day of tracking (i.e. not a full day of tracking)
evening_pts<- trajectory.df %>%
  mutate(date = as.Date(dt), id_day = paste(id, date))%>%
  arrange(id, date)%>%
  group_by(id, date)%>%
  slice_tail()

#select last evening point in the trajectory for every frog to remove 
remove <- evening_pts%>%
  group_by(id)%>%
  slice_tail()

#remove last day pts from evening pts
evening_pts <- anti_join(evening_pts, remove)

#Frogs do not move at night
#Create hypothetical morning points for each date of tracking with coordinates from the last pt of the previous day
#This implies no movement between the last point of the previous day and the next day
morning_pts <- evening_pts%>%
  mutate(date = date + 1, time = "05:30:00", dt = paste(date, time))%>%
  dplyr::select(id, dt, x_utm, y_utm, sex)

#Merge morning points with the full trajectory df and sort it 
morning_pts$dt  <- ymd_hms(morning_pts$dt)
trajectory.df <- bind_rows(trajectory.df, morning_pts)
trajectory.df <- trajectory.df%>%
  arrange(id, dt)%>%
  dplyr::select(id, dt, x_utm, y_utm, sex)

#### CLEAN UP AND TERRITORY CENTERS ####

#Recaculate the dela time and dist because new points added to the df
trajectory.df <- trajectory.df %>%
  group_by(id) %>%
  mutate(time_diff = difftime(dt, lag(dt, n = 1L), units = "min"),
         delta_x = x_utm - lag(x_utm, n = 1L),
         delta_y = y_utm - lag(y_utm, n = 1L),
         dist = sqrt(delta_x^2+delta_y^2))%>%
  dplyr::select(id, sex, dt, x_utm, y_utm, time_diff, dist)

#Calculate territory center point (mean of territory points) for angle rotations 
territory_centers <- territory.df %>%
  group_by(id) %>%
  arrange(id, dt) %>%
  summarize(x_home = mean(x_utm), y_home = mean(y_utm))


#### SPATIOTEMPORAL VARIABLES AND HOMING SUCCESS ####

#Join territory center pts to the dataframe to calculate translocation distance and homing success
trajectory <- trajectory.df %>%
  left_join(territory_centers, by = "id")%>%
  mutate(trans_dist = sqrt((x_home - first(x_utm))^2+(y_home - first(y_utm))^2),
         trans_group = ifelse(trans_dist > 100, "200m", "50m"))%>%
  dplyr::select(id, sex, trans_group, trans_dist, dt,
                x_utm, y_utm, x_home, y_home, time_lag_min = time_diff, dist)

#Truncate trajectories to 4x24h for 50m and 6x24h for 200m translocations
trajectory <- trajectory%>%
  group_by(id)%>%
  mutate(duration = difftime(dt, first(dt), units = "hours"))%>% #time elsapsed in hours since release
  filter(if_else(trans_group == "50m", duration <= 96, duration <= 144))

#Calculate movement and temporal variables
trajectory  <- trajectory%>%
  group_by(id) %>%
  mutate(dist2home = sqrt((x_utm - x_home)^2+(y_utm - y_home)^2), #distance to home center
         dist2last= sqrt((x_utm - last(x_utm))^2+(y_utm - last(y_utm))^2), #distance to end of trajectory
         dist2release = sqrt((x_utm - first(x_utm))^2+(y_utm - first(y_utm))^2), #distance to release site
         homing_cont = last(dist2home), #use distance of last trajectory point to home as continious "homing" variable
         homing_prop = 1 - dist2home/trans_dist, #proportion of the translocation distance (i.e. dist from home center to release site) completed 
         homing_bin = if_else(last(homing_prop) > 0.6, 1, 0), #categorize as "homing frog" if last trajectory point if the frog completed at least 70% of homing 
         homing_cont = last(dist2home), #use distance of last trajectory point to home as continious "homing" variable
         delta_home = lag(dist2home,1) - dist2home, #delta distance moved in relation to home
         home_speed = delta_home/as.numeric(time_lag_min), #movement speed in relation to home (negative when moving away from home, positive when moving towards)
         time_norm = dt - as.POSIXct(date(first(dt)), tz="UTC"), #time elsapsed in hours since since the start of day the release day for plotting in relation to daytime
         nights = as.Date(last(dt)) - as.Date(first(dt))) #number of nights that frog was tracked

## Summarize homing variables per individual
homing_variables <- trajectory %>%
  group_by(id)%>%
  summarize(id = first(id), sex = first(sex), trans_group = first(trans_group),
            relocs = n(), #nomber of relocations per individual
            trans_dist = first(trans_dist), homing_prop = last(homing_prop), 
            homing_bin = last(homing_bin), homing_cont = last(homing_cont),
            dist_cumul = sum(dist, na.rm = TRUE), #cumulative distance of the trajectory
            delta_dist2home = first(dist2home) - last(dist2home), #straight line distance approached towards home between start and end of trajectory
            sc = delta_dist2home/dist_cumul, #straightness coefficient
            orientation_index = sc * abs(delta_dist2home/trans_dist), #Orientation index as SC relative to home corrected by proportion of the translocation distance moved
            track.time = last(duration), #total tracking/homing time
            nights = last(nights), 
            track.daytime = track.time - (nights * 12), #tracking/homing time in hours since release minus the night hours when frogs are immobile
            start_time = first(dt), end_time = last(dt))%>% 
  ungroup()%>%
  mutate(sampl_fq = track.daytime/relocs)%>% #sampling intensity
  dplyr::select(id, sex, trans_group, trans_dist, start_time, end_time,
                relocs, sampl_fq, nights, homing_cont, homing_prop, homing_bin, track.time,
                                  track.daytime, delta_dist2home, dist_cumul, sc, orientation_index)

#### DAILY MOVEMENT ####

#clean up the dataframe
trajectory.df <- trajectory.df %>%
  mutate(date = as.Date(dt), time = str_sub(dt, start = -8))%>%
  dplyr::select(id, sex, date, time, dt, x_utm, y_utm, dist, time_lag_min = time_diff)

#Tally the number relocations for each individual per day
trajectory.df %>%
  group_by(id, date) %>%
  count(id) -> relocs.perday

#Caculate daily cumulative distances and time lags
#Calculate speed
#Join with relocation info
#Add consecutive number per individual to indicate tracking day
daily <- trajectory.df %>%
  group_by(id, date) %>%
  filter(!dist == "NA") %>%
  dplyr::summarize(daily_dist = sum(dist), sex = first(sex), time_lag_min = sum(time_lag_min)) %>%
  mutate(id_day = paste(id, date), daily_speed = (daily_dist/as.numeric(time_lag_min) * 60)) %>%
  left_join(relocs.perday) %>%
  group_by(id) %>%
  mutate(trans_day = row_number())

#### TRAJECTORY ROTATION ####

#ggplot full trajectory
ggplot(trajectory.df, aes(x = x_utm, y = y_utm, color = id)) +
  geom_path(aes(group = id), arrow = arrow(angle = 30, length = unit(0.1, "inches"), ends = "last", type ="closed")) +
  theme_classic() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "none")

#normalize start position
trajectory.df <- trajectory.df %>%
  group_by(id) %>%
  mutate(x_norm = x_utm - x_utm[1], y_norm = y_utm - y_utm[1]) 

#Calculate territory center point (mean of territory points) for angle rotations 
territory_centers <- territory.df %>%
  group_by(id) %>%
  arrange(id, dt) %>%
  summarize(dt = last(dt), x_utm = mean(x_utm), y_utm = mean(y_utm))

### Translocation angle (from territory center to release site) and axis rotation ###
#selecting the earliest date time for the back-home point of each frog
release_pts <-
  trajectory.df %>%
  group_by(id) %>%
  arrange(dt) %>%
  slice(1) %>%
  dplyr::select(id,dt,x_utm,y_utm)

#combining territory center pts with release sites for angle calc
trans_vectors <- territory_centers %>%
  bind_rows(release_pts) %>%
  arrange(id, dt) %>%
  dplyr::select(id,dt,x_utm,y_utm)

#Calculate axis rotation angles directly from x and y with this convoluted ifelse (it works!)
trans_angles <- trans_vectors %>%
  group_by(id) %>%
  mutate(delta_x = x_utm - lag(x_utm, n = 1L), delta_y = y_utm - lag(y_utm, n = 1L), dist = sqrt(delta_x^2+delta_y^2),trans_angle = atan(abs(delta_x)/abs(delta_y)))%>%
  filter(!is.na(trans_angle)) %>%
  mutate(angle_corr = if_else(delta_x >= 0 & delta_y >= 0, trans_angle, if_else(delta_x > 0 & delta_y < 0, pi - trans_angle, if_else(delta_x <= 0 & delta_y <= 0, trans_angle + pi, 2*pi - trans_angle))))

##Rotate each trajectory by it's (reversed) translocation angle and calculate normalized xy coordinates
#Join translocation angles to the trajectory data frame for rotation 
trajectory.df <- left_join(trajectory.df, trans_angles, by = "id")

#Using directly calculated rotation angles with coef of pi 
trajectory.df$x_new <- trajectory.df$x_norm * cos(pi+trajectory.df$angle_corr) - trajectory.df$y_norm * sin(pi+trajectory.df$angle_corr)
trajectory.df$y_new <- trajectory.df$y_norm * cos(pi+trajectory.df$angle_corr) + trajectory.df$x_norm * sin(pi+trajectory.df$angle_corr)

#Clean up dataframe and create a grouping variable 50m or 200m based on translocation distance
trajectory.df <- trajectory.df %>%
  ungroup %>%
  mutate(trans_group = ifelse(dist.y > 100, "200m", "50m")) %>%
  dplyr::select(id, sex, trans_group, trans_dist = dist.y, dt = dt.x, x_utm = x_utm.x, y_utm = y_utm.x, x_new, y_new, dist = dist.x, time_lag_min) 

#Add rotated coordinates to trajectry df with all other variables
trajectory <- left_join(trajectory, trajectory.df)

#Split up the trajectories into 50m and 200m translocations
trajectory50m <- trajectory.df %>%
  filter(trans_group == "50m") 
trajectory200m <- trajectory.df %>%
  filter(trans_group == "200m") 

#### TRAJECTORY POLYGONS ####
#Convert to sf object
tracks_pts <- trajectory.df %>%
  ungroup()%>%
  na.omit()%>%
  dplyr::select(id, sex, trans_group, dt, x_utm, y_utm, x_new, y_new)%>%
  st_as_sf(coords = c("x_utm", "y_utm"), crs = 24817)

#Points to lines
tracks_lines <- tracks_pts %>%
  group_by(id) %>%
  arrange(dt) %>%
  summarize(id = first(id), sex = first(sex), do_union=FALSE)%>%
  st_cast("LINESTRING")

#Lines to buffered plygons. Change "dist" in meters for buffer size
plot(tracks_lines["id"])
tracks_polygons <- st_buffer(tracks_lines, dist = 5)
plot(tracks_polygons["id"])

#add collumn with buffer area
tracks_polygons$area <- st_area(tracks_polygons)

#Covert sf object back to regular data.frame
st_geometry(tracks_polygons) <- NULL

#Add area data to homing data
homing_variables <- left_join(homing_variables, tracks_polygons)
  
#### WEATHER CORRELATIONS ####

#Import weather data
temp <- read.csv(file="ecuador2019_daytime_temp.csv", header = T, sep =",")
rain <- read.csv(file="ecuador2019_rainfall_summaries.csv", header = T, sep =",")

rain$date <- as.Date(rain$date, format = "%Y-%m-%d")
temp$date <- as.Date(temp$date, format = "%Y-%m-%d")

#trans_group info
trans_group <- trajectory.df %>%
  group_by(id)%>%
  summarize(trans_group =first(trans_group))
            
daily <- left_join(daily, trans_group)
daily <- left_join(daily, temp, by = "date")
daily <- left_join(daily, rain, by = "date")

#Plot daily movement by daily temp & sex
ggplot(data = daily, aes(x = temp, y = daily_dist, group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()

#Plot daily movement by daytime rainfall & sex
ggplot(data = daily, aes(x = rain_daytime, y = daily_dist, group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw() + xlim(0,20)

#Plot daily movement by cumulative rainfall & sex
ggplot(data = daily, aes(x = rain_cumul, y = daily_dist, group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()

#### MOVEMENT LATENCIES AND ORIENTATION BY DISTANCE ####
#Angle from start#
#!Note that arguments in atan2(y,x) are flipped to calculate angle from y-axis. Default gives angle from x-axis.
trajectory <- trajectory %>%
  group_by(id)%>%
  mutate(angletohome = atan2(x_new,y_new) *180/pi)

#! Set distance thresholds from start in meters and the wanted orienation distance from release
#minimum distance from start
min_cut <- 5
#maximum distance from start
max_cut <- 15
#wanted orientation distance thresholds
od <- 10
#filter by distance from start
trajectory_cut <- trajectory %>%
  filter(dist2release >= min_cut & dist2release <= max_cut)%>%
  droplevels()

#Get rows closes to the specified orientation dist "od" from the range between min_cut and max_cut
orient1 <- trajectory_cut%>%
  group_by(id)%>%
  slice(which.min(abs(dist2release - od)))%>%
  dplyr::select(sex, trans_group, homing_bin,
                latency1 = duration, orient_dist1 = dist2release,
                orient_angle1 = angletohome)%>%
  mutate(latency1 = latency1, latency_daytime1 = if_else(latency1 < 12, latency1, latency1 - as.integer(latency1/24)*12))

#! Set distance thresholds from start in meters and the wanted orienation distance from release
#minimum distance from start
min_cut <- 20
#maximum distance from start
max_cut <- 60
#wanted orientation distance thresholds
od <- 40
#filter by distance from start
trajectory_cut <- trajectory %>%
  filter(dist2release >= min_cut & dist2release <= max_cut)%>%
  droplevels()

#Get rows closes to the specified orientation dist "od" from the range between min_cut and max_cut
orient2 <- trajectory_cut%>%
  group_by(id)%>%
  slice(which.min(abs(dist2release - od)))%>%
  dplyr::select(sex, trans_group, homing_bin,
                latency2 = duration, orient_dist2 = dist2release,
                orient_angle2 = angletohome)%>%
  mutate(latency2 = latency2, latency_daytime2 = if_else(latency2 < 22, latency2, latency2 - as.integer(latency2/24)*22))

#join both orientation frames
orient <- full_join(orient1, orient2)

#export angle data 
write.csv(orient, file="os_angular.csv")

#change factor order
orient$trans_group <- factor(orient$trans_group, levels = c("50m", "200m"))

#Boxplot latencies to leave the given radius
ggplot(orient, aes(x= sex, y=latency1, fill = sex)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitterdodge(0.2), size = 4, alpha=0.4) +
  theme_bw(20) +
  facet_wrap(~trans_group) +
  stat_summary(geom = "point", fun = "mean", col = "black", size = 3, shape = 24, fill = "red")

#histogram of orientation angles
ggplot(orient, aes(x = orient_angle1)) +
  geom_histogram(aes(y = ..density..), binwidth = 22.5, col = "black", fill = "grey") +
  geom_density(aes(col = sex), size = 1) +
  facet_wrap(~trans_group) +
  theme_bw(20)

#### TRAJECTORY PLOTS ####

### Plot start nomalized and rotated trajectories ###
trajectory.df$trans_group = factor(trajectory.df$trans_group, levels=c('50m','200m'))

### Plot start nomalized and rotated 200m trjectories with circles ###
ggplot() +
  geom_circle(aes(x0 = 0, y0 = 200, r = 15), linetype = "solid", color ="grey", fill = "grey", size = 0.25, alpha = 0.5) +
  geom_circle(aes(x0 = 0, y0 = 0, r = 50), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 100), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 150), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 200), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 250), linetype = "solid", size = 2, color = "#FC4E07") +
  geom_segment(aes(x = -250, xend = 250, y = 0, yend = 0), linetype = "solid", size = 0.25, color = "grey") +
  geom_segment(aes(x = 0, xend = 0, y = 250, yend = -250), linetype = "solid", size = 0.25, color = "grey") +
  geom_path(data = trajectory200m, aes(x = x_new, y = y_new, color = sex, group = id), alpha = 0.8,
            lwd = 1, arrow = arrow(angle = 30, length = unit(0.15, "inches"), ends = "last", type ="closed")) +
  geom_path(data = trajectory200m, aes(x = x_new, y = y_new, group = id), alpha = 0.8,
            lwd = 0.25, linetype = "dashed")  +
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1) +
  annotate("text", x=200, y=60, label= "200m", size = 5, angle = 10) +
  annotate("text", x=100, y=35, label= "100m", size = 5, angle = 10) +
  annotate("text", x=0, y=230, label= "Home", size = 5)

### Plot start nomalized and rotated 50m trjectories with circles ###
ggplot() +
  geom_circle(aes(x0 = 0, y0 = 45, r = 15), linetype = "solid", color ="grey", fill = "grey", size = 0.25, alpha = 0.5) +
  geom_circle(aes(x0 = 0, y0 = 0, r = 15), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 30), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 45), linetype = "solid", size = 0.25, color = "darkgrey") +
  geom_circle(aes(x0 = 0, y0 = 0, r = 60), linetype = "solid", size = 2, color = "#FC4E07") +
  geom_segment(aes(x = -60, xend = 60, y = 0, yend = 0), linetype = "solid", size = 0.25, color = "grey") +
  geom_segment(aes(x = 0, xend = 0, y = 60, yend = -60), linetype = "solid", size = 0.25, color = "grey") +
  geom_path(data = trajectory50m, aes(x = x_new, y = y_new, color = sex, group = id), alpha = 0.8,
            lwd = 1, arrow = arrow(angle = 30, length = unit(0.15, "inches"), ends = "last", type ="closed"))  +
  geom_path(data = trajectory50m, aes(x = x_new, y = y_new, group = id), alpha = 0.8,
            lwd = 0.25, linetype = "dashed")  +
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1) +
  annotate("text", x=44.5, y=14, label= "45m", size = 5, angle = 10) +
  annotate("text", x=30, y=11, label= "30m", size = 5, angle = 10) +
  annotate("text", x=0, y=55, label= "Home", size = 5)

#### TEMPORAL PLOTS ####

### Split up the trajectories into 50m and 200m translocations
trajectory50m <- trajectory %>%
  filter(trans_group == "50m")
trajectory200m <- trajectory %>%
  filter(trans_group == "200m")
trajectory$trans_group = factor(trajectory$trans_group, levels=c('50m','200m'))

##Temporal homing plot
dawn<- 6
dusk <- -6
xstart<-seq(dusk, length.out=10, by= 24)
xend<-seq(dawn, length.out=10, by= 24)
rects<-data.frame(xstart,xend)

#Temporal patterns 50m by homing
temp_plot50 <-
  ggplot(trajectory50m, aes(x=as.numeric(time_norm, units = "hours")/24, y=dist2home, group = id, color=homing_cont)) +
  geom_segment(aes(x = 0.5, xend = 4.5, y = 20, yend = 20), linetype = "dashed", size = 0.5, alpha = 0.7, color = "black") +
  geom_rect(data = rects, aes(xmin = xstart/24, xmax = xend/24,
                              ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes = FALSE) +
  geom_point(shape = 21, size = 3) + geom_line(size=1.5, alpha = 0.8) +
  xlab("Tracking day") + ylab("Distance to home center (m)") +
  theme(text = element_text(size = 20),
        legend.position="right",
        panel.background = element_rect(fill = 'white', colour = "black"),
        aspect.ratio = 1) +
  labs(colour = "End distance\nto home center") +
  scale_color_viridis(option = "E", limits = c(0, 75)) +
  xlim(0.5, 4.7) +
  ylim(0, NA) 
temp_plot50 


#Temporal patterns 200m by homing
temp_plot200 <-
  ggplot(trajectory200m, aes(x=as.numeric(time_norm, units = "hours")/24, y=dist2home, group = id, color=homing_cont)) +
  geom_segment(aes(x = 0.5, xend = 6.5, y = 20, yend = 20), linetype = "dashed", size = 0.5, alpha = 0.7, color = "black") +
  geom_rect(data = rects, aes(xmin = xstart/24, xmax = xend/24,
                              ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes = FALSE) +
  geom_point(shape = 21, size = 3) + geom_line(size=1.5, alpha = 0.8) +
  xlab("Tracking day") + ylab("Distance to home center (m)") +
  theme(text = element_text(size = 20),
        legend.position="right",
        panel.background = element_rect(fill = 'white', colour = "black"),
        aspect.ratio = 1) +
  labs(colour = "End distance\nto home center") +
  scale_color_viridis(option = "E", limits = c(0, 400)) +
  xlim(0.5, 6.5) +
  ylim(0, 400) 
temp_plot200 

#### HOMING SUCCESS BARPLOTS ####

#Homing success percent 50m males
percent_50_m <-homing_variables%>%
filter(trans_group == "50m" & sex == "male")%>%
summarize(percent = sum(homing_bin)/n() * 100)
percent_50_m <- percent_50_m$percent
  
#Homing success percent 50m females
percent_50_f <-homing_variables%>%
  filter(trans_group == "50m" & sex == "female")%>%
  summarize(percent = sum(homing_bin)/n() * 100)
percent_50_f <- percent_50_f$percent

#Homing success percent 200m males
percent_200_m <-homing_variables%>%
  filter(trans_group == "200m" & sex == "male")%>%
  summarize(percent = sum(homing_bin)/n() * 100)
percent_200_m <- percent_200_m$percent

#Homing success percent 200m females
percent_200_f <-homing_variables%>%
  filter(trans_group == "200m" & sex == "female")%>%
  summarize(percent = sum(homing_bin)/n() * 100)
percent_200_f <- percent_200_f$percent

#Barplot of homing success with percents as title
#50m
bar_50m <- homing_variables %>%
  filter(trans_group == "50m") %>%
  ggplot(aes(x = sex, fill = sex, alpha = as.factor(homing_bin))) +
  geom_bar(position = "fill", color = "black", size = 2) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  coord_flip() +
  theme(legend.position="none",
        aspect.ratio = 0.4,
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks.y = element_blank(),
        #axis.ticks.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand =expansion(add = c(0,0.1)), breaks = seq(0, 1, 0.5), 
                     labels = c("0", "50", "100 %")) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(title = paste("males =",percent_50_m,"%",
                     " females =",percent_50_f, "%"), y = "Percent homing")
bar_50m

#Barplot of homing success with percents as title
#200m
bar_200m <- homing_variables %>%
  filter(trans_group == "200m") %>%
  ggplot(aes(x = sex, fill = sex, alpha = as.factor(homing_bin))) +
  geom_bar(position = "fill", color = "black", size = 2) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  coord_flip() +
  theme(legend.position="none",
        aspect.ratio = 0.4,
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks.y = element_blank(),
        #axis.ticks.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand =expansion(add = c(0,0.1)), breaks = seq(0, 1, 0.5), 
                     labels = c("0", "50", "100 %")) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(title = paste("males =",percent_200_m,"%",
                     " females =",percent_200_f, "%"), y = "Percent homing")
bar_200m
  
#### CIRCULAR PLOTS ####
#Circular plot by sex
orient50_m <- filter(orient, trans_group == "50m" & sex == "male")
orient50_m <- filter(orient50_m, !(is.na(orient_angle1)))
corient50_m <- circular(orient50_m$orient_angle1, units = "degrees", zero = pi/2, rotation = "clock")
orient50_f <- filter(orient, trans_group == "50m" & sex == "female")
orient50_f <- filter(orient50_f, !(is.na(orient_angle1)))
corient50_f <- circular(orient50_f$orient_angle1, units = "degrees", zero = pi/2, rotation = "clock")

orient200_m <- filter(orient, trans_group == "200m" & sex == "male")
orient200_m <- filter(orient200_m, !(is.na(orient_angle2)))
corient200_m <- circular(orient200_m$orient_angle2, units = "degrees", zero = pi/2, rotation = "clock")
orient200_f <- filter(orient, trans_group == "200m" & sex == "female")
orient200_f <- filter(orient200_f, !(is.na(orient_angle2)))
corient200_f <- circular(orient200_f$orient_angle2, units = "degrees", zero = pi/2, rotation = "clock")

##Plot full circular info by sex
## start with blank plot, then add group-specific points
dev.new(); 
par(mai = c(1, 1, 0.1,0.1))
par(mar=c(0.5, 0.5, 0.5, 0.5))
pdf("os_circ_50m.pdf")
plot(corient50_m, bg = rgb(0, 0.749, 0.769), pch = 21, cex = 1.5, lwd = 2, stack = TRUE, bin = 60,
     xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), sep = 0.05, shrink = 1, tcl.text = -0.125,
     control.circle=circle.control(lwd = 2))
ticks.circular(circular(seq(0,2*pi,pi/2)), tcl=0.075)
par(new = T)
plot(corient50_f, bg=rgb(0.973, 0.463, 0.427), pch = 21, cex=1.5, lwd = 2, stack=T, bins=60, sep = -0.05, shrink= 1.3, axes = FALSE,
     control.circle=circle.control(lwd = 1))
arrows.circular(mean(corient50_m), y = rho.circular(corient50_m), col = rgb(0, 0.749, 0.769), lwd = 5)
arrows.circular(mean(corient50_f), y = rho.circular(corient50_f), col = rgb(0.973, 0.463, 0.427), lwd = 5)
par(new = T)
plot(corient50_m, col = NA, shrink= 2.5, axes = FALSE, control.circle=circle.control(lty = 2, lwd = 1))
ticks.circular(circular(seq(0,2*pi,pi/8)), tcl=0.2)
lines(density.circular(corient50_m, bw=30), shrink= 1, col = rgb(0, 0.749, 0.769, 0.7), lwd=2, lty=1)
lines(density.circular(corient50_f, bw=30), col = rgb(0.973, 0.463, 0.427, 0.7), lwd=2, lty=1)
#rose.diag(corient50, bins=16, cex=1.5, prop=2, add=TRUE, axes = FALSE)
dev.off()  

##Plot simplified circular info by sex
## start with blank plot, then add group-specific points
dev.new(); par(mai = c(1, 1, 0.1,0.1))
par(mar=c(0.5, 0.5, 0.5, 0.5))
pdf("os_circ_50m_simple.pdf")
plot(corient50_m, bg = rgb(0, 0.749, 0.769), pch = 21, cex = 3, lwd = 2, stack = TRUE, bin = 36,
     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), sep = 0.04, shrink = 1, tcl.text = -0.125,
     control.circle=circle.control(lwd = 2))
ticks.circular(circular(seq(0,2*pi,pi/2)), tcl=0.075)
par(new = T)
plot(corient50_f, bg=rgb(0.973, 0.463, 0.427), pch = 21, cex=3, lwd = 2, stack=T, bins=36, sep = -0.05, shrink= 1.7, axes = FALSE,
     control.circle=circle.control(lwd = 1))
arrows.circular(mean(corient50_m), y = rho.circular(corient50_m), col = rgb(0, 0.749, 0.769), lwd = 7)
arrows.circular(mean(corient50_f), y = rho.circular(corient50_f), col = rgb(0.973, 0.463, 0.427), lwd = 7)
dev.off()  

##Plot full circular info by sex
## start with blank plot, then add group-specific points
dev.new(); par(mai = c(1, 1, 0.1,0.1))
par(mar=c(0.5, 0.5, 0.5, 0.5))
pdf("os_circ_200m.pdf")
plot(corient200_m, bg = rgb(0, 0.749, 0.769), pch = 21, cex = 1.5, lwd = 2, stack = TRUE, bin = 60,
     xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), sep = 0.05, shrink = 1, tcl.text = -0.125,
     control.circle=circle.control(lwd = 2))
ticks.circular(circular(seq(0,2*pi,pi/2)), tcl=0.075)
par(new = T)
plot(corient200_f, bg=rgb(0.973, 0.463, 0.427), pch = 21, cex=1.5, lwd = 2, stack=T, bins=60, sep = -0.05, shrink= 1.3, axes = FALSE,
     control.circle=circle.control(lwd = 1))
arrows.circular(mean(corient200_m), y = rho.circular(corient200_m), col = rgb(0, 0.749, 0.769), lwd = 5)
arrows.circular(mean(corient200_f), y = rho.circular(corient200_f), col = rgb(0.973, 0.463, 0.427), lwd = 5)
par(new = T)
plot(corient200_m, col = NA, shrink= 2.5, axes = FALSE, control.circle=circle.control(lty = 2, lwd = 1))
ticks.circular(circular(seq(0,2*pi,pi/8)), tcl=0.2)
lines(density.circular(corient200_m, bw=30), shrink= 1, col = rgb(0, 0.749, 0.769, 0.7), lwd=2, lty=1)
lines(density.circular(corient200_f, bw=30), col = rgb(0.973, 0.463, 0.427, 0.7), lwd=2, lty=1)
#rose.diag(corient200, bins=16, cex=1.5, prop=2, add=TRUE, axes = FALSE)
dev.off()  

##Plot simplified circular info by sex
## start with blank plot, then add group-specific points
dev.new(); par(mai = c(1, 1, 0.1,0.1))
par(mar=c(0.5, 0.5, 0.5, 0.5))
pdf("os_circ_200m_simple.pdf")
plot(corient200_m, bg = rgb(0, 0.749, 0.769), pch = 21, cex = 3, lwd = 2, stack = TRUE, bin = 36,
     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), sep = 0.04, shrink = 1, tcl.text = -0.125,
     control.circle=circle.control(lwd = 2))
ticks.circular(circular(seq(0,2*pi,pi/2)), tcl=0.075)
par(new = T)
plot(corient200_f, bg=rgb(0.973, 0.463, 0.427), pch = 21, cex=3, lwd = 2, stack=T, bins=36, sep = -0.05, shrink= 1.7, axes = FALSE,
     control.circle=circle.control(lwd = 1))
arrows.circular(mean(corient200_m), y = rho.circular(corient200_m), col = rgb(0, 0.749, 0.769), lwd = 7)
arrows.circular(mean(corient200_f), y = rho.circular(corient200_f), col = rgb(0.973, 0.463, 0.427), lwd = 7)
dev.off()  

#### EXPORT DATA: FULL TRAJECTORIES WITH ATTRIBUTES ####
#Join with rotated trajectories for trajectory plotting
trajectory_export <- trajectory%>%
  ungroup()%>%
  dplyr::select(id, sex, trans_group, trans_dist, homing_prop, homing_bin,
                dt, x_utm, y_utm, x_new, y_new, dist2home, dist2release,
                duration, time_norm)%>%
  mutate(id_coded = str_sub(id, start = 3, end = -2))
write.csv(trajectory_export, file="os_homing_traj_export.csv")

#### EXPORT DATA: HOMING VARIABLES####
#Summarize daily movements and weather parameters for entire tracking period per id
perid <- daily %>%
  group_by(id) %>%
  arrange(id, trans_day)%>%
  summarize(day1_dist = first(daily_dist), day1_time = first(time_lag_min), day1_speed = first(daily_speed),
            day1_temp = first(temp), day1_rain = first(rain_24h),
            mean_temp = mean(temp), cumul_rain = sum(rain_24h), trans_day = first(trans_day))

#add weather variables to homing summary
homing_export <- left_join(homing_variables, perid)

#add orientation angles to homing summary
homing_export <- left_join(homing_export, orient)

# Export homing data for analysis
homing_export <- homing_export  %>%
  dplyr::select(id, sex, trans_group, relocs, homing_bin, homing_prop,
                track.time, track.daytime, dist_cumul, area, sc,
                mean_temp, cumul_rain)

#os19_23f has too few pts for SC.
#Replace respective cells with NA
homing_export[10,11] <- NA

write.csv(homing_export, file="os_homing_export.csv")
