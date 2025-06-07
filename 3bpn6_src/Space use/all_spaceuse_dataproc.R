#################### Description ####################

# Author: Andrius Pašukonis
# Email: andrius.pasukonis@cefe.cnrs.fr
# Script last edited 2022-05-11

# for Pašukonis et al. 2022: sex differences in poison frog spatial behavior
# Data processing and export for poison frog space use

#### IMPORT AND CLEAN UP DATA ####
library(dplyr)
library(maptools)
library(adehabitatHR)
library(ks)
library(ggplot2)
library(sp)

#choose your working directory
setwd("~/")

##read in, display, and summarize 'csv' data
#all_species_space_use contains all points (5448) manually combined from the following cleaned up GIS files
#Contains 3 behavioral categories: "w_op_sex", "w_same_sex" and "pool_move"
tracks <- read.table("all_spaceuse_raw.csv", header = T, sep =",")
head(tracks)

#check sample sizes: total relocs
tracks%>%filter(species == "afemo" & sex == "m")%>%nrow()
tracks%>%filter(species == "afemo" & sex == "f")%>%nrow()
tracks%>%filter(species == "osylv")%>%nrow() 
tracks%>%filter(species == "dtinc")%>%nrow() 
#check sample sizes: ids
tracks%>%filter(species == "afemo" & sex == "m")%>%count(id)
tracks%>%filter(species == "afemo" & sex == "f")%>%count(id)  
tracks%>%filter(species == "osylv")%>%count(id)  
tracks%>%filter(species == "dtinc")%>%count(id)  

#summarize the number of relocations by id
table(tracks$id)

#timestamp
tracks$dt <- as.POSIXct(paste(tracks$date,tracks$time),  format="%d.%m.%y %H:%M")
tracks$date <- as.Date(tracks$dt)
head(tracks)

#Round up the coordinates to 1m to reduce the measurment error
tracks$x_utm <- round(tracks$x_utm, digits = 1)
tracks$y_utm <- round(tracks$y_utm, digits = 1)

#### TRACKING DURATION, DELTA TIME, DELTA DISTANCE ####

#Add collumns with delta distance and delta time 
tracks <- tracks %>%
  group_by(id) %>%
  mutate(time_diff = difftime(dt, lag(dt, n = 1L), units = "min"), delta_x = x_utm - lag(x_utm, n = 1L), delta_y = y_utm - lag(y_utm, n = 1L), dist = sqrt(delta_x^2+delta_y^2))

##Tracking duration
#first pt of trajectory
tracks %>%
  group_by(id) %>%
  slice_head() -> start 
#last pt of trajectory
tracks %>%
  group_by(id) %>%
  slice_tail() -> end 
duration <- left_join(start, end, by = c("id", "sex", "species")) #join first and last pt of trajectory

#calculate durations and clean up dataframe
duration <- duration %>%
  mutate(duration = difftime(dt.y,dt.x, units = "days")) %>%
  dplyr::select(id, sex, species, start_time = dt.x, end_time = dt.y, duration)

#### DOWNSAMPLING ####

#Loop to downsample data based on minimum time and distance between points
#Time parameter is set in "while(diff <= time in minutes | dist > in meters)"
#downsampled pts that were less than 60min apart unless frog moved more than 20m during the time window
tracks_dsmpl <- data.frame()
ids <- unique(tracks$id)
for(i in ids) {
  traj = subset(tracks, tracks$id == i)
  for(i in 1:nrow(traj)) {
    diff <- difftime(traj$dt[i+1], traj$dt[i], units = "mins")
    delta_x <- traj$x_utm[i+1] - traj$x_utm[i]
    delta_y <- traj$y_utm[i+1] - traj$y_utm[i]
    dist <- sqrt(delta_x^2+delta_y^2)
    if(is.na(diff)) {break}
    while(diff <= 60 & dist < 20) {
      traj <- traj[-(i+1),]
      diff <- difftime(traj$dt[i+1], traj$dt[i], units = "mins")
      delta_x <- traj$x_utm[i+1] - traj$x_utm[i]
      delta_y <- traj$y_utm[i+1] - traj$y_utm[i]
      dist <- sqrt(delta_x^2+delta_y^2)
      if(is.na(diff)) {break}
    } 
  }
  tracks_dsmpl <- bind_rows(tracks_dsmpl, traj)
}

tracks_dsmpl <- tracks_dsmpl %>%
  group_by(id) %>%
  mutate(time_diff = difftime(dt, lag(dt, n = 1L), units = "min"), delta_x = x_utm - lag(x_utm, n = 1L), delta_y = y_utm - lag(y_utm, n = 1L), dist = sqrt(delta_x^2+delta_y^2))
tracks <- tracks_dsmpl

#check if and which exceptions were made based on distance threshold
check <- subset(tracks_dsmpl, tracks_dsmpl$time_diff <= 60)
check # 5 instances in Afemo, 3 instances in Dtinc

#check sample sizes: total relocs
tracks%>%filter(species == "afemo" & sex == "m")%>%nrow()
tracks%>%filter(species == "afemo" & sex == "f")%>%nrow()
tracks%>%filter(species == "osylv")%>%nrow() 
tracks%>%filter(species == "dtinc")%>%nrow() 

#### DAILY MOVEMENT ####

#Keep attributes for later joins
info <-  distinct(tracks, date, species, id, .keep_all=TRUE)
info_perid <-  distinct(tracks, id, .keep_all=TRUE)

#Tally the number relocations for each individual per day
tracks %>%
  group_by(date) %>%
  dplyr::count(id) -> relocs.perday

#Aggregate dela distance
daily <- tracks %>%
  group_by(date,id) %>%
  filter(!dist == "NA") %>%
  dplyr::summarize(daily_dist = sum(dist, na.rm=TRUE)) %>%
  arrange(id, date)

#Tally the number dates each frog was tracked
daily %>% ungroup() %>% count(id) %>% rename (days_tracked = n) -> days_tracked

#Join with extra attributes
daily <- dplyr::left_join(daily, info)
daily <- left_join(daily, relocs.perday)

##Daily maximum distances moved
#Filter tracks to remove dist == NA days (1st pt of traj and only 1 reloc per day ) for max dist calculations
no_nas <- filter(tracks, !dist == "NA")
splitted <- split(no_nas[,c("x_utm","y_utm")], paste(no_nas$id, no_nas$date))
distances <- lapply(splitted, dist)
max_dist <-  lapply(distances, max)
daily$max_dist <- unsplit(max_dist, paste(daily$id, daily$date))

#Maximum individual distances moved
splitted_id <- split(no_nas[,c("x_utm","y_utm")], no_nas$id)
distances <- lapply(splitted_id, dist)
max_dist_perid <-  lapply(distances, max)
daily$max_dist_perid <- unsplit(max_dist_perid, daily$id)

  ## Remove first row (tagging date) and last row (untagging day) per individual  ##

#Sort the daily, group by id and add index of id+date as a new variable
daily %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(id_day = paste(id, date)) -> daily_grouped

#select last row of each id group and make new variable
daily_grouped %>%  
  top_n(1, date) -> last_rows

#select first row of each id group and make new variable
daily_grouped %>%  
  top_n(-1, date) -> first_rows  

#filter daily_grouped by removing rows that match "id_day" of "last_rows" and "first_rows"
#filter for days that have less than 2 relocations
#clean up variables
daily_select <- daily_grouped %>%
  filter(!(id_day %in% last_rows$id_day), !(id_day %in% first_rows$id_day), n > 2)%>%
  dplyr::select(id, date, sex, species, relocs = n, daily_dist, max_dist, max_dist_perid)

#check sample sizes: ids
daily_select %>%filter(species == "afemo" & sex == "m")%>%count(id)
daily_select %>%filter(species == "afemo" & sex == "f")%>%count(id)  
daily_select %>%filter(species == "osylv")%>%count(id)  
daily_select %>%filter(species == "dtinc")%>%count(id)  

#summarize by individual
tracks_sum <- daily_select %>%
  dplyr::summarize(sex = first(sex), species = first(species),
                   mean.relocs = mean(relocs), mean_cumul = mean(daily_dist), cumul_sd = sd(daily_dist), 
                   mean_max = mean(max_dist), max_sd = sd(max_dist), max_dist_perid = max(max_dist_perid))
#relocs per frog
relocs <- tracks%>% count(id)%>% rename(relocs = n)

#Join with summary stats
tracks_sum  <- left_join(tracks_sum, duration)
tracks_sum  <- left_join(tracks_sum, relocs)
tracks_sum  <- left_join(tracks_sum, days_tracked)

#### BEHAVIORAL CATEGORIES ####

##Extracting distinct behavior classes for each date
behaviors <- distinct(tracks, date, species, id, beh_simple, .keep_all=TRUE)

beh_pools <- behaviors %>%
  group_by(date, id)%>%
  filter(beh_simple == "pool_move")

beh_w_same_sex <- behaviors %>%
  group_by(date, id) %>%
  filter(beh_simple == "w_same_sex")

beh_w_op_sex <- behaviors %>%
  group_by(date, id) %>%
  filter(beh_simple == "w_op_sex")

##Join daily movements with behavioral data in seperate dfs per behavior
#pool_move
daily_pools <- daily_select%>%
  left_join(beh_pools)%>%
  dplyr::select(id, date, sex,species, relocs, daily_dist, max_dist, behavior = "beh_simple") 

daily_pools$behavior <- as.character(daily_pools$behavior)
daily_pools[is.na(daily_pools$behavior)==TRUE,]$behavior <- "base" #recode NAs as "base" for dall days without "pool_move
daily_pools$behavior <- as.factor(daily_pools$behavior)

#w_op_sex  
daily_sex <- daily_select%>%
  left_join(beh_w_op_sex)%>%
  dplyr::select(id, date, sex, species, relocs, daily_dist, max_dist, behavior = "beh_simple")

daily_sex$behavior <- as.character(daily_sex$behavior)
daily_sex[is.na(daily_sex$behavior)==TRUE,]$behavior <- "base"
daily_sex$behavior <- as.factor(daily_sex$behavior)

# w_same_sex 
daily_same_sex <- daily_select%>%
  left_join(beh_w_same_sex)%>%
  dplyr::select(id, date, sex, species, relocs, daily_dist, max_dist, behavior = "beh_simple")

daily_same_sex$behavior <- as.character(daily_same_sex$behavior)
daily_same_sex[is.na(daily_same_sex$behavior)==TRUE,]$behavior <- "base"
daily_same_sex$behavior <- as.factor(daily_same_sex$behavior)

# Combing w_op_sex and pool_moves in one df  
daily_sex <- dplyr::rename(daily_sex, beh_sex = behavior) 
daily_pools <- dplyr::rename(daily_pools, beh_pool = behavior) 
daily_beh <- left_join(daily_sex, daily_pools)
daily_beh$behavior <- paste(daily_beh$beh_sex, daily_beh$beh_pool)

#Recode combined behavior class so that days with w_op_sex and pool_move only belong to pool_move
#!When both "w_op_sex" and "pool_move" occured, it is classified as "pool_move" because those movements are larger in extent (occured on 12 days)
daily_beh$behavior <- dplyr::recode(daily_beh$behavior, "base base" = "base", "w_op_sex pool_move" = "pool_move", "base pool_move" = "pool_move", "w_op_sex base" = "w_op_sex")

#Calculate average daily movement by behavioral category per individual for plotting 
daily_beh_mean <- daily_beh %>%
  group_by(id, behavior) %>%
  dplyr::summarize(daily_dist = mean(daily_dist), max_dist = mean(max_dist), sex = first(sex),
                   species = first(species))

#### WEATHER DATA ####
#Import all weather data
rain2016 <- read.csv("parare2016_rainfall_summaries.csv", header = T, sep =",")
rain2018 <- read.csv("parare2018_rainfall_summaries.csv", header = T, sep =",")
temp2016_ml <- read.csv("parare2016_daytime_temp_plateau.csv", header = T, sep =",")
temp2016_il <- read.csv("parare2016_daytime_temp_island.csv", header = T, sep =",")
temp_ecuador2017 <- read.csv("florida2017_daytime_temp.csv", header = T, sep =",")
temp2018_il <- read.csv("parare2018_daytime_temp_island.csv", header = T, sep =",")

##A.femoralis weather corelates
rain <- bind_rows(rain2016, rain2018)
temp <- bind_rows(temp2016_il, temp2018_il)
rain$date <- as.Date(rain$date, format = "%Y-%m-%d")
temp$date <- as.Date(temp$date, format = "%Y-%m-%d")
af <- daily_beh %>%
  filter(species == "afemo")%>%
  left_join(temp, by = "date")%>%
  left_join(rain, by = "date")

#Corelation of rain and temp
ggplot(data = af, aes(x = rain_daytime, y = temp)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(af$rain_daytime, af$temp) #cor -0.57

#Plot by temp
ggplot(data = af, aes(x = temp, y = log(daily_dist), group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(log(af$daily_dist), af$temp) #cor -0.4

#Plot by daytime rain (some dates missing rain values)
ggplot(data = af, aes(x = rain_daytime, y = log(daily_dist), group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(log(af$daily_dist), af$rain_daytime) #cor 0.25

#Plot by cumul rain (some dates missing rain values)
ggplot(data = af, aes(x = rain_cumul, y = log(daily_dist), group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(log(af$daily_dist), af$rain_cumul) #cor 0.33
##Temperature is a better and more consistent predictor than rain

## D.tinctorius weather corelates
rain2016$date <- as.Date(rain2016$date, format = "%Y-%m-%d")
temp2016_ml$date <- as.Date(temp2016_ml$date, format = "%Y-%m-%d")
dt <- daily_beh %>%
  filter(species == "dtinc")%>%
  left_join(temp2016_ml, by = "date")%>%
  left_join(rain2016, by = "date")

#Corelation of rain and temp
ggplot(data = dt, aes(x = rain_daytime, y = temp)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(dt$rain_daytime, dt$temp) #cor -0.37

#Plot by temp
ggplot(data = dt, aes(x = temp, y = log(daily_dist), group = sex, color = sex)) +
  geom_point() + geom_smooth(method = "loess") + theme_bw()
cor.test(log(dt$daily_dist), dt$temp) #no correlation

#Plot by daytime rain (some dates missing rain values)
ggplot(data = dt, aes(x = rain_daytime, y = log(daily_dist), group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(log(dt$daily_dist), dt$rain_daytime) #no correlation

#Plot by cumul rain (some dates missing rain values)
ggplot(data = dt, aes(x = rain_cumul, y = log(daily_dist), group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(log(dt$daily_dist), dt$rain_cumul) #cor 0.15
##Only cumul rain slightly correlates with daily_dist

##O. sylvatica weather corelates (no rain data for os17)
temp_ecuador2017$date <- as.Date(temp_ecuador2017$date, format = "%Y-%m-%d")
os <- daily_beh %>%
  filter(species == "osylv")%>%
  left_join(temp_ecuador2017, by = "date")

#Plot by temp
ggplot(data = os, aes(x = temp, y = daily_dist, group = sex, color = sex)) +
  geom_point() + stat_smooth(method=lm) + theme_bw()
cor.test(log(os$daily_dist), os$temp) #no correlation

##Merge all species with weather info and clean up
daily_full <- bind_rows(af,dt,os)
daily_full <- dplyr::select(daily_full, id, species, sex, date, relocs, behavior, 
                            daily_dist, temp)
#Plot by temp
ggplot(data = daily_full, aes(x = temp, y = log(daily_dist), group = species, color = species)) + 
  geom_point() + stat_smooth(method=lm) + theme_bw() 

#### HOME RANGE adehabitatHR ####
##Filter data for frogs tracked > 7 days with at least 30 relocations
exclude <- tracks_sum%>%
  filter(duration < 7 & relocs < 30)
tracks_hr <- tracks%>%
  filter(!id %in% exclude$id)%>%
  droplevels()
#consider removing Jane (t06) as an outlier from the data (because she showed no site fidelity)
#consider removing Tarzan (t08) as an outlier from the data (because she showed no site fidelity)

#check sample sizes: ids
tracks_hr%>%filter(species == "afemo" & sex == "m")%>%count(id)
tracks_hr%>%filter(species == "afemo" & sex == "f")%>%count(id)  
tracks_hr%>%filter(species == "osylv")%>%count(id)  
tracks_hr%>%filter(species == "dtinc")%>%count(id)  

#Covert to spatial data frame 
data_sp<-tracks_hr
coordinates(data_sp) <- c("x_utm", "y_utm")
class(data_sp)

##MCPs
#home range with 95% minimum convex ploygon (MCP95) 
mcp95<- mcp(data_sp["id"], percent = 95, unout = "m2")
#100% MCP for all points
mcp100 <- mcp(data_sp["id"], percent = 100, unout = "m2")
#50% MCP for core area
mcp50<- mcp(data_sp["id"], percent = 50, unout = "m2")

#combine MCP data
mcp100<- as.data.frame(mcp100)%>%rename(mcp100 = area)
mcp95<- as.data.frame(mcp95)%>%rename(mcp95 = area)
mcp50<- as.data.frame(mcp50)%>%rename(mcp50 = area)
mcp_data <- mcp100%>%full_join(mcp95)%>%full_join(mcp50)

##KUD: home range with utilization density and 'href"
kde <- kernelUD(data_sp["id"], h = "href", grid = 200, extent = 1.5, same4all = FALSE)

#get 95% and 50% of utilization density 
khr95 <- getverticeshr(kde, percent = 95, unout = "m2")
khr50 <- getverticeshr(kde, percent = 50, unout = "m2")

khr95<- as.data.frame(khr95)%>%rename(ud95_href = area)
khr50<- as.data.frame(khr50)%>%rename(ud50_href = area)
kud_data <- full_join(khr95, khr50, by = "id")

##Full data summary##
tracks_sum <- left_join(tracks_sum, mcp_data)
tracks_sum <- left_join(tracks_sum, kud_data)
#write.csv(tracks_sum, file="all_frogs_move_summary.csv")

#### HOME RANGE ks ####
locs_hpi <- as.data.frame(dplyr::select(tracks_hr, id, x_utm, y_utm))

### 95% contour with hpi plug-ing method
#Define countour level
cont=c(95)
#Initiate empty SpatialPolygons to store all ids in the loop
allfrogs_polys <- SpatialPolygons(list())
#Initiate empty df to store the metadarta
allfrogs_info <- data.frame(ID = numeric(0), Cont = numeric(0), frog_id = character(0))
#the following loop allows to do the calculation for all individual at once
frogs <- sort(unique(locs_hpi$id))
area.frame <- data.frame(name=frogs)
cnt <- paste0('cont_',cont)
area.frame[,cnt] <- NA
spp_all<-list()

for (kk in 1:length(frogs)) {
  tt <- locs_hpi[locs_hpi$id == frogs[kk],]
  m_hd<-matrix(c(tt$x, tt$y), ncol=2)
  
  hpi_hd<-Hpi(x = m_hd) #instead of the plug-in method (Hpi) others can be used here
  
  kde_hd<-kde(x=m_hd, H=hpi_hd)
  
  hts <- contourLevels(kde_hd, prob=(100-cont)/100, approx=TRUE)
  
  col <- grey(seq(0.9,0,length=length(hts)+1))
  if (length(col)<100) col <- rep(col, length=100)
  
  plot(kde_hd$x[,1], kde_hd$x[,2], type="n",xlab=kde_hd$names[1], ylab=kde_hd$names[2]) 
  areas <- vector()
  polys <- data.frame(ID= numeric(0), cont= integer(0))
  tt <- 1
  
  ## calculate contours
  for (i in 1:length(hts))
  {
    scale <- 1
    
    if (hts[i]>0 | !is.null(kde_hd$deriv.order))
      contour(kde_hd$eval.points[[1]], kde_hd$eval.points[[2]], kde_hd$estimate*scale,
              level=hts[i]*scale, add=TRUE, drawlabels=TRUE, labcex=1, col=col[i], lwd=1)
    clines <- contourLines(kde_hd$eval.points[[1]], kde_hd$eval.points[[2]], kde_hd$estimate*scale,
                           level=hts[i]*scale)
    
    A <- 0
    for (jj in 1:length(clines)){
      x <- clines[[jj]][["x"]]
      y <- clines[[jj]][["y"]]
      A = A + 0.5* abs( sum( x[1:(length(x)-1)]*y[2:length(x)] - y[1:(length(x)-1)]*x[2:length(x)] ) )
      
      pts <- cbind(clines[[jj]]$x,clines[[jj]]$y)
      
      # you need to make sure that the last point and first coincide
      pts = rbind(pts, pts[1,])
      
      sp = SpatialPolygons( list(  Polygons(list(Polygon(pts)), 1)))
      
      sp <- spChFIDs(sp, as.character(tt))
      
      if (tt==1) {
        sp.comb <- sp
      } else {
        sp.comb <- spRbind(sp.comb,sp)
      }
      polys <- rbind(polys,c(tt,cont[i]))
      #    polys[tt] <- cont[i] 
      tt <- tt + 1
    }
    # data frame: rows = names; columns = cont
    areas[i] <- A
    area.frame[kk,i+1] <- areas[i]
  }
  
  
  names(polys) <- c("ID","Cont")
  polys$frog_id <- frogs[kk]
  spp<- SpatialPolygonsDataFrame(sp.comb,data=polys)
  
  allfrogs_polys <- raster::bind(allfrogs_polys,sp.comb)
  allfrogs_info <- rbind(allfrogs_info, polys)
  #Output as shape file with the individual's names and e.g. the corresponding contour level can be changed
  #writeOGR(spp,dsn="./shp_allfrogs_50HR",sprintf("Hpi_50%s",frogs[kk]), driver="ESRI Shapefile",overwrite_layer=T) 
  
  spp_all[[kk]]<-spp
  
} # end kk loop
#Output is area.frame
kud95_hpi <- area.frame%>%dplyr::rename(id = name) 

### 50% contour with hpi plug-ing method
#Define countour level
cont=c(50)
#Initiate empty SpatialPolygons to store all ids in the loop
allfrogs_polys <- SpatialPolygons(list())
#Initiate empty df to store the metadarta
allfrogs_info <- data.frame(ID = numeric(0), Cont = numeric(0), frog_id = character(0))
#the following loop allows to do the calculation for all individual at once
frogs <- sort(unique(locs_hpi$id))
area.frame <- data.frame(name=frogs)
cnt <- paste0('cont_',cont)
area.frame[,cnt] <- NA
spp_all<-list()

for (kk in 1:length(frogs)) {
  tt <- locs_hpi[locs_hpi$id == frogs[kk],]
  m_hd<-matrix(c(tt$x, tt$y), ncol=2)
  
  hpi_hd<-Hpi(x = m_hd) #instead of the plug-in method (Hpi) others can be used here
  
  kde_hd<-kde(x=m_hd, H=hpi_hd)
  
  hts <- contourLevels(kde_hd, prob=(100-cont)/100, approx=TRUE)
  
  col <- grey(seq(0.9,0,length=length(hts)+1))
  if (length(col)<100) col <- rep(col, length=100)
  
  plot(kde_hd$x[,1], kde_hd$x[,2], type="n",xlab=kde_hd$names[1], ylab=kde_hd$names[2]) 
  areas <- vector()
  polys <- data.frame(ID= numeric(0), cont= integer(0))
  tt <- 1
  
  ## calculate contours
  for (i in 1:length(hts))
  {
    scale <- 1
    
    if (hts[i]>0 | !is.null(kde_hd$deriv.order))
      contour(kde_hd$eval.points[[1]], kde_hd$eval.points[[2]], kde_hd$estimate*scale,
              level=hts[i]*scale, add=TRUE, drawlabels=TRUE, labcex=1, col=col[i], lwd=1)
    clines <- contourLines(kde_hd$eval.points[[1]], kde_hd$eval.points[[2]], kde_hd$estimate*scale,
                           level=hts[i]*scale)
    
    A <- 0
    for (jj in 1:length(clines)){
      x <- clines[[jj]][["x"]]
      y <- clines[[jj]][["y"]]
      A = A + 0.5* abs( sum( x[1:(length(x)-1)]*y[2:length(x)] - y[1:(length(x)-1)]*x[2:length(x)] ) )
      
      pts <- cbind(clines[[jj]]$x,clines[[jj]]$y)
      
      # you need to make sure that the last point and first coincide
      pts = rbind(pts, pts[1,])
      
      sp = SpatialPolygons( list(  Polygons(list(Polygon(pts)), 1)))
      
      sp <- spChFIDs(sp, as.character(tt))
      
      if (tt==1) {
        sp.comb <- sp
      } else {
        sp.comb <- spRbind(sp.comb,sp)
      }
      polys <- rbind(polys,c(tt,cont[i]))
      #    polys[tt] <- cont[i] 
      tt <- tt + 1
    }
    # data frame: rows = names; columns = cont
    areas[i] <- A
    area.frame[kk,i+1] <- areas[i]
  }
  
  
  names(polys) <- c("ID","Cont")
  polys$frog_id <- frogs[kk]
  spp<- SpatialPolygonsDataFrame(sp.comb,data=polys)
  
  allfrogs_polys <- raster::bind(allfrogs_polys,sp.comb)
  allfrogs_info <- rbind(allfrogs_info, polys)
  #Output as shape file with the individual's names and e.g. the corresponding contour level can be changed
  #writeOGR(spp,dsn="./shp_allfrogs_50HR",sprintf("Hpi_50%s",frogs[kk]), driver="ESRI Shapefile",overwrite_layer=T) 
  
  spp_all[[kk]]<-spp
  
} # end kk loop
#Output is area.frame
kud50_hpi <- area.frame%>%dplyr::rename(id = name) 

#Merge contour levels
kud_hpi <- kud95_hpi%>%
  full_join(kud50_hpi)%>%
  rename(ud95_hpi = cont_95, ud50_hpi = cont_50)
#write.csv(kud_hpi,"all_frogs_hr_kud_hpi.csv",row.names=F)

#Merge with other variables
tracks_sum <- full_join(tracks_sum, kud_hpi)

#### BOXPLOTS: DAILY MOVEMENT ####

#FIGURE EXPORT: daily movement by sex
daily_plot <-tracks_sum%>%
  ggplot(aes(x= sex, y=log(mean_cumul), fill = species)) + 
  theme_bw(20) +
  geom_boxplot(width= 0.6, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, color = "black", size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = "ln Daily movement (m)") +
  facet_wrap(~species, labeller = labeller(species=c("afemo" = "A. femoralis",
                                                     "dtinc" = "D. tinctorius",
                                                     "osylv" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) 
daily_plot

#Daily movement by sex and species, log transformed violin plots
ggplot(daily_select, aes(x= sex, y=log(daily_dist))) + 
  theme_bw(20) +
  geom_violin(aes(fill = species)) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + labs(y = "ln Daily movement (m)") +
  geom_jitter(aes(shape=sex, fill = species), position=position_jitterdodge(0.2), size = 4, alpha=0.2) + 
  facet_wrap(~species, labeller = labeller(
    species=c(afemo ="A. femoralis", dtinc = "D. tinctorius", osylv = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18), strip.text = element_text(face = "italic")) + 
  scale_x_discrete(labels= c("F", "M"))

#### BOXPLOTS: DAILY MOVEMENT BY BEHAVIOR ####

# prepare data for plotting
daily_beh_mean<- daily_beh_mean%>%
  ungroup()%>%
  filter(behavior != "pool_move" | sex != "m" | species != "osylv") %>% # remove 4 o.sylvatica males seen next to pools from parental care days
  dplyr::add_row(id = "os999", behavior = "pool_move", daily_dist = 999, max_dist = NA, sex = "m", species = "osylv")%>% #add dummy rows to keep boxplot cols same width
  dplyr::add_row(id = "dt999", behavior = "pool_move", daily_dist = 999, max_dist = NA, sex = "f", species = "dtinc")%>%
  dplyr::mutate(behavior = factor(behavior, levels=c("base", "w_op_sex", "pool_move")))

#FIGURE EXPORT: Daily movement by sex and behavior boxplot
beh_plot <- 
  ggplot(data = daily_beh_mean, aes(x= sex, y=log(daily_dist), fill = behavior)) + 
  theme_bw(20) +
  geom_boxplot(aes(fill = behavior), outlier.shape = NA,  width = 0.7, alpha = 0.7) +
  geom_jitter(aes(shape = behavior),
              position=position_jitterdodge(0.2), stroke = 1, color = "black", size = 3, alpha = 0.6) +
  scale_fill_manual(values=c("#440154FF", "#E7B800FF","#35B779FF"), labels=c("Other days", "Mating days", "Parental days")) +
  scale_shape_manual(values=c(21, 24, 23), labels=c("Other days", "Mating days", "Parental days")) +
  labs(y = "ln Daily travel (m)") +
  facet_wrap(~species, labeller = labeller(species=c("afemo" = "A. femoralis",
                                                     "dtinc" = "D. tinctorius",
                                                     "osylv" = "O. sylvatica"))) +
  theme(legend.position= c(0.84,0.85), legend.background = element_rect(size=0.5, linetype="solid",  colour ="black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic", size = 18),
        aspect.ratio = 3) +
  scale_x_discrete(labels= c("F", "M")) +
  guides(fill = guide_legend(reverse=TRUE), shape = guide_legend(reverse=TRUE), color = guide_legend(reverse=TRUE)) +
  labs(fill = "Behavior", shape = "Behavior", color = "Behavior") +
  ylim(0,5.5)
beh_plot

#### BOXPLOTS: HR ####

# #remove Jane (t06) as an outlier from the data (because she showed no site fidelity)
# mcp <- mcp[!mcp$id == "t06",]
# #remove Tarzan (t08) as an outlier from the data (because she showed no site fidelity)
# mcp <- mcp[!mcp$id == "t08",]
# head(mcp)

#FIGURE EXPORT: home range plot
hr_plot <-tracks_sum%>%
  ggplot(aes(x= sex, y=log(ud95_hpi), fill = species)) + 
  theme_bw(18) +
  geom_boxplot(width= 0.7, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, color = "black", size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = expression(ln~Home~range~(m^2))) +
  facet_wrap(~species, labeller = labeller(species=c("afemo" = "A. femoralis",
                                                     "dtinc" = "D. tinctorius",
                                                     "osylv" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(0, 8)
hr_plot 


#FIGURE EXPORT: mcp per species boxplot 
mcp_plot <-tracks_sum%>%
  ggplot(aes(x= sex, y=log(mcp100), fill = species)) + 
  theme_bw(18) +
  geom_boxplot(width= 0.7, outlier.shape = NA) +
  geom_jitter(aes(group = sex),
              position=position_jitterdodge(0.4), shape = 21, stroke = 1, color = "black", size = 2, alpha = 0.6) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + 
  labs(y = expression(Movement~extent~(Log(m^2)))) +
  facet_wrap(~species, labeller = labeller(species=c("afemo" = "A. femoralis",
                                                     "dtinc" = "D. tinctorius",
                                                     "osylv" = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18),
        strip.text = element_text(face = "italic"), aspect.ratio = 4) +
  scale_x_discrete(labels= c("F", "M"), expand = expansion(add = 1)) +
  ylim(0, NA)
mcp_plot

#Log mcp CI plot sex and species  
ggplot(tracks_sum, aes(x= sex, y=log(mcp100))) + 
  theme_bw(20) +
  stat_summary(fun.data="mean_cl_boot", geom="crossbar", aes(fill=species), size=.5, alpha=.9)  +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + labs(y = expression("100%mcp"~m^2)) +
  geom_jitter(aes(shape=sex), position=position_jitterdodge(0.2), size = 4, alpha=0.4) + 
  facet_wrap(~species, labeller = labeller(
    species=c(afemo ="A. femoralis", dtinc = "D. tinctorius", osylv = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18), strip.text = element_text(face = "italic")) + 
  scale_x_discrete(labels= c("F", "M"))

#Log KDE CI plot sex and species
ggplot(tracks_sum, aes(x= sex, y=log(ud95_hpi))) + 
  theme_bw(20) +
  stat_summary(fun.data="mean_cl_boot", geom="crossbar", aes(fill=species), size=.5, alpha=.9)  +
  scale_fill_manual(values=c("#E7B800", "#0072B2","#FC4E07")) + 
  theme(legend.position="none") + labs(y = expression("Log (95%UD ("~m^2~"))")) +
  geom_jitter(aes(shape=sex), position=position_jitterdodge(0.2), size = 4, alpha=0.4) + 
  facet_wrap(~species, labeller = labeller(
    species=c(afemo ="A. femoralis", dtinc = "D. tinctorius", osylv = "O. sylvatica"))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 18), strip.text = element_text(face = "italic")) + 
  scale_x_discrete(labels= c("F", "M"))

#### EXPORT DATA ####

#Daily movements
daily_export <- daily_full %>%
  ungroup() %>%
  dplyr::select(id, species, sex, date, relocs, behavior, temp, daily_dist) 
write.csv(daily_export, file="all_spaceuse_daily.csv")

# Export HR data for analysis
summary_export <- tracks_sum %>%
  dplyr::select(id, species, sex, duration, days_tracked, relocs, mean_cumul_dist = mean_cumul, max_dist_perid,
                mcp50, mcp95, mcp100, ud50_hpi, ud95_hpi)
write.csv(summary_export, file="all_spaceuse_hr.csv")
