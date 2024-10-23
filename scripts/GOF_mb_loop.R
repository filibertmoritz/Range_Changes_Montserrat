
##### Montserrat Forest Bird Counts - GOF test as loop over all species #####
##### written in Sep 2024 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

# install.packages('scales')
library(scales)
# install.packages('data.table')
library(data.table)
library(reshape)
library(lubridate)
library(unmarked)
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename

##### 2: load the prepared data and make last preparation #####

# set working directory and load prepared data 
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
# setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")

load(file = 'data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData') # change to the current year (most recent year with prepared data)

###### 2.1: set YEAR and SPECIES and rm_point that should be analysed ####
# YEAR <- 2024 # set the most recent year
SPECIES # all species prepared data is available for 
rm_point<-c(99,76) # remove two points which are not independent

###### 2.2: check the data and remove unneeded stuff #####

head(countdata) # table with all birdcounts
head(obsCov) # all observation level Covs
rm(obsCov) # they are not needed since all information is available in countdata
head(siteCov) # site level Covs
(species_names <- species) # species codes and their English names
rm(species)

##### 3: prepare data for unmarkedMultFrame ####

###### 3.1: prepare siteCovs ####

head(siteCov)
siteCov <- siteCov %>% 
  rename_with(tolower) %>% 
  rename(canopy = canopy_cover, alt = elevation) %>% 
  filter(!point %in% rm_point) %>%  # remove all 2 points that are not independent
  arrange(point)

###### 3.2: prepare obsCov #####

obsCov <- countdata %>% 
  select(year, Point, Count, Rain, Wind, day, time, activity) %>% 
  rename_with(tolower) %>% 
  filter(!point %in% rm_point) %>% 
  arrange(point, year, count)

###### 3.3: prepare yearlySiteCovs ####

yearlySiteCovYear_fact <- obsCov %>% # these are the yearlysiteCovs year factor for random effects
  group_by(point,year) %>%
  summarise(season = as.factor(mean(year) - 2010)) %>% # code year as.numeric
  spread(key = year, value = season) %>% # spread a key-value pair across multiple columns
  arrange(point)
yearlySiteCovYear_num <- obsCov %>% # these are the yearlysiteCovs year as.numeric for trend hypothesis 
  group_by(point,year) %>%
  summarise(season = as.numeric(mean(year) - 2010)) %>% # code year as.numeric
  spread(key = year, value = season) %>% # spread a key-value pair across multiple columns
  arrange(point)
yearlySiteCov <- list(year_fact=yearlySiteCovYear_fact[,2:ncol(yearlySiteCovYear_fact)], year_num = yearlySiteCovYear_num[,2:ncol(yearlySiteCovYear_num)])

###### 3.3: prepare numPrimary (number of survey years) #### 

numPrimary <- length(unique(countdata$year)) # calculates number of years

###### 3.4: prepare countdata ####

head(countdata)
countdata <- countdata %>% select(year, Point, Count, all_of(SPECIES)) %>% # select all columns needed and all SPECIES
  rename(point = Point, count = Count) %>% 
  filter(!point %in% rm_point) # remove the points that should not be analysed

##### 4: start loop over all species ####

# set number of simulations for goodness of fit test and create input objects
nsim <- 3 # number of simulations should be very high (>1000)
gof <- data.frame()
gof_objects <- list()

for(i in 1:length(SPECIES)){
  
  # prepare occupancy data for each species 
  occdata <- countdata %>% 
    select(year, point, count, SPECIES[i]) %>% 
    rename(n = SPECIES[i]) %>% 
    mutate(occupancy = ifelse(n > 0, 1, 0)) %>% # convert in detection/non-detection
    mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
    select(-n, -year, -count) %>% # remove unneeded columns 
    spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
    arrange(point) # sort rows in order of point number
  
  # create unmarkedMultFrame 
  umf <- unmarkedMultFrame(y = occdata[,2:ncol(occdata)], siteCovs = siteCov, yearlySiteCovs = yearlySiteCov, 
                           obsCovs = obsCov, numPrimary = numPrimary)
  
  # scale numeric variables to fitting problems 
  siteCovs(umf)[c(2,4:6)] <- scale(siteCovs(umf)[c(2,4:6)]) # scale elevation, dbh and teeheight
  obsCovs(umf)[c(1,4,6:8)] <- scale(obsCovs(umf)[c(1,4,6:8)]) # scale activity, rain, day and time (all numeric variables)
  yearlySiteCovs(umf)[2] <- scale(yearlySiteCovs(umf)[2])
  
  # built global model 
  global_model <- colext(~alt:treeheight+dbh+canopy, ~year_fact+alt, ~year_fact+alt, ~day+time+I(time^2)+rain+wind+activity+location, data = umf, se = T) # global model which is a year_fact corrected shift model
  
  # perform MacKenzie-Bailey goodness of fit test by calling the function mb.gof.test
  SPECIES[i]
  (mb <- mb.gof.test(global_model, nsim = nsim, parallel = T)) # perform gof, increase nsim
  c_hat <- mb$c.hat.est # save c-hat estimates
  p_value = mb$p.value # save p_value for modsel table

  # save results in data.frame and RDS
  result <- data.frame(species = SPECIES[i], c_hat = c_hat , p_value = p_value)
  gof <- rbind(gof, result)
  gof_objects[[i]] <- mb
  names(gof_objects)[i] <- SPECIES[i]
  saveRDS(mb, file = sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Range_Changes_Montserrat/output/data/GOF/%s_gof_mb.rds', SPECIES)) # save mb gof as RDS
}


##### 5: check results #### 

# save data in a data.frame 
gof_mb_list <- list()
parameter <- c('species', 'c_hat', 'p_value') 
gof_mb_df <- data.frame(array(data = NA, dim = c(length(SPECIES), length(parameter)))) # prepare df 
names(gof_mb_df) <- parameter

for(i in 1:length(SPECIES)){
  gof_mb_list[[i]] <- readRDS(file = sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/GOF/%s_gof_mb.rds', SPECIES[i]))
  names(gof_mb_list)[i] <- SPECIES[i]
  c_hat <- gof_mb_list[[i]]$c.hat.est # extract c-hat
  p_value <- gof_mb_list[[i]]$p.value # extract p-value
  gof_mb_df[i,] <- c(SPECIES[i], c_hat, p_value)
}


