
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
nsim <- 1000 # number of simulations should be very high (>1000)
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
  
  # define function for Goodness-of-fit test, taken from https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html#occupancy-model
  fitstats <- function(global_model) { # first create function
    observed <- getY(global_model@data)
    expected <- fitted(global_model)
    resids <- residuals(global_model)
    sse <- sum(resids^2,na.rm=TRUE)
    chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
  }
  
  # perform goodness of fit by calling the function within parboot
  SPECIES[i]
  pb <- parboot(global_model, fitstats, nsim=nsim, report=1)  ### increase nsim
  pb # check result
  c_hat <- pb@t0[2]/mean(pb@t.star[,2]) # calculate c-hat
  
  # save results in data.frame and RDS
  result <- data.frame(species = SPECIES[i], SSE = pb@t0[1] , SSE_p = pb@t0[4], Chisq = pb@t0[2], Chisq_p = pb@t0[4], freemanTukey = pb@t0[3], freemanTukey_p = pb@t0[4], c_hat = c_hat)
  gof <- rbind(gof, result)
  gof_objects[[i]] <- pb
  names(gof_objects)[i] <- SPECIES[i]
  #saveRDS(pb, file = sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/GOF/%s_gof_pb.rds', SPECIES[i])) # save the gof as .rds data 
}


##### 5: check results #### 

# save data in a data.frame 
gof_list <- list()
parameter <- c('species', 'SSE', 'SSE_p', 'Chisq', 'Chisq_p', 'freemanTukey', 'freemanTukey_p', 'c_hat') 
gof_df <- data.frame(array(data = NA, dim = c(length(SPECIES), length(parameter)))) # prepare df 
names(gof_df) <- parameter

for(i in 1:length(SPECIES)){
  gof_list[[i]] <- readRDS(file = sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/GOF/%s_gof_pb.rds', SPECIES[i]))
  names(gof_list)[i] <- SPECIES[i]
  c_hat <- gof_list[[i]]@t0[2]/mean(gof_list[[i]]@t.star[,2]) # calculate c-hat
  gof_df[i,] <- c(SPECIES[i], as.numeric(gof_list[[i]]@t0[1]), NA, as.numeric(gof_list[[i]]@t0[2]), NA, as.numeric(gof_list[[i]]@t0[3]), NA, as.numeric(c_hat))
}


# manually add p-value, format results and export
gof_list$MTOR
gof_list$FOTH
gof_list$BRQD
gof_list$TREM
gof_list$ACHU
gof_list$PTCA
gof_list$PETH
gof_list$GTCA
gof_list$SBTH
gof_list$SNPI
gof_list$CAEL
gof_list$BANA


print(SPECIES)
gof_df$SSE_p <- c(MTOR = 0.962, FOTH = 0.651, BRQD = 0.992, TREM = 0.660, ACHU =  0.934, PTCA = 0.864, PETH = 0.583, GTCA = 0.617, SBTH = 0.707, SNPI = 0.971, CAEL = 0.457, BANA = 0.579)
gof_df$Chisq_p <- c(MTOR = 0.993, FOTH = 0.928, BRQD = 0.999, TREM = 0.987, ACHU =  0.794, PTCA = 0.941, PETH = 0.646, GTCA = 0.806, SBTH = 0.916, SNPI = 0.993, CAEL = 0.428, BANA = 0.569)
gof_df$freemanTukey_p <- c(MTOR = 0.973, FOTH = 0.633, BRQD = 0.984, TREM = 0.825, ACHU =  0.857, PTCA = 0.829, PETH = 0.504, GTCA = 0.517, SBTH = 0.578, SNPI = 0.922, CAEL = 0.376, BANA = 0.538)
gof_df[, 'SSE'] <- round(as.numeric(gof_df[, 'SSE']), digits = 3)
gof_df[, 'Chisq'] <- round(as.numeric(gof_df[, 'Chisq']), digits = 3)
gof_df[, 'freemanTukey'] <- round(as.numeric(gof_df[, 'freemanTukey']), digits = 3)
gof_df$c_hat <- round(as.numeric(gof_df$c_hat), digits = 3)

gof_df %>% format(scientific = F)

fwrite(gof_df %>% format(scientific = F), file = 'C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Range_Changes_Montserrat/output/data/GOF/GOF_pb_table.csv')
