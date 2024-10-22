
##### Montserrat Forest Bird Counts - data import and analysis with 'ubms' #####
##### written in Sep 2024 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

# install.packages('scales')
library(scales)
# install.packages('data.table')
library(data.table)
library(reshape)
library(lubridate)
library('ubms')
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename
# install.packages('ubms')
library('ubms')

##### 2: load the prepared data and make last preparation #####

# set working directory and load prepared data 
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
load(file = 'data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData') # change to the current year (most recent year with prepared data)

###### 2.1: set YEAR and SPECIES and rm_point that should be analysed ####
# YEAR <- 2024 # set the most recent year
SPECIES # all species prepared data is available for 
SPECIES <- c('ACHU') # fill in SPECIES the analysis should be made for
rm_point<-c(99,76) # remove two points which are not independent

###### 2.2: check the data and remove unneeded stuff #####

head(countdata) # table with all birdcounts
head(obsCov) # all observation level Covs
rm(obsCov) # they are not needed since all information is available in countdata
head(siteCov) # site level Covs
species # species codes and their English names

##### 3: prepare data for unmarkedMultFrame ####

###### 3.1: prepare siteCovs ####

head(siteCov)
siteCov <- siteCov %>% 
  rename_with(tolower) %>% 
  rename(canopy = canopy_cover, alt = elevation) %>% 
  filter(!point %in% rm_point) %>%  # remove all 2 points that are not independent
  arrange(point) %>% mutate(point = as.factor(point))

###### 3.2: prepare obsCov #####

obsCov <- countdata %>% 
  select(year, Point, Count, Rain, Wind, day, time, activity) %>% 
  rename_with(tolower) %>% 
  filter(!point %in% rm_point) %>% 
  arrange(point, year, count) %>% 
  mutate(point = as.factor(point))

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

###### 3.3: prepare countdata ####

head(countdata)
countdata <- countdata %>% select(year, Point, Count, all_of(SPECIES)) %>% # select all the columns, including the correct species for analysis
  rename(n = all_of(SPECIES)) %>% # this filters for the one species that
  rename_with(tolower) %>% 
  filter(!point %in% rm_point) # remove the points that should not be analysed

occdata <- countdata %>% 
  mutate(occupancy = ifelse(n > 0, 1, 0)) %>% # convert in detection/non-detection
  mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
  select(-n, -year, -count) %>% # remove unneeded columns 
  spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
  arrange(point) # sort rows in order of point number

###### 3.4: prepare numPrimary (number of survey years) #### 

numPrimary <- length(unique(countdata$year)) # calculates number of years

###### 3.6: check dimensions and input data.frames ####

head(siteCov)
head(obsCov)
head(yearlySiteCov)
head(occdata)

dim(occdata)
dim(siteCov)
nrow(obsCov)/(3*numPrimary)
dim(yearlySiteCovYear)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdata[,2:ncol(occdata)], siteCovs = siteCov, yearlySiteCovs = yearlySiteCov, 
                         obsCovs = obsCov, numPrimary = numPrimary)
summary(umf) # looks good
str(umf) # looks good, but pay attention: not all variables from this data set are coded in an appropriate way (as.factor())

###### 3.4: scale numeric variables to fitting problems ####

siteCovs(umf)[c(2,4:6)] <- scale(siteCovs(umf)[c(2,4:6)]) # scale elevation, dbh and teeheight
obsCovs(umf)[c(4,6:8)] <- scale(obsCovs(umf)[c(4,6:8)]) # scale activity, rain, day and time (all numeric variables)
yearlySiteCovs(umf)[2] <- scale(yearlySiteCovs(umf)[2])
str(umf)
summary(umf)

##### 4: multi-season, dynamic single species occupancy model in 'ubms' ####

fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se = T) # the 'unmarked' model using frequentist stats 
fm1_b <- stan_colext(~1, ~1, ~1, ~1, data = umf, se = T, chains = 3, iter = 1000) # the 'ubms' model using the bayesian approach, chains specifies the MCMC and iter the number of iterations within each MCMC
fm2 <- colext(~alt+treeheight, ~1, ~1, ~activity, data = umf, se = T)
fm2_b <- stan_colext(~alt+treeheight, ~1, ~1, ~activity, data = umf, se = T, chains = 3, iter = 1000) 

# check MCMC convergence 
summary(fm1_b, submodel = 'state') # MCMC convergence seems to be okay (Rhat ~ 1)
traceplot(fm1_b, inc_warmup = T) # visually check for convergence failure in MCMC chains, grey is burn-in/warmup

# check the model fit 
resid <- residuals(fm2_b, submodel = 'state')
ubms::plot_residuals(fm2_b, submodel = 'state') # hm, that does not look perfect

# calculate MacKenzie-Bailey goodness-of-fit
?ubms::gof
gof_fm1_b <- ubms::gof(fm1_b, draws = 100) # that's somehow not possible because of NAs in the occupancy data 

# produce a modsel table with waic
fitList_b <- fitList(fm1_b, fm2_b) # first create the fitList
modSel(fitList_b) %>% mutate(waic = -2*elpd, delta_waic = -2*elpd_diff) # create modSel table and calculate waic by hand 

# calculate waic and loo
waic(fm2_b) # recommends loo 
loo(fm2_b) # the estimates are equal



# try to troubleshoot the problems with the gof by removing year 2020 and filling up all NAs by 0 in the occdata 

umf_NA <- unmarkedMultFrame(y = occdata[,2:ncol(occdata)] %>% 
                           select(-`2020_1`, -`2020_2`, -`2020_3`) %>% 
                           mutate(across(everything(), ~ replace_na(., 0))), 
                         siteCovs = siteCov, 
                         yearlySiteCovs = lapply(yearlySiteCov, function(yearlySiteCovYear) yearlySiteCovYear %>% select(-`2020`)), 
                         obsCovs = obsCov %>% filter(!year == 2020), 
                         numPrimary = numPrimary-1)
str(umf_NA)
siteCovs(umf_NA)[c(2,4:6)] <- scale(siteCovs(umf_NA)[c(2,4:6)]) # scale elevation, dbh, teeheight and canopy
obsCovs(umf_NA)[c(4,6:8)] <- scale(obsCovs(umf_NA)[c(4,6:8)]) # scale activity, rain, day and time (all numeric variables)
yearlySiteCovs(umf_NA)[2] <- scale(yearlySiteCovs(umf_NA)[2])
str(umf_NA)
summary(umf_NA)

fm1_NA <- colext(~1, ~1, ~1, ~1, data = umf_NA, se = T) # the 'unmarked' model using frequentist stats 
fm1_b_NA <- stan_colext(~1, ~1, ~1, ~1, data = umf_NA, se = T, chains = 3, iter = 100) # the 'ubms' model using the bayesian approach, chains specifies the MCMC and iter the number of iterations within each MCMC
fm2_NA <- colext(~alt+treeheight, ~1, ~1, ~activity, data = umf_NA, se = T)
fm2_b_NA <- stan_colext(~alt+treeheight+(1|point), ~1, ~1, ~activity, data = umf_NA, se = T, chains = 3, iter = 1000) 
fm3_b_NA <- stan_colext(~alt:treeheight+dbh, ~year+alt, ~year+alt, ~day+time+rain+wind+activity+location, data = umf_NA, se = T, chains = 3, iter = 1000)
fm3_b_NA_random <- stan_colext(~alt:treeheight+dbh+(1|point), ~year_num+alt, ~year_num+alt, ~day+time+rain+wind+activity+location+(1|point+year_fact), data = umf_NA, se = T, chains = 3, iter = 1000) # this is a nested random effect 
fm3_b_NA_random2 <- stan_colext(~alt:treeheight+dbh+(1|point), ~year_num+alt, ~year_num+alt, ~day+time+rain+wind+activity+location+(1|point)+(1|year_fact), data = umf_NA, se = T, chains = 3, iter = 20000) #  
fm3_b_NA_random3 <- stan_colext(~alt:treeheight+dbh, ~year_num+alt, ~year_num+alt, ~day+time+rain+wind+activity+location+(1|point)+(1|year_fact), data = umf_NA, se = T, chains = 3, iter = 1000) # random effect only for det with point and year

# calculate goodness-of-fit
(fm3_b_NA_random_gof <- gof(fm3_b_NA_random, draws = 1000))
(fm3_b_NA_random2_gof <- gof(fm3_b_NA_random2, draws = 1000)) # ACHU psi with random effect for point, and 20 000 iterations
(fm3_b_NA_random3_gof <- gof(fm3_b_NA_random3, draws = 1000)) # ACHU psi without random effect and 1 000 iterations
fm3_b_NA_gof # the ppp is 0 which indicates lack of fit









