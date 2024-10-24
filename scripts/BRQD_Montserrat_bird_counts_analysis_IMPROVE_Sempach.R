
##### Montserrat Forest Bird Counts - data import and analysis with 'unmarked' #####
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
# setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
# setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\UKOT\\Montserrat\\Analysis\\Population_status_assessment\\AnnualMonitoring\\Montserrat")

load(file = 'data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData') # change to the current year (most recent year with prepared data)

###### 2.1: set YEAR and SPECIES and rm_point that should be analysed ####
# YEAR <- 2024 # set the most recent year
SPECIES # all species prepared data is available for 
SPECIES <- c('BRQD') # fill in SPECIES the analysis should be made for
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
dim(yearlySiteCov$year_fact)
dim(yearlySiteCov$year_num)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdata[,2:ncol(occdata)], siteCovs = siteCov, yearlySiteCovs = yearlySiteCov, 
                         obsCovs = obsCov, numPrimary = numPrimary)
summary(umf) # looks good
str(umf) # looks good, but pay attention: not all variables from this data set are coded in an appropriate way (as.factor())

###### 3.4: scale numeric variables to fitting problems ####

siteCovs(umf)[c(2,4:6)] <- scale(siteCovs(umf)[c(2,4:6)]) # scale elevation, dbh and teeheight
obsCovs(umf)[c(1,4,6:8)] <- scale(obsCovs(umf)[c(1,4,6:8)]) # scale activity, rain, day and time (all numeric variables)
yearlySiteCovs(umf)[2] <- scale(yearlySiteCovs(umf)[2])
str(umf)
summary(umf)

##### 4: multi-season, dynamic single species occupancy model ####

###### 4.1: assess fit of the global model using a GOF run in another Script, ATTENTION: DECIDE FOR ONE OF THEM ####

# load pb gof test performed in another script 'GOF_pb_loop' on global model
(pb <- readRDS(file = sprintf('output/data/GOF/%s_gof_pb.rds', SPECIES))) # check results to decide
(c_hat_pb <- pb@t0[2]/mean(pb@t.star[,2])) # calculate c-hat 
c_hat_pb <- ifelse(c_hat_pb < 1, yes = 1, no = c_hat_pb) # this sets c_hat to 1 if c_hat <1
# the p-values have to be added manually from the pb object
print(pb)
p_value_SSE <- 0.992
p_value_Chisq <- 0.999
p_value_freemanTukey <- 0.984

###### 4.2: fit models for detection probability p() first for modSel ####

# fit models for detection probability p() manually, go on with fitList(), modSel()
fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se = T)
fm2 <- colext(~1, ~1, ~1, ~day, data = umf, se = T)
# add time and I(time^2) and decide which one is better
fm3 <- colext(~1, ~1, ~1, ~day+time, data = umf, se = T)
fm4 <- colext(~1, ~1, ~1, ~day+I(time^2), data = umf, se = T)
aictab(list(fm3, fm4), modnames = c('time', 'I(time^2)'), second.ord = T, c.hat = c_hat_pb) # time has lower AICc/OAICc, continue with time 
# continue with the better time predictor, either time or I(time^2)
fm5 <- colext(~1, ~1, ~1, ~day+time+rain, data = umf, se = T)
fm6 <- colext(~1, ~1, ~1, ~day+time+rain+wind, data = umf, se = T)
fm7 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm8 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)

p_fitList <- list(fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm8)
names(p_fitList) <- lapply(p_fitList, function(x) formula(x)) # set formulas as model names 
(p_modSel_df <- aictab(cand.set = p_fitList, c.hat = c_hat_pb) %>% # create a model comparison table with QAICc or AICc depending on c-hat from gof
  mutate(step = 'p'))
# best submodel for p(): day+time+rain+wind+activity+location, AICc difference to second best 0.6 - go on with this fm8 best one

###### 4.2: fit models for initial occupancy psi() first for modSel ####

fm9 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add treeheight
fm10 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm11 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add dbh
fm12 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm13 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm14 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm15 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add canopy
fm16 <- colext(~canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm17 <- colext(~dbh+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm18 <- colext(~alt+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm19 <- colext(~treeheight+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm20 <- colext(~alt+dbh+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm21 <- colext(~alt+treeheight+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm22 <- colext(~treeheight+dbh+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm23 <- colext(~alt+treeheight+dbh+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
# add interaction between alt and treeheight
fm24 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm25 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm26 <- colext(~alt:treeheight+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm27 <- colext(~alt:treeheight+dbh+canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)

# put the fitted models in a fitList() and rank them by AICc or QAICc in modSel()
psi_fitList <- list(fm7, fm9, fm10, fm11, fm12, fm13, fm14, fm15, fm16, fm17, fm18, fm19, fm20, fm21, fm22, fm23, fm24, fm25, fm26, fm27) # don't forget to include the best model from the last modeling step!
names(psi_fitList) <- lapply(psi_fitList, function(x) formula(x)) # set formulas as model names 
(psi_modSel_df <- aictab(cand.set = psi_fitList, c.hat = c_hat_pb) %>% 
    mutate(step = 'psi'))
# best sub-model for psi(): ~treeheight+canopy, AICc difference to second best of 2.37 (~treeheight+dbh+canopy)

###### 4.3: fit models for extinction and colonisation probability for modSel ####

fm28 <- colext(~treeheight + canopy, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T) # constant model
fm29 <- colext(~treeheight + canopy, ~alt, ~1, ~day+time+rain+wind+activity, data = umf, se = T) # expansion model
fm30 <- colext(~treeheight + canopy, ~1, ~alt, ~day+time+rain+wind+activity, data = umf, se = T) # contraction model 
fm31 <- colext(~treeheight + canopy, ~alt, ~alt, ~day+time+rain+wind+activity, data = umf, se = T) # shift model
fm32 <- colext(~treeheight + canopy, ~year_num, ~year_num, ~day+time+rain+wind+activity, data = umf, se = T)  # year_num model (trend), this model will exclude the possibility that observed changes are just annual changes
fm33 <- colext(~treeheight + canopy, ~year_fact, ~year_fact, ~day+time+rain+wind+activity, data = umf, se = T)  # year_fact model, this model will exclude the possibility that observed changes are just annual changes
# correct for year_num (trend) and explore alt effects
fm34 <- colext(~treeheight + canopy, ~year_num+alt, ~year_num, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year_num - expansion
fm35 <- colext(~treeheight + canopy, ~year_num, ~year_num+alt, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year_num - contraction
fm36 <- colext(~treeheight + canopy, ~year_num+alt, ~year_num+alt, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year_num - shift
# correct for year_fact and explore alt effects
fm37 <- colext(~treeheight + canopy, ~year_fact+alt, ~year_fact, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year_fact - expansion
fm38 <- colext(~treeheight + canopy, ~year_fact, ~year_fact+alt, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year_fact - contraction
fm39 <- colext(~treeheight + canopy, ~year_fact+alt, ~year_fact+alt, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year_fact - shift, also  global model

# put the fitted models in a list and rank them by QAIC in aictab
g_e_fitList <- list(constant = fm28, expansion = fm29, contraction = fm30,
                    shift = fm31, year_num = fm32, year_fact = fm33,  
                    year_num_expansion = fm34, year_num_contraction = fm35, year_num_shift = fm36, 
                    year_fact_expansion = fm37, year_fact_contraction = fm38, year_fact_shift = fm39, nullmodel = fm1)
names(g_e_fitList) <- lapply(g_e_fitList, function(x) formula(x)) # set formulas as model names
model_names <- aictab(cand.set = list(constant = fm28, expansion = fm29, contraction = fm30, shift = fm31, year_num = fm32, year_fact = fm33,  
                                      year_num_expansion = fm34, year_num_contraction = fm35, year_num_shift = fm36, year_fact_expansion = fm37, year_fact_contraction = fm38, year_fact_shift = fm39, 
                                      nullmodel = fm1), 
                      c.hat = c_hat_pb)[, 1] # get model names in the correct order 
(g_e_modSel_df <- aictab(cand.set = g_e_fitList, c.hat = c_hat_pb) %>%
    mutate(step = 'g_e',  # Add the step indicator
           model = model_names))  # include short model names
# best sub-model for g_e(): contraction ~1, ~alt, AICc difference to second best shift is 1.3

##### 5: Explore best model and export the first things ####

###### 5.1: Best model ####
best_model <- fm31 # save best model as best_model
saveRDS(best_model, file = sprintf('output/data/best_model/%s_best_model.rds', SPECIES)) # save model on local storage
summaryOD(best_model, c.hat = c_hat_pb) # adjusted summary statistics with c-hat, if c-hat = 1, there is no difference to the normal summary() function 
names(best_model) # get names from the submodels

###### 5.2: Model selection tables ####
modSel_export <- rbind(as.data.frame(g_e_modSel_df) %>% rename(formula = Modnames) %>% select(model, formula, step, everything()), 
                       as.data.frame(psi_modSel_df) %>% rename(formula = Modnames) %>% mutate(model = NA) %>% select(model, formula, step, everything()), 
                       as.data.frame(p_modSel_df) %>% rename(formula = Modnames) %>% mutate(model = NA) %>% select(model, formula, step, everything()))
modSel_export <- modSel_export %>% mutate(species = SPECIES, c_hat = c_hat_pb, p_value_SSE = p_value_SSE, p_value_Chisq = p_value_Chisq, p_value_freemanTukey = p_value_freemanTukey) %>%  # add species name and p_value from gof saved earlier
  select(species, everything()) %>% # change order of columns
  arrange(match(step, c('g_e','psi','p')), AICc) # check for correct order and for the correct information criterion (AICc or QAICc)
fwrite(modSel_export, file = sprintf('output/data/model_selection/%s_modSel.csv', SPECIES)) # export as .csv file

##### 6: draw inference from the best model ####

###### 6.1: Extract Occupancy values for population trajectory ####

# using empirical bayes estimates of occupancy averaged across all points
occupancy <- ranef(best_model)

# manipulate data structure of occupancy estimates, summarize and store data in data.frame
occupancy_data <- as.data.frame(bup(occupancy, stat = 'mean')) %>% # posterior mean by stat = ''
  gather(key = 'Year',value = 'occu') %>% # transfer from wide to long format 
  mutate(Year = as.numeric(str_replace(Year,'V',''))) %>% # delete V before number of season
  mutate(Year = Year + 2010) %>% # calculate year 
  group_by(Year) %>%
  summarise(Occupancy = mean(occu))

# calculate confidence intervals and store them in occupancy_data data.frame for plotting
occupancy_confint <- confint(occupancy, level = 0.95) # 95% CI
occupancy_data$lower_cl <- apply(occupancy_confint[,1,],2, mean)
occupancy_data$upper_cl <- apply(occupancy_confint[,2,],2, mean)

# remove data from 2020 because in this year surveys didn't take place 
occupancy_data <- occupancy_data %>% filter(!Year == '2020')

# make a quick plot out of curiosity
occupancy_data %>% ggplot() + 
  geom_point(mapping = aes(x = Year, y = Occupancy), size = 3, col = 'firebrick') + 
  geom_errorbar(aes(ymin = lower_cl,ymax = upper_cl, x = Year), width = 0.2) + 
  scale_y_continuous(limits = c(0, 1.05), breaks=seq(0,1,0.2),labels=seq(0,1,0.2))  +
  scale_x_continuous(limits = c(2011, 2024), breaks=seq(2011,2023,2),labels=seq(2011,2023,2)) + 
  labs(y="Mean occupancy", title= paste0(SPECIES, ' Occupancy Trajectory')) + 
  theme_bw()

# save occupancy_data for plotting it later
fwrite(occupancy_data, file = sprintf('output/data/occupancy_data_ranef/%s_occupancy_data.csv', SPECIES))

###### 6.2: Make predictions for colonisation and extinction if elevation is included as predictor in the best model ####

# only predict ext as only this one seems to change with elevation

# create input df with for prediction, 
nd <- data.frame(day = 0, time = 0, rain = 0, wind = 1, activity = max(umf@obsCovs$activity, na.rm = T), # use maximum bird activity, lowest wind speed, mean of time and day = 0 
                 location = 'midslope', treeheight = 0, # location with valley or midslope used, mean scaled treeheight used, should be 0 
                 dbh = 0, canopy = 0, # mean scaled dbh and canopy used, should be 0 
                 year_num = 0, year_fact = as.factor(round(numPrimary/2, 0)), # mean of scaled year, should be 0, otherwise calculated mean year and took it as factor 
                 alt = seq(from = min(umf@siteCovs$alt), # alt: scaled altitude for prediction,
                           to = max(umf@siteCovs$alt), by = 0.02),
                 elevation = rescale(seq(from = min(umf@siteCovs$alt), to = max(umf@siteCovs$alt), by = 0.02), # elevation: rescaled altitude for plotting and an exact ecological meaning
                                    to = c(min(siteCov$alt), max(siteCov$alt)), # to = output range as c()
                                    from = c(min(umf@siteCovs$alt), max(umf@siteCovs$alt)))) # from = input range as c()

# predict values for col and ext using modavgPred() on a list that just contains the best_model, c-hat, input data from nd (pred are the same if c-hat = 1)
#pred_col <- as.data.frame(modavgPred(cand.set = list(best_model), newdata = nd, c.hat = c_hat_pb, parm.type = 'gamma')) %>%  # for parm.type choose one of 'psi', 'gamma', 'epsilon', 'detect'
 mutate(Type = 'Colonisation', Elevation = nd$elevation, Species = SPECIES) %>% # add Elevation and Species for easier plotting
 rename(Predicted = mod.avg.pred, SE = uncond.se, lower = lower.CL, upper = upper.CL) %>% 
 select(Predicted, SE, lower, upper, Type, Elevation, Species) # select only the needed columns 

pred_ext <- as.data.frame(modavgPred(cand.set = list(best_model), newdata = nd, c.hat = c_hat_pb, parm.type = 'epsilon')) %>%  # for parm.type choose one of 'psi', 'gamma', 'epsilon', 'detect'
  mutate(Type = 'Extinction', Elevation = nd$elevation, Species = SPECIES) %>% # add Elevation and Species for easier plotting 
  rename(Predicted = mod.avg.pred, SE = uncond.se, lower = lower.CL, upper = upper.CL) %>% 
  select(Predicted, SE, lower, upper, Type, Elevation, Species) # select only the needed columns 

pred_colext <- bind_rows(pred_ext) # connect tables 
fwrite(pred_colext, file = sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Range_Changes_Montserrat/output/data/pred_col_ext/%s_pred_colext.csv', SPECIES)) # no predictions because alt is not included in the best model for either col or ext

# quickly plot col-ext dynamics against elevation 
pred_colext %>%
  ggplot(aes(x = Elevation, y = Predicted, colour = Type, fill = Type)) +
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.2) +
  labs(x="Elevation", y="Predicted probability", title=paste0(SPECIES, ' Elevational Range Dynamics')) +
  theme_bw()

###### 6.3: Extract all coefficients ####

# this calculates all estimates with inflated SE and CI using c-hat 
summary(best_model) # uninflated effect sizes for comparison
best_model_summary <- summaryOD(best_model, c.hat = c_hat_pb) # summary statistics corrected for overdispersion by c-hat on logit-scale
(estimates_best_model <- as.data.frame(best_model_summary$outMat) %>% 
    mutate(species = SPECIES, component = names(coef(best_model))) %>% 
    select(component, everything()) %>% 
    rename(SE = se, lower = lowlim, upper = upplim))

# export coefficient table 
fwrite(estimates_best_model, sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Range_Changes_Montserrat/output/data/coefficients_best_model/%s_best_model_coefficients.csv', SPECIES))

##### 8: export some prepared data #### 

saveRDS(umf, file = sprintf('C:/Users/filib/Documents/Praktika/Sempach/Montserrat/Range_Changes_Montserrat/output/data/prepared_data_umf/%s_unmarkedMultFrame.rds', SPECIES))

##### 9: notes ####

