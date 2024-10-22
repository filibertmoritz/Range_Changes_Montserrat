
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
  arrange(point)

###### 3.2: prepare obsCov #####

obsCov <- countdata %>% 
  select(year, Point, Count, Rain, Wind, day, time, activity) %>% 
  rename_with(tolower) %>% 
  filter(!point %in% rm_point) %>% 
  arrange(point, year, count)

###### 3.3: prepare yearlySiteCovs ####

yearlySiteCovYear <- obsCov %>% 
  group_by(point,year) %>%
  summarise(season = as.numeric(mean(year) - 2010)) %>% # code year as.numeric
  spread(key = year, value = season) %>% # spread a key-value pair across multiple columns
  arrange(point)
yearlySiteCov <- list(year=yearlySiteCovYear[,2:ncol(yearlySiteCovYear)])

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
yearlySiteCovs(umf)[1] <- scale(yearlySiteCovs(umf)[1])
str(umf)
summary(umf)

##### 4: multi-season, dynamic single species occupancy model ####

###### 4.1: fit models for detection probability p() first for modSel ####

# fit models for detection probability p() manually, go on with fitList(), modSel()
fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se = T)
fm2 <- colext(~1, ~1, ~1, ~day, data = umf, se = T)
fm3 <- colext(~1, ~1, ~1, ~day+time, data = umf, se = T)
fm4 <- colext(~1, ~1, ~1, ~day+time+rain, data = umf, se = T)
fm5 <- colext(~1, ~1, ~1, ~day+time+rain+wind, data = umf, se = T)
fm6 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm7 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)

# organise models in fitList before printing the modSel-table
p_fitList <- fitList('psi(.)g(.)e(.)p(.)' = fm1, 
                     'psi(.)g(.)e(.)p(~day)' = fm2, 
                     'psi(.)g(.)e(.)p(~day+time)' = fm3, 
                     'psi(.)g(.)e(.)p(~day+time+rain)' = fm4,
                     'psi(.)g(.)e(.)p(~day+time+rain+wind)' = fm5,
                     'psi(.)g(.)e(.)p(~day+time+rain+wind+activity)' = fm6,
                     'psi(.)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm7)
(p_modSel <- modSel(p_fitList))
# best submodel for p(): ~day+time+rain+wind+activity+location, AIC difference to second best 3.78 - go on with this fm7 best one  

# create data.frame from modSel table
p_modSel_df <- as(p_modSel, Class = 'data.frame') %>% 
  mutate(step = 'p')

###### 4.2: fit models for initial occupancy psi() first for modSel ####

fm8 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
# add treeheight
fm9 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm10 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
# add dbh
fm11 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm12 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm13 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm14 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
# add interaction between alt and treeheight
fm15 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm16 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)

# put the fitted models in a fitList() and rank them by AIC in modSel()
psi_fitList <- fitList('psi(.)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm7, 
                       'psi(~alt)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm8, 
                       'psi(~treeheight)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm9, 
                       'psi(~alt+treeheight)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm10,
                       'psi(~dbh)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm11, 
                       'psi(alt+dbh)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm12, 
                       'psi(treeheight+dbh)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm13, 
                       'psi(~alt+treeheight+dbh)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm14, 
                       'psi(alt:treeheight)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm15, 
                       'psi(~alt:treeheight+dbh)g(.)e(.)p(~day+time+rain+wind+activity+location)' = fm16)
print(psi_modSel <- modSel(psi_fitList)) 
# best sub-model for psi(): ~alt+treeheight, AIC difference to second best is 0.22 (~treeheight+dbh) - go on with this best one (fm10)

# create data.frame from modSel table
psi_modSel_df <- as(psi_modSel, Class = 'data.frame') %>% 
  mutate(step = 'psi')

###### 4.3: fit models for extinction and colonization probability for modSel ####

fm17 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm18 <- colext(~alt+treeheight, ~alt, ~1, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm19 <- colext(~alt+treeheight, ~1, ~alt, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm20 <- colext(~alt+treeheight, ~alt, ~alt, ~day+time+rain+wind+activity+location, data = umf, se = T)
fm21 <- colext(~alt+treeheight, ~year, ~year, ~day+time+rain+wind+activity+location, data = umf, se = T)  ## this model will exclude the possibility that observed changes are just annual changes

# put the fitted models in a fitList() and rank them by AIC in modSel()
g_e_fitList <- fitList(constant = fm17, expansion = fm18, contraction = fm19,
                            shift = fm20, year = fm21, nullmodel = fm1)
print(g_e_modSel <- modSel(g_e_fitList))
# best model: year model (fm21) has lowest AIC, the constant model is 2.73 AIC units apart

###### 4.4: Digression to check whether there is an effect of altitude after correcting for year

fm22 <- colext(~alt+treeheight, ~year+alt, ~year, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year - expansion
fm23 <- colext(~alt+treeheight, ~year, ~year+alt, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year - contraction
fm24 <- colext(~alt+treeheight, ~year+alt, ~year+alt, ~day+time+rain+wind+activity, data = umf, se = T)  # corrected year - shift

# put the fitted models again in a fitList() and rank them by AIC in modSel() - just add the new models to old fitList!
g_e_fitList <- fitList(constant = fm17, expansion = fm18, contraction = fm19,
                       shift = fm20, year = fm21, nullmodel = fm1, 
                       year_expansion = fm22, year_contraction = fm23, year_shift = fm24)
print(g_e_modSel <- modSel(g_e_fitList))
# best model: still year model (fm21) with lowest AIC, second best constant model with 2.73 AIC difference, expansion model with 4.04 AIC difference - continue with year model (fm21)

###### 4.4: export all required model selection tables ####

# create data.frame from modSel table
g_e_modSel_df <- as(g_e_modSel, Class = 'data.frame') %>% 
  mutate(step = 'g_e')

# merge data.frames from ACHU modSel tables to one data frame and export it as .csv
modSel_export <- rbind(g_e_modSel_df %>% select(model, formula, step, nPars, AIC, delta, AICwt, cumltvWt), 
                       psi_modSel_df %>% select(model, formula, step, nPars, AIC, delta, AICwt, cumltvWt), 
                       p_modSel_df %>% select(model, formula, step, nPars, AIC, delta, AICwt, cumltvWt))
modSel_export <- modSel_export %>% mutate(species = SPECIES) %>% # add species name
  select(species, everything()) %>% # change order of columns
  arrange(match(step, c('p','psi','g_e')), AIC)
write.csv(modSel_export, file = sprintf('output/data/modSel/%s_modSel_full.csv', SPECIES)) # export as .csv file

# save best model as best_model
best_model <- fm21 # fill in best model
saveRDS(best_model, file = sprintf('output/data/best_models/%s_best_model.rds', SPECIES))

##### 5: Goodness-of-Fit Test: MacKenzie-Bailey Goodness-of-Fit Test by function mb.gof.test()####

# calculate GoF-Test
(best_model_gof <- mb.gof.test(best_model, nsim = 1000)) 

# save output as Rdata
saveRDS(best_model_gof, file = sprintf('output/data/GOF/%s_gof_mb.rds', SPECIES))

##### 6: view results and extract data from best single model ####

# view models
summary(best_model) # best model based on AIC

# get names of values from model
names(best_model)

# extract everything (psi, col, ext, det) - on the logit-scale!
toExportpsi <- summary(best_model@estimates@estimates$psi) %>%
  mutate(Parameter=names(best_model@estimates@estimates$psi@estimates)) %>%
  mutate(Component="Initial Occupancy") 

toExportdet <- summary(best_model@estimates@estimates$det) %>%
  mutate(Parameter=names(best_model@estimates@estimates$det@estimates)) %>%
  mutate(Component="Detection probability")

toExportcol <- summary(best_model@estimates@estimates$col) %>%
  mutate(Parameter=names(best_model@estimates@estimates$col@estimates)) %>%
  mutate(Component="Colonization")

toExportext <- summary(best_model@estimates@estimates$ext) %>%
  mutate(Parameter=names(best_model@estimates@estimates$ext@estimates)) %>%
  mutate(Component="Extinction")

toExport <- bind_rows(toExportpsi, toExportcol, toExportext, toExportdet)

# export everything backTransformed 
toExport_transformed <- cbind(Estimates = plogis(coef(best_model)), 
                              SE_estimates_backTransformed = plogis(SE(best_model)),
                              rbind(plogis(confint(best_model, type = 'psi')), 
                                    plogis(confint(best_model, type = 'col')), 
                                    plogis(confint(best_model, type = 'ext')), 
                                    plogis(confint(best_model, type = 'det'))))

# merge data frames and export
toExport_final <- cbind(toExport, toExport_transformed) %>% 
  rename(Estimate_backTransformed = Estimates, 
         LCI_backTransformed = `0.025`, 
         UCI_backTransformed = `0.975`) %>% 
  select(Component, Parameter, everything())
toExport_final[,3:ncol(toExport_final)] <- round(toExport_final[,3:ncol(toExport_final)], digits = 4)
fwrite(toExport_final, sprintf('output/data/best_model_output_estimates/%s_best_model_output_estimates_all.csv', SPECIES))

##### 7: extraction of values and plotting ####

###### 7.1 plot relationship between colonization and extinction with altitude #####
# code mainly taken from: https://cran.r-project.org/web/packages/unmarked/vignettes/colext.html and steffen

# create input data.frame with values for the prediction, 
nd <- data.frame(day = 0, time = 0, rain = 0, wind = 1, activity = max(umf@obsCovs$activity, na.rm = T), # use maximum bird activity, lowest wind speed, mean of time and day = 0 used
                 location = 0, treeheight = 0, # location with valley or midslope used, mean scaled treeheight used, should be 0 (steffen used 15, why?)
                 dbh = 0, # mean scaled dbh used, should be 0 
                 alt = seq(from = min(umf@siteCovs$alt), 
                           to = max(umf@siteCovs$alt), by = 0.02),
                 altitude = rescale(seq(from = min(umf@siteCovs$alt), to = max(umf@siteCovs$alt), by = 0.02), 
                                    to = c(min(siteCov$alt), max(siteCov$alt)), # to = output range as c()
                                    from = c(min(umf@siteCovs$alt), max(umf@siteCovs$alt)))) # from = input range as c()
# alt: scaled altitude for prediction, altitude: rescaled altitude for plotting

# predict values by input data from nd, add column with re-scaled altitude
#pred_ext <- predict(best_model, type = 'ext', newdata = nd) %>% # don't predict ext because its constant 
#  mutate(Type = 'Extinction', Altitude = nd$altitude)
pred_col <- predict(best_model, type = 'col', newdata = nd) %>%
  mutate(Type = 'Colonization', Altitude = nd$altitude)

# connect tables and save as one .csv for plotting
# pred_ext_col <- bind_rows(pred_ext, pred_col)
fwrite(pred_col, file = sprintf('output/data/pred_col_ext/%s_pred_ext_col.csv', SPECIES)) # only export pred_col

# plot colonization~altitude and extinction~altitude together
pred_ext_col_plot <- pred_col %>%
  ggplot(aes(x = Altitude, y = Predicted, colour = Type, fill = Type)) +
  geom_line(linewidth = 3) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.2) +
  
  labs(x="Elevation (m above sea level)", y="Predicted probability", title='Antillean Crested Hummingbird') + # change Name
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) + # set begin/end to 0.98 for yellow/ext // 0 for purple/col
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  scale_y_continuous(limits = c(0, 1)) + 
  
  theme(panel.background=element_rect(fill="white", colour="black"),
        plot.margin = unit(c(50, 70, 40, 50), "pt"),
        plot.title=element_text(size=20, color='black', margin = margin(b=20)),
        axis.text=element_text(size=15, color="black"), 
        axis.title.y=element_text(size=18, margin = margin(r=15)),
        axis.title.x=element_text(size=18, margin = margin(t=15)),
        legend.text=element_text(size=15),
        legend.title = element_text(size=18),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(linewidth =.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))
pred_ext_col_plot

# save as image 
ggsave(filename = sprintf('output/plot/%s_ext_col_plot.jpg', SPECIES), 
       plot = last_plot(), scale = 1, dpi = 'retina', units = 'px', 
       width = 3000, height = 2200)
saveRDS(pred_ext_col_plot, file = sprintf('output/plot/data/%s_ext_col_plot.rds', SPECIES))

###### 7.2 plot estimated occupancy for time series ####

# using empirical bayes estimates of occupancy averaged across all points
occupancy <- ranef(best_model)

# manipulate data structure, summarize and store data in data.frame
occupancy_data <- as.data.frame(bup(occupancy, stat = 'mean')) %>% # posterior mean by stat = ''
  gather(key = 'Year',value = 'occu') %>% # transfer from wide to long format 
  mutate(Year = as.numeric(str_replace(Year,'V',''))) %>% # delete V before number of season
  mutate(Year = Year + 2010) %>% # calculate year 
  group_by(Year) %>%
  summarise(Occupancy = mean(occu)) # calculate the mean occupancy about all points in one year 

# calculate confidence intervals and store them in occupancy_data data.frame for plotting
occupancy_confint <- confint(occupancy, level = 0.95) # 95% CI
occupancy_data$lower_cl <- apply(occupancy_confint[,1,],2, mean)
occupancy_data$upper_cl <- apply(occupancy_confint[,2,],2, mean)

# remove data from 2020 because in this year surveys didn't take place 
occupancy_data <- occupancy_data %>% filter(!Year == '2020')

# save occupancy_data for later plotting
fwrite(occupancy_data, file = sprintf('output/data/ranef_occupancy_data/%s_occupancy_data.csv', SPECIES))

# plot occupancy over time from occupancy_data data.frame 
occupancy_plot <- ggplot(data = occupancy_data, aes(x = Year, y = Occupancy)) +
  geom_point(size = 3, pch = 16, colour="firebrick") +
  geom_errorbar(aes(ymin = lower_cl,ymax = upper_cl), width = 0.2) +
  scale_y_continuous(limits = c(0, 1.1), breaks=seq(0,1,0.2),labels=seq(0,1,0.2))  +
  scale_x_continuous(limits = c(2010, 2024), breaks=seq(2011,2023,2),labels=seq(2011,2023,2)) +
  labs(x="Year", y="Mean occupancy probability", title='Antillean Crested Hummingbird') +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"),
        plot.margin = unit(c(50, 70, 40, 50), "pt"),
        plot.title=element_text(size=20, color='black', margin = margin(b=20)),
        axis.text=element_text(size=15, color="black"), 
        axis.title.y=element_text(size=18, margin = margin(r=15)),
        axis.title.x=element_text(size=18, margin = margin(t=15)),
        legend.text=element_text(size=15),
        legend.title = element_text(size=18),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))
occupancy_plot

# save as image 
ggsave(filename = sprintf('output/plot/%s_occupancy_plot.jpg', SPECIES), 
       plot = occupancy_plot, scale = 1, dpi = 'retina', units = 'px', 
       width = 3000, height = 2200)
saveRDS(occupancy_plot, file = sprintf('output/plot/data/%s_occupancy_plot.rds', SPECIES))

##### 8: export some prepared data #### 

saveRDS(umf, file = sprintf('output/data/prepared_data/%s_unmarked_mult_frame.rds', SPECIES))

##### 9: notes ####

