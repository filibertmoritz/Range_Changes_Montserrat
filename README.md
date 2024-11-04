# Range_Shifts_Montserrat
Study range-shifts in Montserrat forests birds using dynamic occupancy models in 'unmarked'

This is a description of which steps and scripts were applied in which order: 

1. prepare data using the .... script that extracts raw data from data/Montserrat_Birds_2024.accdb and export data as data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData (ATTENTION: data prep script still has to be added - steffen has the most recent version!)
2. conduct Goodness of Fit Test using scripts/GOF_mb_loop.R (for MacKenzie-Bailey Goodness-of-Fit test) or scripts/GOF_mb_loop.R (for parboot test)
3. run analysis for all species using the scripts ... all the output is stored in the output file, here only data as .rds or .csv files are exported
4. make pretty plots and tables using the script ... by loading in the output from the analysis 
