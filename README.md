# Range Changes Montserrat
Study range-shifts in Montserrat forests birds by employing dynamic occupancy models fitted with colext in 'unmarked' using point count data from the Forest Bird Monitoring in the Centre Hills Forest Reserve, Montserrat, Lesser Antilles. 

## Bird monitoring
The data was surveyed from 2011 to 2024 (except of 2020) with three visits at 86 sampling locations yearly which were at least 200m apart from each other. Each visit of a sampling location started with a 3-minute settling time to afterwards count all birds for 10 minutues. Environmental variables were collected on site level, and on observation level. More detailed information on the monitoring protocol can be found elsewhere (Bambini et al., 2017; Oppel et al., 2014; Parashuram et al., 2015). 

## Data analysis
Description of analysis steps:  
1. Prepared data is available [here](data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData). 
2. Conduct Goodness of Fit Test on the global model with the full set of predictor variables using the [GOF loop script](scripts/GOF_pb_loop.R) for χ² Goodness-of-Fit test. 
3. Run analysis for all species using the scripts as for example [this script](scripts/MTOR_Montserrat_bird_counts_analysis_IMPROVE_Sempach.R) with the four-letter species code in the front of the script name. All the output is stored in the output file, here only data as .rds or .csv files are uploaded.
4. Make pretty plots and tables using the [script](scripts/Results_section.qmd) by loading in the output from the analysis.
