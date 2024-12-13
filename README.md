# Range_Shifts_Montserrat
Study range-shifts in Montserrat forests birds using dynamic occupancy models fitted using colext in 'unmarked'

This is a description of which steps and scripts were applied in which order: 

1. Prepared data is available at [data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData](data/MONTSERRAT_ANNUAL_DATA_INPUT2024.RData). 
2. Conduct Goodness of Fit Test on the global model with the full set of predictor variables using [file](scripts/GOF_pb_loop.R) for χ² Goodness-of-Fit test. 
3. Run analysis for all species using the scripts as for example [file](scripts/MTOR_Montserrat_bird_counts_analysis_IMPROVE_Sempach.R) with the four-letter species code in the front of the script name. All the output is stored in the output file, here only data as .rds or .csv files are exported.
4. Make pretty plots and tables using the script [file](scripts/Results_section.qmd) by loading in the output from the analysis.
