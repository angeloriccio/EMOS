# EMOS

In this directory you will find two R files:

- "run_main.R"
- "run_single_instance.R"

run_main.R is the main source file. The simulation will start with the selected date and pollutant. run_single_instance.R is called by run_main and performs all the computations.

Every computation is based on several input files:

- "coords_est.csv"
- "coords_val.csv"
- "Italy_\<poll\>data_\<date\>_est.csv"
- "Italy_\<poll\>data_\<date\>_val.csv"

"coords_est.csv" stores the coordinates (UTM32) of each monitoring station selected for estimation purposes. "coords_val.csv" stores the same information, but for stations used for validation purposes.
For each selected date, the auxiliary variables used for the spatio-temporal regression are stored in the "Italy_\<poll\>data_\<date\>est.csv" file, where date is the date in dd-mmm-yyyy format.
"Italy\<poll\>data_\<date\>_val.csv" corresponds to the same information, but for validation stations.

At the end of each simulation, four files are saved:

1) Italy_\<poll\>data_\<date\>_EstimationResults.csv
2) Italy_\<poll\>data_\<date\>_ValidationResults.csv
3) Italy_\<poll\>data_\<date\>_PredictedResults.csv
4) Italy_\<poll\>data_\<date\>_PredictionGrid.csv

"Italy_\<poll\>data_\<date\>EstimationResults.csv" stores the results for stations used for estimation for the previous three days. In particular,
the "Estimation" column refers to the estimation made the INLA-SPDE model. "Italy\<poll\>data_\<date\>ValidationResults.csv" corresponds to the same information, but for validation stations. 
"Italy\<poll\>data_\<date\>PredictedResults.csv" corresponds to the information predicted for the same stations as in "Italy\<poll\>data_\<date\>EstimationResults.csv", but for the day N+1.
"Italy\<poll\>data_\<date\>_PredictionGrid.csv" refers to the prediction made on a regular grid (4x4 km resolution) covering the whole Italian country.

In the estimation file the predictors are: "WS" (wind speed), "T2" (temperature at 2 meters), "PBL00" (pbl height at 00UTC), "PBL12" (pbl height at 12UTC), "TP" (total precipitation), "POP" (total population), "ALT" (altitude), "Imperviousness" (imperviousness percentage), "Builtup" (built-up percentage), "roads" (total length of raods segments), "highways" (total length of main raods segments). "HUD", "LUD", "Agriculture", "Forests" are derived from the CORINE land cover 2018 database and correspond to the area of continuous urban fabric, discontinuous urban fabric, agricultural and foreast areas, respectively.

 
