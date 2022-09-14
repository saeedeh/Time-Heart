# Time-Heart #

This repository contains the code for a study on the relationship between heart rate and subsecond time perception. Experiment data is available at: https://osf.io/u4586/
<br/>	

## Code structure ##


### psychometric
Routines for psychometric fitting the obtain bias and threshold for each subject in the bisection task using the Palamedes toolbox
<br/>

### hddm
Routines for DDM analysis using the hddm package.
fitted models are saved in "04-data"
<br/>


### useful_functions.py and peakdet_modules.py
These files only contain functions that are used in other files
<br/>

### 01-check_data.py
Some routines to inspect and check the quality of the ECG and behavioral data.
<br/>

### 00-setup.py
Initial processing of the ECG and trial data are done at
The Result is 3 data farmes (all_trials, all_subjects, and df_blocks) saved as pickle files and then loaded and further processed at "00-setup-load.py"\
<br/>

### 00-setup-load.py
load the initially processed trial and subject data and do more processing, e.g. filter outliers, estimate features from each trial's RR time-series, etc. 
Resutls are added as new columns to the same data structures (all_trials and all_subjects)
This file also contains code for exporting the resulting dataframe to rpy2 to use R functions in python (for regression analysis).
<br/>

### 00-setup-baselineHRV.py
Extract the baseline HRV measures (recodred during rest), and analyze the result

### "00-setup-interoception.py"
Analyze the heartbeat counting task
<br/>

### 02-basic analysis.py
Regression analyezes
<br/>


### 04-timing_analysis
Plotting RR signal across trials overall, as well as in different conditions
<br/>
