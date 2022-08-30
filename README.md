# CETSA
#### Michael Ronzetti, NIH/NCATS/UMD 

____________________________________________________________
## Packages Needed:
-library(tidyverse)
-library(drc)
-library(stringr)
-library(readxl)
-library(stringr)
-library(ggthemes)
-library(cowplot)
-library(hrbrthemes)
-library(ggpubr)
____________________________________________________________
## Functions List:

  - construct_grid(row_num = 16, col_num = 24, pad_num = FALSE)
    - Construct grid with the option for headers compatible with moltenprot
    
  - prepMatLabforMolt(file_loc = './data/rtcetsa_raw.xlsx', sheet = 'Sheet1', col_names = FALSE,start_temp = 37, end_temp = 90)
    - Read in MatLab file and construct back to appropriate plate format.
    
  - retrieveMoltenData(model = 'standard', plate_format = 384)
    - Takes in moltenprot files and reconstructs according to input model. Need to expand for other types.

  - retrieve_FittedCurves(model = 'baseline-fit', start_temp = 37, end_temp = 90)
    - Gather base-line corrected fit curves for the 384-well plate and pivot plate
    - model = c('fit_curves','baseline-fit')
      - fit_curves: Curve fitting to cleaned raw data
      - baseline-fit: Baseline correction and normalization to [0-1].
  
  - bind_fulldf(param_df, curve_df)
    - Construct full data frame with curve fit and parameters for analysis
  
  - kelToCel(df)
    - Convert any columns containing Kelvin values from MoltenProt to Celsius
  
  - add_tempheaders(df, start_temp = 37,end_temp = 90)
  
  - add_rowcol(df, well_num)
  
  - well_assignment()
    - Process temperature headers coming from moltenprot.
  
  - plate_assignment()
    - Bring in platemaps to assign concentration and id to wells. Plate maps need to follow a specific format found in the template file.
  
  - calculate_auc(df)
    - Calculates AUC from curve fit data at each temperature for each well.
    
  - control_grouping(df, control = 'DMSO')
    - Creates a df using input control wells
    
  - control_variability(df)
    - Reads out control group variability
  
  - control_thermogram()
    - Output control group thermogram and %CV.
    
  - dr_fit()
  
  - dr_analysis()
  
  - plate_heatmap()
  
  - export_heatmaps()
  
  - dmso_rss()
  
  - compare_models()
  
  - calculate_meltingparams()
    - Takes full_df and returns the melting curve parameters derived from MoltenProt after subtracting the vehicle
    control from each well.
    
  - plot_volcanos()
  
  - parameter_doseresponse()
  
  - calculate_zscore()
  
  - convert_zscore()
  
  - fit_nullmodel()
  
  - fit_altmodel()
  
  - compute.rss.models()
  
  - compute_parameter.rssmodel()
  
  - dr.thermogram()
____________________________________________________________  
## To Do List (long..)

### Current

  [ ] Fix heatmaps to handle temperatures not integers.
  [X] Kelvin to Celsius
  [X] Compute AUC for each model. 
  [X] Dose-response thermogram
  [X] Control analysis also looks at positive control.
  [X] RSS Diff vs. p-val
  [X] Max Param Diff vs. RSS Diff
  [X] Plot of compounds thermogram parameters vs RSS diff, colored by p-val in one, either, or none.
  [X] EC50 @ Max RSS Difference in thermogram
  [X] EC50 and p-value heatmaps for parameters
  [X] Control Comparison

### MoltenProt Handling

  [ ] Ability to handle different models coming in from moltenprot, or pivot to irreversible model
  [ ] Can we do QC on some of the parameters, slopes, baselines?
  [X] Analysis of different baseline handling options
  
### DRC

  [ ] Detect and handle outliers in data
  [ ] Curve quality calling
  [ ] Ability to fit difference classes of curves using drc
  [ ] Inter and intra plate replicates
