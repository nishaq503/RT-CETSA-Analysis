# #  @@@@@@@   @@@@@@@              @@@@@@@  @@@@@@@@  @@@@@@@   @@@@@@    @@@@@@   
# #  @@@@@@@@  @@@@@@@             @@@@@@@@  @@@@@@@@  @@@@@@@  @@@@@@@   @@@@@@@@  
# #  @@!  @@@    @@!               !@@       @@!         @@!    !@@       @@!  @@@  
# #  !@!  @!@    !@!               !@!       !@!         !@!    !@!       !@!  @!@  
# #  @!@!!@!     @!!    @!@!@!@!@  !@!       @!!!:!      @!!    !!@@!!    @!@!@!@!  
# #  !!@!@!      !!!    !!!@!@!!!  !!!       !!!!!:      !!!     !!@!!!   !!!@!!!!  
# #  !!: :!!     !!:               :!!       !!:         !!:         !:!  !!:  !!!  
# #  :!:  !:!    :!:               :!:       :!:         :!:        !:!   :!:  !:!  
# #  ::   :::     ::                ::: :::   :: ::::     ::    :::: ::   ::   :::  
# # NonParametric Multiparameter Analysis of CETSA/RT-CETSA Experimental Sets
# #
# # Written by: Michael Ronzetti {NIH/NCATS/UMD}
# # Patents: PCT/US21/45184, HHS E-022-2022-0-US-01
# # Main Analysis

library(tidyverse)
source('./functions.R')
# EXPERIMENTAL PARAMETERS AND SETUP
#
# Input experiment parameters here

startTemp <- 37
endTemp <- 95
plate_format <- 384
control <- 'vehicle'
pc <- 'control'

# NOTE A-  that the code to use moltenprot so we don't need that 
# Prepare the MatLab file for MoltenProt processing
# raw_df <-
#   prepMatLabforMolt(file_loc = './data/example_plate.xlsx',
#                     start_temp = startTemp,
#                     end_temp = endTemp)

# Read in the processed MoltProt data and prepare the data frames.
full_param <- retrieveMoltenData(model = 'standard')
curve_df <-
  retrieve_FittedCurves(model = 'baseline-fit',
                        start_temp = startTemp,
                        end_temp = endTemp)
full_df <- full_param %>%
  plate_assignment(., './data/platemap.xlsx')
full_df <- bind_fulldf(full_df, curve_df) %>%
  kelToCel(.)

full_df <- calculate_auc(full_df)

# Perform some preliminary control group analysis of variability
control_df <-
  control_grouping(full_df, control, pc) # Pull out control compound datapoints
control_var <-
  control_variability(control_df) # Read out the control group variability
controlPlot <-
  control_analysis(
    full_df,
    nc = 'vehicle',
    pc = 'control',
    output = 'plot',
    controlDF = control_var
  )
print(controlPlot)

# #Calculate melting parameter difference for each well from MoltenProt
# # full_df <- calculate_meltingparams(full_df) %>%
# #   calculate_zscore() %>%
# #   convert_zscore

# #Derive RSS values for null and alternate model for each compound from full_df
# rss <- compute.rss.models(full_df, rssPlot = TRUE, drPlot = TRUE, plotModel = FALSE)

# #Perform dose-response for each thermogram parameter
# parameters <- compute_parameter.rssmodel(full_df, plotModel = TRUE)

# #Merge these plots for further analysis
# signif.df <- merge(rss, parameters)
# colnames(signif.df)[9] <- 'mannwhit.pval'
# signif.df <- determineSig(signif.df)
# signif.df <- rankOrder(signif.df)

# # Volcano plots comparing the different parameters of analysis against the NPARC RSS Difference
# # Colored by significance test and whether the compound passes any.
# plot_volcanos(signif.df)

# # Plot of RSS Differences vs. p-values for NPARC
# rss.pval.plot(signif.df, savePlot = TRUE)

# #Heatmap of compounds vs. different measurement styles.
# parameter_heatmaps(signif.df, plotHeat = TRUE)

# #Write out signif.df and full_df
# write.csv(x = full_df, file = './data/full_df.csv')
# write.csv(x = signif.df, file = './data/signif_df.csv')

