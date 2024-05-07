## So this is the first part of the R pipeline
## which massage back the data from moltenprot
## batches of 96 wells are stitched back 
## only useful data is kept

## NOTE This only test the `standard` model, should the same for `irrev`
## (filter on different columns)

.First<-function() {
    require("ctest", quietly = TRUE)
    options(editor="vscode")
}
.First()

# install.packages('tidyverse')
# install.packages('readxl')
# install.packages('stringr')
# install.packages('drc')
# install.packages('ggthemes')
# install.packages('cowplot')
# install.packages('hrbrthemes')
# install.packages('ggpubr')
# install.packages('MESS')
# install.packages('devtools')


library(tidyverse)
library(readxl)
library(stringr)
library(drc)
library(ggthemes)
library(cowplot)
library(hrbrthemes)
library(ggpubr)
library(MESS)
library(devtools)

full_df <- read_csv("/Users/antoinegerardin/RT-CETSA-Analysis/test_merged_molten_final.csv",
show_col_types = FALSE
)

# Calculate AUC for each well
calculate_auc2 <- function(df) {
  #Retrieve temperatures to be used for AUC determination.
  auc.df <- df %>%
    dplyr::select(matches('t_\\d'))
  
  #Initialize the AUC column
  df$auc <- NA
  
  # Pivot and clean each row for AUC model
  for (i in 1:nrow(auc.df)) {
    curve_vals <- auc.df[i,] %>%
      pivot_longer(cols = everything(),
                   names_to = 'temp',
                   values_to = 'response')
    curve_vals$temp <- curve_vals$temp %>%
      sub('t_', '', .)
    curve_vals$temp <- as.numeric(curve_vals$temp)
    df$auc[i] <- auc(x = curve_vals$temp, y = curve_vals$response)
  }
  message('AUC Values calculated for ', nrow(auc.df), ' wells.')
  return(df)
}

full_df <- calculate_auc2(full_df)

write.csv(full_df, "test_auc.csv")