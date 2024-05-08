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

# NOTE this process creates spurious columns that should be removed

curve_df <- read_csv("/Users/antoinegerardin/RT-CETSA-Analysis/test_exp_curve_all.csv",
show_col_types = FALSE
)

full_param <- read_csv("/Users/antoinegerardin/RT-CETSA-Analysis/test_exp_param_full.csv",
show_col_types = FALSE
)
platemap_filepath = './data/platemap.xlsx'

# Assign compound ids and concentration from platemap
plate_assignment <- function(df, platemap_file) {
  # read sample sheet from plate file
  id_df <- read_excel(platemap_file, sheet = 'sample') %>%
    # remove first column
    dplyr::select(-1) %>%
    # pivot to get row, col coordinates as columns
    pivot_longer(cols = 1:ncol(.))  %>%
    # rename the columns with all ids
    rename(ncgc_id = value) %>%
    # remove name column
    dplyr::select(-c('name'))
    # NOTE `EMPTY` are considered as vehicle
    id_df$ncgc_id <- gsub('empty', 'vehicle', id_df$ncgc_id)

    # read the concentration from the file
    conc_df <- read_excel(platemap_file, sheet = 'conc') %>%
    # remove first colum
    dplyr::select(-1) %>%
    # pivot
    pivot_longer(., cols = 1:ncol(.)) %>%
    rename(conc = value) %>%
    dplyr::select(-c('name'))

    # add the columns to the datset
    df <- cbind(id_df, conc_df, df)
  message('Plate assignment attached to dataframe.')

  # make sure we have numeric value? (unecessary?)
  df$row <- as.numeric(df$row)
  df$col <- as.numeric(df$col)
  return(df)
}


full_df <- full_param

full_df <- plate_assignment(full_df, platemap_filepath)

# Construct full data frame with curve fit and parameters for analysis
bind_fulldf <- function(param_df, curve_df) {
  df <- cbind(param_df, curve_df)
  return(df)
}

# Concat dataframes.
full_df <- bind_fulldf(full_df, curve_df)

#Convert any columns containing Kelvin values from MoltenProt to Celsius
kelToCel <- function(df) {
  df <- df %>%
    mutate(Tm_fit = Tm_fit - 273.15) %>%
    mutate(T_onset = T_onset - 273.15)
}

# TODO move that before for each dataset
full_df <- full_df %>% dplyr::select(-c('...1')) %>% dplyr::select(-c('...1')) 
print(names(full_df))

full_df <- kelToCel(full_df)

write.csv(full_df, "test_merged_molten_final.csv")