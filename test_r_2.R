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


# same, row, col grid.
col_by_row <-
  expand.grid(row = sprintf('%.2d', 1:16), col = sprintf('%.2d', 1:24)) %>%
  arrange(., row)

# read out another moltenprot result : the Baseline-corrected sheet
exp1_curve <-
  read_excel('./data/cleaned_expt1/Signal_resources/Signal_results.xlsx',
              sheet = 'Baseline-corrected')
exp2_curve <-
  read_excel('./data/cleaned_expt2/Signal_resources/Signal_results.xlsx',
              sheet = 'Baseline-corrected') %>% 
  # TODO remove Temperature otherwise it stays in the data
  dplyr::select(-c('Temperature'))
exp3_curve <-
  read_excel('./data/cleaned_expt3/Signal_resources/Signal_results.xlsx',
              sheet = 'Baseline-corrected') %>%
  dplyr::select(-c('Temperature'))
exp4_curve <-
  read_excel('./data/cleaned_expt4/Signal_resources/Signal_results.xlsx',
              sheet = 'Baseline-corrected') %>%
  dplyr::select(-c('Temperature'))

exp_curve_all <-
# append all wells observations as new cols
# `=` preprend all columns name imported with the lvalue (xp1.A1)
cbind(
  xp1 = exp1_curve,
  xp2 = exp2_curve,
  xp3 = exp3_curve,
  xp4 = exp4_curve
)  %>%
# rename to Temperature
rename(., Temperature = xp1.Temperature) %>%
# add prefix
mutate(., Temperature = paste('val_t_', Temperature, sep = ''))

exp_curve_all <- exp_curve_all %>%
  # pivot transform columns into row combinations (vy creating a column name)
  pivot_longer(cols = 2:ncol(exp_curve_all)) %>%
  # pivoting again to get temperature as columns
  pivot_wider(names_from = Temperature) %>%
  # create a id column called well 
  rownames_to_column() %>% rename('well' = 'rowname') %>%
  # add the grid coordinates
  bind_cols(col_by_row) %>%
  # remove all unused cols
  dplyr::select(-c('name', 'well', 'row', 'col'))


# Add temperature headers to df
# TODO REVIEW this is sketchy 
add_tempheaders <- function(df,
                            start_temp = 37,
                            end_temp = 90) {
  # generate temperature intervals                       
  temperature_df <-
    seq(start_temp, end_temp, by = ((end_temp - start_temp) / (ncol(df) - 1))) %>%
    round(., digits = 1)
  # rewrite all temperatures!
  # TODO CHECK that, that's quite sketchy. Should convert existing temp
  for (i in 1:ncol(df)) {
    colnames(df)[i] <- paste('t_', temperature_df[i], sep = '')
  }
  message('Temperature assignments changed for ',
          ncol(df),
          ' points.')
  return(df)
}

start_temp = 37
end_temp = 90
exp_curve_all <- add_tempheaders(exp_curve_all, start_temp, end_temp)
message('Fit curves retrieved.')

write.csv(exp_curve_all, "test_exp_curve_all.csv")