## Calculate thermogram for each compound for each concentration

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

# NOTE HELPER FUNCTION
# Mutates a binary variable testing each analysis method for significance
# 0 if insignificant
# 1 if significant
determineSig <- function(df, alpha = 0.05) {
  analysisMethods <- dplyr::select(df, contains('pval'))
  analysisMethods <- colnames(analysisMethods)
  analysisMethodsNames <- sub('pval', 'pval.sig', analysisMethods)
  sigVal <- alpha / nrow(df)
  for (i in 1:length(analysisMethods)) {
    df[, analysisMethodsNames[i]] <-
      ifelse(df[, analysisMethods[i]] < sigVal, 1, 0)
    df[is.na(df)] <- 0
  }
  return(df)
}

# Use after determineSig from above
rankOrder <- function(df) {
  # look at the significance columns
  methodSig <- dplyr::select(df, contains('pval.sig'))
  methodSig <- colnames(methodSig)
  methodRank <- sub('pval.sig', 'rankOrder', methodSig)
  methods <- sub('.rankOrder', '', methodRank)
  methodsEC <- paste(methods, '.ec50', sep = '')
  
  # NOTE compute some ranking
  for (i in 1:length(methods)) {
    rank.df <- filter(df, df[, (methodSig[i])] == 1)
    rank.df[, methodRank[i]] <-
      as.integer(rank(rank.df[, methodsEC[i]]))
    df <- left_join(df, rank.df)
  }
  return(df)
}

full_df <- read_csv("test_auc.csv",
show_col_types = FALSE
)

rss <- read_csv("test_rss.csv",
show_col_types = FALSE
)

parameters <- read_csv("test_parameters.csv",
show_col_types = FALSE
)

#Merge these plots for further analysis
# basically a outer join but there is no overlap between columns
signif.df <- merge(rss, parameters)



# NOTE WRONG should be colnames(signif.df)[10] but anyway too hacky
# to rely on index
# colnames(signif.df)[9] <- 'mannwhit.pval'

# check p value for parameter to determine significance
# and also on the whitemann test
signif.df <- determineSig(signif.df)


signif.df <- rankOrder(signif.df)

write.csv(signif.df , "test_signif.csv")