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
# Dose-response curve fit with LL.4 log-logistic fit
# otrace = TRUE; output from optim method is displayed. Good for diagnostic
# trace = TRUE; trace from optim displayed
# robust fitting
#   robust = 'lms': doesn't handle outlier/noisy data well
dr_fit <- function(df) {
  try(expr = {
    # TODO CHHECK what this is
    # dose response model. Comes from drc package
    # https://www.rdocumentation.org/packages/drc/versions/2.5-12/topics/drm
    drm(
      resp ~ conc,
      data = df,
      type = 'continuous',
      # four-parameter log-logistic function
      fct = LL.4(),
      control = drmc(
        errorm = FALSE,
        maxIt = 10000,
        noMessage = TRUE,
      )
    )
  })
}

# NOTE HELPER FUNCTION
# Fit the alternate model log logistic to a set of conc ~ resp values.
# Returns: RSS values
# Requires: same as null_model
fit_altmodel <- function(df,
                         plot.model,
                         graphTitle = '') {
  # model fit with logistic function                        
  alt.model <- dr_fit(df)
  if (plot.model == TRUE) {
    try(jpeg(filename = paste('./data/models/', graphTitle, '.jpg', sep =
                                '')))
    try(plot(
      alt.model,
      main = graphTitle,
      pch = 21,
      cex = 3,
      col = 'black',
      bg = 'cyan',
      lwd = 3
    ))
    try(dev.off()
    )
  }
  message('Alternate Model RSS: ',
          (sum(residuals(alt.model) ^ 2)))
  return(sum(residuals(alt.model) ^ 2))
}

# NOTE helper function
# Fit the null model to a set of conc ~ resp values.
# Returns: RSS values
# Requires; df with concentration and response values
# df[1]: concentration values
# df[2]: response values
fit_nullmodel <- function(df,
                          plot.model,
                          graphTitle = '') {
  # fit resp to linear model where resp = 1                       
  null.model <- lm(resp ~ 1, data = df)

  #Plot if TRUE. For diagnostic use mostly...
  if (plot.model == TRUE) {
    try(jpeg(filename = paste('./data/models/', graphTitle, '.jpg', sep =
                                '')))
    try(plot(
      df$conc,
      df$resp,
      main = graphTitle,
      pch = 21,
      cex = 3,
      col = 'black',
      bg = 'orange',
      lwd = 3
    ))
    try(abline(null.model, col = 'black', lwd = 3))
    try(dev.off()
    )
  }

  #Return squared residuals for null model
  message('NUll Model RSS: ',
          (sum(residuals(null.model) ^ 2)))
  return(sum(residuals(null.model) ^ 2))
}

# NOTE helper function.
# Creates a thermogram with all concentrations of a target for plotting
# Must match exact ncgc_id in well assignment..
dr.thermogram <- function(df, target = '') {
  # Create df with the compound, conc, and temperature columns
  
  df <- df %>%
  # pivot on the temperature datapoints 
  #(create a temp column with temp label and resp column with values)
    pivot_longer(cols = 3:ncol(df),
                 names_to = 'temp',
                 values_to = 'resp')
  # to allow plotting
  df$temp <- as.numeric(sub('t_', '', df$temp))

  # write.csv(df, "test_rss.csv")
  
  # for each compound, draw thermogram for each concentration
  #(dose response thermogram)
  # fill allows to associated a color label automatically
  dr.plot <- ggplot(df, aes(
    y = resp,
    x = temp,
    fill = as.factor(signif(conc)),
    group_by(signif(conc))
  )) +
    geom_line(color = 'black',
              alpha = 0.8,
              size = 1) +
    geom_point(shape = 21, size = 3) +
    theme_minimal() +
    scale_color_viridis_d() +
    labs(
      title = paste('Dose-Response Thermogram for ', target, sep = ''),
      x = 'Temperature [C]',
      y = 'Response',
      fill = 'Concentration'
    ) +
    theme()
  print(dr.plot)
  ggsave(
    filename = paste('./data/models/dr_', target, '.png', sep = ''),
    scale = 1.25,
    dpi = 'retina'
  )
  return(dr.plot)
}

 # Derive RSS values for null and alternate model for each compound from full_df
compute.rss.models <-
  function(df,
           control = 'DMSO',
           plotModel = TRUE,
           rssPlot = TRUE,
           drPlot = TRUE) {
    #Construct tibble of unique compounds names
    rss.df <- tibble(compound = (unique(
      filter(df, ncgc_id != control | ncgc_id != 'vehicle')$ncgc_id
    ))) %>%
      filter(compound != 'control') %>%
      filter(compound != 'vehicle')

    # add empty columns for what we want to calculate
    rss.df$null.model.n <- NA
    rss.df$alt.model.n <- NA
    rss.df$null.model.sum <- NA
    rss.df$alt.model.sum <- NA
    rss.df$null.model.sd <- NA
    rss.df$alt.model.sd <- NA
    rss.df$rss.diff <- NA
    rss.df$mannwhit.pval <- NA
    rss.df$mannwhit.ec50 <- NA
    
    # iterating over compounds
    for (i in 1:nrow(rss.df)) {
      #Construct df for current compound
      # filter rows pertaining to this compounds 
      # (we may have any number of exp at different concentration)
      fit.df <- df %>% filter(ncgc_id == toString(rss.df[i, 1])) %>%
        # select id, concentration, and intensities
        dplyr::select(ncgc_id, conc, starts_with('t_')) %>%
        # not quite sure how those are generated byt filter them out
        dplyr::select(!contains('onset'))

      #Plot out dose-response thermogram here?
      # dr stands for dose-response
      if (drPlot == TRUE) {
        dr.thermogram(fit.df, target = rss.df$compound[i])
      }

      
      
      #Construct a df to hold the rss values until final calculations of mean,sd,N
      # keep track of all temperature columns for each compounds and concentration
      cmpnd.fit.df <- fit.df %>%
        dplyr::select(starts_with('t_'))
      # create a new "tibble" with temp column collecting all temperature as rows
      # and add to placeholder for values to calculate
      cmpnd.fit.df <- tibble(temp = colnames(cmpnd.fit.df))
      cmpnd.fit.df$null <- NA
      cmpnd.fit.df$alt <- NA

      # write.csv(cmpnd.fit.df, "test_rss.csv")
      
      #Iterate through each temperature, construct df, perform rss analysis, and add to cmpnd.fit.df
      for (t in 3:ncol(fit.df)) {
        # for each temperature, get the response at each concentration/value
        current.fit.df <- fit.df %>%
          dplyr::select(1:2, colnames(fit.df)[t])
        # rename temp label to concentration.
        colnames(current.fit.df)[3] <- 'resp'
        # add a column

        # compute residual standard square (rss) error for null model
        cmpnd.fit.df$null[t - 2] <-
          fit_nullmodel(current.fit.df,
                        plot.model = plotModel,
                        graphTitle = as.character(paste(
                          current.fit.df[1, 1], ' Null Model at ', colnames(fit.df)[t], sep = ''
                        )))

        # compute rss for alt model       
        cmpnd.fit.df$alt[t - 2] <-
          fit_altmodel(current.fit.df,
                       plot.model = plotModel,
                       graphTitle = as.character(paste(
                         current.fit.df[1, 1], ' Alternate Model at ', colnames(fit.df)[t], sep = ''
                        )))
        
      }

      # RSS0-RSS1
      # create a new column with the rss diff
      cmpnd.fit.df <- cmpnd.fit.df %>%
        mutate(diff = null - alt)
      
      # update rss.df with proper stats
      #Now, we calculate and assign rss values for both models in the rss.df for this compound.
      rss.df$null.model.n[i] <- length(na.omit(cmpnd.fit.df$null))
      rss.df$alt.model.n[i] <- length(na.omit(cmpnd.fit.df$alt))
      rss.df$null.model.sum[i] <- sum(cmpnd.fit.df$null)
      rss.df$alt.model.sum[i] <- sum(cmpnd.fit.df$alt)
      rss.df$null.model.sd[i] <- sd(cmpnd.fit.df$null)
      rss.df$alt.model.sd[i] <- sd(cmpnd.fit.df$alt)
      rss.df$rss.diff[i] <-
        sum(cmpnd.fit.df$null) - sum(cmpnd.fit.df$alt)
      
      #Perform Mann-Whitney iU test on alternative vs. null model dataframe for compound.
      mann.whit <-
        wilcox.test(x = cmpnd.fit.df$null,
                    y = cmpnd.fit.df$alt,
                    exact = TRUE)
      rss.df$mannwhit.pval[i] <- mann.whit$p.value
      
      #Message out RSS0-RSS1 and p value
      message('RSS Difference for ',
              rss.df[i, 1],
              ': ',
              rss.df$rss.diff[i])
      message('Mann-Whitney U Test p-val: ',
              rss.df$mannwhit.pval[i])
      
      # Construct drc model and derive ec50 if p-val is significant
      if (rss.df$mannwhit.pval[i] <= 0.05) {
        #Find what temperature is the max point
        rss.max.temp <-
          cmpnd.fit.df$temp[cmpnd.fit.df$diff == max(cmpnd.fit.df$diff)]
        #Construct df ready for drc at max temperature
        rss.drc.df <- fit.df %>%
          dplyr::select(conc, one_of(rss.max.temp))
        colnames(rss.drc.df)[2] <- 'resp'
        
        # same dose response fit  for the temperature of max response
        rss.drc.model <- dr_fit(rss.drc.df)

        # collect ec50 from the dr package
        ec50.temp <- rss.drc.model$coefficients[4]
        if (length(ec50.temp != 0)) {
          rss.df$mannwhit.ec50[i] <- signif(ec50.temp, 3)
        }
      }
      
      #Plot the RSS values across the temperature range if true
      if (rssPlot == TRUE) {
        # First, clean up the temperatures
        cmpnd.fit.df$temp <- sub('t_', '', cmpnd.fit.df$temp)
        cmpnd.fit.df$temp <- as.numeric(cmpnd.fit.df$temp)
        
        #Plot RSS as
        rss.plot <- ggplot(cmpnd.fit.df, aes(x = temp, y = diff)) +
          geom_point(shape = 21,
                     size = 4,
                     fill = '#AA4499') +
          theme_minimal() +
          labs(
            title = paste(current.fit.df[1, 1], ' RSS Difference', sep = ''),
            subtitle = paste('Mann-Whitney U pval: ', signif(rss.df$mannwhit.pval[i])),
            sep = '',
            y = 'RSS0-RSS1',
            x = 'Temperature [C]'
          )
        print(rss.plot)
        ggsave(
          filename = paste('./data/models/', current.fit.df[1, 1], '_rss.png', sep =
                             ''),
          scale = 1.25,
          dpi = 'retina'
        )
      }
    }
    return(rss.df)
  }

full_df <- read_csv("test_auc.csv",
show_col_types = FALSE
)

rss <- compute.rss.models(full_df, rssPlot = TRUE, drPlot = TRUE, plotModel = TRUE)
write.csv(rss, "test_rss.csv")