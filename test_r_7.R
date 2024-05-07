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

# NOTE: 2nd occurence of this HELPER FUNCTION
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
# Controls analysis and z' output for groups
# Possible outputs:
# output = 'plot': Cowplot of controls
# output = 'df': Control dataframe
control_analysis <-
  function(df,
           nc = 'vehicle',
           pc = 'control',
           output = '',
           controlDF) {
    controls.df <- df %>%
      filter(ncgc_id == nc | ncgc_id == pc)
    
    #Calculate Z' from controls for each parameter
    # TODO CHECK not done in previous step already?
    test_params <-
      c('Tm_fit',
        'auc')
    Tm.nc.mean <-
      mean(controls.df$Tm_fit[controls.df$ncgc_id == nc])
    Tm.nc.sd <- sd(controls.df$Tm_fit[controls.df$ncgc_id == nc])
    Tm.pc.mean <-
      mean(controls.df$Tm_fit[controls.df$ncgc_id == pc])
    Tm.pc.sd <- sd(controls.df$Tm_fit[controls.df$ncgc_id == pc])
    Tm.z <-
      1 - (((3 * Tm.pc.sd) + (3 * Tm.nc.sd)) / abs(Tm.pc.mean - Tm.nc.mean))
    
    message('Z\' for Tm: ', signif(Tm.z))
    auc.nc.mean <- mean(controls.df$auc[controls.df$ncgc_id == nc])
    auc.nc.sd <- sd(controls.df$auc[controls.df$ncgc_id == nc])
    auc.pc.mean <- mean(controls.df$auc[controls.df$ncgc_id == pc])
    auc.pc.sd <- sd(controls.df$auc[controls.df$ncgc_id == pc])
    auc.z <-
      1 - (((3 * auc.pc.sd) + (3 * auc.nc.sd)) / abs(auc.pc.mean - auc.nc.mean))
    message('Z\' for AUC: ', signif(auc.z))
    
    if (output == 'plot') {
      Tm.plot <-
        ggplot(controls.df, aes(x = ncgc_id, y = Tm_fit, fill = ncgc_id)) +
        geom_boxplot(outlier.alpha = 0, size = 0.75) +
        geom_jitter(shape = 21, size = 3) +
        theme_minimal() +
        scale_fill_hue() +
        labs(title = 'Controls | Tagg',
             subtitle = paste('Z\': ', signif(Tm.z), sep = '')) +
        theme(
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, face = 'bold'),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 12, face = 'bold')
        )
      auc.plot <-
        ggplot(controls.df, aes(x = ncgc_id, y = auc, fill = ncgc_id)) +
        geom_boxplot(outlier.alpha = 0, size = 0.75) +
        geom_jitter(shape = 21, size = 3) +
        theme_minimal() +
        scale_fill_hue() +
        labs(title = 'Controls | AUC',
             subtitle = paste('Z\': ', signif(auc.z), sep = '')) +
        theme(
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, face = 'bold'),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 12, face = 'bold')
        )
      right.grid <-
        plot_grid(Tm.plot, auc.plot, ncol = 1)
      control.grid <-
        plot_grid(
          control_thermogram(controlDF, ncTm = Tm.nc.mean, pcTm = Tm.pc.mean),
          right.grid,
          ncol = 2,
          nrow = 1
        )
      ggsave('./data/controls.png', dpi = 'retina', scale = 1.5)
      return(control.grid)
    }
    if (output == 'df') {
      means <-
        c(Tm.nc.mean,
          auc.nc.mean)
      parameters <- c('Tm_fit', 'auc')
      output.df <- tibble(parameters, means)
      return(output.df)
    }
  }

# Compute the rss difference and significance for each of the parameters
compute_parameter.rssmodel <- function(df, plotModel = FALSE) {
  #Test parameters for standard model
  test_params <-
    c('Tm_fit',
      'auc')
  
  #Construct df of unique compounds and initialize parameter readouts.
  param.rss.df <- tibble(compound = (unique(
    filter(df, ncgc_id != 'control' & ncgc_id != 'vehicle')$ncgc_id
  )))
  param.rss.df$Tm_fit.ec50 <- as.numeric(NA)
  param.rss.df$Tm_fit.pval <- as.numeric(NA)
  param.rss.df$Tm_fit.maxDiff <- as.numeric(NA)
  param.rss.df$auc.ec50 <- as.numeric(NA)
  param.rss.df$auc.pval <- as.numeric(NA)
  param.rss.df$auc.maxDiff <- as.numeric(NA)
  
  # NOTE seems already done before. Remove duplicates
  control.means <- control_analysis(df, output = 'df')
  
  for (i in 1:nrow(param.rss.df)) {
    # filter row by compounds
    cmpnd.fit.df <- df %>%
      filter(ncgc_id == param.rss.df$compound[i])
    #Now iterate through columns in test_params
    for (p in 1:length(test_params)) {
      current_param <- test_params[p]
      current.fit.df <- cmpnd.fit.df %>%
        # note not sure I() is needed here
        dplyr::select(ncgc_id, conc, I(test_params[p]))
      colnames(current.fit.df)[3] <- 'resp'
      current.model <- dr_fit(current.fit.df)
      
      #Workaround to avoid drm that can't converge
      if (class(current.model) != 'list') {
        param.rss.df[i, paste(current_param, '.pval', sep = '')] <-
          # https://www.rdocumentation.org/packages/drc/versions/3.0-1/topics/noEffect
          noEffect(current.model)[3]
        
        # get rss from the fit
        param.rss.df[i, paste(current_param, '.ec50', sep = '')] <-
          summary(current.model)$coefficients[4]
        
        #Calculate the maximum difference in param and subtract negative control mean from it.
        current.fit.df$absDiff <-
          abs(current.fit.df$resp - control.means$means[control.means$parameters ==
                                                          current_param])
        param.rss.df[i, paste(current_param, '.maxDiff', sep = '')] <-
          current.fit.df$resp[current.fit.df$absDiff == max(current.fit.df$absDiff)] - control.means$means[control.means$parameters ==
                                                                                                             current_param]
        
        message('Analyzing Compound ', param.rss.df[i, 1], '...')
        # Tm_fit or AUC
        message(current_param)
        message('EC50: ', param.rss.df[i, paste(current_param, '.ec50', sep =
                                                  '')])
        message('No Effect ANOVA p-val: ', signif(param.rss.df[i, paste(current_param, '.pval', sep =
                                                                          '')]), 1)
        if (plotModel == TRUE) {
          png(
            filename = paste(
              './data/models/',
              param.rss.df[i, 1],
              '_',
              current_param,
              '.png',
              sep = ''
            ),
            bg = 'transparent'
          )
          plot(
            current.model,
            main = paste(
              param.rss.df[i, 1],
              '\n',
              ' NoEffect pval: ',
              signif(param.rss.df[i, paste(current_param, '.pval', sep =
                                             '')]),
              '\n',
              'EC50: ',
              signif(param.rss.df[i, paste(current_param, '.ec50', sep =
                                             '')]),
              '\n',
              current_param
            )
          )
          dev.off()
        }
      }
    }
  }
  return(param.rss.df)
}


full_df <- read_csv("test_auc.csv",
show_col_types = FALSE
)

#Perform dose-response for each thermogram parameter
parameters <- compute_parameter.rssmodel(full_df, plotModel = TRUE)

write.csv(parameters, "test_parameters.csv")