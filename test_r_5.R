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

full_df <- read_csv("test_auc.csv",
show_col_types = FALSE
)


control_grouping <- function(df, control = 'DMSO', pc = 'control') {
  # return the controls, minus the concentration column with is meaningless
  control_df <- filter(df, ncgc_id == control | ncgc_id == pc)
  if (nrow(control_df) == 0) {
    message('No control wells found. Review control input to function.')
  } else
    if (nrow(control_df) > 0) {
      control_df <- control_df %>%
        dplyr::select(-'conc')
      return(control_df)
    }
}

# certain wells are annotated as such
control <- 'vehicle'
pc <- 'control' 

control_df <- control_grouping(full_df, control, pc) # Pull out control compound datapoints

control_variability <-
  function(df, nc = 'vehicle', pc = 'control') {
    #Filter out positive and negative controls into their own df

    write.csv(df, "test_control.csv")

    # vehicle df
    nc.controls.df <- df %>%
      filter(ncgc_id == nc) %>%
      dplyr::select(-c('ncgc_id', 'well', 'row', 'col'))

    # control df
    pc.controls.df <- df %>%
      filter(ncgc_id == pc) %>%
      dplyr::select(-c('ncgc_id', 'well', 'row', 'col'))
    
    #Calculate means, sd, and %CV
    nc.mean.df <-
      # compute mean on 
      apply(nc.controls.df[1:ncol(nc.controls.df)], 2, mean)
    
    nc.sd.df <- apply(nc.controls.df[1:ncol(nc.controls.df)], 2, sd)
    pc.mean.df <-
      apply(pc.controls.df[1:ncol(pc.controls.df)], 2, mean)
    pc.sd.df <- apply(pc.controls.df[1:ncol(pc.controls.df)], 2, sd)
    
    #Calculate %CV
    nc.var.df <- tibble(nc.mean = nc.mean.df, nc.sd = nc.sd.df) %>%
      mutate(nc.cv = (nc.sd / nc.mean) * 100)
    pc.var.df <- tibble(pc.mean = pc.mean.df, pc.sd = pc.sd.df) %>%
      mutate(pc.cv = (pc.sd / pc.mean) * 100)
    analysis_method <- colnames(nc.controls.df)
    var_df <- cbind(analysis_method, nc.var.df, pc.var.df)
    message('Control group variability analyzed.')
    return(var_df)
  }

control_var <-
  control_variability(control_df) # Read out the control group variability

write.csv(control_var, "test_variability.csv")

## After that is analysis of those results

# This is an helper function used later
# Returns thermogram with mean/sd of DMSO curve across temps
control_thermogram <- function(df, pcTm, ncTm) {
  subset_df <- subset(df, grepl('t_', analysis_method)) %>%
    mutate(temp = as.numeric(gsub('t_', '', analysis_method))) %>%
    dplyr::select(-'analysis_method')
  therm_plot <- ggplot(subset_df, aes(x = temp)) +
    geom_line(aes(y = nc.mean),
              size = 1.5,
              alpha = 0.75,
              color = '#88CCEE') +
    geom_errorbar(aes(ymin = nc.mean - nc.sd, ymax = nc.mean + nc.sd),
                  size = 0.5,
                  width = 1) +
    geom_point(
      aes(y = nc.mean),
      size = 3.25,
      shape = 21,
      color = 'black',
      fill = '#88CCEE'
    ) +
    geom_line(aes(y = pc.mean),
              size = 1.5,
              alpha = 0.75,
              color = '#882255') +
    geom_errorbar(aes(ymin = pc.mean - pc.sd, ymax = pc.mean + pc.sd),
                  size = 0.5,
                  width = 1) +
    geom_point(
      aes(y = pc.mean),
      size = 3.25,
      shape = 21,
      color = 'black',
      fill = '#EE3377'
    ) +
    theme_minimal() +
    labs(title = 'Control Thermograms',
         x = 'Temperature [C]',
         y = 'Fraction Unfolded')
  print(therm_plot)
  return(therm_plot)
}


# Controls analysis and z' output for groups
# Possible outputs:
# output = 'plot': Cowplot of controls
# output = 'df': Control dataframe
control_analysis <- function(
  df,
  nc = 'vehicle',
  pc = 'control',
  output = '',
  controlDF) {

    # recreate the subset of vehicle and control observations
    controls.df <- df %>%
      filter(ncgc_id == nc | ncgc_id == pc)
    
    #Calculate Z' from controls for each parameter
    # test params are Tm_fit and auc
    test_params <-
      c('Tm_fit',
        'auc')

    # specifically recalculate the mean, sd, z' for Tm_fit
    Tm.nc.mean <-
      mean(controls.df$Tm_fit[controls.df$ncgc_id == nc])
    Tm.nc.sd <- sd(controls.df$Tm_fit[controls.df$ncgc_id == nc])
    Tm.pc.mean <-
      mean(controls.df$Tm_fit[controls.df$ncgc_id == pc])
    Tm.pc.sd <- sd(controls.df$Tm_fit[controls.df$ncgc_id == pc])
    Tm.z <-
      1 - (((3 * Tm.pc.sd) + (3 * Tm.nc.sd)) / abs(Tm.pc.mean - Tm.nc.mean))
    
    message('Z\' for Tm: ', signif(Tm.z))

    # specifically recalculate the mean, sd, z' for AUC
    auc.nc.mean <- mean(controls.df$auc[controls.df$ncgc_id == nc])
    auc.nc.sd <- sd(controls.df$auc[controls.df$ncgc_id == nc])
    auc.pc.mean <- mean(controls.df$auc[controls.df$ncgc_id == pc])
    auc.pc.sd <- sd(controls.df$auc[controls.df$ncgc_id == pc])
    auc.z <-
      1 - (((3 * auc.pc.sd) + (3 * auc.nc.sd)) / abs(auc.pc.mean - auc.nc.mean))
    message('Z\' for AUC: ', signif(auc.z))
    
    # Generate graph
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

controlPlot <-
  control_analysis(
    full_df,
    nc = 'vehicle',
    pc = 'control',
    output = 'plot',
    controlDF = control_var
  )

print(controlPlot)