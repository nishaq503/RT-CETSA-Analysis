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

signif.df <- read_csv("test_signif.csv",
show_col_types = FALSE
)

# Export heatmaps of EC50 and P-values across analysis parameters
# Pass in parameters df
parameter_heatmaps <- function(df, plotHeat = FALSE) {
  ec50.heat.df <- df %>%
    dplyr::select(compound, contains('ec50')) %>%
    pivot_longer(cols = !compound,
                 names_to = 'parameter',
                 values_to = 'ec50') %>%
    mutate(ec50 = log10(ec50))
  ec50.heat.df$parameter <- ec50.heat.df$parameter %>%
    sub('.ec50', '', .)
  
  ec50.heat.plot <-
    ggplot(ec50.heat.df,
           aes(
             x = parameter,
             y = compound,
             fill = ec50,
             label = signif(ec50)
           )) +
    geom_tile(color = 'black') +
    geom_text(alpha = 0.85, size = 2.5) +
    theme_minimal() +
    scale_fill_gradientn(colors = c('#EE3377', '#DDCC77', '#88CCEE'), ) +
    labs(title = 'EC50 Parameter Comparison',
         fill = 'Log EC50') +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, face = 'bold')
    )
  
  pval.heat.df <- df %>%
    dplyr::select(compound, contains('pval')) %>%
    pivot_longer(cols = !compound,
                 names_to = 'parameter',
                 values_to = 'pval') %>%
    mutate(sigVal = ifelse(pval < (0.05 / length(unique(
      df
    ))), 'Significant', 'Insignificant'))
  pval.heat.df$parameter <- pval.heat.df$parameter %>%
    sub('.pval', '', .)
  pval.heat.plot <-
    ggplot(pval.heat.df,
           aes(
             x = parameter,
             y = compound,
             fill = sigVal,
             label = signif(pval)
           )) +
    geom_tile(color = 'black') +
    geom_text(alpha = 0.85, size = 2.5) +
    theme_minimal() +
    labs(title = 'P-Value Parameter Comparison',
         fill = 'P-Value') +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, face = 'bold'),
      
    )
  if (plotHeat == TRUE) {
    print(pval.heat.plot)
    ggsave('./data/models/pval_heatmap.png',
           dpi = 'retina',
           scale = 1.25)
    print(ec50.heat.plot)
    ggsave('./data/models/ec50_heatmap.png',
           dpi = 'retina',
           scale = 1.25)
  }
}

# Plot out RSS Difference by p-value for the MannWhitney
rss.pval.plot <- function (df, savePlot = FALSE) {
  plot.df <- df %>%
    dplyr::select(compound, rss.diff, mannwhit.pval, mannwhit.ec50)
  plot.df$mannwhit.pval <- log2(plot.df$mannwhit.pval)
  
  rss.plot <-
    ggplot(plot.df,
           aes(x = rss.diff, y = mannwhit.pval, fill = mannwhit.ec50)) +
    geom_point(shape = 21, size = 3.5) +
    theme_minimal() +
    scale_fill_gradient(low = '#EE3377',
                        high = '#88CCEE',
                        na.value = 'grey20') +
    labs(
      title = 'RSS vs. Mann Whitney P-val',
      x = 'RSS0 - RSS1',
      y = 'Log2 Mann Whitney P-val',
      fill = 'NPARC EC50'
    )
  print(rss.plot)
  if (savePlot == TRUE) {
    ggsave('./data/models/rssPvalcomp.png',
      dpi = 'retina',
      scale = 1.25
    )
  }
}

# Print out the volcano plots for each parameter and RSS vs. p-val
plot_volcanos <- function(df, save = TRUE) {
  test_params <-
    c('Tm_fit.maxDiff',
      'auc.maxDiff')
  test_pval <-
    c('Tm_fit.maxDiff',
      'auc.maxDiff')
  
  # Plot out RSS Difference(x) vs. Parameter Difference(y)
  # Conditional fill: grey/alpha if not significant in either
  #   grey/alpha if not significant in either #DDDDDD
  #   teal if by parameter only #009988
  #   orange if by NPARC only #EE7733
  #   wine if by both #882255
  # NEED TO CODE THIS BETTER WTF
  for (i in 1:length(test_params)) {
    current_param <- test_params[i]
    current_pval <- test_pval[i]
    plot.df <- df %>%
      dplyr::select(compound,
                    rss.diff,
                    mannwhit.pval,
                    one_of(current_param),
                    one_of(current_pval))
    # Assign significance testing outcomes
    plot.df$sigVal <-
      case_when((plot.df$mannwhit.pval < 0.05 &
                   plot.df[, current_pval] < 0.05) ~ 'Both',
                (plot.df$mannwhit.pval < 0.05 &
                   plot.df[, current_pval] >= 0.05) ~ 'RSS NPARC',
                (plot.df$mannwhit.pval >= 0.05 &
                   plot.df[, current_pval] < 0.05) ~ 'Parameter',
                (plot.df$mannwhit.pval >= 0.05 &
                   plot.df[, current_pval] >= 0.05) ~ 'Insignificant'
      )
    
    fillvalues <-
      c('Both', 'RSS NPARC', 'Parameter', 'Insignificant')
    colors <- c('#882255', '#EE7733', '#009988', '#DDDDDD')
    volcano_plot <-
      ggplot(plot.df,
             aes(x = rss.diff,
                 y = plot.df[, current_param],
                 label = compound)) +
      geom_point(shape = 21,
                 aes(fill = sigVal),
                 size = 5) +
      theme_minimal() +
      labs(
        title = paste('Residual Variance vs. ', current_param, sep = ''),
        y = paste(current_param, ' Experimental - Vehicle Mean', sep = ''),
        x = 'RSS0 - RSS1 NPARC',
        fill = 'Significance Detected'
      ) +
      scale_fill_manual(breaks = fillvalues, values = colors) +
      theme(legend.position = 'bottom')
    print(volcano_plot)
    ggsave(
      paste('./data/', current_param, '_volcano.png', sep = ''),
      dpi = 'retina',
      scale = 1.25
    )
  }
}

# NOTE DOES NOT WORK
# plot_volcanos(signif.df)

# Plot of RSS Differences vs. p-values for NPARC
rss.pval.plot(signif.df, savePlot = TRUE)

#Heatmap of compounds vs. different measurement styles.
parameter_heatmaps(signif.df, plotHeat = TRUE)
