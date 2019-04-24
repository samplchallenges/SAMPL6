# Analysis of log *P* predictions

General analysis of log *P* predictions include calculated vs predicted log *P* correlation plots and 6 performance statistics (RMSE, MAE, ME, R^2, linear regression slope(m), and error slope(ES)) for all the submissions. 
95%-percentile bootstrap confidence intervals of all the statistics were reported. Error slope (ES) statistic is calculated as the slope of the line fit the the QQ plot of model uncertainty predictions.

## Manifest
- `logP_analysis.py` - Python script that parses submissions and performs the analysis.
- `logP_predictions/` - This directory includes SAMPL6 type III pKa submission files.

- `analysis_outputs/` - All analysis outputs are organized under this directory.
   Also includes submission IDs assigned to each submission.
  - `error_for_each_logP.pdf` - Violin plots that show error distribution of predictions related to each experimental logP. 
  - `logPCorrelationPlots/` - This directory contains plots of predicted vs. experimental logP values with linear regression line (blue) for each method. Files are named by submission ID of each method, which can be found in `statistics_table.pdf`. In correlation plots, dashed black line has slope of 1. Dark and light green shaded areas indicate +-0.5 and +-1.0 logP unit error regions, respectively.
  - `pKaCorrelationPlotsWithSEM/` - This directory contains similar plots to the `pKaCorrelationPlots/` directory with error bars added for Standard Error of the Mean(SEM) of experimental and predicted values for submissions that reported these values. Since experimental logP SEM values are small horizontal error bars are mostly not visible.
  - `AbsoluteErrorPlots/` - This directory contains a bar plot for each method showing the absolute error for each logP prediction compared to experimental value.
  - `StatisticsTables/` - This directory contains machine readable copies of Statistics Table, bootstrap distributions of performance statistics, and overall performance comparison plots based on RMSE and MAE values.
    - `statistics_table.pdf` - A table of performance statistics (RMSE, MAE, ME, R^2, linear regression slope(m), and error slope(ES)) for all the submissions.
    - `RMSE_vs_method_plot.pdf`
    - `RMSE_vs_method_plot_colored_by_method_category.pdf`
    - `RMSE_vs_method_plot_for_Physical_category.pdf`
    - `RMSE_vs_method_plot_for_Empirical_category.pdf`
    - `RMSE_vs_method_plot_for_Mixed_category.pdf
    - `RMSE_vs_method_plot_for_Other_category.pdf` 
    - `MAE_vs_method_plot.pdf`
    - `MAE_vs_method_plot_colored_by_method_category.pdf`
    - `MAE_vs_method_plot_for_Physical_category.pdf`
    - `MAE_vs_method_plot_for_Empirical_category.pdf`
    - `MAE_vs_method_plot_for_Mixed_category.pdf`
    - `MAE_vs_method_plot_for_Other_category.pdf` 
    - `statistics_bootstrap_distributions.pdf` - Violin plots showing bootstrap distributions of performance statistics of each method. Each method is labelled by submission ID.
  - `QQPlots/` - Quantile-Quantile plots for the analysis of model uncertainty predictions. 
    
    
    
