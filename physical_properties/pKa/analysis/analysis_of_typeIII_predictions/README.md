# Analysis of the typeIII macroscopic pKa predictions

General analysis of typeIII macroscopic pKa predictions include calculated vs predicted pKa correlation plots and 5 performance statistics (RMSE, MAE, ME, R^2 and linear regression slope(m) ) for all the submissions. 
95%-percentile bootstrap confidence intervals of all the statistics were reported. 

One-to-one matching of predicted pKas to experimental pKas was performed based on minimum absolute error principle.
Predicted pKas were matched to experimental pKa values that minimize the absolute error.
When more than one predicted pKa were matched to the same experimental pKa, only the predicted pka that has the lowest absolute error was kept. 
Predicted pKas that were not matched to experimental pKas were excluded from this analysis.

## Manifest
- `typeIII_analysis.py` - Python script that parses submissions and performs the analysis.
- `typeIII_predictions/` - This directory includes SAMPL6 type III pKa submission files.
- `analysis_outputs/` - All analysis outputs are organized under this directory.
  - `statistics_table.pdf` - A table of performance statistics (RMSE, MAE, ME, R^2 and linear regression slope(m)) for all the submissions. Also includes submission IDs assigned to each submission.
  - `error_for_each_macroscopic_pKa.pdf` - Violin plots that show error distribution of predictions related to each experimental pKa. 
  - `pKaCorrelationPlots/` - This directory contains plots of predicted vs. experimental pKa values with linear regression line (blue) for each method. Files are named by submission ID of each method, which can be found in `statistics_table.pdf`. In correlation plots, dashed black line has slope of 1. Dark and light green shaded areas indicate +-0.5 and +-1.0 pKa unit error regions, respectively.
  - `StatisticsTables/` - This directory contains machine readable copies of Statistics Table and bootstrap distributions of performance statistics.  
    - `statistics_bootstrap_distributions.pdf` - Violin plots showing bootstrap distributions of performance statistics of each method. Each method is labelled by submission ID.

## Remarks
- Submissions with submission IDs nb001, nb002, nb003, nb004, nb005 and nb005 include non-blind corrections to pKa predictions of only SM22 molecule. pKas of the rest of the molecules in these submissions were blindly predicted before experimental data was released.

- pKa predictions of Epik-sequencial method (submission ID: nb007) were not blind. They were submitted after the submission deadline to be used as a reference method.
