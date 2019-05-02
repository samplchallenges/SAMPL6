# Molecular Statistics Analysis

This directory contains tables and barplots of molecular statistics analysis. 
This analysis was performed to indicate logP values of which molecules of SAMPL6 logP Challenge set were more difficult for methods to predict accurately.
Error statistics (MAE and RMSE) were calculated for each molecule averaging across all methods or for all methods within a method category.
95% confidence intervals were calculated via bootstrapping (10000 samples).

## Manifest
  
- `MAE_vs_molecule_ID_plot.pdf` - Barplot of MAE calculated for each molecule averaging over all prediction methods.
- `RMSE_vs_molecule_ID_plot.pdf` - Barplot of RMSE calculated for each molecule averaging over all prediction methods.
- `molecular_error_statistics.csv` - MAE and RMSE statistics calculated for each molecule averaging over all prediction methods. 
- `molecular_MAE_comparison_between_method_categories.pdf` - Barplot of MAE calculated for each method category for each molecule averaging over all predictions in that method category. Colors of bars indicate method categories.
- `Empirical/` - This directory contains table and barplots of molecular statistics analysis calculated only for methods in Empirical method category. 
  - `molecular_error_statistics_for_Empirical_methods.csv`
  - `MAE_vs_molecule_ID_plot.pdf`
  - `RMSE_vs_molecule_ID_plot.pdf`
- `Mixed/` - This directory contains table and barplots of molecular statistics analysis calculated only for methods in Mixed method category. 
  - `molecular_error_statistics_for_Mixed_methods.csv`
  - `MAE_vs_molecule_ID_plot.pdf`
  - `RMSE_vs_molecule_ID_plot.pdf`
- `Other/` - This directory contains table and barplots of molecular statistics analysis calculated only for methods in Other method category. 
  - `molecular_error_statistics_for_Other_methods.csv`
  - `MAE_vs_molecule_ID_plot.pdf`
  - `RMSE_vs_molecule_ID_plot.pdf`
- `Physical/` - This directory contains table and barplots of molecular statistics analysis calculated only for methods in Physical method category. 
  - `molecular_error_statistics_for_Physical_methods.csv`
  - `MAE_vs_molecule_ID_plot.pdf`
  - `RMSE_vs_molecule_ID_plot.pdf`    
