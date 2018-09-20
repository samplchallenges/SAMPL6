# Analysis by molecule

Statistics across methods grouped by host-guest system. The statistics were generated using only the 10 methods scoring
the lowest RMSE either in the combined OA/TEMOA or CB8-NOBONUS guest set.

## Manifest
- `StatisticsTables/`: Tables reporting RMSE, mean absolute error, mean error, Pearson coefficient of determination,
linear regression slope, and Kendall's tau coefficient for all submissions with 95%-percentile bootstrap confidence
interval in CSV, JSON and PDF format.
- `molecules_error.pdf`: Violin plots of free energy predictions by host-guest system.
