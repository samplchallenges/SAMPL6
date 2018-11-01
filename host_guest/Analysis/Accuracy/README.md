# Analysis of the submissions for the host-guest challenge.

The analysis reports correlation plots and 5 performance statistics (RMSE, MAE, ME, R^2 and linear
regression slope) for all the submissions. We report 95%-percentile bootstrap confidence interval of all the statistics.

## Manifest

- `CB8/`: Statistics and correlation plots for CB8 predictions.
- `OA/`: Statistics and correlation plots for OA predictions.
- `TEMOA/`: Statistics and correlation plots for TEMOA predictions.
- `OA-TEMOA/`: Merged statistics and correlation plots for methods making predictions of both OA and TEMOA guests.
- `CB8-OA-TEMOA/`: Merged statistics and correlation plots for methods making predictions of CB8, OA, and TEMOA guests.
- `MoleculesStatistics/`: Statistics across methods grouped by host-guest system.
- `PaperImages`: Images generated for the overview paper. Most of these have been created automatically by
`../Scripts/analyze_hostguest.py`.
