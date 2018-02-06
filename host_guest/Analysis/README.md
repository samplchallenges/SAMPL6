# Analysis of the submissions for the host-guest challenge.

The current version of the analysis report correlation plots and 5 performance statistics (RMSE, MAE, ME, R^2 and linear
regression slope) for all the submissions. We report 95%-percentile bootstrap confidence interval of all the statistics.
We don't currently include the reported free energy errors into our bootstrap analysis. This will be added in a future
version.

## Manifest

- `ExperimentalMeasurements/`: Experimental data used for the analysis.
- 'CB8/': First analysis of the free energy calculations for the CB8 challenge.
- 'OA-TEMOA/': First analysis of the free energy calculations for the OA/TEMOA challenge.
- 'Scripts/': Python scripts used to analyze the submissions.
