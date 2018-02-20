# Analysis of the submissions for the host-guest challenge.

The current version of the analysis report correlation plots and 5 performance statistics (RMSE, MAE, ME, R^2 and linear
regression slope) for all the submissions. We report 95%-percentile bootstrap confidence interval of all the statistics.
We don't currently include the reported free energy errors into our bootstrap analysis. This will be added in a future
version.

## Manifest

- `ExperimentalMeasurements/`: Experimental data used for the analysis.
- `CB8/`: Statistics and correlation plots for CB8 predictions.
- `OA/`: Statistics and correlation plots for OA predictions.
- `TEMOA/`: Statistics and correlation plots for TEMOA predictions.
- `OA-TEMOA/`: Merged statistics and correlation plots for methods making predictions of both OA and TEMOA guests.
- `CB8-OA-TEMOA/`: Merged statistics and correlation plots for methods making predictions of CB8, OA, and TEMOA guests.
- `MoleculesStatistics/`: Statistics across methods grouped by host-guest system.
- `SAMPLing/`: Preliminary analysis and results of the SAMPLing component of the host-guest challenge.
- `Scripts/`: Python scripts used to analyze the submissions.
- `Submissions`: Participant submissions, including detailed method descriptions.
