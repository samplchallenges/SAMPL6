## Manifest

- `experimental_measurements.pdf`: Summary table of the experimental data after error propagation.
- `experimental_measurements.csv`: Summary table of the experimental data after error propagation in CSV format.
- `experimental_measurements.json`: Summary table of the experimental data after error propagation in JSON format.
- `SAMPL6_Gibb_data.pdf`: Data provided by the Gibb group.
- `SAMPL6_Isaacs_data.pdf`: Data provided by the Isaacs group.
- `generate_tables.py`: Script used to perform error propagation and create the `experimental_measurements.X` files
based on the data provided by the Gibb and Isaacs groups.

## Notes on error propagation

`SAMPL6_Gibb_data.pdf` and `SAMPL6_Isaacs_data.pdf` report the relative uncertainties obtained from the ITC fitting procedure.
The upper bound (1%) was used for errors that were reported to be < 1%. We also included a 3% relative uncertainty in the
titrant concentration assuming the stoichiometry coefficient to be fitted to the ITC data [1]. This is exact only for the
OA/TEMOA sets (with the exception of OA-G5, TEMOA-G5, and TEMOA G7). For the other guests, we may expand the error
analysis to include also the effect of the uncertainties in titrand concentration and cell volume.

## References

[1] Boyce SE, Tellinghuisen J, Chodera JD. Avoiding accuracy-limiting pitfalls in the study of protein-ligand
interactions with isothermal titration calorimetry. bioRxiv. 2015 Jan (doi:10.1101/023796).
