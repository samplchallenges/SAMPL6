# SAMPL6 host-guest files


## Manifest
- `Gibb_SAMPL6_guests.cdxml`: ChemDraw file provided by Bruce Gibb Aug. 18, 2017
- `Gibb_SAMPL6_guests.smi`: Isomeric SMILES format file for the same guests, created by D. Mobley Aug. 18, 2017. Format is "smiles name" on each line. Created by using selecting each structure in ChemDraw, choosing "copy as -> smiles", and then pasting into (a) Picto to check, and (b) this text file to save. Names typed by hand.
- `Isaacs_SAMPL6_guests.cdx`: ChemDraw file provided by Lyle Isaacs Aug. 18, 2017
- `Isaacs_SAMPL6_guests.smi`: Isomeric SMILES format file for the same guests, created by D. Mobley Aug. 18, 2017; the SMILES for palonosetron was updated 8/22/17 after careful treatment of stereochemistry by Isaacs and Mobley. Format is "smiles name" on each line. Created by selecting each structure in ChemDraw, choosing "copy as -> smiles", and then pasting into (a) Picto to check, and (b) this text file to save. Names typed by hand. NOTE: Salts were not copied, only the guest molecules (tartrate and oxalate have been verified not to bind). N2: Gallamine triethiodate SMILES was taken from Wikipedia and manually verified to match the intended structure. N3: Chirality of bicyclic ring system bridgehead carbon did not propagate properly from ChemDraw to SMILES; manually re-added in MarvinSketch and verified using Picto before pasting SMILES here.
- `OctaAcidsAndGuests/`: Directory containing octa acid (OA and TEMOA) input structures and guest structure files
- `CB8AndGuests/`: Directory containing CB8 input structure and guest structure files
- `GenerateInputs.ipynb`: Jupyter notebook using the OpenEye toolkits to generate molecule structure files from the other inputs noted above. NOTE: This is provided for informational purposes; the output files are already available here so it is unnecessary for you to use this notebook.
- `Reference/`: System files used to run the reference calculations for the main SAMPL challenge.
- `SAMPLing/`: Equilibrated system files for the SAMPLing challenge. See [`SAMPLing_instructions.md`](../SAMPLing_instructions.md#files-description) for a detailed description of the files.
