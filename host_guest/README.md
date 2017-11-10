# SAMPL6 host-guest files


## Manifest
- `Gibb_SAMPL6_guests.cdxml`: ChemDraw file provided by Bruce Gibb Aug. 18, 2017
- `Gibb_SAMPL6_guests.smi`: Isomeric SMILES format file for the same guests, created by D. Mobley Aug. 18, 2017. Format is "smiles name" on each line. Created by using selecting each structure in ChemDraw, choosing "copy as -> smiles", and then pasting into (a) Picto to check, and (b) this text file to save. Names typed by hand.
- `Isaacs_SAMPL6_guests.cdx`: ChemDraw file provided by Lyle Isaacs Aug. 18, 2017
- `Isaacs_SAMPL6_guests.smi`: Isomeric SMILES format file for the same guests, created by D. Mobley Aug. 18, 2017; the SMILES for palonosetron was updated 8/22/17 after careful treatment of stereochemistry by Isaacs and Mobley. Format is "smiles name" on each line. Created by selecting each structure in ChemDraw, choosing "copy as -> smiles", and then pasting into (a) Picto to check, and (b) this text file to save. Names typed by hand. NOTE: Salts were not copied, only the guest molecules (tartrate and oxalate have been verified not to bind). N2: Gallamine triethiodate SMILES was taken from Wikipedia and manually verified to match the intended structure. N3: Chirality of bicyclic ring system bridgehead carbon did not propagate properly from ChemDraw to SMILES; manually re-added in MarvinSketch and verified using Picto before pasting SMILES here.
- `OctaAcidsAndGuests/`: Directory containing octa acid (OA and TEMOA) input structures and guest structure files
- `CB8AndGuests/`: Directory containing CB8 input structure and guest structure files
- `GenerateInputs.ipynb`: Jupyter notebook using the OpenEye toolkits to generate molecule structure files from the other inputs noted above. NOTE: This is provided for informational purposes; the output files are already available here so it is unnecessary for you to use this notebook.
- `SAMPLing/`: Equilibrated system files for the SAMPLing challenge (see below for further details).

## SAMPLing challenge files

Equilibrated systems are provided for OA-G3 (5-hexenoic acid), OA-G6 (4-methylpentanoic acid) and CB8-G3 (quinine). 5 different initial configurations are given for each complex. Files are available in Amber (`prmtop`/`rst7`), Gromacs (`top`/`gro`), OpenMM (`xml`) and PDB formats. Solvated system files are available both for the host-guest complex (e.g. `complex`) and the guest alone (e.g. `solvent`).

The `SAMPLing/` folder includes an atom map in JSON format for relative free energy calculations. The ligand atoms of OA-G3 that match the ligand atoms in OA-G6 are given for the systems in complex and in solvent. The file has this format
```json
"complex":
   "unique_atoms_G3": [184, 185, 187, 192, 193, 194, 195, 196]
   "unique_atoms_G6": [185, 186, 189, 192, 193, 194, 195, 196, 197, 202]
   "atom_map_G3_to_G6":
      "197": 198
      "198": 199
      "199": 201
      ...
"solvent":
   ...
```
where `unique_atoms_G3` is a list of atom indices that do not match any G6 atom, and `atom_map_G3_to_G6` matches atoms of G3 to those of G6 by atom index. All indices are 0-based. This map can be used with any of the 5 replicates `OA-G3-X` and `OA-G6-X`.

### Preparation
All the host-guest system files in the `SAMPLing/` directory were prepared using the protocol below.
- Protonation states and initial starting configurations were taken as given by the original `mol2` files in the `OctaAcidsAndGuests/` and `CB8AndGuests/` directories.
- 5 docked complexes were generated with OpenEye `2017.6.1`.
- Hosts and guests were both parametrized with GAFF v1.8 and antechamber. AM1-BCC charges were generated using OpenEye's QUACPAC toolkit through `openmoltools 0.8.1`.
- The systems were solvated in a 12A buffer of TIP3P water molecules using tleap. ParmEd `2.7.3` was used to remove some of the water molecules from the OA complexes to reduce them to have the same number of waters.
- The systems' net charge was neutralized with Na+ and Cl- ions. More Na+ and Cl- ions were added to reach the ionic strength of 60mM for OA/TEMOA systems and 150mM for CB8 to simulate effect of the 10mM and 25mM sodium phosphate buffer used in their respective experiments.
- The system was minimized with the L-BFGS optimization algorithm and equilibrated by running 1ns of Langevin dynamics (BAOAB splitting, 1fs time step) at 298.15K with a Monte Carlo barostat set at 1atm using `OpenMM 7.1.1`. PME was used for long-range electrostatic interactions with a cutoff of 10A. VdW interactions used the same 10A cutoff and a switching distance of 9A.
- After the equilibration, the `System` was serialized into the OpenMM `xml` format. The `rst7` file was generated during the equilibration using the `RestartReporter` in the `parmed.openmm` module. The AMBER `prmtop` and `rst7` files were then converted to GROMACS `top`/`gro` and PDB formats by ParmEd and MDTraj `1.9.1` respectively.
