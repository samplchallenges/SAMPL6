# SAMPL6 Detailed Host-Guest Instructions

## Due date
Your predictions must be uploaded on the D3R SAMPL5 web-page by January 12, 2018.
The experimental results will be released immediately after the challenge closes.
You must use the provided template to upload your predictions to the [SAMPL website](https://drugdesigndata.org/about/sampl6).
A draft template (example) is available for each of the host-guest systems; it is placed in the relevant directory and labeled by the host name; ([host_guest/CB8AndGuests/CB8.txt](host_guest/CB8AndGuests/CB8.txt), [host_guest/OctaAcidsAndGuests/OA.txt](host_guest/OctaAcidsAndGuests/OA.txt), and [host_guest/OctaAcidsAndGuests/TEMOA.txt](host_guest/OctaAcidsAndGuests/TEMOA.txt)).
Additional information on using these templates is provided below.

## Anonymous versus public participation
When you upload your submission, you will have the option of having it treated anonymously.
Anonymous submission means that we may report on your predictions and methods, but not your identity.
Public participations means we may also say who you are. 
Please note that, although we will work to protect the identity of anonymous participants, we cannot make any guarantees.
You may use the D3R website to change your submission’s anonymous/public status until the challenge has closed.
However, after the challenge has closed, you may not change its anonymous/public status.  

## Submission of multiple predictions
Some participants use SAMPL to help evaluate various computational methods.
To accommodate this, multiple prediction set from a single research group or company are allowed. 

## SAMPL6 workshop February 22-23, 2018
Participants are invited to share and discuss their results, as well as the D3R and SAMPL projects more broadly, at the second in-person D3R and SAMPL workshop, which is scheduled for February 22-23, 2018, at UC San Diego, La Jolla, CA.
Note that the workshop immediately follows the Biophysical Society National Meeting in San Francisco.

## Molecular systems
Three aqueous host-guest series are provided. 
The first two comprise hosts OA and TEMOA, each with the same set of 8 guests. (These hosts also appeared in SAMPL5, where they were named OctaAcidH (OAH) and OctaAcidMe (OAMe), respectively.
The name change here is to improve consistency with the experimental literature.)
Measured binding free energies are available for all 16 cases, and binding enthalpies are also expected to be available. 
The third set comprises host cucurbit[7]uril (CB7), with a different set of 14 guests; of these, 11 form the main challenge, and the remaining 3 are regarded as a bonus challenge, due to the special complications they pose.
(See [detailed host-guest description](host_guest_description.md).)

## Computational methods
You may use any method(s) you like to generate your predictions; e.g., MD with implicit or explicit solvent; quantum methods; docking; etc. 

## Files provided
All hosts and guests are provided in mol2, PDB and SDfile formats, prepared by David Mobley (UCI).
A Jupyter notebook which was used in generating these 3D structures from SMILES strings is provided in the data download package if you would like to examine exactly what was done.
A README.md file in `host_guest` gives more details about the files provided.
Disclaimers apply: There is no guarantee or representation that the protonation and tautomer states provided are optimal.
It is also possible that the protonation state of a guest or host will change on binding.
Guest files are provided using the same frame of reference as the host files after some attempt at docking guests into the hosts (code available from Mobley), but there is no expectation that the provided structures show the actual binding mode of the guests.

**Some of the provided files deserve special caveats/particular attention**:
- *CB8 bonus challenge*:
   - For oxaliplatin, CB8-G13, the provided 3D geometry is likely incorrect due to issues with the platinum. If you choose to predict this compound, you will likely need to give this compound particular attention.
   - For aricept, CB8-G12, the nitrogen was protonated by our tools; the proton was arbitrarily placed so that the substituent is equatorial. It could also have been axial. Participants wishing to predict this compound will need to determine whether a proton is expected to be present and, if so, where. 

The organizers plan to carry out a complete set of binding free energy and enthalpy calculations; these will form the basis of a separate host-guest ``SAMPLing'' challenge where the goal is to reproduce computed reference values as efficiently as possible.
Details and input files (converted to formats for a wide variety of simulation packages) will be posted separately as soon as they are available, likely around mid-September.  

The `host_guest` directory/data is organized as follows (with more detail in the README.md there):
- `README.md`: README file describing contents of directory and their organization/how they were generated
- CB8AndGuests: CB8 challenge files, including:
  - mol2 and sdf files for all of the guests (nomenclature: CB8-GN.sdf and .mol2, where N runs from 0 to 13), using the same frame of reference as the host 
  - mol2, sdf and PDB files for CB8
  - README.md briefly describing provenance of files
  - CB8.txt: Template (example) file for submissions.
- OctaAcidsAndGuests: OctaAcid challenge files, including:
  - mol2 and sdf files for all of the guests, with two for each guest (one docked to OA, the other docked to TEMOA). Nomenclature: OA-GN.sdf and .mol2, TEMOA-GN.sdf and .mol2, where N runs from 0 through 7
  - mol2, sdf and PDB files for OA and TEMOA
  - README.md briefly describing provenance of files
  - OA.txt, TEMOA.txt: Template (example) files for submissions.
- GenerateInputs.ipynb: Jupyter notebook using OpenEye toolkits to generate molecules from SMILES strings and dock to hosts
- Gibb_SAMPL6_guests.cdxml: Gibb’s source ChemDraw file
- Gibb_SAMPL6_guests.smi: SMILES of guests from ChemDraw
- Isaacs_SAMPL6_guests.cdxml: Isaacs’ source ChemDraw file
- Isaacs_SAMPL6_guests.smi: SMILES of guests from Isaacs, after some additional curation as described in README.md


## Uploading your predictions
D3R is currently outfitting the SAMPL6 page with the ability to accept your uploaded predictions.
As soon as this is ready, you may upload your predictions.
 If you want to upload more than one set of predictions, generated by different methods, each set must be uploaded as a separate file.
Please use the template provided, as the predictions will be parsed and analyzed with automated scripts.
A complete set of predictions constitutes predicted binding free energies for all host-guest pairs, with predicted numerical uncertainties.
We also encourage predictions of the binding enthalpies, and of the binding free energies and enthalpies of the CB8 bonus cases. 
Incomplete submissions - such as for a subset of compounds - will also be accepted, but will not necessarily be evaluated together with the rest of the submissions.
However, we would emphasize that omission of enthalpies and/or bonus cases will not cause a submission to be regarded as incomplete. 

Names of the prediction files must begin with the name of the host molecule for which it contains predictions (i.e., OA, TEMOA or CB8, case-independently), and must end with an integer indicating which of your predictions for this host it contains.
For example, if you want to submit one prediction file for CB8, you might name it CB8-myname-1.csv, where myname is arbitrary text of your choice. If you submit two prediction files for CB8, you might name them CB8-myname-1.txt and CB8-myname-2.txt

The file will be machine parsed, so correct formatting is essential.

Lines beginning with a hash-tag (#) may be included as comments. These and blank lines will be ignored.

The file must contain the following four components in the following order: your predictions, a name for your computational protocol, a list of the major software packages used, and a long-form methods description. Each of these components must begin with a line containing only the corresponding keyword: Predictions:, Name:, Software:, and Method:, as illustrated in the example files. 
More detailed instructions will be provided on the challenge submission site on the D3R website.


## Reference calculations

Similarly to previous editions of the SAMPL challenge, we will run reference free energy calculations for the host-guest systems included in the primary challenge (i.e., not for the 3 CB8 ligands forming the "bonus" challenge). These are solely intended to provide a point of comparison across all targets for other methods. The simulations will be run using YANK (more details about the method will be added to this section in the near future). We have shared the input files that will be used for the calculations in [`host_guest/Reference/`](host_guest/Reference/) (see below for a description of the setup protocol). Note that, contrary to the SAMPLing challenge, we do _not_ expect participants to use these files for their simulations. The files are available for transparency or in case anyone would like to run the same systems to compare performance.

## Reference files description
Files are available in Amber (`prmtop`/`rst7`), Gromacs (`top`/`gro`), OpenMM (`XML`) and PDB formats. Solvated system files are available both for the host-guest complex (e.g. `complex`) and the guest alone in water (e.g. `solvent`).

All the host-guest system files in [`host_guest/Reference/`](host_guest/Reference/) were prepared using the following protocol:
- We used the most likely protonation states as predicted by Epik `4.0013` from the Schrodinger toolkit at experimental pH. These are identical to those given in the `mol2` files in `host_guest/OctaAcidsAndGuests/` and `host_guest/CB8AndGuests/`.
- The initial starting configurations are taken from the `mol2` files in the `host_guest/OctaAcidsAndGuests/` and `host_guest/CB8AndGuests/` directories. These were determined by docking through the OpenEye toolkit.
- Hosts and guests were both parametrized with GAFF v1.8 and antechamber. AM1-BCC charges were generated using OpenEye's QUACPAC toolkit through `openmoltools 0.8.1`.
- The systems were solvated in a 12A buffer of TIP4P-Ew water molecules using tleap.
- The systems' net charge was neutralized with Na+ and Cl- ions. More Na+ and Cl- ions were added to reach the ionic strength of 60mM for OA/TEMOA systems and 150mM for CB8 to simulate the effect of the 10mM and 25mM sodium phosphate buffer used in their respective experiments.
- `prmtop`, `inpcrd` and `pdb` files for each system were generated by tleap. The `prmtop`/`inpcrd` files were converted to GROMACS (`top`/`gro`) files using ParmEd `2.7.3`.
- The serialized OpenMM system (`XML`) was generated by OpenMM `7.1.1`. It includes a `MonteCarloBarostat` set at 298.15K and 1atm, and it is configured to use PME, a nonbonded cutoff of 11A, a switching distance of 10A, and constraints on all hydrogen bonds lengths.


## Problems

If you notice any issues with any of these files, please contact us via the GitHub issue tracker.
You are also strongly advised to both sign up for the SAMPL e-mail list via the D3R site and sign up for notifications on this GitHub repository in case we have updates to any of the files.

## Pending items, error reports, questions

We will attempt to notify you (and update this repository) during the challenge about the following pending items:
- Updates regarding completion of the experimental studies
- Template for uploading your predictions
- Workshop details
- Availability of files for the SAMPLing challenge

Please feel free to contact us if you notice any errors in the information provided or have questions about SAMPL6; please use the issue tracker connected with this repository, or use our e-mail: samplchallenge@gmail.com

