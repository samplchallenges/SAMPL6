# SAMPL6 Detailed Host-Guest Instructions

## Due date
Your predictions must be uploaded on the D3R SAMPL5 web-page by January 12, 2018.
The experimental results will be released immediately after the challenge closes.
You must use the provided template to upload your predictions to the [SAMPL website](https://drugdesigndata.org/about/sampl6).
(Currently the template is not yet available, but we expect to provide it by mid-September; please contact us if it is not available when you want it.)

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
- OctaAcidsAndGuests: OctaAcid challenge files, including:
  - mol2 and sdf files for all of the guests, with two for each guest (one docked to OA, the other docked to TEMOA). Nomenclature: OA-GN.sdf and .mol2, TEMOA-GN.sdf and .mol2, where N runs from 0 through 7
  - mol2, sdf and PDB files for OA and TEMOA
  - README.md briefly describing provenance of files
- GenerateInputs.ipynb: Jupyter notebook using OpenEye toolkits to generate molecules from SMILES strings and dock to hosts
- Gibb_SAMPL6_guests.cdxml: Gibb’s source ChemDraw file
- Gibb_SAMPL6_guests.smi: SMILES of guests from ChemDraw
- Isaacs_SAMPL6_guests.cdxml: Isaacs’ source ChemDraw file
- Isaacs_SAMPL6_guests.smi: SMILES of guests from Isaacs, after some additional curation as described in README.md


## Uploading your predictions
D3R is currently outfitting the SAMPL6 page with the ability to accept your uploaded predictions, and we are also preparing template files to be used in submitting predictions.
As soon as these are ready, you may upload your predictions.
 If you want to upload more than one set of predictions, generated by different methods, each set must be uploaded as a separate file.
Please use the template provided, as the predictions will be parsed and analyzed with automated scripts.
A complete set of predictions constitutes predicted binding free energies for all host-guest pairs, with predicted numerical uncertainties.
We also encourage predictions of the binding enthalpies, and of the binding free energies and enthalpies of the CB8 bonus cases. 
Incomplete submissions - such as for a subset of compounds - will also be accepted, but will not necessarily be evaluated together with the rest of the submissions.
However, we would emphasize that omission of enthalpies and/or bonus cases will not cause a submission to be regarded as incomplete. 

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

