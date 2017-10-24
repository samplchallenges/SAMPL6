# SAMPL6 pKa Challenge Instructions

Challenge timeframe: Oct 24, 2017 to Jan 10, 2018  

This challenge consists of predicting microscopic and macroscopic pKas of 24 small organic molecules. 
These fragment-like small molecules are selected for their similarity to kinase inhibitors and for experimental tractability. 
Our aim is to evaluate how well current pKa prediction methods perform with drug fragment-like molecules through blind predictions.

Three formats of pKa prediction results will be evaluated:
1. microscopic pKa values and related microstates
2. microstate populations as a function of pH
3. macroscopic pKa values

The following subsections describe the molecules included in this challenge, the experimental conditions and measurements, the quantities to be predicted, and how prediction results must be submitted.

## 24 small molecules included in the pKa challenge
![pKa_challenge_small_molecules](images/pKa_challenge_small_molecules.jpg)
A list of SAMPL6 pKa challenge small molecules canonical isomeric SMILES and molecule IDs can be found here: [/physical_properties/pKa/molecule_ID_and_SMILES.csv](physical_properties/pKa/molecule_ID_and_SMILES.csv). 
Counterions, where present in solid formulations (see experimental details below), were included in canonical isomeric SMILES for the sake of completeness, although no significant effect is expected from the presence of chloride counterions as experiments were conducted in ionic-strength adjusted medium with KCl.

## Experimental details
pKa measurements were collected using spectrophotometric pKa measurements with Sirius T3 instrument by Mehtap Isik from the Chodera Lab at MSKCC with the support of Merck Rahway, Preformulation Department, especially Dorothy Levorse, Timothy Rhodes and Brad Sherborne. 

Small molecules were purchased in powder form. 
10 mg/ml DMSO solutions were prepared and used as stock solutions for preparation of samples, where 1-5 uL of 10 mg/ml DMSO stock solution is diluted in 1.5 mL ionic-strength adjusted water (0.15 M KCl). 
pH titrations with acid (0.5 M HCl, 0.15 M KCl) and base (0.5 M KOH, 0.15 M KCl) and cosolvent addition (80% methanol, 0.15 M KCl) were performed in an automated fashion with with a [Sirius T3 instrument (Sirius Analytical)](http://www.sirius-analytical.com/products/t3).

[The UV-metric pKa measurement protocol of Sirius T3](http://www.sirius-analytical.com/science/pka) measures the change in multiwavelength absorbance in the 250-450 nm UV region of the spectrum while the pH is titrated between pH 1.8 and 12.2 to evaluate pKas. 
Protonation state change of titratable sites near chromophores will modulate the UV absorbance spectra of these chromophores, allowing populations of distinct UV-active species to be resolved as a function of pH. 
To do this, basis spectra are identified and populations extracted via analysis of the pH-dependent multi-wavelength absorbance.
The number of pKas is determined based on the quality of fit between experimental and modeled microstate pH-dependent populations. 

This method is capable of measuring pKas between 2 and 12, when protonatable groups are at a maximum distance of 4-5 atoms away from chromophores, so that a change of protonation state affects the absorbance of UV-chromophore. 
We have selected compounds where titratable groups are close to potential chromophores (conjugated rings), but it is possible to completely miss detection of a pKa value if a titratable group is not in interaction with a UV-chromophore.  

pKa measurements of soluble compounds were performed in ionic-strength adjusted water with 0.15 M KCl. 
Indications for precipitation was observed by visual inspection of samples and also by inspection of UV-spectra plots as a function of pH at 500-600 nm region of the spectra. 
For compounds with insufficient solubility, cosolvent protocol is used where 3 UV-metric pKa measurements were done at different cosolvent:water ratios (typically 30%, 40% and 50% methanol) and  Yasuda-Shedlovsky extrapolation method is used to estimate pKa value at 0% cosolvent.

Three replicate pKa measurements were made for all compounds at room temperature (25°C). 
Replicate measurements were set up from the same compound stock solutions (~10 mg/ml in DMSO) and independent aliquotes were taken to prepare samples for Sirius T3 titration.
Multiwavelength absorbance analysis of Sirius T3 allows for very good resolution of pKas but essentially this method measures macroscopic pKas. 
Microscopic pKas with very close pKa values and overlapping changes in absorbance spectra could be measured as one macroscopic pKa value.

## Due Date

Your predictions must be uploaded on the D3R SAMPL6 web-page by January 10, 2018. 
The experimental results will be released immediately after the challenge closes. 
You must use the provided templates to upload your predictions to the [SAMPL website](https://drugdesigndata.org/about/sampl6). Additional information on using these templates is provided below.

## Computational prediction methods

You may use any method(s) you like to generate your predictions; e.g., molecular mechanics or quantum mechanics based methods, QSPR, empirical pKa prediction tools etc.

## Instructions and Submission Templates
Three types of predictions will be accepted. 
Participants are encouraged to submit their results in all or multiple submission types if their prediction methods allow it.

#### Prediction Type I - microscopic pKas and related microstates
Predicting microscopic pKas and related microstate structures. 
Different protonation states and tautomers constitute different microstates. 
- Fill the `typeI_microscopic_pKas_and_microstates.csv` template for all molecules.
- For each molecule, report as many microscopic pKas as your method predicts.
- Record the pair of microstates IDs of microstate structures pairs (protonated HA and deprotonated A) associated with each microscopic pKa. To determine the microstate ID for your predicted structure, check the csv files and spreadsheets in [SAMPL6/physical_properties/pKa/microstates](SAMPL6/physical_properties/pKa/microstates) that list microscopic species.
- If your predicted structure is not included in the list, contact us to [make a request for new microstate](mehtap.isik@choderalab.org). See more details in the section below ("A warning about enumerated microstates and requesting the missing microstates").
- Report microscopic pKa values to two decimal places (e.g. 10.71).
- Reporting the standard error of the mean (SEM) is optional, but if reported, two decimal places should be provided (e.g. 1.02).

#### Prediction Type II - microstate populations as a function of pH
Predicting the fractional microstate populations between pH interval 2 to 12 with 0.1 pH increments.

- Fill the `typeII_microstate_fractional_populations.csv` template file for all molecules and microstates you have predictions for.
- For each molecule, report as many microstates as your method predicts.
- To determine the microstate ID for your predicted microstate populations, check the csv files and spreadsheets in [SAMPL6/physical_properties/pKa/microstates](SAMPL6/physical_properties/pKa/microstates) that list microscopic species.
- If your predicted structure is not included in the list, contact us to [make a request for new microstate](mehtap.isik@choderalab.org). See more details in the section below ("A warning about enumerated microstates and requesting the missing microstates").
- For each pH, report the *natural logarithm* of the fractional micrstate populations in scientific notation with three decimals of precision (e.g., 1.02e-4).
e.g. For a molecule with only two possible microstates A and B `ln(fractional microstate population) = ln(N_A/(N_A+N_B))` where `N_A` and `N_B` represent percentage of microstate populations of A and B.   
At a pH where 90.0% of the molecules are in microstate B and 10.0% of molecules are in state A  `ln(fractional microstate A population) = ln(0.100/(0.100+0.900)) = -2.30E+00`.  
- If your estimate of `fractional microstate population` is 0, thus `ln(fractional microstate population) = ln(0)`, report as `-infinity`.
- Do not report SEM in this submission type.
- For pH value or microstates which you don't have any estimates for, leave the csv table cell empty. 

#### Prediction Type III - macroscopic pKas
Predicting the value of  macroscopic pKas based between 2 and 12.
- Fill one `typeIII_macroscopic_pKas.csv` template file for all predicted molecules.
- For each molecule, report as many macroscopic pKas as your method predicts.
- Report pKa values to two decimal places (e.g. 10.71).
- Reporting the standard error of the mean (SEM) is optional, but if reported, two decimal places should be provided (e.g. 1.02).

## A warning about enumerated microstates and requesting the missing microstates
A list of microstates and microstate IDs were generated for each molecule to aid parsing of submissions. 
Enumerated list of microstates was not created with the intend to guide computational predictions. 
It is possible that some relevant microscopic species are missing from these lists. 
If your predicted structure is not already included in the microstates list, contact Mehtap Isik (mehtap.isik@choderalab.org). 
Please send us 2D structure depiction and canonical isomeric SMILES of your proposed microstate. 
We will assign a unique microstate ID and include it in the analysis.  

Please do not create a microstate ID yourself. 
It is important that challenge organizers assign unique microstate IDs and keep track.  

## Submission of multiple predictions

Some participants use SAMPL to help evaluate various computational methods. 
To accommodate this, multiple prediction sets from a single research group or company are allowed, even for same type of predictions if they are made by different methods.

## Uploading your predictions

D3R is currently outfitting the SAMPL6 page with the ability to accept your uploaded predictions. 
As soon as this is ready, you may upload your predictions. 
If you want to upload more than one type predictions (see description of three pKa prediction types above) or different set of predictions (same type but generated by different methods), each must be uploaded as a separate file. 
Please use the template provided, as the predictions will be parsed and analyzed with automated scripts. 
Please include predictions related all molecules with same method and same submission type in one file. 

We encourage submitting predictions in all three formats and for all 24 molecules when possible. 
Incomplete submissions - such as for a subset of compounds - will also be accepted, but will not necessarily be evaluated together with the rest of the submissions. 
However, we would emphasize that omission of SEM estimates will not cause a submission to be regarded as incomplete, though we highly encourage including such estimates.

Names of the prediction files must begin with the name of the submission type (i.e., typeIII), and your name and must end with an integer indicating set of prediction. 
For example, if you want to submit one prediction file for type III prediction (macroscopic pKas), you would name it `typeIII-myname-1.csv`, where myname is arbitrary text of your choice. 
If you submit two prediction files of the same submission type, you would name them `typeIII-myname-1.txt` and `typeIII-myname-2.txt`.

The file will be machine parsed, so correct formatting is essential.  Files with the wrong format will not be accepted.

Lines beginning with a hash-tag (#) may be included as comments. 
These and blank lines will be ignored.

The file must contain the following four components in the following order: your predictions, a name for your computational protocol, a list of the major software packages used, and a long-form methods description. 
Each of these components must begin with a line containing only the corresponding keyword: `Predictions:`, `Name:`, `Software:`, and `Method:`, as illustrated in the example files. 
More detailed instructions will be provided on the challenge submission site on the D3R website.  

## Evaluation strategy for computational pKa predictions

Macroscopic pKa predictions will be directly compared to experimental pKa measurements. 
Microscopic pKas and predicted microstates will be evaluated by comparison to each other and how well they recapitulate experimental macroscopic pKas, as well as pKa predictions of ACD which uses a data-trained pKa predictions method.

## Anonymous versus public participation

When you upload your submission, you will have the option of having it treated anonymously. 
Anonymous submission means that we may report on your predictions and methods, but not your identity. 
Public participations means we may also say who you are. 
Please note that, although we will work to protect the identity of anonymous participants, we cannot make any guarantees. 
You may use the D3R website to change your submission’s anonymous/public status until the challenge has closed. 
However, after the challenge has closed, you may not change its anonymous/public status.

## SAMPL6 workshop February 22-23, 2018

Participants are invited to share and discuss their results, as well as the D3R and SAMPL projects more broadly, at the second in-person D3R and SAMPL workshop, which is scheduled for February 22-23, 2018, at UC San Diego, La Jolla, CA. 
Note that the workshop immediately follows the Biophysical Society National Meeting in San Francisco.

## Files provided
- `/physical_properties/pKa/molecule_ID_and_SMILES.csv` - CSV file that indicates SAMPL6 pKa challenge molecule IDs and canonical isomeric SMILES.
- `/physical_properties/pKa/microstates/` - This directory contains files that list of microstate IDs and canonical isomeric SMILES of microstates. Files are separated by molecule ID.
- `/physical_properties/pKa/submission_templates/` - Empty prediction submission template files are located here.
- `/physical_properties/pKa/example_submission_files/` - This directory contains example submission files filled with random values to illustrate expected format.

## Problems, questions and contact

If you notice any issues with any of these files, please contact us via the GitHub issue tracker. 
You are also strongly advised to both sign up for the SAMPL e-mail list via the D3R site and sign up for notifications on this GitHub repository in case we have updates.

Please feel free to contact us if you notice any errors in the information provided or have questions about SAMPL6; please use the issue tracker connected with this repository, or use our e-mail: samplchallenge@gmail.com. 
For specific questions about pKa challenge, use the following email: mehtap.isik@choderalab.org.
