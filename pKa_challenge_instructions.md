# SAMPL6 pKa Challenge Instructions

Challenge timeframe: Oct 24, 2017 to Jan 10, 2018  

This challenge consists of predicting microscopic and macroscopic pKas of 24 small organic molecules. These fragment-like small molecules are selected for their similarity to kinase inhibitors and for experimental tractability. Our aim is to evaluate how well current pKa prediction methods perform with drug fragment-like molecules through blind predictions.

Three formats of pKa prediction results will be evaluated:
1. microscopic pKa values and related microstates
2. microstate populations as a function of pH
3. macroscopic pKa values

The following subsections describe the molecules included in this challenge, the experimental conditions and measurements, the quantities to be predicted, and how prediction results must be submitted.

## Due Date

Your predictions must be uploaded on the D3R SAMPL6 web-page by January 10, 2018. The experimental results will be released immediately after the challenge closes. You must use the provided templates to upload your predictions to the [SAMPL website](https://drugdesigndata.org/about/sampl6). Additional information on using these templates is provided below.

## Small molecules included in the pKa challenge
![pKa_challenge_small_molecules](images/pKa_challenge_small_molecules.jpg)
A list of SAMPL6 pKa challenge small molecules canonical isomeric SMILES and molecule IDs can be found here: [/physical_properties/pKa/molecule_ID_and_SMILES.csv](physical_properties/pKa/molecule_ID_and_SMILES.csv). Counterions, where present in solid formulations, were included in canonical isomeric SMILES for the sake of completeness, although no significal effect is expected from the presence of chloride counterions as experiments were conducted in ionic-strength adjusted medium with KCl.

## Experimental details
pKa measurements were collected using spectrophotometric pKa measurements with Sirius T3 instrument by Mehtap Isik from the Chodera Lab at MSKCC with the support of Merck Rahway, Preformulation Department. 

Small molecules were purchased in powder form. 10 mg/ml DMSO solutions were prepared and used as stock solutions for preparation of samples, where 1-5 uL of 10 mg/ml DMSO stock solution is diluted in 1.5 mL ionic-strength adjusted water (0.15 M KCl). pH titrations with acid (0.5 M HCl, 0.15 M KCl) and base (0.5 M KOH, 0.15 M KCl) and cosolvent addition (80% Methanol, 0.15 M KCl) were performed in an automated fashion with Sirius T3.

[UV-metric pKa measurement protocol of Sirius T3](http://www.sirius-analytical.com/science/pka) evaluates pKas based on change of absorbance in the UV region of the spectrum while a pH titration is being performed between pH 1.8 and 12.2. UV-absorbance spectra of small molecule chromophores can change as a function of pH,  based on the protonation state change of groups close to the chromophore. Change in multiple wavelengths as a function of pH is fitted to a model to determine pKas. Number of pKas in the model is selected based on better fit between experimental and modeled microstates curves. 

This method is capable of measuring pKas between 2 and 12, when protonatable groups are at a maximum distance of 4-5 atoms away from chromophores, so that a change of protonation state affects the absorbance of UV-chromophore. We have selected compounds where titratable groups are close to potential chromophores (conjugated rings), but it is possible to completely miss detection of a pKa value if a titratable group is not in interaction with a UV-chromophore.  

pKa measurements of soluble compounds were performed in ionic-strength adjusted water with 0.15 M KCl. Indications for precipitation was observed by visual inspection of samples and also by inspection of UV-spectra plots as a function of pH at 500-600 nm region of the spectra. For compounds with insufficient solubility, cosolvent protocol is used where 3 UV-metric pKa measurements were done at different cosolvent:water ratios (typically 30%, 40% and 50% methanol) and  Yasuda-Shedlovsky extrapolation method is used to estimate pKa value at 0% cosolvent.

Three replicate pKa measurements were made for all compounds at room temperature (25°C). Multiwavelength absorbance analysis of Sirius T3 allows very good resolution pKas but essentially this method measures macroscopic pKas. Microscopic pKas with with very close pKa values and overlapping changes in absorbance spectra could be measured as one macroscopic pKa value.

## Instructions and Submission Templates
Three types of predictions will be accepted. Participants are encouraged to submit their results in all or multiple submission types as it fits to their prediction methods.

#### Type I - microscopic pKas and related microstates
Predicting microscopic pKa values and related microstates structures. Different protonation states and tautomers constitute different microstates. 
- Fill the `microscopic_pKa_and_microstates.csv` template for all molecules.
- For each molecule, report as many microscopic pKas based on your best guest.
- Record the pair of microstates IDs of microstate structures pairs (protonated HA and deprotonated A) associated for each microscopic pKa. To determine the microstate ID for your predicted structure check the spreadsheets in [SAMPL6/physical_properties/pKa/microstates](SAMPL6/physical_properties/pKa/microstates).
- If your predicted structure is not included in the list, contact us to [make a request for new microstate](mehtap.isik@choderalab.org). See more details are in missing microstates section.
- Report microscopic pKa values with 2 decimals (e.g. 10.71).
- Reporting standard error of the mean (SEM) is optional and should be reported with 2 decimals. (e.g. 1.02).

#### Type II - microstate populations as a function of pH
Predicting fractional population of microstates between pH interval 2 to 12 with 0.1 pH increments, based on microscopic pKa predictions.
- Fill the `microstate_fractional_populations.csv` template file for all molecules and microstates you have predictions for.
- To determine the microstate ID for your predicted microstate populations check the spreadsheets in [SAMPL6/physical_properties/pKa/microstates](SAMPL6/physical_properties/pKa/microstates).
- If your predicted structure is not included in the list, contact us to [make a request for new microstate](mehtap.isik@choderalab.org). See more details are in missing microstates section.
- Report natural logaritm of fractional microstate populations with scientific notation and 3 decimals of precision, for each pH.  
e.g. For a molecule with only two possible microstates A and B `ln(fractional microstate population) = ln(N_A/(N_A+N_B))` where `N_A` and `N_B` represent number of molecules at microstate A and B.   
At a pH where 90.0% of the molecules are in microstate B and 10.0% of molecules are in state A  `ln(fractional microstate A population) = ln(0.100/(0.100+0.900)) = -2.30E+00`.  
- Do not report SEM in this submission type.
- For pH value or microstates which you don't have any estimates for, leave the csv table cell empty. 

#### Type III - macroscopic pKas
Predicting the value of  macroscopic pKas based between 2 and 12.
- Fill one `macroscopic_pKa.csv` template file for all predicted molecules.
- Report as many macroscopic pKas based on your predictions.
- Report pKa values with 2 decimals (e.g. 10.71).
- Reporting standard error of the mean (SEM) is optional and should be reproted with 2 decimals. (e.g. 1.02).

## A warning about enumerated microstates and requesting missing microstates
A list of microstates and microstate IDs were generated for each molecule to aid parsing of submissions. Enumerated list of microstates was not created with the intend to guide computational predictions. It is possible that some relevant microscopic species are missing from these lists. If your predicted structure is not already included in the microstates list, contact Mehtap Isik (mehtap.isik@choderalab.org). Please send us 2D structure depiction and canonical isomeric SMILES of your proposed microstate. We will assign a unique microstate ID and include it in the analysis.  

Please do not create a microstate ID yourself. It is important that challenge organizers assign unique microstate IDs and keep track. 

## Evaluation strategy for computational pKa predictions
Macroscopic pKa predictions will be directly compared to experimental pKa measurements. Microscopic pKas and predicted microstates will be evaluated by comparison to each other and how well they recapitulate experimental macroscopic pKas, as well as pKa predictions of ACD which uses a data-trained pKa predictions method.


