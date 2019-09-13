# SAMPL6 Physicochemical Propery Prediction Challenges

## pKa Challenge

This challenge consists of predicting microscopic and macroscopic pKas of 24 small organic molecules.
These fragment-like small molecules are selected for their similarity to kinase inhibitors and for experimental tractability.
Our aim is to evaluate how well current pKa prediction methods perform with drug fragment-like molecules through blind predictions.

For detailed instructions for pKa challenge: [SAMPL6/pKa_challenge_instructions.md](https://github.com/MobleyLab/SAMPL6/blob/pKa/pKa_challenge_instructions.md)

Challenge start date: Oct 25, 2017   
Challenge submission due: Jan 23, 2018  

Experimental pKa measurements were made available in this repository after the pKa challenge deadline and can be found here: [pKa/experimental_data/](pKa/experimental_data/)

Performance evaluation of pKa challenge can be found here: [/pKa/analysis/](/pKa/analysis/)

#### Reference calculations for pKa challenge

We collected reference pKa predictions with several established techniques (Epik, Jaguar, Chemicalize, and MoKa), though we were late completing these calculations and they were only added well after the challenge closed. 
These calculations are added with submission IDs beginning with "nb" to clearly distinguish them from formal submissions as "non-blind".
Reference calculations were submitted by Bas Rustenburg, Mehtap Isik, and Thomas Fox.
- `nb007`: TypeIII predictions with Epik Scan
- `nb008`: TypeI prediction with Epik Microscopic
- `nb009`: TypeII predictions with Epic Microscopic
- `nb010`: TypeIII predictions with Epik Microscopic
- `nb011`: TypeI predictions with Jaguar
- `nb012`: TypeII predictions with Jaguar
- `nb013`: TypeIII predictions with Jaguar
- `nb015`: TypeIII predictions with Chemicalize (ChemAxon)
- `nb016`: TypeI predictions with MoKa
- `nb017`: TypeIII predictions with MoKa


## log *P* Challenge

This challenge consists of predicting the octanol-water partition coefficients (log *P*) of 11 small molecules that resemble fragments of small molecule protein kinase inhibitors. Our aim is to evaluate how well current models can capture the transfer free energy of small molecules between different solvent environments through blind predictions.

For detailed instructions for log *P* challenge: [SAMPL6/logP_challenge_instructions.md](/logP_challenge_instructions.md)

Challenge start date: Nov 1, 2018  
Challenge submission due: Mar 22, 2019  

Experimental log *P* measurements were added to this repository after the log *P* challenge deadline and can be found here: [logP/experimental_data/](logP/experimental_data/)

Performance evaluation of log P challenge can be found here: [/logP/analysis/](/logP/analysis/)

#### Reference calculations for log *P* challenge

As in many previous SAMPL challenges, the Mobley group ran reference log *P* calculations with several established techniques (in this cased based on alchemical free energy calculations), though we were late completing these calculations and they were only added well after the challenge closed.
These calculations are added with submission IDs beginning with "REF" to clearly distinguish them from formal submissions, and analysis of statistics/performance overall with and without the reference calculations are being provided.
All reference calculations used a comparable protocol with Yank for solvation free energies, but varied the force field, water model and/or water content of octanol.
- `REF01`: YANK, GAFF force field, TIP3P-FB water, wet octanol
- `REF02`: YANK, GAFF force field, TIP3P water, wet octanol
- `REF03`: YANK, GAFF force field, OPC water, wet octanol
- `REF04`: YANK, SMIRNOFF force field, TIP3P-FB water, wet octanol
- `REF05`: YANK, SMIRNOFF force field, TIP3P water, wet octanol
- `REF06`: YANK, SMIRNOFF force field, OPC water, wet octanol
- `REF07`: YANK, GAFF force field, TIP3P water, dry octanol
- `REF08`: YANK, SMIRNOFF force field, TIP3P water, dry octanol
