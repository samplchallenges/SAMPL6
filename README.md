# The SAMPL6 Blind Prediction Challenge for Computational Chemistry

This repository gives challenge details and inputs for the SAMPL6 challenge. 
This cycle we have migrated the data download package to GitHub so it will be version controlled and more broadly acccessible.
**Because these files are available publicly, we have no record of who downloads them. Therefore, you should sign up for notifications**.
Specifically, if you want to receive updates if we uncover any problems, it is imperative that you either (a) sign up for the SAMPL e-mail list via the D3R site, or (b) sign up for notifications of changes to this GitHub repository (the ``Watch'' button, above); ideally you would do both. 

## What's Here
- [Challenge Overview](#challenge-overview)
- `host_guest`: Directory containing inputs for the host-guest challenges, as well as supporting files and a README detailing their organization
- [Host-guest challenge instructions](host_guest_instructions.md): Detailed instructions on the host-guest component of the challenge.
- [Detailed host-guest description](host_guest_description.md): Detailed description of the hosts, guests, and background information.

## What's forthcoming
- Physical property challenge files (see description, below) 
- SAMPLing challenge files (see description, below) for both host-guest and physical properties
 
## Challenge Overview 
(This is reproduced from the [SAMPL6 Website](https://drugdesigndata.org/about/sampl6))

SAMPL6 includes challenges based on aqueous host-guest binding data (binding free energies and, optionally, binding enthalpies) for three different host molecules, and on physical properties (specifically, distribution coefficients and possibly solubilities), for a set of fragment-like molecules.
Host-guest systems test simulation methods, force fields, and solvent models, in the context of binding, without posing the setup issues and computational burden of protein simulations.
The physical properties offer efficient tests of force field accuracy when detailed simulations are used, and can also be used to test continuum solvation models and knowledge-based prediction methods. 

Machine-readable structure files for the hosts, guests, and physical property solutes will be provided when each challenge component opens.
Note that the physical property component of SAMPL6 is expected to have a later opening date than the host-guest part, because the measurements are being done on a later schedule. 

SAMPL6, will also introduce a new challenge component, the ``SAMPLing challenge'', where participants will be provided with input files for prepared systems that the organizers will be simulating, and evaluated on how efficiently their calculations; i.e., with the least simulation time and/or energy calls,  to within a specified tolerance of the organizers’ well-converged “gold standard” answers for the provided systems, using cutoffs and long-range treatments recommended by the organizers.
Participants will still be responsible for preparing the provided systems for use with their chosen free energy method, but provided files will include fully solvated systems including all force field parameters and other details.

We thank Drs. Bruce Gibb (Tulane U.) and Lyle Isaacs (U. Maryland) for providing the host-guest data, and John Chodera, Mehtap Isik, and Merck for the distribution coefficient data.

Further information on the host-guest and physical property components of SAMPL6 follow.
Tentatively, the host-guest components of the challenge opens Aug. 23 with submissions due Jan. 12.

### Gibb Deep Cavity Cavitand (Octa Acids) binding of guests

One host-guest series is based on the Gibb Deep Cavity Cavitands (GDCCs), or octa-acids, previously used in SAMPL4 and SAMPL5.
The two hosts, OA and TEMOA (previously OAH and OAME) are identical, except that TEMOA has four additional methyl groups, which alter the shape and depth of the hydrophobic cavity.
Both were developed in the laboratory of Dr. Bruce Gibb (Tulane U), who will provide binding free energies and enthalpies, measured by ITC, for eight guest molecules interacting with each host.
The measurements are done in 10 mM sodium phosphate buffer at pH 11.7 ± 0.1, and T = 298 K.
Host OA is described here: doi:10.1021/ja200633d; and host TEMOA is described here doi:10.1007/s10822-013-9690-2.
There are also a number of papers from SAMPL4 and SAMPL5 which discuss calculations for these systems, as summarized, respectively, in doi:10.1007/s10822-014-9735-1 and doi:10.1007/s10822-016-9974-4.
Existing benchmark datasets based on these hosts also may be of interest for those preparing to tackle these new complexes: https://github.com/MobleyLab/benchmarksets; this ``perpetual'' review paper also provides a good introduction to the sampling and experimental issues which are known to be relevant in these systems. 

### Cucubit[8]uril (CB8) binding of guests

This host-guest series is based on the host cucurbit[8]uril (CB8), which was used in SAMPL3, as previously summarized (DOI 10.1007/s10822-012-9554-1).
CB8 is the eight-membered relative of cucurbit[7uril, which was used in several other prior SAMPL challenges.
Data will be provided for ~14 guests, including several FDA approved drugs.
Background information on CB8 may be found in a number of publications, including DOI 10.1021/jp2110067, 10.1002/chem.201403405, and 10.1021/ja055013x.

### Physical properties
The SAMPL6 physical property challenge will center on predicting physicochemical properties for 25-50  fragment- and drug-like small molecules that small molecule protein kinase inhibitors (or fragments thereof).
Because the SAMPL5 logD challenge highlighted the difficulty in correctly predicting transfer free energies involving protonation states, we will provide participants with experimental pKa values for these compounds.
We will ask participants to predict distribution coefficients (logD) at a single pH and (as a separate challenge), provided the measurements can be completed in time, pH-dependent solubilities for these compounds.

The experimental data being measured include:

- pKa values, measured by electrochemical and/or UV-metric titration
- pH-dependent distribution coefficients (logD) of one or both of the following types:
  - water and cyclohexane (as in SAMPL5) 
  - water and octanol (new in SAMPL6) 
- pH-dependent solubility measurements performed using CheqSol (tentatively)

These measurements will be performed on Sirius T3 instruments from Sirius Analytical at Merck’s Rahway site.
The exact size of the dataset will depend on practical data collection throughput.
An initial batch of ~25 fragment-like compounds is currently being assayed, with the prospect for additional measurements performed subsequently.
Post-challenge follow-up experiments are possible and will be conducted as needed.

Distribution coefficients were included in the SAMPL5 challenge (overview doi:10.1007/s10822-016-9954-8 and experiment doi:10.1007/s10822-016-9971-7; JCAMD special issue https://link.springer.com/journal/10822/30/11/page/1); in many cases, they were predicted as if they were partition coefficients, using solvation free energies in the relevant solvents.
The difference between distribution coefficients (logD, which reflects the transfer free energy at a given pH including the effects of accessing all equilibrium  protonation states of the solute in each phase) and partition coefficients (logP, which reflects the free energy of transfer for the neutral form only) proved particularly important.
In some cases, other effects like the presence of small amount of water in cyclohexane may also have played a role.

Tentatively, this portion of the challenge will begin October 15 and close January 15; this page will be updated as additional details become clear.
We hope to provide a preliminary list of compounds by September 15 to allow participants to begin planning.
