# Manifest

For each of the 5 configurations (0 through 4) of the three
compounds (CB8-G3, OA-G3, OA-G6), the following files are included:

AMBER/
	complex.prmtop
	complex.rst7
	solvent.prmtop
	solvent.rst7
CHARMM/
	complex.inp
	complex.psf
	complex.crd
	complex.rtf
	complex.prm
	solvent.inp
	solvent.psf
	solvent.crd
	solvent.rtf
	solvent.prm
DESMOND/
	complex.cms
	solvent.cms
GROMACS/
	complex.top
	complex.gro
	solvent.top
	solvent.gro
GROMACS_old/
	complex.top
	complex.gro
	solvent.top
	solvent.gro
OPENMM/
	complex.xml
	solvent.xml
LAMMPS/
	complex.lmp
	complex.input
	solvent.lmp
	solvent.input
PDB/
	complex.pdb
	solvent.pdb

SAMPL6_energyoutput_conversion_longcutoff.txt
SAMPL6_energyoutput_conversion_shortcutoff.txt

runfiles/
	min_SAMPL6.in
	grompp_SAMPL6.mdp
	onepoint_SAMPL6.cfg	 
	min_default.in
	grompp_default.mdp
	onepoint_default.cfg  

# Explanation

Files in `OPENMM/`, `AMBER/`, `PDB/` and `GROMACS_old` were generated as describe
in `SAMPLing_instructions.md`. Other simulation input files are generated as described here.

The SAMPL6_energyoutput_conversion files contain the outputs of
InterMol called with two different cutoff schemes, and the runfiles/
directories contain the parameter files used to generated the energies.

Conversion was carried out from the AMBER files listed above to
the GROMACS, LAMMPS, CHARMM, and DESMOND files using InterMol
(https://github.com/shirtsgroup/InterMol, Git hash f691465,
May 24,2017), using the ParmEd library
(https://github.com/ParmEd/ParmEd, Git hash 0bab490, Dec 11, 2017,
with PR #935 merged in for increased energy precision).

The command used was `python convert.py --odir $PROGRAMDIR --gromacs --charmm --lammps --desmond --amber --energy --amb_in AMBER/$SYSTEMTYPE.prmtop AMBER/$SYSTEMTYPE.rst --gropath $GROPATH --amberpath $AMBERPATH --lmppath $LMPSPATH --charmmpath $CHARMMPATH --inefile runfiles/min_default.in -gs runfiles/grompp_default.mdp -ds runfiles/onepoint_default.cfg -as runfiles/min_default.in -ls pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1.0e-5\n\n -cs nbond inbfrq -1 imgfrq -1 -\nelec ewald pmew fftx 48 ffty 48 fftz 48 kappa 0.34 order 4 -\nvdw vips cutnb 12. cutim 12. ctofnb 10. ctonnb 9.`

Where $PROGRAMDIR is the program (one of 'GROMACS', 'DESMOND', 'LAMMPS', 'CHARMM'), $SYSTEMTYPE is one of 'complex' or 'solvent', $GROPATH is the path of the `bin` GROMACS binary used to calculate energies, $AMBERPATH is the path of the `sander` binary used to calculate energies, $LMPPATH is the path of the LAMMPS binary used to calculate energies, and $CHARMMPATH is the path of CHARMM binary used to calculate energies.  DESMOND files use the `desmond` binary that is in the user-defined $SCHRODINGER path.

Programs used were:

* AMBER: Amber 16 release, `sander` binary
* CHARMM: Developmental Version 40b2   February 15, 2016  
* DESMOND: Part of the academic release of Schrodinger2016-1.
* LAMMPS: Downloaded on 16 Feb 2016
* GROMACS: gmx mdrun, VERSION 5.1.2 (double precision)

For the energy comparisons in `SAMPL6_energyoutput_conversion_longcutoff.txt`, the same command was used except for the following differences in the options: `--inefile runfiles/min_SAMPL6.in -gs runfiles/grompp_SAMPL6.mdp -ds runfiles/onepoint_SAMPL6.cfg -as runfiles/min_SAMPL6.in -ls pair_style lj/cut/coul/long 14.0 14.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n -cs nbond inbfrq -1 imgfrq -1 -\nelec ewald pmew fftx 48 ffty 48 fftz 48 kappa 0.22310095 order 4 -\nvdw vips cutnb 14. cutim 14. ctofnb 14. ctonnb 14.`

The AMBER `prmtop` and `rst7` files were converted to GROMACS
`top`/`gro` and PDB formats by ParmEd version 2.7.3. These were the
first originally converted, but used a less precise coversion metric,
and are stored in `GROMACS_old/`.  The parameters in the `GROMACS/`
folder are more accurate in energies, the the parameters in the
`GROMACS_old/` folder, but the difference should be negligible
compared to the statistical error in the free energy calculations.

The energies in the `amber` entries were the energies obtained by
converting the AMBER files to GROMACS, and then back to AMBER using
ParmEd routines. They are always in agreement within 0.0001 to 0.0002
kcal/mol (0.0004184 kJ/mol = 0.0001 kcal/mol is the limit of
resolution in `sander`) for total energies of magnitude 10^5 kcal/mol,
which means that essentially all of the rest of the difference of
energies is due to the ways that the programs calculate energies,
particularly the method by which the long range component of nonbonded
energies are calculated (see Shirts et al. below for further
discussion of energies), such as small differences in tapering
functions for Lennard-Jones terms and in Ewald summation algorithms.

A known issue in comparing energies between different programs is the
difference in the Coulomb's law constant, as analyzed in Shirts et
al. "Lessons learned from comparing molecular dynamics engines on the
SAMPL5 dataset", J. Comput.-Aided Mol. Design, 31: 1, pp 147-161
(2017) (DOI:10.1007/s10822-016-9977-1.  This is particularly
noticeable for AMBER and CHARMM, which for historical reasons have the
furthest choices from each other. The difference in Coulomb's law
constants directly affects the Coulombic energies, and explains 99% of
the difference in Coulomb 1-4 interactions, and approximately 70% of
the difference in the total Coulomb interactions for standard cutoff
schemes. This can be seen in the fact that the GROMACS, LAMMPS and
DESMOND energies are generally always within 1 kJ/mol from each other,
whereas the AMBER and CHARMM energies are further off (2-4 kJ/mol from
the others.)

We have not corrected the topology files for this difference in
Coulomb's law constant, as it is not expected that these differences
with significantly affect free energies calculated. The difference is
an overall scaling to the Coulomb energy, and will mostly cancel when
taking free energy differences, resulting in an error that is of
order of 1 part in 20,000.

