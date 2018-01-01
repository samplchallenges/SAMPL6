# Manifest

For each of the 5 configurations (0 through 4) of the three
compounds(CB8-G3, OA-G3, OA-G6), the following files are included:

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

In addition, at the top level, the directory contains the files:


runfiles/
	min_SAMPL6.in    
	grompp_SAMPL6.mdp
	onepoint_SAMPL6.cfg	 
	min_default.in
	grompp_default.mdp
	onepoint_default.cfg  

Conversion was carried out from the AMBER files listed above to the
GROMACS, LAMMPS, CHARMM, and DESMOND using InterMol
(https://github.com/shirtsgroup/InterMol, git hash git hash  f691465, May 24,2017), using the
ParmED library (https://github.com/ParmEd/Parmed, git hash  0bab490, Dec 11, 2017, with PR #935 merged in for increased energy precision).

The command used was `python convert.py --odir $PROGRAMDIR --gromacs --charmm --lammps --desmond --amber --energy --amb_in AMBER/$SYSTEMTYPE.prmtop AMBER/$SYSTEMTYPE.rst --gropath $GROPATH --amberpath $AMBERPATH --lmppath $LMPSPATH --charmmpath $CHARMMPATH --inefile runfiles/min_default.in -gs runfiles/grompp_default.mdp -ds runfiles/onepoint_default.cfg -as runfiles/min_default.in -ls pair_style lj/cut/coul/long 9.0 9.0\npair_modify tail yes\nkspace_style pppm 1.0e-5\n\n -cs nbond inbfrq -1 imgfrq -1 -\nelec ewald pmew fftx 48 ffty 48 fftz 48 kappa 0.34 order 4 -\nvdw vips cutnb 12. cutim 12. ctofnb 10. ctonnb 9.`

Where $PROGRAMDIR is the program (one of 'GROMACS', 'DESMOND', 'LAMMPS', 'CHARMM'), $SYSTEMTYPE is one of 'complex' or 'solvent', $GROPATH is the path of the `bin` GROMACS binary used for , $AMBERPATH is the path of the `sander` binary used to calculate energies, $LMPPATH is the path of the LAMMPS binary used to calculate energies, and $CHARMMPATH is the path of CHARMM binary used to calculate energies.  Desmond usage uses the `desmond` binary that is in the user-defined $SCHRODINGER path.

Programs used were:

* AMBER: Amber 16 release, `sander` binary
* CHARMM: Developmental Version 40b2   February 15, 2016  
* DESMOND: Part of the academic release of Schrodinger2016-1.
* LAMMPS: Downloaded on 16 Feb 2016
* GROMACS: gmx mdrun, VERSION 5.1.2 (double precision)

For the energy comparisons in `SAMPL6_energyoutput_conversion_longcutoff.txt`, the same command was used except differences in the options `--inefile runfiles/min_SAMPL6.in -gs runfiles/grompp_SAMPL6.mdp -ds runfiles/onepoint_SAMPL6.cfg -as runfiles/min_SAMPL6.in -ls pair_style lj/cut/coul/long 14.0 14.0\npair_modify tail yes\nkspace_style pppm 1e-8\n\n -cs nbond inbfrq -1 imgfrq -1 -\nelec ewald pmew fftx 48 ffty 48 fftz 48 kappa 0.22310095 order 4 -\nvdw vips cutnb 14. cutim 14. ctofnb 14. ctonnb 14.`

The OPENMM, GROAMCS_old, and PDB files were converted using the
default options of ParmEd (Git hash ???). The parameters for GROMACS are
therefore slightly more accurate the the parameters in the GROMACS_old
file.

The energies in the `amber` entries were the energies obtained by
converting the AMBER files to GROMACS, and then back to AMBER using
ParmEd routines. They are always in agreement within 0.0001 to 0.0002
kcal/mol (0.0004184 kJ/mol = 0.0001 kcal/mol is the limit of
resolution in `sander`) for total energies of magnitude 10^5 kcal/mol, which means
that essentially all of the rest of the difference of energies is due
to the ways that the programs calculate energies, particularly the
method by which the long range component of nonbonded energies are
calculated (see Shirts et al. below for further discussion of energies).

A known issue in comparing energies between different programs is the
difference in the permittivity constant between the programs, as
analyzed in Shirts et al. "Lessons learned from comparing molecular
dynamics engines on the SAMPL5 dataset", J. Comput.-Aided Mol. Design,
31: 1, pp 147-161 (2017) (DOI:10.1007/s10822-016-9977-1.  This is
particularly noticible for AMBER and CHARMM, which for historical
reasons have the furthest choices from each other. The difference in
permittivity constants affects the Coulombic energies, and explains
99% of the difference in Coulomb 1-4 interactions, and approximately
70% of the difference in the total Coulomb interactions for standard
cutoff schemes. This can be seen in the fact that the GROMACS,LAMMPS
and DESMOND energies are generally always within 1 kJ/mol from each
other, whereas the AMBER and CHARMM energies are further off. 

We have not corrected the topology files for this difference, as it is
not expected that these differences with significantly affect free
energies calculated, since it is an overall scaling to the Coulomb
energy, and will mostly cancel, resulting in an error that will be of
order of 1 part in 20,000.


