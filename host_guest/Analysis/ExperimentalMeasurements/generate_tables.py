#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import math
import csv
import json
from collections import OrderedDict

import numpy as np
from simtk import unit as u


# =============================================================================
# CONSTANTS
# =============================================================================

T = 298 * u.kelvin
R = u.MOLAR_GAS_CONSTANT_R
RELATIVE_TITRANT_CONC_ERROR = 0.03

CB8_GUESTS_SMILES_PATH = '../../Isaacs_SAMPL6_guests.smi'
OA_GUESTS_SMILES_PATH = '../../Gibb_SAMPL6_guests.smi'

# Experimental results as provided by the Gibb and Isaacs groups.
# The error is relative. None means that the error is <1%.
EXPERIMENTAL_DATA = OrderedDict([
    ('OA-G0', OrderedDict([
        ('Ka', 1.47e+4 / u.molar), ('dKa', 0.02),
        ('DH', -4.84 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 0.84 * u.kilocalories_per_mole), ('dTDS', 0.04),
        ('n', 1), ('DG', -5.68 * u.kilocalories_per_mole)
    ])),
    ('OA-G1', OrderedDict([
        ('Ka', 2.57e+3 / u.molar), ('dKa', 0.01),
        ('DH', -5.52 * u.kilocalories_per_mole), ('dDH', 0.01),
        ('TDS', -0.86 * u.kilocalories_per_mole), ('dTDS', 0.1),
        ('n', 1), ('DG', -4.65 * u.kilocalories_per_mole)
    ])),
    ('OA-G2', OrderedDict([
        ('Ka', 1.4e+6 / u.molar), ('dKa', 0.01),
        ('DH', -12.07 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -3.69 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -8.38 * u.kilocalories_per_mole)
    ])),
    ('OA-G3', OrderedDict([
        ('Ka', 6.21e+3 / u.molar), ('dKa', None),
        ('DH', -7.53 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -2.35 * u.kilocalories_per_mole), ('dTDS', 0.01),
        ('n', 1), ('DG', -5.18 * u.kilocalories_per_mole)
    ])),
    ('OA-G4', OrderedDict([
        ('Ka', 1.64e+5 / u.molar), ('dKa', 0.01),
        ('DH', -6.92 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 0.19 * u.kilocalories_per_mole), ('dTDS', 0.13),
        ('n', 1), ('DG', -7.11 * u.kilocalories_per_mole)
    ])),
    ('OA-G5', OrderedDict([
        ('Ka', 2.33e+3 / u.molar), ('dKa', 0.01),
        ('DH', -5.31 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -0.71 * u.kilocalories_per_mole), ('dTDS', 0.02),
        ('n', 1), ('DG', -4.59 * u.kilocalories_per_mole)
    ])),
    ('OA-G6', OrderedDict([
        ('Ka', 4.37e+3 / u.molar), ('dKa', None),
        ('DH', -5.29 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -0.33 * u.kilocalories_per_mole), ('dTDS', 0.11),
        ('n', 1), ('DG', -4.97 * u.kilocalories_per_mole)
    ])),
    ('OA-G7', OrderedDict([
        ('Ka', 3.6e+4 / u.molar), ('dKa', None),
        ('DH', -7.44 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -1.23 * u.kilocalories_per_mole), ('dTDS', 0.02),
        ('n', 1), ('DG', -6.22 * u.kilocalories_per_mole)
    ])),

    ('TEMOA-G0', OrderedDict([
        ('Ka', 2.81e+4 / u.molar), ('dKa', None),
        ('DH', -7.85 * u.kilocalories_per_mole), ('dDH', 0.02),
        ('TDS', -1.77 * u.kilocalories_per_mole), ('dTDS', 0.08),
        ('n', 1), ('DG', -6.06 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G1', OrderedDict([
        ('Ka', 2.4e+4 / u.molar), ('dKa', 0.04),
        ('DH', -8.25 * u.kilocalories_per_mole), ('dDH', 0.04),
        ('TDS', -2.27 * u.kilocalories_per_mole), ('dTDS', 0.14),
        ('n', 1), ('DG', -5.97 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G2', OrderedDict([
        ('Ka', 9.82e+4 / u.molar), ('dKa', None),
        ('DH', -9.27 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -2.46 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -6.81 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G3', OrderedDict([
        ('Ka', 1.28e+4 / u.molar), ('dKa', 0.04),
        ('DH', -8.86 * u.kilocalories_per_mole), ('dDH', 0.02),
        ('TDS', -3.25 * u.kilocalories_per_mole), ('dTDS', 0.04),
        ('n', 1), ('DG', -5.6 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G4', OrderedDict([
        ('Ka', 5.12e+5 / u.molar), ('dKa', None),
        ('DH', -8.87 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -1.08 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -7.79 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G5', OrderedDict([
        ('Ka', 1.13e+3 / u.molar), ('dKa', None),
        ('DH', -7.96 * u.kilocalories_per_mole), ('dDH', 0.01),
        ('TDS', -3.8 * u.kilocalories_per_mole), ('dTDS', 0.03),
        ('n', 1), ('DG', -4.16 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G6', OrderedDict([
        ('Ka', 9.12e+3 / u.molar), ('dKa', 0.02),
        ('DH', -6.19 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -0.79 * u.kilocalories_per_mole), ('dTDS', 0.07),
        ('n', 1), ('DG', -5.4 * u.kilocalories_per_mole)
    ])),
    ('TEMOA-G7', OrderedDict([
        ('Ka', 1.07e+3 / u.molar), ('dKa', None),
        ('DH', -8.33 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -4.2 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -4.13 * u.kilocalories_per_mole)
    ])),

    ('CB8-G0', OrderedDict([
        ('Ka', 8.06e+4 / u.molar), ('dKa', 0.05),
        ('DH', -4.22 * u.kilocalories_per_mole), ('dDH', 0.02),
        ('TDS', 2.48 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -6.69 * u.kilocalories_per_mole)
    ])),
    ('CB8-G1', OrderedDict([
        ('Ka', 4.03e+5 / u.molar**2), ('dKa', 0.04),
        ('DH', -5.05 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 2.6 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 0.5), ('DG', -7.65 * u.kilocalories_per_mole)
    ])),
    ('CB8-G2', OrderedDict([
        ('Ka', 4.08e+5 / u.molar), ('dKa', 0.06),
        ('DH', -6.5 * u.kilocalories_per_mole), ('dDH', 0.01),
        ('TDS', 1.16 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -7.66 * u.kilocalories_per_mole)
    ])),
    ('CB8-G3', OrderedDict([
        ('Ka', 5.34e+4 / u.molar), ('dKa', 0.07),
        ('DH', -2.46 * u.kilocalories_per_mole), ('dDH', 0.03),
        ('TDS', 3.99 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -6.45 * u.kilocalories_per_mole)
    ])),
    ('CB8-G4', OrderedDict([
        ('Ka', 5.13e+5 / u.molar**3), ('dKa', 0.04),
        ('DH', -9.83 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -2.03 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 0.33), ('DG', -7.8 * u.kilocalories_per_mole)
    ])),
    ('CB8-G5', OrderedDict([
        ('Ka', 9.9e+5 / u.molar), ('dKa', 0.06),
        ('DH', -3.18 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 5.0 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -8.18 * u.kilocalories_per_mole)
    ])),
    ('CB8-G6', OrderedDict([
        ('Ka', 1.3e+6 / u.molar), ('dKa', 0.06),
        ('DH', -5.69 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 2.65 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -8.34 * u.kilocalories_per_mole)
    ])),
    ('CB8-G7', OrderedDict([
        ('Ka', 2.08e+7 / u.molar), ('dKa', 0.14),
        ('DH', -6.48 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 3.5 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -9.98 * u.kilocalories_per_mole)
    ])),
    ('CB8-G8', OrderedDict([
        ('Ka', 8.26e+9 / u.molar), ('dKa', 0.04),
        ('DH', -14.4 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', -0.88 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -13.5 * u.kilocalories_per_mole)
    ])),
    ('CB8-G9', OrderedDict([
        ('Ka', 2.29e+6 / u.molar), ('dKa', 0.1),
        ('DH', -4.63 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 4.05 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -8.68 * u.kilocalories_per_mole)
    ])),
    ('CB8-G10', OrderedDict([
        ('Ka', 1.05e+6 / u.molar), ('dKa', 0.09),
        ('DH', -2.0 * u.kilocalories_per_mole), ('dDH', 0.01),
        ('TDS', 6.22 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -8.22 * u.kilocalories_per_mole)
    ])),
    ('CB8-G11', OrderedDict([
        ('Ka', 4.98e+5 / u.molar), ('dKa', 0.06),
        ('DH', -2.11 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 5.67 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -7.77 * u.kilocalories_per_mole)
    ])),
    ('CB8-G12a', OrderedDict([
        ('Ka', 1.67e+7 / u.molar), ('dKa', 0.025),
        ('DH', -9.16 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 0.697 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -9.86 * u.kilocalories_per_mole)
    ])),
    ('CB8-G12b', OrderedDict([
        ('Ka', 1.46e+5 / u.molar**2), ('dKa', 0.01),
        ('DH', -4.83 * u.kilocalories_per_mole), ('dDH', None),
        ('TDS', 2.23 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 2), ('DG', -7.05 * u.kilocalories_per_mole)
    ])),
    ('CB8-G13', OrderedDict([
        ('Ka', 1.61e+5 / u.molar), ('dKa', 0.02),
        ('DH', -6.8 * u.kilocalories_per_mole), ('dDH', 0.01),
        ('TDS', 0.31 * u.kilocalories_per_mole), ('dTDS', None),
        ('n', 1), ('DG', -7.11 * u.kilocalories_per_mole)
    ])),
])


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_smiles(file_path):
    """Return the list of guests names and SMILES."""
    guests = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            smiles, name = line.split(' ', 1)
            guests.append([smiles, name])
    return guests


def compute_DG(Ka, dKa):
    """Compute the free energy from the association constant.

    Parameters
    ----------
    Ka : simtk.Quantity
        Association constant.
    dKa : simtk.Quantity
        Association constant uncertainty.

    Returns
    -------
    DG : simtk.Quantity
        Binding free energy.
    dDG : simtk.Quantity
        Binding free energy uncertainty.

    """
    concentration_unit = 1 / Ka.unit
    DG = -R * T * np.log(Ka*concentration_unit)
    # Propagate error.
    if dKa is None:
        dDG = None
    else:
        dDGdKa = -R * T / Ka  # Derivative dDG(Ka)/dKa.
        dDG = np.sqrt(dDGdKa**2 * dKa**2)
    return DG, dDG


def compute_TDS(DG, dDG, DH, dDH):
    """Compute the entropy from free energy and enthalpy.

    Parameters
    ----------
    DG : simtk.Quantity
        Free energy.
    dDG : simtk.Quantity
        Free energy uncertainty.
    DH : simtk.Quantity
        Enthalpy.
    dDH : simtk.Quantity
        Enthalpy uncertainty.

    Returns
    -------
    TDS : simtk.Quantity
        Entrop.
    dTDS : simtk.Quantity
        Binding free energy uncertainty.

    """
    TDS = DH - DG
    dTDS = np.sqrt(dDH**2 + dDG**2)
    return TDS, dTDS


def strip_units(quantities):
    for k, v in quantities.items():
        if isinstance(v, u.Quantity):
            # We only have energies and association constants.
            if 'Ka' in k:
                quantities[k] = v.value_in_unit(v.unit)
            else:
                quantities[k] = v.value_in_unit(u.kilocalories_per_mole)


def reduce_to_first_significant_digit(quantity, uncertainty):
    first_significant_digit = math.floor(math.log10(abs(uncertainty)))
    quantity = round(quantity, -first_significant_digit)
    uncertainty = round(uncertainty, -first_significant_digit)
    return quantity, uncertainty


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # Load names and SMILES of guests.
    molecule_names = {
        'CB8': load_smiles(CB8_GUESTS_SMILES_PATH),
        'OA' : load_smiles(OA_GUESTS_SMILES_PATH),
        'TEMOA' : load_smiles(OA_GUESTS_SMILES_PATH),
    }

    output_dict = OrderedDict()
    upper_bound_molecules = dict(Ka=set(), DH=set(), TDS=set())

    for system_name, system_data in EXPERIMENTAL_DATA.items():
        host_name, guest_name = system_name.split('-')
        guest_idx = int(guest_name[1:3])

        # Load SMILES and common name of the molecule.
        molecule_smiles, molecule_name = molecule_names[host_name][guest_idx]

        # Create entry in the output dictionary.
        output_dict[system_name] = OrderedDict([
            ('name', molecule_name),
            ('SMILES', molecule_smiles),
        ])
        output_dict[system_name].update(system_data)
        system_data = output_dict[system_name]  # Shortcut.

        # Incorporate the relative concentration uncertainties into quantities.
        for k in ['Ka', 'DH']:
            quantity = system_data[k]
            relative_uncertainty = system_data['d' + k]
            # Use upper-bound of 1% if <1% is reported. Keep track of these molecules.
            if relative_uncertainty is None:
                upper_bound_molecules[k].add(system_name)
                relative_uncertainty = 0.01
            # Incorporate the relative concentration uncertainties into quantities.
            relative_uncertainty += RELATIVE_TITRANT_CONC_ERROR
            # Convert relative to absolute errors.
            system_data['d' + k] = abs(quantity * relative_uncertainty)

        # Propagate Ka and DH error into DG and TDS.
        DG, dDG = compute_DG(system_data['Ka'], system_data['dKa'])
        system_data['dDG'] = dDG
        TDS, dTDS = compute_TDS(system_data['DG'], system_data['dDG'],
                                system_data['DH'], system_data['dDH'])
        system_data['dTDS'] = dTDS

        # Strip units.
        strip_units(system_data)

        # Consistency checks.
        computed_DG = DG.value_in_unit(u.kilocalories_per_mole)
        computed_TDS = TDS.value_in_unit(u.kilocalories_per_mole)
        assert np.isclose(system_data['DG'], system_data['DH'] - system_data['TDS'], atol=0.020000000000001, rtol=0.0)
        assert np.isclose(np.around(computed_TDS, decimals=2), system_data['TDS'], atol=0.0200000000000001, rtol=0.0)
        assert np.isclose(np.around(computed_DG, decimals=2), system_data['DG'], atol=0.0200000000000001, rtol=0.0)

        # Report only error most significant digit.
        for k in ['Ka', 'DH', 'TDS', 'DG']:
            quantity, uncertainty = system_data[k], system_data['d' + k]
            if uncertainty is not None:
                system_data[k], system_data['d' + k] = reduce_to_first_significant_digit(quantity, uncertainty)

    # Create output JSON file.
    with open('experimental_measurements.json', 'w') as f:
        json.dump(output_dict, f)

    # Create output CSV file.
    # Convert single dict to list of dicts.
    csv_dicts = []
    for system_id, system_data in output_dict.items():
        csv_dict = OrderedDict([('ID', system_id)])
        csv_dict.update(system_data)
        csv_dicts.append(csv_dict)
    with open('experimental_measurements.csv', 'w') as f:
        writer = csv.DictWriter(f, csv_dicts[0].keys(), delimiter=';')
        writer.writeheader()
        writer.writerows(csv_dicts)

    # Create a LaTex table.
    os.makedirs('PDFTable', exist_ok=True)
    old_host = ''
    with open('PDFTable/experimental_measurements.tex', 'w', encoding='utf-8') as f:
        f.write('\\documentclass{article}\n'
                '\\usepackage[a4paper,margin=0.4in,tmargin=0.5in,landscape]{geometry}\n'
                '\\usepackage{tabu}\n'
                '\\pagenumbering{gobble}\n'
                '\\begin{document}\n'
                '\\begin{center}\n'
                '\\begin{tabu}')

        # Cell alignment.
        field_names = ['ID', 'name', '$K_a$ (M$^{-1}$)', '$\\Delta G$ (kcal/mol) $^{(a)}$', '$\\Delta H$ (kcal/mol)', '$T\\Delta S$ (kcal/mol) $^{(b)}$', '$n$']
        f.write('{| ' + ' | '.join(['c' for _ in range(len(field_names))]) + ' |}\n')

        # Table header.
        f.write('\\hline\n')
        f.write('\\rowfont{\\bfseries} ' + ' & '.join(field_names) + ' \\\\\n')
        f.write('\\hline\n')

        # Print lines.
        for csv_dict in csv_dicts:
            # Separate hosts with a double horizontal line.
            host_name = csv_dict['ID'].split('-')[0]
            if host_name != old_host:
                f.write('\\hline\n')
                old_host = host_name

            row = '{ID} & {name}'
            for k in ['Ka', 'DG', 'DH', 'TDS']:
                row += ' & '

                # Report Ka in scientific notation.
                if k == 'Ka':
                    first_significant_digit = math.floor(math.log10(abs(csv_dict['d' + k])))
                    csv_dict['d' + k] /= 10**first_significant_digit
                    csv_dict[k] /= 10**first_significant_digit
                    row += '('
                row += '{' + k + '} +- {d' + k + '}'
                if k == 'Ka':
                    row += ') $\\times$ 10'
                    if first_significant_digit != 1:
                        row += '$^{{{{{}}}}}$'.format(first_significant_digit)

                # Check if we used the upperbound.
                superscript = ''
                # if k != 'DG' and csv_dict['ID'] in upper_bound_molecules[k]:
                #     superscript += 'a'
                if k == 'Ka':
                    if csv_dict['n'] == 0.33:
                        superscript += 'd'
                    elif csv_dict['n'] == 0.5 or csv_dict['n'] == 2:
                        superscript += 'c'
                if superscript != '':
                    row += ' $^{{(' + superscript + ')}}$'

            row += (' & {n} \\\\\n'
                    '\\hline\n')
            f.write(row.format(**csv_dict))

        f.write('\end{tabu}\end{center}\\vspace{5mm}\n'
                'All quantities are reported as point estimate +- statistical error from the ITC data fitting procedure. '
                'The upper bound ($1\%$) was used for errors reported to be $<1\%$. We also included a 3\% relative '
                'uncertainty in the titrant concentration assuming the stoichiometry coefficient to be fitted to the ITC '
                'data [1]. This is exact only for the OA/TEMOA sets (with the exception of OA-G5, TEMOA-G5, and TEMOA G7). '
                'For the other guests, we may expand the error analysis to include also the effect of the uncertainties '
                'in titrand concentration and cell volume. \\\\\n'
                '($^a$) Statistical errors were propagated from the $K_a$ measurements. \\\\\n'
                '($^b$) All experiments were performed at 298 K. \\\\\n'
                '($^c$) Units of M$^{-2}$. \\\\\n'
                '($^d$) Units of M$^{-3}$.\n'
                '\end{document}\n')
