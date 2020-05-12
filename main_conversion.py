#!/usr/bin/env python3
"""
    Main function for computing chi3 values from dalton calculated
    gamma values using the rotational averaging formule of
        Kwak 2015 Rigorous theory of molecular orientational nonlinear optics
    and Lorentz-Lorenz local field corrections
    
    Author: Richard Thurston 6/13/2019
"""

from sys import argv
from copy import deepcopy
import numpy as np
import sympy as sp
import scipy.constants as constants
import itertools

import parse_gammas
import gamma_to_chi3_conversion

def display_chi3_elements(chi3_symbols, chi3_values):
    coordinates = [0, 1, 2]
    for i in coordinates:
        for j in coordinates:
            for k in coordinates:
                for l in coordinates:
                    display(Math("{0} = {1} \\text{{ }} \\text{{m}}^2 / \\text{{V}}^2".format(sp.latex(chi3_symbols[i][j][k][l]), chi3_values[i][j][k][l])))

def print_chi3_elements(chi3_symbols, chi3_values):
    """
    when formating with latex
    {0} = {1} \\text{{ }} \\text{{m}}^2 / \\text{{V}}^2
    """
    coordinates = [0, 1, 2]
    for i in coordinates:
        for j in coordinates:
            for k in coordinates:
                for l in coordinates:
                    print("{0} = {1} m^2 / V^2".format(sp.latex(chi3_symbols[i][j][k][l]), chi3_values[i][j][k][l]))

def number_density_of_liquid_nitrobenzene():
    #calculate number density of liquid nitrobenzene
    NB_mass_density = 1.199 #g/cm^3
    NB_molar_mass = 123.11 # g/mol
    NB_numb_density_per_cm3 = np.multiply(np.true_divide(NB_mass_density, NB_molar_mass), constants.N_A) # numb/cm^3
    NB_numb_density_per_m3 =  np.multiply(NB_numb_density_per_cm3, 100**3)
    return NB_numb_density_per_m3

def refractive_index_sq_of_liq_NB(lambda_um):
    """
    temp 20 deg C
    from https://refractiveindex.info/?shelf=organic&book=nitrobenzene&page=Kedenburg
    wavelen of light must be in microns
    """
    if np.isclose(lambda_um, 0.0):
        return 1 + 1.30628 + 0.00502
    else:
        return 1 + ((1.30628*np.square(lambda_um))/(np.square(lambda_um) - 0.02268)) + ((0.00502*np.square(lambda_um))/(np.square(lambda_um) - 0.18487))

def convert_hartree_freq_to_microns(freq):
    if np.isclose(freq, 0.0):
        return 0.0
    else:
        E_hartee_conv = constants.physical_constants['Hartree energy'][0]
        freq_SI = constants.h*constants.c/(freq*E_hartee_conv)
        freq_micron = freq_SI/constants.micron
        return freq_micron

def determine_drive_and_probe(freqs):
    freqs_keys = list(deepcopy(freqs).keys())
    drive_keys = []
    freq_combos = [(freq1_key, freq2_key, np.isclose(freqs[freq1_key], -freqs[freq2_key]))
                    for freq1_key, freq2_key in itertools.combinations(freqs.keys(), 2)]
    for freq1_key, freq2_key, is_drive in freq_combos:
        if is_drive:
            drive_keys.append(freq1_key)
            drive_keys.append(freq2_key)
            freqs_keys.remove(freq1_key)
            freqs_keys.remove(freq2_key)
            break
    try:
        assert(len(freqs_keys) == 1)
    except AssertionError as err:
        print('Printing frequency combinations: (freq1, freq2, np.isclose(freq1, -freq2))')
        print(freq_combos)
        raise AssertionError("Drive Frequency not found in file")
    probe_key = freqs_keys[0]
    return probe_key, drive_keys

def main_conversion(file_path, should_return=False):
    gamma_tuples, freqs, warning_flag = parse_gammas.extract_gammas_and_freq_from_file(file_path)
    if len(gamma_tuples) > 0: # check for valid dalton run
        probe_key, drive_keys = determine_drive_and_probe(freqs)
        freqs['A'] = -sum(freqs.values()) # add signal field
        freqs_micron = {key: convert_hartree_freq_to_microns(freq) for key, freq in freqs.items()}
        gamma_rot_ave, chi3_rot_ave, chi3_sym = gamma_to_chi3_conversion.compute_and_display_chi3_from_raw_gamma(
            gamma_tuples,
            number_density_of_liquid_nitrobenzene(),
            refractive_index_sq_of_liq_NB,
            freqs_micron['A'],
            freqs_micron['B'],
            freqs_micron['C'],
            freqs_micron['D']
        )
        
        # choose terms for the effective chi3
        indices = {'signal': 1, 'probe': 0, 'drive1': 0, 'drive2': 1}
        map1 = {field_key: indices[index_key] for index_key, field_key in zip(['signal', 'probe', 'drive1', 'drive2'], ['A', probe_key, drive_keys[0], drive_keys[1]])}
        map2 = {field_key: indices[index_key] for index_key, field_key in zip(['signal', 'probe', 'drive1', 'drive2'], ['A', probe_key, drive_keys[1], drive_keys[0]])}
        
        chi3_eff_sym = sp.symbols("chi^(3)_eff")
        chi3_eff_expr = chi3_sym[map1['A']][map1['B']][map1['C']][map1['D']] + chi3_sym[map2['A']][map2['B']][map2['C']][map2['D']]
        chi3_eff_value = chi3_rot_ave[map1['A']][map1['B']][map1['C']][map1['D']] + chi3_rot_ave[map2['A']][map2['B']][map2['C']][map2['D']]
        gamma_eff_value = gamma_rot_ave[map1['A']][map1['B']][map1['C']][map1['D']] + gamma_rot_ave[map2['A']][map2['B']][map2['C']][map2['D']]
        if should_return:
            return (chi3_eff_sym, 
            chi3_eff_expr, 
            chi3_eff_value, 
            gamma_eff_value, 
            gamma_rot_ave, 
            chi3_rot_ave,
            chi3_sym,
            {"micron": freqs_micron, "hartree": freqs, 'probe_key': probe_key, 'drive_keys': drive_keys},
            warning_flag)
        else:
            print_chi3_elements(chi3_sym, chi3_rot_ave)
            print("\n{0} = {1}".format(sp.latex(chi3_eff_sym), sp.latex(chi3_eff_expr)))
            print("{0} = {1} m^2 / V^2".format(sp.latex(chi3_eff_sym), chi3_eff_value))
            print("Warning Flag Raised = {}".format(warning_flag))
    return []

if __name__ == "__main__":
	main_conversion(argv[1])