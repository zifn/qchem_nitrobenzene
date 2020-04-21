#!/usr/bin/env python3
"""
    Main function for computing chi3 values from dalton calculated
    gamma values using the rotational averaging formule of
        Kwak 2015 Rigorous theory of molecular orientational nonlinear optics
    and Lorentz-Lorenz local field corrections
    
    Author: Richard Thurston 6/13/2019
"""

from sys import argv
import numpy as np
import sympy as sp
import scipy.constants as constants

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

def main_conversion(file_path, should_return=False):
    gamma_tuples, freqs, warning_flag = parse_gammas.extract_gammas_and_freq_from_file(file_path)
    if len(gamma_tuples) > 0: # check for valid dalton run
        freqs.insert(0, -sum(freqs))
        lambda_out, lambda_1, lambda_2, lambda_3 = [convert_hartree_freq_to_microns(freq) for freq in freqs]
        gamma_rot_ave, chi3_rot_ave, chi3_sym = gamma_to_chi3_conversion.compute_and_display_chi3_from_raw_gamma(
            gamma_tuples,
            number_density_of_liquid_nitrobenzene(),
            refractive_index_sq_of_liq_NB,
            lambda_out, 
            lambda_1, 
            lambda_2, 
            lambda_3
        )
        chi3_eff_sym = sp.symbols("chi^(3)_eff")
        chi3_eff_expr = chi3_sym[1][0][0][1] + chi3_sym[1][0][1][0]
        chi3_eff_value = chi3_rot_ave[1][0][0][1] + chi3_rot_ave[1][0][1][0]
        gamma_eff_value = gamma_rot_ave[1][0][0][1] + gamma_rot_ave[1][0][1][0]
        if should_return:
            return (chi3_eff_sym, chi3_eff_expr, chi3_eff_value, gamma_eff_value, gamma_rot_ave, chi3_rot_ave, chi3_sym, {"micron": [lambda_out, lambda_1, lambda_2, lambda_3], "hartree": freqs}, warning_flag)
        else:
            print_chi3_elements(chi3_sym, chi3_rot_ave)
            print("\n{0} = {1}".format(sp.latex(chi3_eff_sym), sp.latex(chi3_eff_expr)))
            print("{0} = {1} m^2 / V^2".format(sp.latex(chi3_eff_sym), chi3_eff_value))
            print("Warning Flag Raised = {}".format(warning_flag))
    return []

if __name__ == "__main__":
	main_conversion(argv[1])