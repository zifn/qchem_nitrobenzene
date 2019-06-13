#!/usr/bin/env python3
"""
    Main function for computing chi3 values from dalton calculated
    gamma values using the rotational averaging formule of
        Kwak 2015 Rigorous theory of molecular orientational nonlinear optics
    and Lorentz-Lorenz local field corrections
    
    Author: Richard Thurston 6/13/2019
"""

from sys import argv
import scipy.constants as constants

import parse_gammas
import gamma_to_chi3_conversion

def number_density_of_liquid_nitrobenzene():
    #calculate number density of liquid nitrobenzene
	NB_mass_density = 1.199 #g/cm^3
    NB_molar_mass = 123.11 # g/mol
    NB_numb_density_per_cm3 = np.multiply(np.true_divide(NB_mass_density, NB_molar_mass), constants.N_A) # numb/cm^3
    NB_numb_density_per_m3 =  np.multiply(NB_numb_density_g_per_cm3, 100**3)
    return NB_numb_density_per_m3

def main():
    lambda_out = .8 # 800nm OKE pulses
    lambda_1 = .8
    lambda_2 = .8
    lambda_3 = .8
    
    gamma_tuples = parse_gammas.extract_gammas_from_file(argv[1])
    if len(gamma_tuples) > 0:
        gamma_rot_ave, chi3_rot_ave = gamma_to_chi3_conversion.compute_and_display_chi3_from_raw_gamma(gamma_tuples,
                                                                                                       number_density_of_liquid_nitrobenzene(), 
                                                                                                       lambda_out, 
                                                                                                       lambda_1, 
                                                                                                       lambda_2, 
                                                                                                       lambda_3)

if __name__ == "__main__":
	main()