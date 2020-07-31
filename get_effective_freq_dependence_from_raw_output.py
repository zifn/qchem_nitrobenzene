from sys import argv
import os
import itertools

from scipy import constants as consts
import pandas as pd
import numpy as np
import sympy as smp

import main_conversion

def is_dalton_output_file(file):
    root, ext = os.path.splitext(file)
    return ext == ".out"

def find_parameter_used(file_name, tag='freq-'):
    if tag in file_name:
        _, temp = file_name.split(tag)
        value = temp.split('_')[0]
        return float(value)
    else:
        return np.nan

def convert_hartree_freq_to_eV(freq_hartee):
    E_hartree_conversion = consts.physical_constants['Hartree energy in eV'][0]
    freq_eV = freq_hartee*E_hartree_conversion
    return freq_eV

def main(dir_path, output_file_name, should_return=False, should_display_chi3=False):
    folder_location = os.path.normpath(dir_path)
    output_file_path = os.path.join(folder_location, output_file_name)
    output_data = []
    for file in os.listdir(folder_location):
        if is_dalton_output_file(file):
            print('found file  = {}'.format(file))
            freq = find_parameter_used(file)
            freq_perm = find_parameter_used(file, 'perm-')
            file_output = main_conversion.main_conversion(os.path.join(folder_location, file), True)
            if file_output:
                chi3_eff_sym, chi3_eff_expr, chi3_eff_value, gamma_eff_value, gamma_rot_ave, chi3_rot_ave, chi3_sym, freqs, warning_flag = file_output
                probe_key = freqs['probe_key']
                drive_keys = freqs['drive_keys']
                #print("freqs[hartree] = " ,freqs['hartree'])
                output_data.append({'freq probe (Hartree)': freqs['hartree'][probe_key],
                                    'freq probe (eV)': abs(convert_hartree_freq_to_eV(freqs['hartree'][probe_key])),
                                    'freq probe (nm)': freqs['micron'][probe_key]*1000,
                                    'freq drive (Hartree)': abs(freqs['hartree'][drive_keys[0]]),
                                    'freq drive (eV)': abs(convert_hartree_freq_to_eV(freqs['hartree'][drive_keys[0]])),
                                    'freq drive (nm)': abs(freqs['micron'][drive_keys[0]])*1000,
                                    'freq_permutation': freq_perm,
                                    'gamma effective': gamma_eff_value,
                                    'chi3 effective': chi3_eff_value,
                                    'warning_flag': warning_flag})
    output_df = pd.DataFrame(output_data)
    #print('saving data to -> {}'.format(output_file_path))
    output_df.to_csv(output_file_path)
    if should_display_chi3:
        from IPython.display import Math, Latex
        display(Math("{0} = {1}".format(smp.latex(chi3_eff_sym), smp.latex(chi3_eff_expr))))
    if should_return:
        return output_df

if __name__ == "__main__":
    """
    expects two input parameters
    the path to the directory containing the data
    a name for the output csv style data
    """
    main(argv[1], argv[2])