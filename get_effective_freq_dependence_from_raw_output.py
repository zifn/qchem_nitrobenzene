from sys import argv
import os
import itertools

from scipy import constants as consts
import pandas as pd
import numpy as np

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
    
def determine_drive_and_probe(lambda_1, lambda_2, lambda_3):
    freqs = [lambda_1, lambda_2, lambda_3]
    freq_drive = 0
    freq_probe = 0
    freq_combos = [(freq1, freq2, np.isclose(freq1, -freq2)) for freq1, freq2 in itertools.combinations(freqs, 2)]
    for freq1, freq2, is_drive in freq_combos:
        if is_drive:
            freq_drive = abs(freq1)
            freqs.remove(freq1)
            freqs.remove(freq2)
            break
    try:
        assert(len(freqs) == 1)
    except AssertionError as err:
        print('Printing frequency combinations: (freq1, freq2, np.isclose(freq1, -freq2))')
        print(freq_combos)
        raise AssertionError("Drive Frequency not found in file")
    freq_probe = freqs[0]
    return freq_probe, freq_drive

def main(dir_path, output_file_name, should_return=False):
    folder_location = os.path.normpath(dir_path)
    output_file_path = os.path.join(folder_location, output_file_name)
    output_data = []
    for file in os.listdir(folder_location):
        if is_dalton_output_file(file):
            #print('found file  = {}'.format(file))
            freq = find_parameter_used(file)
            freq_perm = find_parameter_used(file, 'perm-')
            file_output = main_conversion.main_conversion(os.path.join(folder_location, file), True)
            if file_output:
                chi3_eff_sym, chi3_eff_expr, chi3_eff_value, gamma_eff_value, gamma_rot_ave, chi3_rot_ave, chi3_sym, lambda_out, lambda_1, lambda_2, lambda_3, warning_flag = file_output
                freq_probe, freq_drive = determine_drive_and_probe(lambda_1, lambda_2, lambda_3)
                output_data.append({'freq probe (Hartree)': freq_probe,
                                    'freq probe (eV)': abs(convert_hartree_freq_to_eV(freq_probe)),
                                    'freq drive (Hartree)': freq_drive,
                                    'freq drive (eV)': abs(convert_hartree_freq_to_eV(freq_drive)),
                                    'freq_permutation': freq_perm,
                                    'gamma effective': gamma_eff_value,
                                    'chi3 effective': chi3_eff_value,
                                    'warning_flag': warning_flag})
    output_df = pd.DataFrame(output_data)
    print('saving data to -> {}'.format(output_file_path))
    output_df.to_csv(output_file_path)
    if should_return:
        return output_df

if __name__ == "__main__":
    """
    expects two input parameters
    the path to the directory containing the data
    a name for the output csv style data
    """
    main(argv[1], argv[2])