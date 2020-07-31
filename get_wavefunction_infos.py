from sys import argv
import os
import itertools

import pandas as pd
from scipy import constants as consts

import parse_wavefunction_info

def is_dalton_output_file(file):
    root, ext = os.path.splitext(file)
    return ext == ".out"

def convert_hartree_freq_to_eV(freq_hartee):
    E_hartree_conversion = consts.physical_constants['Hartree energy in eV'][0]
    freq_eV = freq_hartee*E_hartree_conversion
    return freq_eV

def main(dir_path, output_file_name):
    folder_location = os.path.normpath(dir_path)
    output_file_path = os.path.join(folder_location, output_file_name)
    output_data = []
    for file in os.listdir(folder_location):
        if is_dalton_output_file(file):
            print('found file  = {}'.format(file))
            output = parse_wavefunction_info.extract_state_info_from_file(os.path.join(folder_location, file))
            if output["State"] != None:
                output_data.append(output)

    output_df = pd.DataFrame(output_data)
    min_energy = min(output_df["Energy (a.u.)"])
    output_df["Rel Energy (eV)"] = convert_hartree_freq_to_eV(output_df["Energy (a.u.)"] - min_energy)
    
    output_df.to_csv(output_file_path)
    return output_df

if __name__ == "__main__":
    """
    expects two input parameters
    the path to the directory containing the data
    a name for the output csv style data
    """
    df = main(argv[1], argv[2])
    print(df)