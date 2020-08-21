from sys import argv
import os
import itertools

import pandas as pd
from scipy import constants as consts

import parse_ground_to_excite_state_moments
import parse_excited_transition_dipoles

def is_dalton_output_file(file):
    root, ext = os.path.splitext(file)
    return ext == ".out"

def convert_hartree_freq_to_eV(freq_hartee):
    E_hartree_conversion = consts.physical_constants['Hartree energy in eV'][0]
    freq_eV = freq_hartee*E_hartree_conversion
    return freq_eV

def main(es_es_dipole_path, gs_es_dipole_path, output_file_path=None):
    excite_dipoles = parse_excited_transition_dipoles.extract_info_from_file(es_es_dipole_path)
    ground_dipoles = parse_ground_to_excite_state_moments.extract_info_from_file(gs_es_dipole_path)
    dipoles = excite_dipoles + ground_dipoles
    df_dipoles = pd.DataFrame(dipoles)
    if output_file_path != None:
        df_dipoles.to_csv(output_file_path)
    return df_dipoles

if __name__ == "__main__":
    """
    expects two input parameters
    the path to the directory containing the data
    a name for the output csv style data
    """
    _, excited_transition_dipole_path, ground_state_transition_dipole_path, output_file_path = argv
    df = main(excited_transition_dipole_path, ground_state_transition_dipole_path, output_file_path)
    print(df)