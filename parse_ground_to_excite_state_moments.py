#!/usr/bin/env python3
from sys import argv
import traceback
import math
import scipy.constants as constants
"""
take a dalton output file like:
    globus\\nb_trans_dipole_hf-mp2-mcscf_determ_06_GW_dunning.out
and extract the dipole transition moments for the different excited states
"""

def remove_repeated_white_space(line):
    temp = str(line)
    while "  " in temp:
        temp = temp.replace("  ", " ")
    if line[0] == " ":
        temp = temp[1:]
    return temp

def eV_to_hartree(energy_eV):
    converstion_factor = constants.physical_constants["electron volt-hartree relationship"][0]
    return energy_eV*converstion_factor

def hartree_to_eV(energy_hartree):
    converstion_factor = constants.physical_constants["Hartree energy in eV"][0]
    return energy_eV*converstion_factor

def extract_info_from_file(file_path):
    labels = {"ZDIPLEN": 'z',"YDIPLEN": 'y',"XDIPLEN": 'x'}
    
    with open(file_path) as f:
        lines = f.readlines()
        transition_dipole_flag = False
        transition_dipoles = []
        temp_transition_dipoles = []
        temp = {}
        temp_label = ""
        try:
            for index, line in enumerate(lines):
                if "*** RSPCTL MICROITERATIONS CONVERGED" in line:
                    transition_dipole_flag = True
                if  transition_dipole_flag:
                    if  "@ Singlet transition operator label:" in line:
                        _, temp_line = line.split(":")
                        temp_line = remove_repeated_white_space(temp_line)
                        temp_line = temp_line.split(" ")[0]
                        temp_label = labels[temp_line]
                    if "@ STATE NO:" in line:
                        state_line, moment_line, energy_line = line.split("*")
                        _, state_no = state_line.split(":")
                        _, trans_moment = moment_line.split(":")
                        _, energy = energy_line.split(":")
                        temp = {
                                    "mu_label": temp_label,
                                    "mu_sym": None,
                                    "mu_spin": None,
                                    "B_state_no": 0,
                                    "B_sym": None,
                                    "B_spin": None,
                                    "C_state_no": int(state_no),
                                    "C_sym": None,
                                    "C_spin": None,
                                    "B_energy": 0.0,
                                    "C_energy": eV_to_hartree(float(energy)),
                                    "moment": float(trans_moment),
                                    }
                        temp_transition_dipoles.append(temp)
                        temp = {}
                    if "*** @ Excit." in line:
                        _, operator_sym_line, ref_sym_line, excite_sym_line = line.split("sym")
                        operator_sym = remove_repeated_white_space(operator_sym_line).split(" ")[0]
                        ref_sym = remove_repeated_white_space(ref_sym_line).split(" ")[0]
                        excite_sym = remove_repeated_white_space(excite_sym_line).split(" ")[1]
                        for temp_transition_dipole in temp_transition_dipoles:
                            temp_transition_dipole["B_sym"] = int(ref_sym)
                            temp_transition_dipole["C_sym"] = int(excite_sym)
                            temp_transition_dipole["mu_sym"] = int(operator_sym)
                        transition_dipoles += temp_transition_dipoles
                        temp_transition_dipoles = []
                        transition_dipole_flag = False
        except ValueError as err:
            print("error parsing file - {}".format(filename))
            print("errored on line - {}".format(line))
            traceback.print_tb(err.__traceback__)
    
    return transition_dipoles

if __name__ == "__main__":
    import pandas as pd
    transition_dipoles = extract_info_from_file(argv[1])
    df_transition_dipoles = pd.DataFrame(transition_dipoles)
    df_transition_dipoles.to_csv(argv[2])
    print(df_transition_dipoles)

