#!/usr/bin/env python3
from sys import argv
import traceback
import math
"""
take a dalton output file like:
    globus/wavefunctions_0033/state-1_sym-2_freqd-0.1_freqp-0.1_spin-1_CN_disp-0.0_ONO_rot-0.0_NBopt_dunningZ-2.out
and extract basic state information
"""

def remove_repeated_white_space(line):
    temp = str(line)
    while "  " in temp:
        temp = temp.replace("  ", " ")
    if line[0] == " ":
        temp = temp[1:]
    return temp

def extract_state_info_from_file(file_path):
    state, spin_mult, sym, charge , energy, dipole_mag = None, None, None, None, None, None
    with open(file_path) as f:
        lines = f.readlines()
        dipole_flag = False
        ABACUS_results_flag = False
        flagged_index = math.inf
        try:
            for index, line in enumerate(lines):
                if "@    Spin multiplicity:" in line:
                    _, value = line.split(":")
                    spin_mult = int(value)
                elif "@    Spatial symmetry:" in line:
                    temp = remove_repeated_white_space(line)
                    temp = temp.split(" ")
                    sym = int(temp[3])
                elif "@    Total charge of molecule" in line:
                    _, value = line.split(":")
                    charge = float(value)
                elif "@    State number:" in line:
                    _, value = line.split(":")
                    state = int(value)
                elif "@    Final MCSCF energy:" in line:
                    _, value = line.split(":")
                    energy = float(value)
                elif "FINAL RESULTS from ABACUS" in line:
                    ABACUS_results_flag = True
                elif ABACUS_results_flag and "Dipole moment" in line and "components" not in line:
                    dipole_flag = True
                    flagged_index = index + 4
                elif dipole_flag and flagged_index == index:
                    temp_line = remove_repeated_white_space(line)
                    temp = temp_line.split(" ")
                    dipole_mag = float(temp[1])
        except ValueError as err:
            print("error parsing file - {}".format(filename))
            print("errored on line - {}".format(line))
            traceback.print_tb(err.__traceback__)

    output = {
    "State": state,
    "Spin multiplicity": spin_mult,
    "Symmetry": sym,
    "Charge": charge,
    "Energy (a.u.)": energy,
    "dipole (Debye)": dipole_mag
    }
    
    return output

if __name__ == "__main__":
    state_info = extract_state_info_from_file(argv[1])
    print(state_info)
        

