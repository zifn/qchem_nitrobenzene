#!/usr/bin/env python3
from sys import argv
import traceback
import math
"""
take a dalton output file like:
    output_files\\nb_dft_opt_static_quad_dipole_hf-mp2-mcscf_determ_05_state1_roots3333.out
and extract transition dipole information
"""

def remove_repeated_white_space(line):
    temp = str(line)
    while "  " in temp:
        temp = temp.replace("  ", " ")
    if line[0] == " ":
        temp = temp[1:]
    return temp

def extract_info_from_file(file_path):
    labels = {"ZDIPLEN": 'z',"YDIPLEN": 'y',"XDIPLEN": 'x'}
    reference_dipole = {"z": 0, "y": 0, "x": 0}
    
    with open(file_path) as f:
        lines = f.readlines()
        transition_dipole_flag = False
        ABACUS_flag = False
        dipole_flag = False
        transition_dipoles = []
        temp = {}
        try:
            for index, line in enumerate(lines):
                if "@ Transition moment <B | A - <A> | C> in a.u. for" in line:
                    transition_dipole_flag = True
                if "FINAL RESULTS from ABACUS" in line:
                    ABACUS_flag = True
                if ABACUS_flag:
                    if "Dipole moment components" in line:
                        dipole_flag = True
                    if dipole_flag:
                        if "      z " in line or "      y " in line or "      x " in line:
                            temp_line = remove_repeated_white_space(line)
                            key, mu_au, mu_debye, mu_Cm = temp_line.split(" ")
                            reference_dipole[key] = float(mu_au)
                        if "Units:   " in line:
                            dipole_flag = False
                if transition_dipole_flag:
                    if "@ A operator label,    symmetry, spin:" in line:
                        _, temp_line = line.split(":")
                        temp_line = remove_repeated_white_space(temp_line)
                        A_label, A_sym, A_spin = temp_line.split(" ")
                        temp["A_label"] = labels[A_label]
                        temp["A_sym"] = int(A_sym)
                        temp["A_spin"] = int(A_spin)
                    if "@ B excited state no., symmetry, spin:" in line:
                        _, temp_line = line.split(":")
                        temp_line = remove_repeated_white_space(temp_line)
                        B_state_no, B_sym, B_spin = temp_line.split(" ")
                        temp["B_state_no"] = int(B_state_no)
                        temp["B_sym"] = int(B_sym)
                        temp["B_spin"] = int(B_spin)
                    if "@ C excited state no., symmetry, spin:" in line:
                        _, temp_line = line.split(":")
                        temp_line = remove_repeated_white_space(temp_line)
                        C_state_no, C_sym, C_spin = temp_line.split(" ")
                        temp["C_state_no"] = int(C_state_no)
                        temp["C_sym"] = int(C_sym)
                        temp["C_spin"] = int(C_spin)
                    if "@ B and C excitation energies, moment:" in line:
                        _, temp_line = line.split(":")
                        temp_line = remove_repeated_white_space(temp_line)
                        B_energy, C_energy, transition_moment = temp_line.split(" ")
                        temp["B_energy"] = float(B_energy)
                        temp["C_energy"] = float(C_energy)
                        temp["transition_moment"] = float(transition_moment)
                        transition_dipoles.append(temp)
                        temp = {}
                        transition_dipole_flag = False
        except ValueError as err:
            print("error parsing file - {}".format(filename))
            print("errored on line - {}".format(line))
            traceback.print_tb(err.__traceback__)

    for transition_dipole in transition_dipoles:
        ref_dipole_componet = reference_dipole[transition_dipole["A_label"]]
        transition_dipole["transition_moment"] = transition_dipole["transition_moment"] + ref_dipole_componet
    
    return transition_dipoles

if __name__ == "__main__":
    import pandas as pd
    transition_dipoles = extract_info_from_file(argv[1])
    df_transition_dipoles = pd.DataFrame(transition_dipoles)
    df_transition_dipoles(argv[2])
    print(df_transition_dipoles)
        

