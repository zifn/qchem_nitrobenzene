#!/usr/bin/env python3
from sys import argv
import traceback
"""
take a dalton output file like:
    output_files/nb_dft_opt_static_quad_dipole_hf-mp2-mcscf_determ_05_state1.out 
and extract the dipole transition moments for the different excited states
"""

def is_start_of_new_transition_moment(line):
    key = "@ Transition moment <B | A - <A> | C> in a.u."
    return key in line
    
def parse_dipole_label(line):
    axes = {"ZDIPLEN": 'z', "YDIPLEN": 'y', "XDIPLEN": 'x'}
    _, values = line.split(":")
    values = remove_repeated_white_space(values)
    raw_axis, symmetry, spin = values.split(" ")
    axis = axes[raw_axis]
    return axis, int(symmetry), int(spin)
    
def parse_state_label(line):
    _, values = line.split(":")
    values = remove_repeated_white_space(values)
    state, symmetry, spin = values.split(" ")
    return int(state), int(symmetry), int(spin)
    
def parse_energies_and_moment(line):
    _, values = line.split(":")
    values = remove_repeated_white_space(values)
    b_energy, c_energy, moment = values.split(" ")
    return float(b_energy), float(c_energy), float(moment)

def remove_repeated_white_space(line):
    temp = str(line)
    while "  " in temp:
        temp = temp.replace("  ", " ")
    if line[0] == " ":
        temp = temp[1:]
    return temp
    
def parse_transition_moment_grouping(lines):
    dipole_info = parse_dipole_label(lines[0])
    state_B_info = parse_state_label(lines[1])
    state_C_info = parse_state_label(lines[2])
    energy_and_moment = parse_energies_and_moment(lines[4])
    output = {
        "B state level": state_B_info[0],
        "B state symmetry": state_B_info[1],
        "B State Spin": state_B_info[2],
        "B State energy": energy_and_moment[0],
        "C state level": state_C_info[0],
        "C state symmetry": state_C_info[1],
        "C State Spin": state_C_info[2],
        "C State energy": energy_and_moment[1],
        "dipole axis": dipole_info[0],
        "dipole symmetry": dipole_info[1],
        "dipole spin": dipole_info[2],
        "transition moment": energy_and_moment[2],
    }
    
    return output 

def extract_dipole_transition_moments_from_file(filename):
    with open(filename) as f:
        lines = f.readlines()
        grouped_lines = []
        temp = []
        lines_to_store = -1 # 5
        warning_flag = False
        for line in lines:
            if is_start_of_new_transition_moment(line):
                temp = []
                lines_to_store = 5
                continue
            if lines_to_store > 0:
                lines_to_store -= 1
                temp.append(line)
            elif lines_to_store == 0:
                grouped_lines.append(temp)
                lines_to_store -= 1
        
        processed_transition_moments = []
        for group in grouped_lines:
            processed_transition_moments.append(parse_transition_moment_grouping(group))
            
        return processed_transition_moments

if __name__ == "__main__":
    import pandas as pd
    transition_moments = extract_dipole_transition_moments_from_file(argv[1])
    pd_transition_moments = pd.DataFrame(transition_moments)
    pd_transition_moments .to_csv(argv[2])
    print(pd_transition_moments)
        

