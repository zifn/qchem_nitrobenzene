from sys import argv
import os
import pandas as pd

import main_conversion

def is_dalton_output_file(file):
    root, ext = os.path.splitext(file)
    return ext == ".out"

def find_frequency_used(file):
    _, temp = file.split('freq-')
    value = temp.split('_')[0]
    return value

def main(dir_path, output_file_name):
    folder_location = os.path.normpath(dir_path)
    output_file_path = os.path.join(folder_location, output_file_name)
    output_data = []
    for file in os.listdir(folder_location):
        if is_dalton_output_file(file):
            print('found file  = {}'.format(file))
            freq = find_frequency_used(file)
            file_output = main_conversion.main_conversion(os.path.join(folder_location, file), True)
            if file_output:
                chi3_eff_sym, chi3_eff_expr, chi3_eff_value, gamma_eff_value, gamma_rot_ave, chi3_rot_ave, chi3_sym, lambda_out, lambda_1, lambda_2, lambda_3 = file_output
                output_data.append({'freq (Hartree)': freq, 'freq (nm)': abs(lambda_1), 'gamma effective': gamma_eff_value, 'chi3 effective': chi3_eff_value})
    output_df = pd.DataFrame(output_data)
    print('saving data to -> {}'.format(output_file_path))
    output_df.to_csv(output_file_path)

if __name__ == "__main__":
    """
    expects two input parameters
    the path to the directory containing the data
    a name for the output csv style data
    """
    main(argv[1], argv[2])