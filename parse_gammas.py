#!/usr/bin/env python3
from sys import argv
import traceback
"""
take a dalton output file like output_files/NB_rotated_static_gamma_STO-3G_hf.out and extract the
gamma info
"""

def convert(k):
    """convert "Y;Y,Z,Z" to (1,1,2,2)
    output is a tuple
    """
    convert_axes = {
        "X": 0,
        "Y": 1,
        "Z": 2
    }
    k = k.replace(";", ",")
    parts = k.split(',')
    converted = [convert_axes[part] for part in parts]
    return tuple(converted)


def parse_gamma(gamma_line):
    """convert line from dalton to
    (
        (int, int, int, int),
        float
    )
    """
    _, line = gamma_line.split('(')
    k, v = line.split(')')
    k = convert(k)
    v = float(v)
    return (k, v)

def parse_freq(freq_line):
    """convert line from dalton to
    (
        (int, int, int, int),
        float
    )
    """
    _, _, v = freq_line.split()
    v = float(v)
    return v

def extract_gammas_and_freq_from_file(filename):
    with open(filename) as f:
        lines = f.readlines()
        gammas = []
        freqs = []
        warning_flag = False
        for line in lines:
            try:
                if "@ gamma(" in line:
                    gammas.append(parse_gamma(line))
                if "@ B-freq:" in line or "@ C-freq:" in line or "@ D-freq:" in line:
                    freqs.append(parse_freq(line))
                if "@WARNING" in line:
                    warning_flag = True
            except ValueError as err:
                print("error parsing file - {}".format(filename))
                print("errored on line - {}".format(line))
                traceback.print_tb(err.__traceback__)
                return [], [], True
        print("freqs = ",freqs)
        return gammas, freqs, warning_flag

if __name__ == "__main__":
    gamma_output, freqs = extract_gammas_and_freq_from_file(argv[1])
    print(gamma_output)
    print(freqs)
        

