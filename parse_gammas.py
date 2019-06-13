#!/usr/bin/env python3
from sys import argv
"""
take a dalton output file like output_files/NB_rotated_static_gamma_STO-3G_hf.out and extract the
gamma info
"""

def strip_gamma(k):
    return k.replace("gamma(", "").replace(")", "")


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
    _, k, v = gamma_line.split()
    k = convert(strip_gamma(k))
    v = float(v)
    return (k, v)


def extract_gammas_from_file(filename):
    with open(filename) as f:
        lines = f.readlines()
        gammas = []
        for line in lines:
            if "@ gamma(" in line:
                gammas.append(parse_gamma(line))
        return gammas

if __name__ == "__main__":
    # jank tests because I wanted to use pytest but then remembered my system python setup is borked
    gamma_output = extract_gammas_from_file(argv[1])
    print(gamma_output)
        

