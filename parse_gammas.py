#!/usr/bin/env python3

"""
take a dalton output file like output_files/NB_rotated_static_gamma_STO-3G_hf.out and extract the
gamma info
"""

def strip_gamma(k):
    return k.replace("gamma(", "").replace(")", "")


def convert(k):
    """convert "Y;Y,Z,Z" to "(1,1,2,2)"
    output is a string because json keys need to be strings
    """
    convert_axes = {
        "X": 0,
        "Y": 1,
        "Z": 2
    }
    k = k.replace(";", ",")
    parts = k.split(',')
    converted = [convert_axes[part] for part in parts]
    return str(tuple(converted))


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

    
def test_parse_gamma():
    test_data = [
        ("@ gamma(Y;Y,Z,Z)       -190.6259", ("(1, 1, 2, 2)", -190.6259)),
        ("@ gamma(Z;Z,Z,Z)       6205.8187", ("(2, 2, 2, 2)", 6205.8187)),
        ("@ gamma(X;X,X,X)        -18.3970", ("(0, 0, 0, 0)", -18.3970))
    ]
    for test_input, expected in test_data:
        assert parse_gamma(test_input) == expected, "got {0}, expected {1}".format(test_input, expected)


def test_strip_gamma():
    test_data = [
        ("gamma(Y;Y,Z,Z)", "Y;Y,Z,Z")
    ]
    for test_input, expected in test_data:
        test_output = strip_gamma(test_input)
        assert test_output == expected, test_output


def test_convert_to_ints():
    test_data = [
        ("Y;Y,Z,Z", "(1, 1, 2, 2)")
    ]
    for test_input, expected in test_data:
        test_output = convert(test_input)
        assert test_output == expected, "output: {0}, type: {1}, expected: {2}".format(test_output, type(test_output), expected)


def test():
    test_strip_gamma()
    test_convert_to_ints()
    test_parse_gamma()
    print("\nHooray, it worked!\n")


if __name__ == "__main__":
    # jank tests because I wanted to use pytest but then remembered my system python setup is borked
    test()

