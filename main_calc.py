import os
from sys import argv
import numpy as np
import util_calc


def make_dal_file(file_path, freq_hartree):
    dal_input = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
**WAVE FUNCTIONS
.HF
.MP2
.MCSCF
*CONFIGURATION INPUT
.SPIN MULTIPLICITY
1
.SYMMETRY
1
.INACTIVE ORBITALS
15 10 2 0
.CAS SPACE
1 1 2 3 
.ELECTRONS
10
*ORBITAL INPUT
.MOSTART
.HUCKEL
*OPTIMIZATION
.DETERMINANTS
*CI VECTOR
.PLUS COMBINATIONS
**RESPONSE
*CUBIC
.DIPLEN
.BFREQ
1
{0}
.CFREQ
1
-{0}
.DFREQ
1
{0}
**END OF DALTON INPUT""".format(freq_hartree)
    with open(file_path, 'w') as dal_file:
        dal_file.write(dal_input)
    
    return file_path

def make_mol_file(file_path):
    mol_input=  """BASIS
aug-cc-pCVDZ
 nitrobenzene
 Generated by the Dalton Input File Plugin for Avogadro
Atomtypes=4 Angstrom Generators=2 X Y
Charge=6.0 Atoms=4
C        0.0000000000            0.0000000000           -0.0265759623
C       -1.2190534920            0.0000000000           -0.6979856948
C       -1.2105095870            0.0000000000           -2.0897394457
C        0.0000000000            0.0000000000           -2.7841484631
Charge=7.0 Atoms=1
N        0.0000000000            0.0000000000            1.4549291696
Charge=8.0 Atoms=1
O       -1.0847503743            0.0000000000            2.0235779798
Charge=1.0 Atoms=3
H       -2.1410602955            0.0000000000           -0.1327964319
H       -2.1492827751            0.0000000000           -2.6311811056
H        0.0000000000            0.0000000000           -3.8683953478    
"""
    with open(file_path, 'w') as mol_file:
        mol_file.write(mol_input)
    
    return file_path

def main():
    cwd = os.getcwd()
    dal_file_path = os.path.join(cwd, "temp_dal_file.dal")
    mol_file_path = os.path.join(cwd, "temp_mol_file.dal")
    output_dir = os.path.join(cwd, "output_files")

    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)

    for freq_hartree in np.linspace(0, 0.06, 21):
        output_file_path = os.path.join(output_dir, "freq-{}_cubic_response_NBopt_dunningZ-2.out".format(freq_hartree))
        stdout_output_file_path = os.path.join(output_dir, "freq-{}_cubic_response_NBopt_dunningZ-2.stdout".format(freq_hartree))
        cmd_to_run = ['./dalton', '-mb', '8000', '-o', str(output_file_path), str(make_dal_file(dal_file_path, freq_hartree)), str(make_mol_file(mol_file_path))]
        try:
            print("running next calculation freq = {}".format(freq_hartree))
            print("\t running command - {}".format(cmd_to_run))
            exit_code, stdout, stderr = util_calc.run_cmd(cmd_to_run)
            with open(stdout_output_file_path, 'w') as file:
                file.write(stdout)
        except:
            print("error with output_file_path tyring next calculation")
            pass

if __name__ == "__main__":
    main()