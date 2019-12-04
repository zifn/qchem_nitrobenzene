import os
import traceback
from sys import argv
import numpy as np
import util_calc

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def make_dal_file(file_path, freq_hartree, state):
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
.STATE
{1}
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
**END OF DALTON INPUT""".format(freq_hartree, int(state))
    with open(file_path, 'w') as dal_file:
        dal_file.write(dal_input)
    
    return file_path

def make_mol_file(file_path, CN_displacement=0, ONO_rotation=0):
    nitrogen_eq = np.array((0.0000000000, 0.0000000000, 1.4549291696))
    oxygen_eq = np.array((-1.0847503743, 0.0000000000, 2.0235779798))
    ONO_half_angle = angle_between(nitrogen_eq, oxygen_eq)
    NO_distance = np.linalg.norm(oxygen_eq - nitrogen_eq)
    N_new_location = nitrogen_eq + np.array((0.0000000000, 0.0000000000, CN_displacement))
    O_new_location = np.array((NO_distance*np.sin(ONO_half_angle+ONO_rotation/2), 0.0000000000, N_new_location[2] + NO_distance*np.cos(ONO_half_angle+ONO_rotation/2)))
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
N        0.0000000000            0.0000000000            {0}
Charge=8.0 Atoms=1
O       {1}            {2}            {3}
Charge=1.0 Atoms=3
H       -2.1410602955            0.0000000000           -0.1327964319
H       -2.1492827751            0.0000000000           -2.6311811056
H        0.0000000000            0.0000000000           -3.8683953478    
""".format(N_new_location[2], O_new_location[0], O_new_location[1], O_new_location[2])
    with open(file_path, 'w') as mol_file:
        mol_file.write(mol_input)
    
    return file_path

def main():
    cwd = "qchem_nitrobenzene"#os.getcwd()
    dal_file_path = os.path.join(cwd, "temp_dal_file.dal")
    mol_file_path = os.path.join(cwd, "temp_mol_file.mol")
    output_dir = os.path.join(cwd, "main_calc_output_files")

    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)

    state = 2
    for freq_hartree in np.linspace(0, 0.06, 21):
        output_file_path = os.path.join(output_dir, "state-{}_freq-{}_cubic_response_NBopt_dunningZ-2.out".format(state, freq_hartree))
        stdout_output_file_path = os.path.join(output_dir, "state-{}_freq-{}_cubic_response_NBopt_dunningZ-2.stdout".format(state, freq_hartree))
        next_dal_file_path = make_dal_file(dal_file_path, freq_hartree, state)
        next_mol_file_path = make_mol_file(mol_file_path)
        cmd_to_run = ['./dalton', '-mb', '8000', '-o', str(output_file_path), str(next_dal_file_path), str(next_mol_file_path)]
        try:
            print("running next calculation freq = {}".format(freq_hartree))
            print("\t running command - {}".format(cmd_to_run))
            exit_code, stdout, stderr = util_calc.run_cmd(cmd_to_run)
            with open(stdout_output_file_path, 'w') as file:
                file.write(str(stdout))
        except Exception as err:
            print("An error occured trying next calculation")
            traceback.print_tb(err.__traceback__)
            pass

if __name__ == "__main__":
    main()
