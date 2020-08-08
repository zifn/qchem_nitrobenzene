import os
import traceback
from sys import argv
import numpy as np
import json
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

def make_dal_file(file_path, freq_hartree_probe, freq_hartree_drive, state, spin_mult, max_ittr, symmetry, run_response=True):
    
    dal_wave = """**WAVE FUNCTIONS
.HF
.MP2
.MCSCF
*CONFIGURATION INPUT
.SPIN MULTIPLICITY
{0}
.SYMMETRY
{1}
.INACTIVE ORBITALS
15 10 0 0
.CAS SPACE
1 1 6 3 
.ELECTRONS
14
*ORBITAL INPUT
.MOSTART
.HUCKEL
*OPTIMIZATION
.DETERMINANTS
.STATE
{2}
*CI VECTOR
.PLUS COMBINATIONS
""".format(int(spin_mult), int(symmetry),  int(state))

    if run_response:
        fd = freq_hartree_drive
        fp = freq_hartree_probe
        freq_1, freq_2, freq_3 = -fd, fd, fp
    
        dal_input = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
"""

        dal_response = """**RESPONSE
*CUBIC
.MAX IT
{0}
.DIPLEN
.BFREQ
1
{1}
.CFREQ
1
{2}
.DFREQ
1
{3}
**END OF DALTON INPUT""".format(max_ittr, str(freq_1), str(freq_2), str(freq_3))

        dalton_dal_text = dal_input + dal_wave + dal_response
    else:
        dal_input = """**DALTON INPUT
.RUN PROPERTIES
.DIRECT
"""
        dal_properties = """**PROPERTIES
"""
        dal_end = "**END OF DALTON INPUT"
        dalton_dal_text = dal_input + dal_wave + dal_properties + dal_end

    with open(file_path, 'w') as dal_file:
        dal_file.write(dalton_dal_text)
    
    return file_path

def make_mol_file(file_path, CN_displacement=0, ONO_rotation=0):
    """
    ONO_rotation: Positive ONO_rotation value increases the ONO angle in radians from equilibrium
    CN_displacement: Positive CN_displacement increases the CN bond length
    """

    nitrogen_eq = np.array((0.0000000000, 0.0000000000, 1.4549291696))
    oxygen_eq = np.array((-1.0847503743, 0.0000000000, 2.0235779798))
    ONO_half_angle = angle_between(nitrogen_eq, oxygen_eq - nitrogen_eq)
    NO_distance = np.linalg.norm(oxygen_eq - nitrogen_eq)
    N_new_location = nitrogen_eq + np.array((0.0000000000, 0.0000000000, CN_displacement))
    O_new_location = np.array((NO_distance*np.sin(ONO_half_angle+ONO_rotation/2),
                               0.0000000000,
                               N_new_location[2] + NO_distance*np.cos(ONO_half_angle+ONO_rotation/2)))
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

def parse_config(json_path):
    file_path = os.path.normpath(json_path)
    with open(file_path) as config_file:
        raw_json = json.load(config_file)
    output_dir = os.path.join(raw_json["output_file_dir"])
    mpi_process_numb = raw_json["mpi_processes_numb"]
    memory_mb = raw_json['memory_mb']
    max_itter = raw_json["max_itter"]
    states = raw_json["states"]
    run_response =  bool(raw_json["run_response"])
    symmetries = raw_json["symmetries"]
    spin_mults = raw_json["spin_multiplicities"]
    hartree_freqs_probe = np.linspace(raw_json["hartree_freqs_probe"]['start'], raw_json["hartree_freqs_probe"]['end'], raw_json["hartree_freqs"]['points'])
    hartree_freqs_drive = np.linspace(raw_json["hartree_freqs_drive"]['start'], raw_json["hartree_freqs_drive"]['end'], raw_json["hartree_freqs"]['points'])
    CN_displacements = np.linspace(raw_json["CN_displacements"]['start'], raw_json["CN_displacements"]['end'], raw_json["CN_displacements"]['points'])
    ONO_rotations = np.linspace(raw_json["ONO_rotations"]['start'], raw_json["ONO_rotations"]['end'], raw_json["ONO_rotations"]['points'])
    return output_dir, states, spin_mults, hartree_freqs_probe, hartree_freqs_drive, CN_displacements, ONO_rotations, mpi_process_numb, memory_mb, max_itter, symmetries, run_response


def main(json_config_path):
    output_dir, states, spin_mults, hartree_freqs_probe, hartree_freqs_drive, CN_displacements, ONO_rotations, mpi_process_numb, memory_mb, max_itter, symmetries, run_response = parse_config(json_config_path)
    
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    
    temp_dir = os.path.join(output_dir, "temp")
    if os.path.isdir(temp_dir) == False:
        os.mkdir(temp_dir)
    
    for state in states:
        for symmetry in symmetries:
            for spin_mult in spin_mults:
                for CN_displacement in CN_displacements:
                    for ONO_rotation in ONO_rotations:
                        for hartree_freq_probe, hartree_freq_drive in zip(hartree_freqs_probe, hartree_freqs_drive):
                            root_name = "state-{}_sym-{}_freqd-{}_freqp-{}_spin-{}_CN_disp-{}_ONO_rot-{}_NBopt_dunningZ-2".format(state, symmetry, hartree_freq_probe, hartree_freq_drive, spin_mult, CN_displacement, ONO_rotation)
                            output_file_path = os.path.join(output_dir, root_name + ".out")
                            stdout_output_file_path = os.path.join(output_dir, root_name + ".stdout")
                            dal_file_path = os.path.join(temp_dir, root_name + ".dal")
                            mol_file_path = os.path.join(temp_dir, root_name + ".mol")
                            
                            next_dal_file_path = make_dal_file(dal_file_path, hartree_freq_probe, hartree_freq_drive, state, spin_mult, max_itter, symmetry, run_response)
                            next_mol_file_path = make_mol_file(mol_file_path, CN_displacement, ONO_rotation)
                            cmd_to_run = ['./dalton', '-mb', str(memory_mb), '-N', str(mpi_process_numb), '-o', str(output_file_path), str(next_dal_file_path), str(next_mol_file_path)]
                            if not os.path.isfile(stdout_output_file_path):
                                try:
                                    print("running next calculation: freq probe {}, freq drive {}, state {}, sym {}, spin {}, CN_disp {}, ONO_rot {}".format(hartree_freq_probe, hartree_freq_drive, state, symmetry, spin_mult, CN_displacement, ONO_rotation))
                                    print("\t running command - {}".format(cmd_to_run))
                                    exit_code, stdout, stderr = util_calc.run_cmd(cmd_to_run)
                                    with open(stdout_output_file_path, 'w') as file:
                                        file.write(str(stdout))
                                except Exception as err:
                                    print("An error occured trying next calculation")
                                    print(err.args)
                                    traceback.print_tb(err.__traceback__)
                                    pass

if __name__ == "__main__":
    main(argv[1])
