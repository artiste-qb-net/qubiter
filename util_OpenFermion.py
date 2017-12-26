from SEO_writer import *


def openfermion_qasm_to_qbtr_files(file_prefix, num_bits, qasm_gates):
    """
    The OpenFermion library at github performs certain tasks that arise when
    applying quantum computing to quantum chemistry. The file
    openfermion/_utils/_trotter_exp_to_qgates.py contains various functions
    such as pauli_exp_to_qasm() that return a list of strings (or an
    iterable of such strings) called here `qasm_gates'. Each item in the
    list `qasm_gates` denotes a gate in IBM qasm. This function takes a
    `qasm_gates` string list as input and creates English and Picture files
    that are translations of `qasm_gates'.

    References
    ----------
    1. OpenFermion repo at github
    2. https://arxiv.org/abs/1710.07629,
       "OpenFermion: The Electronic Structure Package for Quantum Computers"


    Parameters
    ----------
    file_prefix : str
    num_bits :  int
    qasm_gates : list[str]

    Returns
    -------
    None

    """
    emb = CktEmbedder(num_bits, num_bits)
    wr = SEO_writer(file_prefix, emb)
    for gate in qasm_gates:
        gate_words = gate.split()
        if gate_words[0] == 'H':
            tar_bit_pos = int(gate_words[1])
            wr.write_one_bit_gate(tar_bit_pos, OneBitGates.had2)
        elif gate_words[0] == 'CNOT':
            control_bit = int(gate_words[1])
            target_bit = int(gate_words[2])
            wr.write_cnot(control_bit, target_bit)
        elif gate_words[0][0] == 'R':
            xyz = gate_words[0][1]
            tar_bit_pos = int(gate_words[2])
            rad_ang = float(gate_words[1])/2
            if xyz == 'x':
                axis = 1
            elif xyz == 'y':
                axis = 2
            elif xyz == 'z':
                axis = 3
            else:
                assert False
            wr.write_one_bit_gate(tar_bit_pos, OneBitGates.rot_ax,
                                  [rad_ang, axis])
        else:
            assert False
    wr.close_files()

if __name__ == "__main__":
    file_prefix = 'io_folder/openfermion_qasm_to_qbtr'
    num_bits = 4
    qasm_gates = ['H 0',
                  'CNOT 1 3',
                  'Rx 1.5707963267948966 3',
                  'Rx -1.5707963267948966 2',
                  'Rz  .6  2']
    openfermion_qasm_to_qbtr_files(file_prefix, num_bits, qasm_gates)


