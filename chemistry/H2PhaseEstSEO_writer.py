from PhaseEstSEO_writer import *
from chemistry.H2EvolOpData import *
from chemistry.MolEvolOpSEO_writer import *

"""

This script calls the constructor PhaseEstSEO_writer() to write a phase
estimation circuit for estimating the ground state of an H2 molecule.
Various other support classes must be called before PhaseEstSEO_writer can
be called.

"""


if __name__ == "__main__":
    data = H2EvolOpData()
    data.set_arg_list()

    num_probe_bits = 3
    num_orbitals = 4
    num_bits = num_probe_bits + num_orbitals + 1
    atom_writer_cls = MolEvolOpSEO_writer
    cls_arg_list = data.get_arg_list()
    file_prefix = "chem_io_folder//H2_ground_state"
    emb = CktEmbedder(num_bits, num_bits)
    PhaseEstSEO_writer(True, num_probe_bits,
        atom_writer_cls, cls_arg_list, file_prefix, emb)
