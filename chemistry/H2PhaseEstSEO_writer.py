from PhaseEstSEO_writer import *
from chemistry.H2EvolOpData import *
from chemistry.MolEvolOpSEO_writer import *

"""

This script calls the constructor PhaseEstSEO_writer() to write a phase
estimation circuit for estimating the ground state energy of an H2 molecule.
Various other support classes must be called before PhaseEstSEO_writer can
be called.

"""


if __name__ == "__main__":

    num_probe_bits = 3
    num_orbitals = 4
    num_bits = num_probe_bits + num_orbitals + 1

    file_prefix = "chem_io_folder/H2_ground_state"
    emb = CktEmbedder(num_bits, num_bits)
    atom_wr = MolEvolOpSEO_writer(data=H2EvolOpData(),
                            do_write=False)
    wr = PhaseEstSEO_writer(do_write=True,
                            num_probe_bits=num_probe_bits,
                            atom_writer=atom_wr,
                            file_prefix=file_prefix,
                            emb=emb)
