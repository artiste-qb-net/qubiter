"""
    This script is meant to illustrate how to use Qubiter to simulate ( i.e.,
    predict) the outcome of an IBM Quantum Experience experiment
"""

from SEO_writer import *
from SEO_simulator import *
import numpy as np

# Number of qubits is 5.
# Note that we use "bit" for both qbits and cbits.
num_bits = 5

# Use a trivial circuit embedder that embeds 5 qubits into same 5 qubits.
emb = CktEmbedder(num_bits, num_bits)

# Open a writer, tell it where to write to.
# We will use zero bit last (ZL) convention
file_prefix = 'io_folder/ibm_qe_test'
wr = SEO_writer(file_prefix, emb)

# Write Pauli matrices sigx, sigy, sigz at position 2
wr.write_one_bit_gate(2, OneBitGates.sigx)
wr.write_one_bit_gate(2, OneBitGates.sigy)
wr.write_one_bit_gate(2, OneBitGates.sigz)

# write 1 qubit Hadamard matrix at position 3
wr.write_one_bit_gate(3, OneBitGates.had2)

# $$S = diag[1, i] = exp(i\pi/4) diag[exp(-i\pi/4), exp(i\pi/4)]
#     = exp(i\pi/4)exp(-i \pi/4 sigz)$$
# $$S^\dagger = exp(-i\pi/4)exp(+i\pi/4 sigz)$$
#
# $$T = \sqrt{S}= diag[1, exp(i\pi/4)] = exp(i\pi/8)exp(-i\pi/8sigz)$$
# $$T^\dagger = exp(-i\pi/8)exp(+i\pi/8 sigz)$$
#
# Write S, S^\dagger, T, T^dagger at position=2, up to a global phase factor
z_axis = 3
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [-np.pi/4, z_axis]) # S(2)
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [+np.pi/4, z_axis])
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [-np.pi/8, z_axis]) # T(2)
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [+np.pi/8, z_axis])

# write CNOT = sigx(target_pos)^n(control_pos)
control_pos = 3
target_pos = 0
trols = Controls.new_knob(num_bits, control_pos, kind=True)
wr.write_controlled_one_bit_gate(
    target_pos, trols, OneBitGates.sigx)

# Close English and Picture files
wr.close_files()

# Specify initial state vector for simulation.
# This example corresponds to |0>|0>|1>|1>|0>.
# In ZL convention, last ket corresponds to bit 0.
init_st_vec = SEO_simulator.get_standard_basis_st([0, 0, 1, 1, 0])

# Open a simulator. This automatically
# multiplies quantum circuit in given file.
sim = SEO_simulator(file_prefix, num_bits, init_st_vec)

# Print description of final state vector
sim.describe_fin_st(print_st_vec=False)
