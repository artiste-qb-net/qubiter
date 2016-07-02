from SEO_writer import *
from SEO_simulator import *
import pprint as pp
import numpy as np

"""This script is meant to illustrate how to use Qubiter to simulate ( i.e.,
predict) the outcome of an IBM Quantum Experience experiment """

# number of qubits is 5
# note that we use "bit" for both qbits and cbits
num_bits = 5
# use a trivial circuit embedder that embeds 5 qubits into same 5 qubits
emb = CktEmbedder(num_bits, num_bits)
# open a writer, tell it where to write to.
# We will use zero bit last (ZL) convention
wr = SEO_writer('io_folder//ibm_qe_test', emb)

# write Pauli matrices sigx, sigy, sigz at position 2
wr.write_one_bit_gate(2, OneBitGates.sigx)
wr.write_one_bit_gate(2, OneBitGates.sigy)
wr.write_one_bit_gate(2, OneBitGates.sigz)
# write 1 qubit Hadamard matrix at position 3
wr.write_one_bit_gate(3, OneBitGates.had2)

z_axis = 3
# 'pos' means position
# write exp(i*ang_rads*sigz) at pos=2 with ang_rads = pi/2, -pi/2, i.e.
# sqrt(sigz), what they call S and S^\dagger
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [np.pi/2, z_axis])
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [-np.pi/2, z_axis])
# rite exp(i*ang_rads*sigz) at pos=2 with ang_rads = pi/4, -pi/4, i.e.
# sqrt sqrt (sigz), what they call T and T^\dagger
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [np.pi/4, z_axis])
wr.write_one_bit_gate(2, OneBitGates.rot_ax, [-np.pi/4, z_axis])

# write CNOT = sigx(target_pos)^n(control_pos)
control_pos = 3
target_pos = 0
trols = Controls.new_knob(num_bits, control_pos, kind=True)
wr.write_controlled_one_bit_gate(
    target_pos, trols, OneBitGates.sigx)

wr.close_files()

# specify initial state vector for simulation
# this example corresponds to |0>|0>|1>|1>|0>
init_st_vec = SEO_simulator.get_standard_basis_st([0, 0, 1, 1, 0])
# open a simulator. This automatically
# multiplies quantum circuit in given file.
sim = SEO_simulator('io_folder//ibm_qe_test', num_bits, init_st_vec)

# the final state vector is the first in a list
fin_st_vec = sim.cur_st_vec_list[0]

print('total probability of final state vector (should be one)=',
    sim.get_total_prob(fin_st_vec))

print('dictionary with key=qubit, value=final (P(0), P(1))')
pp.pprint(sim.get_bit_probs(fin_st_vec))
