from SEO_writer import *
from SEO_simulator import *
import numpy as np

num_bits = 3
emb = CktEmbedder(num_bits, num_bits)
file_prefix = 'io_folder/teleportation-fedor'
wr = SEO_writer(file_prefix, emb)

z_axis = 3
def init_ckt():
    wr.write_one_bit_gate(0, OneBitGates.had2)  # H(0)
    wr.write_one_bit_gate(0, OneBitGates.rot_ax, [-np.pi/8, z_axis])  # T(0)
    wr.write_one_bit_gate(0, OneBitGates.had2)  # H(0)
    wr.write_one_bit_gate(0, OneBitGates.rot_ax, [-np.pi/4, z_axis])  # S(0)
init_ckt()
wr.close_files()

init_st_vec = SEO_simulator.get_standard_basis_st_vec([0, 0, 0])
sim = SEO_simulator(file_prefix, num_bits, init_st_vec)
sim.describe_fin_st(print_st_vec=True)

wr = SEO_writer(file_prefix, emb)
init_ckt()
wr.write_one_bit_gate(2, OneBitGates.had2)  # H(2)

control_pos = 2
target_pos = 1
trols = Controls.new_knob(num_bits, control_pos, kind=True)
wr.write_controlled_one_bit_gate(target_pos, trols, OneBitGates.sigx)

control_pos = 0
target_pos = 1
trols = Controls.new_knob(num_bits, control_pos, kind=True)
wr.write_controlled_one_bit_gate(target_pos, trols, OneBitGates.sigx)

wr.write_one_bit_gate(0, OneBitGates.had2)  # H(0)

wr.write_MEAS(tar_bit_pos=0, kind=0)
wr.write_MEAS(tar_bit_pos=1, kind=0)
wr.close_files()

init_st_vec = SEO_simulator.get_standard_basis_st_vec([1, 1, 0])
sim = SEO_simulator(file_prefix, num_bits, init_st_vec)
sim.describe_fin_st(print_st_vec=True)




