from SEO_writer import *
from SEO_simulator import *
from StateVec import *
from Plotter import *
import numpy as np

num_bits = 3
emb = CktEmbedder(num_bits, num_bits)
file_prefix = '../io_folder/teleportation-with-ifs'
wr = SEO_writer(file_prefix, emb)

ckt = 'ibm'
if ckt == 'ibm':
    #from IBM
    wr.write_one_bit_gate(0,
                          OneBitGates.rot,
                          list(np.pi/180*np.array([20, 68, 46])))
    wr.write_PRINT("ALL")

    wr.write_one_bit_gate(1, OneBitGates.had2)
    wr.write_cnot(control_bit=1, target_bit=2)
    #wr.write_one_bit_gate(0, OneBitGates.rot_ax, [-np.pi/8, 2])
    wr.write_cnot(control_bit=0, target_bit=1)
    wr.write_one_bit_gate(0, OneBitGates.had2)
    wr.write_PRINT("ALL")

    wr.write_MEAS(0, kind=2)
    wr.write_MEAS(1, kind=2)

    wr.write_PRINT("ALL")
    wr.write_IF_M_beg(Controls.new_knob(num_bits, 0, True))
    wr.write_one_bit_gate(2, OneBitGates.sigz)
    wr.write_IF_M_end()
    wr.write_PRINT("ALL")
    wr.write_IF_M_beg(Controls.new_knob(num_bits, 1, True))
    wr.write_one_bit_gate(2, OneBitGates.sigx)
    wr.write_IF_M_end()
    wr.write_PRINT("ALL")
elif ckt == 'met':
    #from Metcalf paper
    wr.write_one_bit_gate(0,
                          OneBitGates.rot,
                          list(np.pi/180*np.array([0, 0, 0])))
    wr.write_PRINT("ALL")

    wr.write_one_bit_gate(2, OneBitGates.had2)
    wr.write_cnot(control_bit=2, target_bit=1)
    wr.write_cnot(control_bit=0, target_bit=1)
    wr.write_one_bit_gate(0, OneBitGates.had2)
    wr.write_one_bit_gate(0, OneBitGates.sigz)
    wr.write_one_bit_gate(1, OneBitGates.sigz)
    wr.write_PRINT("ALL")

    wr.write_MEAS(0, kind=2)
    wr.write_MEAS(1, kind=2)

    wr.write_PRINT("ALL")
    wr.write_IF_M_beg(Controls.new_knob(num_bits, 0, True))
    wr.write_one_bit_gate(2, OneBitGates.sigz)
    wr.write_IF_M_end()
    wr.write_PRINT("ALL")
    wr.write_IF_M_beg(Controls.new_knob(num_bits, 1, True))
    wr.write_one_bit_gate(2, OneBitGates.sigx)
    wr.write_IF_M_end()
    wr.write_PRINT("ALL")
else:
    assert False

wr.close_files()

init_st_vec = StateVec.get_standard_basis_st_vec([0, 0, 0])
sim = SEO_simulator(file_prefix, num_bits, init_st_vec)

np.set_printoptions(precision=2, suppress=True)

# print(",,", sim.cached_sts[2]["pure"])
den_mat1 = StateVec.get_den_mat(num_bits, sim.cached_sts[2])
den_mat1_df = Plotter.get_den_mat_df(num_bits, den_mat1)
print("\nden_mat1=\n", den_mat1_df)

# print(",,.", sim.cached_sts[15]["0T1T"])
den_mat2 = StateVec.get_den_mat(num_bits, sim.cur_st_vec_dict)
den_mat2_df = Plotter.get_den_mat_df(num_bits, den_mat2)
print("\nden_mat2=\n", den_mat2_df)
# Plotter.plot_phasors(['den_mat2'], den_mat_df_list=[den_mat2_df])
print("impurity of den_mat2=",
      StateVec.get_impurity(den_mat2))

tr01_den_mat2 = StateVec.get_partial_tr(num_bits, den_mat2, {0, 1})
tr01_den_mat2_df = Plotter.get_den_mat_df(1, tr01_den_mat2)
print("\ntr01_den_mat2=\n", tr01_den_mat2_df)
print("impurity of tr01_den_mat2=",
      StateVec.get_impurity(tr01_den_mat2))
