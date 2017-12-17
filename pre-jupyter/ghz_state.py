from SEO_writer import *
from SEO_simulator import *
import numpy as np


def write_ghz_plus(file_prefix, ghz_only=True, meas=None):
    num_bits = 3
    z_axis = 3
    emb = CktEmbedder(num_bits, num_bits)
    print('-------------------', file_prefix)
    wr = SEO_writer(file_prefix, emb)
    wr.write_one_bit_gate(0, OneBitGates.had2)
    wr.write_one_bit_gate(1, OneBitGates.had2)
    wr.write_one_bit_gate(2, OneBitGates.sigx)

    control_pos = 0
    target_pos = 2
    trols = Controls.new_knob(num_bits, control_pos, kind=True)
    wr.write_controlled_one_bit_gate(
        target_pos, trols, OneBitGates.sigx)

    control_pos = 1
    target_pos = 2
    trols = Controls.new_knob(num_bits, control_pos, kind=True)
    wr.write_controlled_one_bit_gate(
        target_pos, trols, OneBitGates.sigx)

    wr.write_one_bit_gate(0, OneBitGates.had2)
    wr.write_one_bit_gate(1, OneBitGates.had2)
    wr.write_one_bit_gate(2, OneBitGates.had2)

    if not ghz_only:
        for pos in range(3):
            m = meas[pos]
            if m == 1:
                wr.write_one_bit_gate(pos, OneBitGates.had2)
            elif m == 2:
                wr.write_one_bit_gate(pos,
                    OneBitGates.rot_ax, [np.pi/4, z_axis])  # S^\dagger(pos)
                wr.write_one_bit_gate(pos, OneBitGates.had2)
            else:
                assert False
    wr.close_files()
    pic_file = file_prefix + '_' + str(num_bits) + '_ZLpic.txt'
    with open(pic_file) as f:
        print(f.read())
    init_st_vec = SEO_simulator.get_standard_basis_st_vec([0, 0, 0])
    sim = SEO_simulator(file_prefix, num_bits, init_st_vec)
    sim.describe_fin_st(print_st_vec=True, do_pp=True,
                        omit_zero_amps=True)
    fin_st_vec = sim.cur_st_vec_list[0]
    print('Prob(bit0=i, bit1=j, bit2=k) for i, j,k,=0,1:')
    prob_arr = np.abs(fin_st_vec)**2
    print(prob_arr)
    mean = prob_arr[0, 0, 0]  \
            + prob_arr[0, 1, 1] \
            + prob_arr[1, 0, 1] \
            + prob_arr[1, 1, 0] \
            - prob_arr[1, 1, 1] \
            - prob_arr[0, 0, 1] \
            - prob_arr[1, 0, 0] \
            - prob_arr[0, 1, 0]
    print('mean=', mean)
    return mean

# sigz(0)sigz(1)sigz(2) measurement
file_prefix = 'io_folder/ghz_zzz_meas'
mean_zzz = write_ghz_plus(file_prefix, ghz_only=True)

# sigy(0)sigy(1)sigx(2) measurement
file_prefix = 'io_folder/ghz_yyx_meas'
mean_yyx = write_ghz_plus(file_prefix, ghz_only=False, meas=[2, 2, 1])

# sigy(0)sigx(1)sigy(2) measurement
file_prefix = 'io_folder/ghz_yxy_meas'
mean_yxy = write_ghz_plus(file_prefix, ghz_only=False, meas=[2, 1, 2])

# sigx(0)sigy(1)sigy(2) measurement
file_prefix = 'io_folder/ghz_xyy_meas'
mean_xyy = write_ghz_plus(file_prefix, ghz_only=False, meas=[1, 2, 2])

# sigx(0)sigx(1)sigx(2) measurement
file_prefix = 'io_folder/ghz_xxx_meas'
mean_xxx = write_ghz_plus(file_prefix, ghz_only=False, meas=[1, 1, 1])


print('-----------------------')
print('mean_yyx =', mean_yyx)
print('mean_yxy =', mean_yxy)
print('mean_xyy =', mean_xyy)
print('mean_xxx =', mean_xxx)