from SEO_writer import *
from SEO_simulator import *
import numpy as np

smat = np.matrix([[1, 0], [0, 1j]])
tmat = np.matrix([[1, 0], [0, np.exp(1j*np.pi/4)]])
had = np.matrix([[1, 1], [1, -1]])/np.sqrt(2)
sigx = np.matrix([[0, 1], [1, 0]])
sigy = np.matrix([[0, -1j], [1j, 0]])
sigz = np.matrix([[1, 0], [0, -1]])
id2 = np.matrix([[1, 0], [0, 1]])

def exp_mat2(theta, vec4):
    # vec4 is 4 dimensional np.matrix. Zero component not used.
    unit_vec = np.array([0, vec4[1], vec4[2], vec4[3]])
    unit_vec = unit_vec/np.linalg.norm(unit_vec)

    mat = unit_vec[1]*sigx + unit_vec[2]*sigy + unit_vec[3]*sigz
    return np.cos(theta)*id2 + 1j*mat*np.sin(theta)

roty = exp_mat2(np.pi/8, np.array([0, 0, 1, 0]))
sigw = (sigx + sigz)/np.sqrt(2)
sigv = (-sigx + sigz)/np.sqrt(2)

print(np.linalg.norm(sigw - roty.getH()*sigz*roty))
print(np.linalg.norm(sigv - roty*sigz*roty.getH()))

roty1 = np.exp(-1j*np.pi/8)*smat.getH()*had*tmat*had*smat
print(np.linalg.norm(roty - roty1))


def write_zz_plus(file_prefix, zz_only=True, extra_had=False, t_herm=False):
    num_bits = 2
    z_axis = 3
    emb = CktEmbedder(num_bits, num_bits)
    print('-------------------', file_prefix)
    wr = SEO_writer(file_prefix, emb)
    wr.write_one_bit_gate(0, OneBitGates.had2)

    control_pos = 0
    target_pos = 1
    trols = Controls.new_knob(num_bits, control_pos, kind=True)
    wr.write_controlled_one_bit_gate(
        target_pos, trols, OneBitGates.sigx)

    if not zz_only:
        if extra_had:
            wr.write_one_bit_gate(0, OneBitGates.had2)

        wr.write_one_bit_gate(1, OneBitGates.rot_ax, [-np.pi/4, z_axis]) # S(1)
        wr.write_one_bit_gate(1, OneBitGates.had2)   # H(1)
        if t_herm:
            pm_one = -1
        else:
            pm_one = 1
        wr.write_one_bit_gate(1, OneBitGates.rot_ax,
                              [-pm_one*np.pi/8, z_axis]) # T(1) if pm_one=1
        wr.write_one_bit_gate(1, OneBitGates.had2)   # H(1)
    wr.close_files()
    pic_file = file_prefix + '_' + str(num_bits) + '_ZLpic.txt'
    with open(pic_file) as f:
        print(f.read())
    init_st_vec = SEO_simulator.get_standard_basis_st_vec([0, 0])
    sim = SEO_simulator(file_prefix, num_bits, init_st_vec)
    sim.describe_fin_st(print_st_vec=True, do_pp=True,
                        omit_zero_amps=True)
    fin_st_vec = sim.cur_st_vec_list[0]
    print('Prob(bit0=j, bit1=k) for j,k=0,1:')
    prob_arr = np.abs(fin_st_vec)**2
    print(prob_arr)
    mean = prob_arr[0, 0] \
           + prob_arr[1, 1] \
           - prob_arr[0, 1] \
           - prob_arr[1, 0]
    print('<siga(0)sigb(1)>=', mean)
    return mean

# sigz(0)sigz(1) measurement
file_prefix = 'io_folder/bell_zz_meas'
mean_zz = write_zz_plus(file_prefix, zz_only=True)

# sigz(0)sigw(1) measurement
file_prefix = 'io_folder/bell_zw_meas'
mean_zw = write_zz_plus(file_prefix, zz_only=False, extra_had=False, t_herm=False)

# sigz(0)sigv(1) measurement
file_prefix = 'io_folder/bell_zv_meas'
mean_zv = write_zz_plus(file_prefix, zz_only=False, extra_had=False, t_herm=True)

# sigx(0)sigw(1) measurement
file_prefix = 'io_folder/bell_xw_meas'
mean_xw = write_zz_plus(file_prefix, zz_only=False, extra_had=True, t_herm=False)

# sigx(0)sigv(1) measurement
file_prefix = 'io_folder/bell_xv_meas'
mean_xv = write_zz_plus(file_prefix, zz_only=False, extra_had=True, t_herm=True)


print('-----------------------')
print('mean_zw + mean_zv + mean_xw - mean_xv - 2*np.sqrt(2)=',
      mean_zw + mean_zv + mean_xw - mean_xv - 2*np.sqrt(2))