"""
# Quantum Phase Estimation, Iterative
https://arxiv.org/abs/1512.06860
https://arxiv.org/abs/1605.03590



```
|       H
|       exp(i*alpha(k)*sigz)
U^k-----@
|       H
|       measure n(k) here
```


$\alpha(k+1) = \pi n(k) + 2\alpha(k)$

$\alpha(k) = 2\pi  2^{k-1}\sum_{b=0}^{k} \frac{n(b)}{2^{b}}$

$U = e^{i*rads*\sigma_Z}$
"""

from SEO_writer import *
from SEO_simulator import *
import numpy as np

rads = 2*np.pi*(1/16 + 1/8 + 1e-8)
z_axis = 3
num_bits = 2
num_reps = 15
file_prefix = 'io_folder/ph_est_iterative'

emb = CktEmbedder(num_bits, num_bits)

alpha = 0
ptr_state = 0
ptr_st_list = []
for k in range(num_reps):
    print('--------k=', k)
    # refresh angle alpha to twice its previous value plus
    # 2\pi times latest measurement of pointer qubit
    alpha = 2*alpha + np.pi*ptr_state/2
    print('rads, alpha/2^(num_reps)=', rads, alpha/(1 << num_reps-2))

    wr = SEO_writer(file_prefix, emb)

    # write circuit
    wr.write_one_bit_gate(0, OneBitGates.had2)

    wr.write_one_bit_gate(0, OneBitGates.rot_ax, [alpha, z_axis])

    control_pos = 0
    target_pos = 1
    trols = Controls.new_knob(num_bits, control_pos, kind=True)
    wr.write_controlled_one_bit_gate(
        target_pos, trols, OneBitGates.rot_ax, [(1 << k)*rads, z_axis])

    wr.write_one_bit_gate(0, OneBitGates.had2)
    wr.close_files()

    # simulate circuit
    init_st_vec = SEO_simulator.get_standard_basis_st_vec([0, 0])
    sim = SEO_simulator(file_prefix, num_bits, init_st_vec)
    sim.describe_fin_st(print_st_vec=True, do_pp=True, omit_zero_amps=True)

    # find final state of pointer qubit
    fin_st_vec = sim.cur_st_vec_list[0]
    # dictionary with key=qubit, value=final (P(0), P(1))
    bit_to_probs = sim.get_bit_probs_from_st_vec(fin_st_vec)
    p0, p1 = bit_to_probs[0]
    if p0 > p1:
        ptr_state = 0
    else:
        ptr_state = 1
    ptr_st_list.append(ptr_state)
    print('ptr_state=', ptr_state)
print('---------------------')
print('timeline of bit 0 measurements', ptr_st_list)
print("rads, alpha(num_reps-1)/2^(num_reps-2)", rads, alpha/(1 << num_reps-2))
