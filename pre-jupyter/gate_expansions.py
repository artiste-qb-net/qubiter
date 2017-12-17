from CktExpander import *
from SEO_writer import *

file_prefix = 'io_folder/expansions_examples'
num_bits = 3
emb = CktEmbedder(num_bits, num_bits, range(num_bits))
wr = SEO_writer(file_prefix, emb)

trols1 = Controls(num_bits)
trols1.bit_pos_to_kind = {0: True}
trols1.refresh_lists()

trols2 = Controls(num_bits)
trols2.bit_pos_to_kind = {0: True, 1: False}
trols2.refresh_lists()

wr.write_NOTA('simple cnot ( sigx(1)^n(0) )')
wr.write_controlled_one_bit_gate(1, trols1, OneBitGates.sigx)

wr.write_NOTA('controlled sigy ( sigy(1)^n(0) )')
wr.write_controlled_one_bit_gate(1, trols1, OneBitGates.sigy)

# let exp(i*(radx*sigx + rady*sigy + radz*sigz)
wr.write_NOTA('controlled Y,Z rotation ( rot(1)^n(0) with radx = 0 )')
wr.write_controlled_one_bit_gate(1, trols1, OneBitGates.rot,
                                 [0.0, np.pi/3, np.pi/4])

wr.write_NOTA('controlled rotation ( rot(1)^n(0) )')
wr.write_controlled_one_bit_gate(1, trols1, OneBitGates.rot,
                                 [np.pi/5, np.pi/3, np.pi/4])

wr.write_NOTA('2-controlled not (sigx(2)^(nbar(1)n(0))')
wr.write_controlled_one_bit_gate(2, trols2, OneBitGates.sigx)

wr.write_NOTA('swap of 0 and 1')
wr.write_bit_swap(0, 1)

wr.write_NOTA('swap of 1 and 2 controlled by 0')
wr.write_controlled_bit_swap(1, 2, trols1)

wr.close_files()


CktExpander(file_prefix, num_bits)
SEO_reader(file_prefix + '_X1', num_bits)