from FouSEO_writer import *
from quantum_CSD_compiler.Tree import *
from quantum_CSD_compiler.DiagUnitarySEO_writer import *


num_bits = 3
init_unitary_mat = FouSEO_writer.fourier_trans_mat(1 << num_bits)
emb = CktEmbedder(num_bits, num_bits)
file_prefix = './io_folder/csd_test'
t = Tree(True, file_prefix, emb, init_unitary_mat, verbose=False)
t.close_files()
file = file_prefix + '_3_ZLpic.txt'
with open(file) as f:
    print(f.read())

file = file_prefix + '_3_eng.txt'
with open(file) as f:
    print(f.read())

file_prefix = "./io_folder/d_unitary_exact_check"
num_bits = 4
num_angles = (1 << num_bits)
emb = CktEmbedder(num_bits, num_bits)
rad_angles = list(np.random.rand(num_angles)*2*np.pi)
wr = DiagUnitarySEO_writer(file_prefix, emb, 'exact', rad_angles)
wr.write()
wr.close_files()
file = file_prefix + '_4_ZLpic.txt'
with open(file) as f:
    print(f.read())

matpro = SEO_MatrixProduct(file_prefix, num_bits)
exact_mat = DiagUnitarySEO_writer.du_mat(rad_angles)
err = np.linalg.norm(matpro.prod_arr - exact_mat)
print("diag unitary error=", err)

file_prefix = "./io_folder/plexor_exact_check"
num_bits = 4
num_angles = (1 << (num_bits-1))
emb = CktEmbedder(num_bits, num_bits)
rad_angles = list(np.random.rand(num_angles)*2*np.pi)
wr = MultiplexorSEO_writer(file_prefix, emb, 'exact', rad_angles)
wr.write()
wr.close_files()
file = file_prefix + '_4_ZLpic.txt'
with open(file) as f:
    print(f.read())

matpro = SEO_MatrixProduct(file_prefix, num_bits)
exact_mat = MultiplexorSEO_writer.mp_mat(rad_angles)
err = np.linalg.norm(matpro.prod_arr - exact_mat)
print("multiplexor error=", err)

