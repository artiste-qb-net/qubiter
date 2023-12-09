from qubiter.SEO_simulator import *
from qubiter.StateVec import *
# from Controls import *


class SEO_simulator_mp(SEO_simulator):
    """
    This class is a mp (matrix product) child of class SEO_simulator. An
    object of this class is created inside class SEO_MatrixProduct. The
    purpose of this class is to override the use_MEAS() function of its
    parent class SEO_simulator. The new use_MEAS() returns an error message
    if an English file with a MEAS line is being read.
    """

    def use_MEAS(self, tar_bit_pos, kind):
        """
        Returns an error message if an English file with a MEAS line is
        being read.

        Parameters
        ----------
        tar_bit_pos : int
        kind : int|bool

        Returns
        -------
        None

        """
        assert False, "this lightweight class cannot handle MEAS"\
            "or its concomitant branches"

    def use_PRINT(self, style, line_num):
        """
        If circuit has any PRINT statements, skip them.

        Parameters
        ----------
        style : str
        line_num : int

        Returns
        -------
        None

        """
        pass


class SEO_MatrixProduct:
    """
    The class SEO_simulator has an initial state vector as an input, and it
    calculates the evolution of that state vector after each line of an
    English file. This class, on the other hand, has no input (initial) nor
    output (evolved) state vectors. Instead, the class calculates the
    product of the matrices corresponding to each line (gate) of an English
    file.

    In order to accomplish this goal, this class calls SEO_simulalor_mp
    repeatedly using as initial state vector all the standard basis vectors
    (2^num_qbits of them). Then the class assembles the product matrix that
    we seek by stacking on top of each other all the 2^num_qbits final
    evolved state vectors. Admittedly, this is a very slow, inefficient way
    of finding the sought for matrix product. However, it works fine for a
    small number of qubits. It can be used to check that gate expansions
    agree with what they are supposed to be an expansion of.

    Though slow, the strategy used by this class has the advantages that it
    was very easy to code, and that every time it is shown to work,
    SEO_simulator is shown to work too, 2^num_bit times. This strategy is
    also easier to parallelize if the need for that arises at some later date.

    Attributes
    ----------
    prod_arr : np.ndarray
        the product matrix obtained by multiplying each line of the input
        English file.

    """

    def __init__(self, file_prefix, num_qbits):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
            Prefix of English file being read
        num_qbits : int
            number of bits in English file begin read.

        Returns
        -------


        """

        self.prod_arr = None

        num_comps = 1 << num_qbits
        fin_list = []
        for s in range(num_comps):
            spin_dir_list = [(s >> bpos) & 1 for bpos in
                             reversed(range(num_qbits))]
            init_st_vec = StateVec.get_standard_basis_st_vec(
                spin_dir_list, ZL=True)
            # print("---", spin_dir_list)
            # print(init_st_vec)
            sim = SEO_simulator_mp(file_prefix, num_qbits,
                                   init_st_vec=init_st_vec)
            fin_st_vec = sim.cur_st_vec_dict["pure"]
            fin = StateVec.get_traditional_st_vec
            fin_list.append(fin)
            # print(fin_st_vec)
            # print(fin)
        self.prod_arr = np.vstack(fin_list).transpose()
        # print(self.prod_arr)


if __name__ == "__main__":
    from qubiter.FouSEO_writer import *

    def main():
        num_qbits = 3
        emb = CktEmbedder(num_qbits, num_qbits)
        file_prefix = 'matrix_prod_test'
        wr = FouSEO_writer(True, file_prefix, emb)
        wr.close_files()
        mp = SEO_MatrixProduct(file_prefix, num_qbits)
        prod = mp.prod_arr
        exact = FouSEO_writer.fourier_trans_mat(1 << num_qbits)
        # print(prod)
        # print(exact)
        err = np.linalg.norm(prod - exact)
        print("error=", err)
    main()
