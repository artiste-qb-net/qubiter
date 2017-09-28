import Utilities as ut
from SEO_writer import *
from HadamardTransform import *
from SEO_MatrixProduct import *


class DiagUnitarySEO_writer(SEO_writer):
    """
    A diagonal unitary (d-unitary) is square diagonal matrix whose diagonal
    entries are magnitude 1 complex numbers. Any d-unitary gate can be
    represent in an English file by a single line starting with DIAG.

    This class is a child of SEO_writer. It adds to its parent class the
    ability to write a d-unitary gate in various styles.

    When style = 'one_line', this class writes a d-unitary as a single line
    in an English file, the same away the parent class SEO_writer would
    write it.

    When style = 'exact', this class writes a d-unitary as a SEO expansion
    occupying multiple lines of an English file. For this style, the class
    writes an exact expansion described in Refs. 1 and 2 below.

    Actually, this class can write more than a mere d-unitary. It can write
    a controlled d-unitary too, meaning it can attach T or F controls to the
    intrinsic controls of the d-unitary. Those intrinsic controls are
    represented by percent signs in Picture files and by "half moon" nodes
    in the arxiv papers cited below. In class Controls, T controls are of
    kind True, F controls are of kind False, and intrinsic controls are
    assigned an int for kind.

    It is important to note that the bits in the d-unitary being written
    will be in the following order, in order of increasing bit position:

    1. T controls
    2. F controls
    3. intrinsic d-unitary controls
    4. grounded bits, if any

    References
    ----------
    1. R.R. Tucci, A Rudimentary Quantum Compiler(2cnd Ed.)
    https://arxiv.org/abs/quant-ph/9902062

    2. R.R. Tucci, How to Compile Some NAND Formula Evaluators,
    https://arxiv.org/abs/0706.0479

    3. R.R. Tucci, Oracular Approximation of Quantum Multiplexors and
    Diagonal Unitary Matrices, https://arxiv.org/abs/0901.3851

    Attributes
    ----------
    num_of_F_trols : int
        The number of False controls of the controlled d-unitary.
    num_gbits : int
        number of gbits (grounded bits). This is only needed for the
        oracular approximation. Grounded bits are extra ancilla bits that
        have been initialized to the ground state (state |0>).
    num_of_T_trols : int
        The number of True controls of the controlled d-unitary.
    rad_angles : list(float)
        list of angles in radians. These angles are the parameters
        specifying an DIAG gate. If the DIAG has N intrinsic controls,
        there are 2^N angles.
    style : str
        must equal either 'one_line' or exact'.

    """

    def __init__(self, file_prefix, emb, style, rad_angles=None,
                 num_T_trols=0, num_F_trols=0, num_gbits=0, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        emb : CktEmbedder
        style : str
        rad_angles : list(float)
        num_T_trols : int
        num_F_trols : int
        num_gbits : int
        kwargs : dict()

        Returns
        -------
        None

        """
        self.style = style
        self.rad_angles = rad_angles
        if rad_angles:
            self.rad_angles = ut.centered_rads1(rad_angles)
        self.num_T_trols = num_T_trols
        self.num_F_trols = num_F_trols
        self.num_gbits = 0
        if style == 'oracular':
            self.num_gbits = num_gbits

        num_bits = emb.num_bits_bef
        assert num_bits >= 1, "d-unitary must have at least 1 qubit"

        ntf = num_T_trols + num_F_trols
        num_MP_trols = num_bits - ntf - num_gbits
        assert num_MP_trols > 0
        if rad_angles:
            assert len(rad_angles) == (1 << num_MP_trols), \
                "wrong  number of d-unitary angles"

        SEO_writer.__init__(self, file_prefix, emb, **kwargs)

    def write_one_line(self):
        """
        Writes in English file a one line representation of the d-unitary.

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        nt = self.num_T_trols
        nf = self.num_F_trols
        ntf = nt + nf
        num_MP_trols = num_bits - ntf - self.num_gbits
        trols = Controls(num_bits)
        trols.bit_pos_to_kind = dict(enumerate(
                [True]*nt + [False]*nf + list(range(num_MP_trols))
            ))
        trols.refresh_lists()

        self.write_controlled_diag_unitary_gate(trols, self.rad_angles)

    def write_exact(self):
        """
        Writes in English file a multiple line, exact representation of the
        d-unitary.

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        nt = self.num_T_trols
        nf = self.num_F_trols
        ntf = nt + nf
        num_MP_trols = num_bits - ntf - self.num_gbits
        rads_arr = np.array(ut.centered_rads1(self.rad_angles))
        if np.linalg.norm(rads_arr) < 1e-6:
            print("unit d-unitary")
            return

        conj_rads = HadamardTransform.ht(num_MP_trols, rads_arr)
        num_factors = (1 << num_MP_trols)
        f, lazy = BitVector.lazy_advance(0, 0)  # start at f=1
        cur_rot_bpos = 0
        prev_rot_bpos = 0
        cur_bvec = BitVector(num_MP_trols+1, 1)  # start at 1
        prev_bvec = BitVector(num_MP_trols+1, 0)
        diff_bvec = BitVector(num_MP_trols+1, 0)

        TF_dict = dict(enumerate([True]*nt + [False]*nf))
        trols1 = Controls(num_bits)
        trols1.bit_pos_to_kind = TF_dict.copy()
        trols1.refresh_lists()
        trols2 = Controls(num_bits)

        def write_cnots(diff_bvec, init_prev_T_bit):
            prev_T_bit = init_prev_T_bit
            while True:
                cur_T_bit = diff_bvec.find_T_bit_to_left_of(prev_T_bit)
                if cur_T_bit == -1:
                    break
                trols2.bit_pos_to_kind = TF_dict.copy()
                trols2.bit_pos_to_kind[cur_T_bit + ntf] = True
                trols2.refresh_lists()
                self.write_controlled_one_bit_gate(
                    ntf + init_prev_T_bit, trols2, OneBitGates.sigx)
                prev_T_bit = cur_T_bit

        norma = np.power(np.sqrt(2), num_MP_trols)
        # for first A factor, f = 0, just global phase
        # write conditioned global phase
        global_ph = conj_rads[0]*norma/len(conj_rads)
        if abs(global_ph) > 1e-6:
            self.write_controlled_one_bit_gate(ntf, trols1,
                    OneBitGates.phase_fac, [global_ph])

        while f < num_factors:
            cur_bvec.dec_rep = lazy
            # Since we have excluded f=0, f always has at least one T bit.
            cur_rot_bpos = cur_bvec.find_rightmost_T_bit()
            # print(cur_bvec.get_bit_string(), cur_rot_bpos)
            rads = ut.centered_rads(conj_rads[cur_bvec.dec_rep]/norma)
            if abs(rads) < 1e-6:
                pass
            else:
                # If cur_rot_bpos equals (doesn't equal) prev_rot_bpos,
                # then there is (isn't) cancellation between:
                # (1)the c-nots sigma_x(cur_rot_bpos)^n()
                # contributed by the right part of the current A factor
                # and
                # (2)the c-nots sigma_x(prev_rot_bpos)^n()
                # contributed by the left part of the previous A factor.

                if cur_rot_bpos == prev_rot_bpos:
                    diff_bvec = BitVector.new_with_T_on_diff(
                        cur_bvec, prev_bvec)
                    write_cnots(diff_bvec, cur_rot_bpos)
                else:
                    write_cnots(prev_bvec, prev_rot_bpos)
                    write_cnots(cur_bvec, cur_rot_bpos)
                    diff_bvec = BitVector.copy(cur_bvec)

                self.write_controlled_one_bit_gate(
                    ntf + cur_rot_bpos, trols1, OneBitGates.rot_ax, [rads, 3])
                prev_bvec = BitVector.copy(cur_bvec)
                prev_rot_bpos = cur_rot_bpos

            f, lazy = BitVector.lazy_advance(f, lazy)

        # Don't forget the leftmost c-nots
        write_cnots(prev_bvec, prev_rot_bpos)

    def write(self):
        """
        Main write function of this class. All other write functions are
        internal. This function writes a d-unitary in the style specified by
        the parameter self.style.

        Returns
        -------
        None

        """
        if all([abs(ang) < 1e-6 for ang in self.rad_angles]):
            print("unit d-unitary detected")
            return None

        if self.style == 'one_line':
            self.write_one_line()
        elif self.style == 'exact':
            self.write_exact()
        else:
            assert False, "unsupported d-unitary expansion style"

    @staticmethod
    def du_mat(rad_angles, herm_conj=False):
        """
        This function returns a square numpy array whose diagonal is the
        component-wise exp(1j* ) of rad_angles.

        Parameters
        ----------
        rad_angles : list(float)
        herm_conj : bool
            If True, uses exp(-1j*rad_angles).
            If False, uses exp(+1j*rad_angles)

        Returns
        -------
        np.ndarray

        """
        sign = 1
        if herm_conj:
            sign = -1
        return np.diag(np.exp(1j*sign*np.array(rad_angles)))

if __name__ == "__main__":
    nt = 1
    nf = 2
    num_MP_trols = 3
    num_angles = (1 << num_MP_trols)
    rad_angles = list(np.random.rand(num_angles)*2*np.pi)

    for style in ['one_line', 'exact']:
        num_gbits = 0
        if style == 'oracular':
            num_gbits = 3
        num_bits = nt + nf + num_MP_trols + num_gbits
        emb = CktEmbedder(num_bits, num_bits)
        file_prefix = "../io_folder/d_unitary_test_" + style
        wr = DiagUnitarySEO_writer(file_prefix, emb, style, rad_angles,
            num_T_trols=nt, num_F_trols=nf, num_gbits=num_gbits)
        wr.write()
        wr.close_files()

    file_prefix = "../io_folder/d_unitary_exact_check"
    num_bits = 4
    num_angles = (1 << num_bits)
    emb = CktEmbedder(num_bits, num_bits)
    rad_angles = list(np.random.rand(num_angles)*2*np.pi)
    # av = sum(rad_angles)/len(rad_angles)
    # rad_angles = list(np.array(rad_angles)-av)
    wr = DiagUnitarySEO_writer(file_prefix, emb, 'exact', rad_angles)
    wr.write()
    wr.close_files()
    matpro = SEO_MatrixProduct(file_prefix, num_bits)
    exact_mat = DiagUnitarySEO_writer.du_mat(rad_angles)
    print(np.linalg.norm(matpro.prod_arr - exact_mat))
    # print(matpro.prod_arr)
    # print(np.diag(exact_mat))
