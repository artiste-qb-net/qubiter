import Utilities as ut
from SEO_writer import *
from HadamardTransform import *
from SEO_MatrixProduct import *


class MultiplexorSEO_writer(SEO_writer):
    """
    A multiplexor is a matrix of the form

        [ cc, ss]
        [-ss, cc]

    where cc and ss are square diagonal matrices, both of the same
    dimension, satisfying cc^2 + ss^2=1.

    A multiplexor gate can be represent in an English file by a single line
    starting with MP_Y.

    This class is a child of SEO_writer. It adds to its parent class the
    ability to write a multiplexor gate in various styles.

    When style = 'one_line', this class writes a multiplexor as a single
    line in an English file, the same away the parent class SEO_writer would
    write it.

    When style = 'exact', this class writes a multiplexor as a SEO expansion
    occupying multiple lines of an English file. For this style, the class
    writes an exact expansion described in Refs. 1 and 2 below.

    When style = 'oracular', this class writes a multiplexor as a SEO
    expansion occupying multiple lines of an English file. For this style,
    the class writes an approximation (called the "oracular approximation")
    described in Ref.3 below.

    Actually, this class can write more than a mere multiplexor. It can
    write a controlled multiplexor too, meaning it can attach T or F
    controls to the intrinsic controls of the multiplexor. Those intrinsic
    controls are represented by percent signs in Picture files and by "half
    moon" nodes in the arxiv papers cited below. In class Controls,
    T controls are of kind True, F controls are of kind False, and intrinsic
    controls are assigned an int for kind.

    It is important to note that the bits in the multiplexor being written
    will be in the following order, in order of increasing bit position:

    1. T controls
    2. F controls
    3. intrinsic multiplexor controls
    4. target bit of multiplexor
    5. grounded bits, if any

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
        The number of False controls of the controlled multiplexor
    num_gbits : int
        number of gbits (grounded bits). This is only needed for the
        oracular approximation. Grounded bits are extra ancilla bits that
        have been initialized to the ground state (state |0>).
    num_of_T_trols : int
        The number of True controls of the controlled multiplexor
    rad_angles : list(float)
        list of angles in radians. These angles are the parameters
        specifying an MP_Y gate. If the MP_Y has N intrinsic controls,
        there are 2^N angles.
    style : str
        must equal either 'one_line', exact' or 'oracular'.

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
        assert num_bits >= 2, "multiplexor must have at least 2 qubits"

        ntf = num_T_trols + num_F_trols
        num_MP_trols = num_bits - ntf - num_gbits - 1
        assert num_MP_trols > 0
        if rad_angles:
            assert len(rad_angles) == (1 << num_MP_trols), \
                "wrong  number of multiplexor angles"

        SEO_writer.__init__(self, file_prefix, emb, **kwargs)

    def write_one_line(self):
        """
        Writes in English file a one line representation of the multiplexor.

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        nt = self.num_T_trols
        nf = self.num_F_trols
        ntf = nt + nf
        num_MP_trols = num_bits - ntf - self.num_gbits - 1
        trols = Controls(num_bits)
        trols.bit_pos_to_kind = dict(enumerate(
                [True]*nt + [False]*nf + list(range(num_MP_trols))
            ))
        trols.refresh_lists()

        tar_bit_pos = ntf + num_MP_trols
        self.write_controlled_multiplexor_gate(tar_bit_pos,
                trols, self.rad_angles)

    def write_exact(self):
        """
        Writes in English file a multiple line, exact representation of the
        multiplexor.

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        nt = self.num_T_trols
        nf = self.num_F_trols
        ntf = nt + nf
        num_MP_trols = num_bits - ntf - self.num_gbits - 1
        rads_arr = np.array(ut.centered_rads1(self.rad_angles))
        if np.linalg.norm(rads_arr) < 1e-6:
            print("unit multiplexor")
            return

        conj_rads = HadamardTransform.ht(num_MP_trols, rads_arr)
        num_factors = (1 << num_MP_trols)

        cur_bvec = BitVector(num_MP_trols+1, 0)  # start at zero
        prev_bvec = BitVector(num_MP_trols+1, 0)

        TF_dict = dict(enumerate([True]*nt + [False]*nf))
        trols1 = Controls(num_bits)
        trols1.bit_pos_to_kind = TF_dict.copy()
        trols1.refresh_lists()
        trols2 = Controls(num_bits)

        def write_cnots(diff_bvec):
            prev_T_bit = num_MP_trols
            while True:
                cur_T_bit = diff_bvec.find_T_bit_to_right_of(prev_T_bit)
                if cur_T_bit == -1:
                    break
                trols2.bit_pos_to_kind = TF_dict.copy()
                trols2.bit_pos_to_kind[cur_T_bit + ntf] = True
                trols2.refresh_lists()
                self.write_controlled_one_bit_gate(
                    ntf + num_MP_trols, trols2, OneBitGates.sigx)
                prev_T_bit = cur_T_bit

        norma = np.power(np.sqrt(2), num_MP_trols)
        f = 0
        lazy = 0
        while f < num_factors:
            rads = conj_rads[cur_bvec.dec_rep]/norma
            if abs(rads) < 1e-6:
                pass
            else:
                diff_bvec = BitVector.new_with_T_on_diff(cur_bvec, prev_bvec)
                write_cnots(diff_bvec)
                self.write_controlled_one_bit_gate(
                    ntf + num_MP_trols, trols1, OneBitGates.rot_ax, [rads, 2])
                prev_bvec = BitVector.copy(cur_bvec)
            f, lazy = BitVector.lazy_advance(f, lazy)
            cur_bvec.dec_rep = lazy

        # Don't forget the leftmost c-nots:
        diff_bvec = prev_bvec
        write_cnots(diff_bvec)

    def write_oracular(self):
        """
        Writes in English file a multiple line, approximate representation
        of the multiplexor.

        Returns
        -------

        """
        num_bits = self.emb.num_bits_bef
        nt = self.num_T_trols
        nf = self.num_F_trols
        ntf = nt + nf
        num_MP_trols = num_bits - ntf - self.num_gbits - 1
        bit_precision = self.num_gbits
        num_angles = len(self.rad_angles)
        ang_bools = [True]*num_angles

        TF_dict = dict(enumerate([True]*nt + [False]*nf))
        side_trols = Controls(num_bits)
        center_trols = Controls(num_bits)

        def write_omega(tar_bit_pos, ang_bools):
            for b in range(num_angles):
                if ang_bools[b]:
                    side_trols.bit_pos_to_kind = {
                        c + ntf: ((b >> c) & 1) == 1
                        for c in range(num_MP_trols)
                    }
                    side_trols.bit_pos_to_kind.update(TF_dict)
                    side_trols.refresh_lists()
                    self.write_controlled_one_bit_gate(tar_bit_pos, side_trols,
                                                       OneBitGates.sigx)

        bit_pos = ntf + num_MP_trols  # this is the target of Ry, gbits follow
        for k in range(1, bit_precision+1):
            bit_pos += 1
            center_trols.bit_pos_to_kind = {bit_pos: True}
            center_trols.bit_pos_to_kind.update(TF_dict)
            center_trols.refresh_lists()
            for b in range(num_angles):
                # self.rad_angles in [-\pi, \pi] by now
                fraction = (self.rad_angles[b] + np.pi)/(2*np.pi)
                x = int(fraction*(1 << k))  # keep only int part
                ang_bools[b] = (x & 1 == 1)
            write_omega(bit_pos, ang_bools)
            self.write_controlled_one_bit_gate(
                    ntf + num_MP_trols,
                    center_trols,
                    OneBitGates.rot_ax, [2*np.pi/(1 << k), 2])
            write_omega(bit_pos, ang_bools)

    def write(self):
        """
        Main write function of this class. All other write functions are
        internal. This function writes a multiplexor in the style specified
        by the parameter self.style.

        Returns
        -------
        None

        """
        if all([abs(ang) < 1e-6 for ang in self.rad_angles]):
            print("unit multiplexor detected")
            return None

        if self.style == 'one_line':
            self.write_one_line()
        elif self.style == 'exact':
            self.write_exact()
        elif self.style == 'oracular':
            self.write_oracular()
        else:
            assert False, "unsupported multiplexor expansion style"

    @staticmethod
    def mp_mat(rad_angles, herm_conj=False):
        """
        This function returns a numpy array with the multiplexor matrix

            [ cc, ss]
            [-ss, cc]

        in it, where cc (ss) is the component-wise cosine (sine) of rad_angles

        Parameters
        ----------
        rad_angles : list(float)
        herm_conj : bool
            When herm_conj=True, it uses ss = sine(-rad_angles), whereas if
            herm_conj=False, it uses ss = sine(+rad_angles)

        Returns
        -------
        np.ndarray

        """
        num_angles = len(rad_angles)
        num_rows = 2*num_angles
        mat = np.zeros((num_rows, num_rows))
        sign = 1
        if herm_conj:
            sign = -1
        for k in range(num_angles):
            c = np.cos(rad_angles[k])
            s = np.sin(rad_angles[k])
            mat[k, k] = c
            mat[k+num_angles, k+num_angles] = c
            mat[k, k+num_angles] = sign*s
            mat[k+num_angles, k] = -sign*s
        return mat

if __name__ == "__main__":
    nt = 1
    nf = 2
    num_MP_trols = 2
    num_angles = (1 << num_MP_trols)
    rad_angles = list(np.random.rand(num_angles)*2*np.pi)

    for style in ['one_line', 'exact', 'oracular']:
        num_gbits = 0
        if style == 'oracular':
            num_gbits = 3
        num_bits = nt + nf + num_MP_trols + 1 + num_gbits
        emb = CktEmbedder(num_bits, num_bits)
        file_prefix = "../io_folder/plexor_test_" + style
        wr = MultiplexorSEO_writer(file_prefix, emb, style, rad_angles,
            num_T_trols=nt, num_F_trols=nf, num_gbits=num_gbits)
        wr.write()
        wr.close_files()

    file_prefix = "../io_folder/plexor_exact_check"
    num_bits = 4
    num_angles = (1 << (num_bits-1))
    emb = CktEmbedder(num_bits, num_bits)
    rad_angles = list(np.random.rand(num_angles)*2*np.pi)
    wr = MultiplexorSEO_writer(file_prefix, emb, 'exact', rad_angles)
    wr.write()
    wr.close_files()
    matpro = SEO_MatrixProduct(file_prefix, num_bits)
    exact_mat = MultiplexorSEO_writer.mp_mat(rad_angles)
    print(np.linalg.norm(matpro.prod_arr - exact_mat))
