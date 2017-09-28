from SEO_writer import *
from SEO_simulator import *
from CktEmbedder import *
from BitVector import *
import itertools as it
import collections as co


class CGateSEO_writer(SEO_writer):
    """
    This class is a child of SEO_writer. When one_line=True, it writes a
    single line for an c_u2 (controlled U(2) gate).  When one_line=False,
    it writes an expansion of that gate. When expand_1c_u2=True, the c_u2 is
    expanded so that its sub-component 1c_u2 (singly controlled U(2) gates)
    are expanded into cnots and qubit rotations. If expand_1c_u2=False,
    the 1c_u2 are not expanded.

    If we say a gate is controlled, it may have 1 or more controls (it might
    be singly or multiply controlled).

    Global phase factors are ignored, so expansions equal the original line
    up to a phase factor.

    References
    ----------

    1. CGateSEO_writer.pdf, by Robert Tucci, included with Qubiter source
    code.

    Attributes
    ----------
    do_checking : bool
        Does some checking of algebra
    expand_1c_u2 : bool
        When this is True, the c_u2 (controlled U(2) gate) is expanded so
        that its sub-component 1c_u2 (uni controlled U(2) gates) are
        expanded into cnots and qubit rotations. If this is False, the 1c_u2
        are not expanded.
    one_line : bool
        When this is True, it writes a single line for the c_u2 ( controlled
        U(2) gate). When False, it writes an expansion of that gate.

    """

    def __init__(self, file_prefix, emb,
            one_line=True, expand_1c_u2=False, do_checking=False,
                 verbose=False, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        emb : CktEmbedder
        one_line : bool
        expand_1c_u2 : bool
        do_checking : bool
        verbose : bool
        kwargs : dict

        Returns
        -------
        None

        """
        self.one_line = one_line
        self.expand_1c_u2 = expand_1c_u2
        self.do_checking = do_checking
        self.verbose = verbose
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)

    @staticmethod
    def su2_mat_prod(su2_pair1, su2_pair2):
        """
        An SU(2) matrix can be expressed as exp(i*theta*sig_n) where theta
        is a real number and sig_n = n \cdot sigma. Here n is a 3 dim real
        UNIT vector and sigma = [sigx, sigy, sigz], where sigx, sigy and
        sigz are the 3 Pauli matrices. We define the su2_pair of this SU(2)
        matrix as the list [theta, n], with n expressed as a numpy array.
        One can prove by Taylor expansion that
    
        exp(i*theta*sig_n) = c + i*sig_n*s
    
        where c = cos(theta) and s = sin(theta).
    
        This subroutine maps su2_pair1, su2_pair2 --> su2_pair
    
        where
    
        exp(i*theta*sig_n) = exp(i*theta1*sig_n1) exp(i*theta2*sig_n2)
    
        Parameters
        ----------
        su2_pair1 : [float, np.array]
            su2_pair of left matrix in matrix product
    
        su2_pair2 : [float, np.array]
            su2_pair of right matrix in matrix product
    
        Returns
        -------
        [float, np.array]
    
        """
    
        theta1, n1 = su2_pair1
        theta2, n2 = su2_pair2
        s1, c1 = np.sin(theta1), np.cos(theta1)
        s2, c2 = np.sin(theta2), np.cos(theta2)
        n = s1*c2*n1 + s2*c1*n2 - s1*s2*np.cross(n1, n2)
        mag = np.linalg.norm(n)
        n /= mag
        c = c1*c2 - s1*s2*np.dot(n1, n2)
        s = mag
        theta = np.arctan2(s, c)
        return [theta, n]

    def write_1c_u2(self, tar_bit_pos, trol_bit_pos, rads_list, delta=None):
        """
        Writes an expansion of an 1c_u2 (singly controlled U(2) matrix). In
        general, such an expansion will contain 3 cnots, but for special
        cases taken here into account, it's possible to get away with using
        only 2 or 1 cnots.

        Parameters
        ----------
        tar_bit_pos : int
            target bit position
        trol_bit_pos : int
            control bit position

        rads_list : list[float]
            list of 3 angles in radians. If it equals [radx, rady, radz],
            then U(2) gate given by e^{i*delta} exp(i*(radx*sigx + rady*sigy
            + radz*sigz))

        delta : float|None
            U(2) gate being controlled equals e^{i*delta} times SU(2) gate

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        trols = Controls.new_knob(num_bits, trol_bit_pos, True)
        if not self.expand_1c_u2:
            self.write_controlled_one_bit_gate(
                tar_bit_pos, trols, OneBitGates.rot, rads_list)
            return

        def write_delta_rot():
            if delta:
                self.write_one_bit_gate(
                    trol_bit_pos, OneBitGates.rot_ax, [-delta/2, 3])

        def write_cnot():
            self.write_controlled_one_bit_gate(
                tar_bit_pos, trols, OneBitGates.sigx)

        def write_rot(rads_list1, herm_conj=False):
            if not herm_conj:
                rads_list2 = rads_list
            else:
                rads_list2 = list(-np.array(rads_list1))
            self.write_one_bit_gate(
                tar_bit_pos, OneBitGates.rot, rads_list2)

        rads_x, rads_y, rads_z = rads_list
        theta_w = np.sqrt(rads_x**2 + rads_y**2 + rads_z**2)
        cw, sw = np.cos(theta_w), np.sin(theta_w)
        wx, wy, wz = [rads_x/theta_w, rads_y/theta_w, rads_z/theta_w]
        TOL = 1E-6

        if abs(theta_w - np.pi/2) < TOL:  # 1 cnot, 0 or 2 target rots
            if abs(wy) < TOL and abs(wz) < TOL and \
                    abs(delta + np.pi/2) < TOL:  # simple cnot, 0 target rots
                write_cnot()
            else:  # 1 cnot, 2 target rots
                theta_a = np.pi/2
                ax = np.sqrt((wx+1)/2)
                # ax != 0 or else wy=wz=0, which was already considered
                ay = wy/(2*ax)
                az = wz/(2*ax)
                rads_list_a = list(theta_a*np.array([ax, ay, az]))

                write_rot(rads_list_a, herm_conj=True)
                write_cnot()
                write_rot(rads_list_a)
                write_delta_rot()

                if self.do_checking:
                    mat_w = np.matrix(OneBitGates.rot(*rads_list))
                    mat_a = np.matrix(OneBitGates.rot(*rads_list_a))
                    mat_sigx = np.matrix(OneBitGates.sigx())
                    diff = mat_w - mat_a*mat_sigx*mat_a.getH()
                    err = np.linalg.norm(diff)
                    if err > TOL:
                        print("1 cnot, 2 rots")
                        print(diff)
                        assert False
            
        elif abs(wx) < TOL:  # 2 cnots, 2 target rots
            # this is the same as the general case (2 cnots, 3 target rots)
            # with alp=beta
            alp = np.arctan2(wz*sw, cw)/2
            gamma = np.arctan2(sw*wy,
                               np.sqrt(cw**2 + (wz*sw)**2))
            su2_pair_a = CGateSEO_writer.su2_mat_prod(
                [alp, np.array([0, 0, 1])],
                [gamma/2, np.array([0, 1, 0])])
            rads_list_a = list(su2_pair_a[0]*su2_pair_a[1])

            write_cnot()
            write_rot(rads_list_a, herm_conj=True)
            write_cnot()
            write_rot(rads_list_a)
            write_delta_rot()

            if self.do_checking:
                mat_w = np.matrix(OneBitGates.rot(*rads_list))
                mat_a = np.matrix(OneBitGates.rot(*rads_list_a))
                mat_sigx = np.matrix(OneBitGates.sigx())
                diff = mat_w - mat_a*mat_sigx*mat_a.getH()*mat_sigx
                err = np.linalg.norm(diff)
                if err > TOL:
                    print("2 cnot, 2 rots", rads_list)
                    print(diff)
                    assert False

        else:  # 2 cnots, 3 target rots
            theta1 = np.arctan2(wz*sw, cw)
            theta2 = np.arctan2(wx, wy)
            alp = (theta1 + theta2)/2
            beta = (theta1 - theta2)/2
            gamma = np.arctan2(sw*np.sqrt(wx**2 + wy**2),
                               np.sqrt(cw**2 + (wz*sw)**2))
            su2_pair_a = CGateSEO_writer.su2_mat_prod(
                [alp, np.array([0, 0, 1])],
                [gamma/2, np.array([0, 1, 0])])
            su2_pair_b = CGateSEO_writer.su2_mat_prod(
                [-gamma/2, np.array([0, 1, 0])],
                [-(alp+beta)/2, np.array([0, 0, 1])])
            su2_pair_c = [(beta-alp)/2, np.array([0, 0, 1])]
            rads_list_a = list(su2_pair_a[0]*su2_pair_a[1])
            rads_list_b = list(su2_pair_b[0]*su2_pair_b[1])
            rads_list_c = list(su2_pair_c[0]*su2_pair_c[1])

            write_rot(rads_list_c)
            write_cnot()
            write_rot(rads_list_b)
            write_cnot()
            write_rot(rads_list_a)
            write_delta_rot()

            if self.do_checking:
                mat_id = np.matrix(OneBitGates.phase_fac(0.0))
                mat_w = np.matrix(OneBitGates.rot(*rads_list))
                mat_a = np.matrix(OneBitGates.rot(*rads_list_a))
                mat_b = np.matrix(OneBitGates.rot(*rads_list_b))
                mat_c = np.matrix(OneBitGates.rot(*rads_list_c))
                mat_sigx = np.matrix(OneBitGates.sigx())
                diff = mat_w - mat_a*mat_sigx*mat_b*mat_sigx*mat_c
                err = np.linalg.norm(diff)
                if err > TOL:
                    print("2 cnot, 3 rots, identity 1")
                    print(diff)
                    assert False
                diff = mat_id - mat_a*mat_b*mat_c
                err = np.linalg.norm(diff)
                if err > TOL:
                    print("2 cnot, 3 rots, identity 2")
                    print(diff)
                    assert False

    # def write_gen_n_controlled_u2(self, n_index_list, rads_list, delta=None):
    #     """
    #     Writes an expansion for a U(2) matrix W(num_bits-1) that is
    #     controlled by a "generalized n" equal to GN = n(n_index_list). Thus,
    #     the gate written by this function equals W(num_bits-1)^GN.
    #     Generalized n's are defined in the reference CktExpander.pdf
    #
    #     Parameters
    #     ----------
    #     n_index_list : list[int]
    #         indices of the generalized n
    #
    #     rads_list : list[float]
    #         list of 3 angles in radians. If it equals [radx, rady, radz],
    #         then U(2) gate given by e^{i*delta} exp(i*(radx*sigx + rady*sigy
    #         + radz*sigz))
    #
    #     delta : float|None
    #         U(2) gate being controlled equals e^{i*delta} times SU(2) gate
    #
    #     Returns
    #     -------
    #     None
    #
    #     """
    #
    #     num_bits = self.emb.num_bits_bef
    #
    #     for k in range(len(n_index_list)-1):
    #         tar_pos = n_index_list[k+1]
    #         trol_pos = n_index_list[k]
    #         trols = Controls.new_knob(num_bits, trol_pos, True)
    #         self.write_controlled_one_bit_gate(
    #             tar_pos, trols, OneBitGates.sigx)
    #
    #     self.write_1c_u2(num_bits - 1, n_index_list[-1], rads_list, delta)
    #
    #     for k in reversed(range(len(n_index_list)-1)):
    #         tar_pos = n_index_list[k+1]
    #         trol_pos = n_index_list[k]
    #         trols = Controls.new_knob(num_bits, trol_pos, True)
    #         self.write_controlled_one_bit_gate(
    #             tar_pos, trols, OneBitGates.sigx)
    #
    # def write_internal(self, rads_list, delta=None):
    #     """
    #     This internal function is used in write() and is less general than
    #     the latter. It expands an c_u2 into a product of several
    #     "generalized n" controlled U(2) gates.
    #
    #     Parameters
    #     ----------
    #     rads_list : list[float]
    #         list of 3 angles in radians. If it equals [radx, rady, radz],
    #         then U(2) gate given by e^{i*delta} exp(i*(radx*sigx + rady*sigy
    #         + radz*sigz))
    #
    #     delta : float|None
    #         U(2) gate being controlled equals e^{i*delta} times SU(2) gate
    #
    #     Returns
    #     -------
    #     None
    #
    #     """
    #     num_bits = self.emb.num_bits_bef
    #     for num_boxes in range(1, num_bits):
    #         for comb in it.combinations(range(0, num_bits-1), num_boxes):
    #             n_index_list = sorted(list(comb), reverse=True)
    #             sign = 1
    #             if len(comb) % 2 == 0:
    #                 sign = -1
    #             new_rads_list = list(
    #                 sign*np.array(rads_list)/(1 << (num_bits-2)))
    #             if delta:
    #                 new_delta = sign*delta/(1 << (num_bits-2))
    #             else:
    #                 new_delta = None
    #             self.write_gen_n_controlled_u2(
    #                 n_index_list, new_rads_list, new_delta)

    def write_internal(self, rads_list, delta=None):
        """
        This internal function is used in write() and is less general than
        the latter. It expands an c_u2 into a product of 1c_u2 with
        intervening cnots.

        In the CGateSEO_writer.pdf documentation, we show that any c_u2 can
        be expanded into a product of several "generalized n" controlled U(
        2) gates of the form  W(num_bits-1)^GN, wherein U(2) matrix W(
        num_bits-1) is controlled by a "generalized n" equal to GN = n(
        n_index_list)

        Since the factors W(num_bits-1)^GN in the product commute amongst
        themselves, it is possible and convenient to order them in Gray code
        order (Qubiter knows about Gray Code via its class BitVector).
        Ordering them in Gray Code allows this function to cancel some cnots
        from adjacent GN.

        An earlier version of this function, now commented, did not use Gray
        Code and used more cnots than this one.

        Parameters
        ----------
        rads_list : list[float]
            list of 3 angles in radians. If it equals [radx, rady, radz],
            then U(2) gate given by e^{i*delta} exp(i*(radx*sigx + rady*sigy
            + radz*sigz))

        delta : float|None
            U(2) gate being controlled equals e^{i*delta} times SU(2) gate

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        num_trols = num_bits-1
        max_f = (1 << num_trols)-1

        def write_cnot(tar_bpos, trol_bpos):
            trol = Controls.new_knob(num_bits, trol_bpos, True)
            self.write_controlled_one_bit_gate(tar_bpos, trol,
                                               OneBitGates.sigx)
        def write_cnot_stair(bvec):
            tar_bpos = bvec.find_rightmost_T_bit()
            trol_bpos = tar_bpos
            while True:
                trol_bpos = bvec.find_T_bit_to_left_of(trol_bpos)
                if trol_bpos == -1:
                    break
                write_cnot(tar_bpos, trol_bpos)

        cur_bvec = BitVector(num_trols, 0)
        prev_bvec = BitVector(num_trols, 0)
        f, lazy = 0, 0
        f, lazy = BitVector.lazy_advance(f, lazy)
        cur_bvec.dec_rep = lazy
        while f <= max_f:
            if self.verbose:
                print("\nf, lazy", f, lazy)
                self.write_NOTA(str(prev_bvec) + "->" + str(cur_bvec))
            sign = 1
            if cur_bvec.get_num_T_bits() % 2 == 0:
                sign = -1
            new_rads_list = list(
                sign*np.array(rads_list)/(1 << (num_bits-2)))
            if delta:
                new_delta = sign*delta/(1 << (num_bits-2))
            else:
                new_delta = None

            diff_bvec = BitVector.new_with_T_on_diff(cur_bvec, prev_bvec)
            diff_bpos = diff_bvec.find_rightmost_T_bit()
            min_prev_bpos = prev_bvec.find_rightmost_T_bit()
            min_cur_bpos = cur_bvec.find_rightmost_T_bit()

            if f > 1:  # first 1c_u2 has not cnots preceding it
                if min_cur_bpos == min_prev_bpos:
                    write_cnot(min_cur_bpos, diff_bpos)
                else:
                    write_cnot_stair(prev_bvec)
                    write_cnot_stair(cur_bvec)

            u2_trol_bpos = min_cur_bpos
            self.write_1c_u2(
                num_bits - 1, u2_trol_bpos, new_rads_list, new_delta)
            prev_bvec = BitVector.copy(cur_bvec)
            f, lazy = BitVector.lazy_advance(f, lazy)
            cur_bvec.dec_rep = lazy

    def write_hads(self, trol_kinds, herm_conj=False):
        """
        Writes a chain of cnots that are useful when some of the controls of
        the c_u2 being considered are n_bar = P_0 = |0><0| instead of n =
        P_1 = |1><1|. We are using the identity sigx n sigx = nbar to
        convert n's to nbar's.

        Parameters
        ----------
        trol_kinds : list[bool]
            A list of control kinds. True for n=P_1 and False for n_bar=P_0
        herm_conj : bool
            When this is True, writes Hermitian conjugate of expansion.

        Returns
        -------
        None

        """
        num_trols = len(trol_kinds)
        if not herm_conj:
            range1 = range(num_trols)
        else:
            range1 = reversed(range(num_trols))
        for k in range1:
            if not trol_kinds[k]:
                self.write_one_bit_gate(num_trols-k-1, OneBitGates.sigx)

    def write(self, trol_kinds, u2_fun, fun_arg_list=None):
        """
        This is the most general function of this class. All other functions
        of the class are mostly internal and are called by this function.
        This function achieves the main goal of the class, which is to give
        various expansions of an c_u2 (controlled U(2) matrix). For
        one_line=True, this function just calls
        write_controlled_one_bit_gate() of the parent class. For
        one_line=False, it gives an expansion of the c_u2.

        Parameters
        ----------
        trol_kinds : list[bool]
            list of control types. Type is False if nbar=P_0 and True if n=P_1
        u2_fun : function
            One of the functions in class OneBitGates
        fun_arg_list : list[int|float]
            list of arguments of u2_fun

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef
        tar_bit_pos = num_bits-1
        num_trols = num_bits-1
        assert len(trol_kinds) == num_trols
        trols = Controls(num_bits)
        trols.bit_pos_to_kind = {k: trol_kinds[num_bits-k-2]
                                 for k in range(0, num_bits-1)}
        trols.refresh_lists()

        if self.one_line or num_trols == 0:
            self.write_controlled_one_bit_gate(tar_bit_pos,
                    trols, u2_fun, fun_arg_list)
            return

        # insert opening Hadamards for controls equal to n_bar = |0><0|
        self.write_hads(trols.kinds)

        if u2_fun == OneBitGates.P_0_phase_fac:
            rads = fun_arg_list[0]
            self.write_internal([0, 0, rads / 2], rads / 2)
        elif u2_fun == OneBitGates.P_1_phase_fac:
            rads = fun_arg_list[0]
            self.write_internal([0, 0, -rads / 2], rads / 2)
        elif u2_fun == OneBitGates.sigx:
            if num_bits == 2:
                # If it's a CNOT, no expansion necessary
                # Control must be set to True because
                # opening and closing Hadamards take care of False
                trols1 = Controls.new_knob(num_bits, 0, True)
                self.write_controlled_one_bit_gate(
                    tar_bit_pos, trols1, OneBitGates.sigx)
            else:
                self.write_internal([np.pi / 2, 0, 0], -np.pi / 2)
        elif u2_fun == OneBitGates.sigy:
            self.write_internal([0, np.pi / 2, 0], -np.pi / 2)
        elif u2_fun == OneBitGates.sigz:
            self.write_internal([0, 0, np.pi / 2], -np.pi / 2)
        elif u2_fun == OneBitGates.had2:
            rads = np.pi/(2*np.sqrt(2))
            self.write_internal([rads, 0, rads], -np.pi / 2)
        elif u2_fun == OneBitGates.rot_ax:
            rads = fun_arg_list[0]
            axis = fun_arg_list[1]
            if axis == 1:
                self.write_internal([rads, 0, 0])
            elif axis == 2:
                self.write_internal([0, rads, 0])
            elif axis == 3:
                self.write_internal([0, 0, rads])
            else:
                assert False
        elif u2_fun == OneBitGates.rot:
            self.write_internal(fun_arg_list)
        else:
            assert False, "writing an unsupported controlled gate"

        # insert closing Hadamards for controls equal to n_bar = |0><0|
        self.write_hads(trols.kinds, herm_conj=True)

if __name__ == "__main__":

    num_bits_bef = 4
    num_bits_aft = 5
    bit_map = list(range(num_bits_bef))
    emb = CktEmbedder(num_bits_bef, num_bits_aft, bit_map)

    # trol_kinds in ZL convention
    trol_kinds = [True, False, False]

    wr = CGateSEO_writer('io_folder/cgate_expansions', emb,
            do_checking=True, verbose=False)

    u2_fun_to_fun_arg_list = co.OrderedDict((
        (OneBitGates.P_0_phase_fac, [np.pi/3]),
        (OneBitGates.P_1_phase_fac, [np.pi/3]),
        (OneBitGates.sigx, None),
        (OneBitGates.sigy, None),
        (OneBitGates.sigz, None),
        (OneBitGates.had2, None),
        (OneBitGates.rot_ax, [np.pi/3, 2]),
        (OneBitGates.rot, [np.pi/3, np.pi/6, np.pi/3])
    ))

    for u2_fun, fun_arg_list in u2_fun_to_fun_arg_list.items():

        wr.write_NOTA('--------new u2 gate --------------------------')
        for one_line in [True, False]:
            wr.one_line = one_line
            if one_line:
                wr.write(trol_kinds, u2_fun, fun_arg_list)
            else:
                for expand_1c_u2 in [False, True]:
                    wr.expand_1c_u2 = expand_1c_u2
                    wr.write_NOTA('--------expand_1c_u2=' + str(expand_1c_u2))
                    print("\n", u2_fun,
                          "one_line=", one_line, "expand=", expand_1c_u2)
                    wr.write(trol_kinds, u2_fun, fun_arg_list)

    wr.close_files()

    # a check that an expansion multiplies to original
    num_bits = 5
    emb = CktEmbedder(num_bits, num_bits)
    # trol_kinds in ZL convention
    trol_kinds = [True, False, False, False]
    file_prefix = 'io_folder/cgate_expan_mat_prod'

    wr = CGateSEO_writer(file_prefix, emb)

    u2_fun = OneBitGates.rot_ax
    rads = np.pi/3

    wr.one_line = True
    wr.write_NOTA("one line=True-----------")
    wr.write(trol_kinds, u2_fun, [rads, 2])

    wr.write_NOTA("herm. conj, one line=False-----------")
    wr.one_line = False
    wr.write(trol_kinds, u2_fun, [-rads, 2])

    wr.close_files()

    from SEO_MatrixProduct import *

    mp = SEO_MatrixProduct(file_prefix, num_bits)
    id_mat = np.diag(np.ones((1<<num_bits,)))
    err = np.linalg.norm(mp.prod_arr - id_mat)
    print("err=", err)


