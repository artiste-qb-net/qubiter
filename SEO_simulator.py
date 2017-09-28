import numpy as np
import copy as cp
# import pprint as pp
from SEO_reader import *
from OneBitGates import *
from StateVec import *
# import Utilities as ut


class SEO_simulator(SEO_reader):
    """
    This class simulates the evolution of a quantum state vector.

    This class has SEO_reader as a parent. Each line of an English file is
    read by the parent class and handed over to the use_ functions of this
    simulator class. The use_  functions multiply the current state vector by
    the unitary matrix that represents the latest line read.

    An initial state vector can be entered via the constructor or else it is
    set to the ground state automatically by the constructor.

    3 kinds (called 0, 1, 2) of measurements MEAS are allowed. A type 0
    measurement inserts a projector |0><0| = n = P_0 at the target bit. A
    type 1 measurement inserts a projector |1><1| = nbar = P_1 at the target
    bit. A type 2 measurement stores a copy of the state vector after |0><0|
    has been applied, and another copy after |1><1| has been applied.

    cur_st_vec_list is a list of state vectors on num_bits qubits. We will
    refer to each state vec in the list as a branch. Initially, this list
    contains a single branch. A measurement MEAS of kinds 0 or 1 does not
    change the number of branches in the list, but a measurement of kind 2
    doubles their number.

    Note that since projectors are not unitary matrices, the branches of
    cur_st_vec_list are not expected to be normalized except when there is
    only a single branch.

    If cur_st_vec_list consists of a list of branches (numpy arrays) denoted
    by [|br0>, |br1>, |br2>, ...], then one can construct the density matrix
    of that state as \rho = |br0><br0| + |br1><br1| + |br2><br2| + .... In
    other words, cur_st_vec_list is just a particular way of storing the
    density matrix of a state. A state with a single branch is a pure state,
    but a state with more than one branch may not be.

    Attributes
    ----------
    cur_st_vec_dict : dict(str, np.ndarray)
        dictionary with key= branch_key and values= numpy array for state
        vector. If there is a single item in dict, it's a pure state and the
        key is "pure". If there are more than one items , the state is mixed
        and the branch key is a string that uniquely characterizes the
        measured controls. For example, if qubit 2 has been measured True
        and qubit 4 has been measured False, the branch key will be '4F2T'.
        This only applies to type 2 measurements.
    """

    # rrtucci: combines my java classes:
    # LineList, UnitaryMat, SEO_readerMu

    def __init__(self, file_prefix, num_bits,
                 init_st_vec=None, verbose=False):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        init_st_vec : np.array
            get this using the functions get_ground_st_vec() or
            get_standard_basis_st_vec()
        verbose : bool

        Returns
        -------

        """
        if init_st_vec is None:
            init_st_vec = StateVec.get_ground_st_vec(num_bits)
        self.cur_st_vec_dict = {"pure": init_st_vec}
        self.verbose = verbose

        SEO_reader.__init__(self, file_prefix, num_bits, verbose)

    @staticmethod
    def branch_is_part_of_mcase(br_trols, case_trols):
        """
        Returns True iff the controls br_trols defining a branch of
        self.cur_st_vec_dict agree with the controls case_trols defining a
        case measured.

        Parameters
        ----------
        br_trols : Controls
        case_trols : Controls

        Returns
        -------
        bool

        """
        assert br_trols is not None
        assert case_trols is not None
        br_dict = br_trols.bit_pos_to_kind
        case_dict = case_trols.bit_pos_to_kind
        ans = True
        for bit in case_dict.keys():
            if bit in br_dict.keys():
                if case_dict[bit] != br_dict[bit]:
                    ans = False
                    break
            else:
                ans = False
                break
        # print('\n')
        # print('br_trols', br_trols.bit_pos_to_kind)
        # print('case_trols', case_trols.bit_pos_to_kind)
        # print(' is part of case', ans)
        return ans

    def get_controls_from_br_key(self, br_key):
        """
        Returns a Controls object built from br_key. br_key is assumed to be
        a str key for self.cur_st_vec_dict

        Parameters
        ----------
        br_key : str

        Returns
        -------
        Controls

        """
        assert br_key != "pure"
        x = br_key.replace("T", " T ")
        x = x.replace("F", " F ")
        li = x.split()
        # print("li", li)
        bit_pos = []
        kinds = []
        for s in li:
            if s.isdigit():
                bit_pos.append(int(s))
            else:
                if s == 'T':
                    kinds.append(True)
                else:
                    kinds.append(False)
        trols = Controls(self.num_bits)
        trols.bit_pos_to_kind = dict(zip(bit_pos, kinds))
        trols.refresh_lists()
        return trols

    @staticmethod
    def get_br_key_with_new_link(br_key, new_bit_pos, new_kind):
        """
        Say new_bit_pos=2 and new_kind=True. This returns a new branch
        key that adds the link 2T to the old br_key

        Parameters
        ----------
        br_key : str
        new_bit_pos : int
        new_kind : bool

        Returns
        -------
        str

        """
        x = str(new_bit_pos) + ("T" if new_kind else "F")
        if br_key == "pure":
            new_br_key = x
        else:
            new_br_key = br_key + x
        return new_br_key

    def evolve_by_controlled_bit_swap(self, bit1, bit2, controls):
        """
        Evolve each branch of cur_st_vec_dict by controlled bit swap iff the 
        bit swap line is outside of an if_m block, or if it is inside such a 
        block, evolve branch if it satisfies self.mcase_trols. 

        Parameters
        ----------
        bit1 : int
        bit2 : int
            bit1 and bit2 are the positions of the 2 bits being swapped.
        controls : Controls

        Returns
        -------
        None

        """
        assert bit1 != bit2, "swapped bits must be different"
        for bit in [bit1, bit2]:
            assert -1 < bit < self.num_bits
            assert bit not in controls.bit_pos
 
        slicex = [slice(None)]*self.num_bits
        num_controls = len(controls.bit_pos_to_kind)
        for k in range(num_controls):
            assert isinstance(controls.kinds[k], bool)
            if controls.kinds[k]:  # it's True
                slicex[controls.bit_pos[k]] = 1
            else:  # it's False
                slicex[controls.bit_pos[k]] = 0

        # components that are fixed are no longer axes
        scout = 0
        for bit in range(self.num_bits):
            if bit == bit1:
                new1 = scout
            if bit == bit2:
                new2 = scout
            if bit not in controls.bit_pos:
                scout += 1

        # perm = permutation that will use in transpose()
        perm_len = scout
        perm = list(range(perm_len))
        perm[new1], perm[new2] = perm[new2], perm[new1]

        # br = branch
        for br_key in self.cur_st_vec_dict.keys():
            if self.cur_st_vec_dict[br_key] is None:
                continue
            evolve_br = False
            if not self.measured_bits or not self.mcase_trols:
                evolve_br = True
            else:
                br_trols = self.get_controls_from_br_key(br_key)
                if SEO_simulator.branch_is_part_of_mcase(
                        br_trols, self.mcase_trols):
                    evolve_br = True
            if evolve_br:
                vec = self.cur_st_vec_dict[br_key][slicex]
                self.cur_st_vec_dict[br_key][slicex] = \
                    vec.transpose(perm)

    def evolve_by_controlled_one_bit_gate(self,
                tar_bit_pos, controls, one_bit_gate):
        """
        Evolve each branch of cur_st_vec_dict by controlled one bit gate ( 
        from class OneBitGates). Note one_bit_gate is entered as np.ndarray. 
        Evolve each branch of cur_st_vec_dict iff the controlled one bit 
        gate line is outside of an if_m block, or if it is inside such a 
        block, evolve branch if it satisfies self.mcase_trols. 

        Parameters
        ----------
        tar_bit_pos : int
            bit position of target of one bit gate.
        controls : Controls
        one_bit_gate : np.ndarray

        Returns
        -------
        None

        """
        assert -1 < tar_bit_pos < self.num_bits

        vec_slicex = [slice(None)]*self.num_bits
        num_controls = len(controls.bit_pos_to_kind)
        for k in range(num_controls):
            assert isinstance(controls.kinds[k], bool)
            if controls.kinds[k]:  # it's True
                vec_slicex[controls.bit_pos[k]] = 1
            else:  # it's False
                vec_slicex[controls.bit_pos[k]] = 0
        vec_slicex = tuple(vec_slicex)

        # components that are fixed are no longer axes
        scout = 0
        for bit in range(self.num_bits):
            if bit == tar_bit_pos:
                new_tar = scout
            if bit not in controls.bit_pos:
                scout += 1

        perm_len = scout
        # example tar = 2
        # want to map [2, 0, 1, 3, 4] back to [0, 1, 2, 3, 4]

        # this didn't work
        # use perm 0>2, 1>0, 2>1, 3>3, 4>4
        # perm = [new_tar] + list(range(new_tar))
        # perm += list(range(new_tar+1, perm_len))

        # use perm 2>0, 0>1, 1>2, 3>3, 4>4
        perm = list(range(1, new_tar+1)) + [0]
        perm += list(range(new_tar+1, perm_len))

        # br = branch
        assert not(self.mcase_trols and not self.measured_bits)
        for br_key in self.cur_st_vec_dict.keys():
            if self.cur_st_vec_dict[br_key] is None:
                continue
            evolve_br = False
            if not self.measured_bits or not self.mcase_trols:
                evolve_br = True
            else:
                br_trols = self.get_controls_from_br_key(br_key)
                if SEO_simulator.branch_is_part_of_mcase(
                        br_trols, self.mcase_trols):
                    evolve_br = True
            if evolve_br:
                vec = self.cur_st_vec_dict[br_key][vec_slicex]
                # Axes 1 of one_bit_gate and new_tar of vec are summed over.
                #  Axis 0 of one_bit_gate goes to the front of all the axes
                # of new vec. Use transpose() to realign axes.
                vec = np.tensordot(one_bit_gate, vec, ([1], [new_tar]))
                self.cur_st_vec_dict[br_key][vec_slicex] = \
                    np.transpose(vec, axes=perm)

    def finalize_next_line(self):
        """
        Prints running documentary at the end of the reading of each line.

        Returns
        -------
        None

        """
        if self.verbose:
            print('\n')
            print(self.split_line)
            print('line number = ', self.line_count)
            print('operation = ', self.num_ops)
            st_vecs = self.cur_st_vec_dict
            StateVec.describe_st_vec_dict(st_vecs,
                                          # print_st_vec=True,
                                          show_probs=True)

    def use_DIAG(self, trols, rad_angles):
        """
        This method should not be called. If called, it explains why it
        shouldn't be called.

        Parameters
        ----------
        trols : Controls
        rad_angles : list[float]

        Returns
        -------

        """
        assert False, \
            "This class cannot simulate a circuit containing " \
            "raw diagonal unitaries DIAG. Work around: use first our" \
            "DiagUnitaryExpander class to expand " \
            "DIAGs into simpler gates."

    def use_HAD2(self, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for had2.


        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneBitGates.had2()
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_MEAS(self, tar_bit_pos, kind):
        """
        Overrides the parent class use_ function.

        For kind 0 (resp., 1) measurements, it applies |0><0| (resp.,
        |1><1|) to each branch of cur_st_vec_dict.

        For kind 2 measurements, it first creates a list concatenating x
        plus a deep copy of x, where x is the cur_st_vec_dict. Next,
        it applies P_0 to the first half of the list and P_1 to the second
        half.

        Parameters
        ----------
        kind : int
        tar_bit_pos: int

        Returns
        -------
        None

        """
        # slicex = slice index
        slicex = [slice(None)]*self.num_bits
        # br = branch
        if kind in [0, 1]:
            b = 1 if kind == 0 else 0
            for br_key in self.cur_st_vec_dict:
                st_vec = self.cur_st_vec_dict[br_key]
                if st_vec is not None:
                    slicex[tar_bit_pos] = b
                    # set projection |b=0><b=0| to zero for kind=1
                    st_vec[tuple(slicex)] = 0
                    tot_prob = StateVec.get_total_prob(st_vec)
                    if tot_prob < 1e-8:
                        # this didn't work
                        # st_vec = None
                        self.cur_st_vec_dict[br_key] = None
                    slicex[tar_bit_pos] = slice(None)
                # self.cur_st_vec_dict[br_key] = st_vec
        elif kind == 2:
            old_st_vec_dict = cp.deepcopy(self.cur_st_vec_dict)
            self.cur_st_vec_dict = {}
            for br_key in old_st_vec_dict.keys():
                new_T_key = self.get_br_key_with_new_link(br_key,
                        tar_bit_pos, True)
                new_F_key = self.get_br_key_with_new_link(br_key,
                        tar_bit_pos, False)
                # print('new keys=', new_F_key,',', new_T_key)

                self.cur_st_vec_dict[new_T_key] = \
                    cp.deepcopy(old_st_vec_dict[br_key])
                self.cur_st_vec_dict[new_F_key] = \
                    cp.deepcopy(old_st_vec_dict[br_key])

                for b, new_key in enumerate([new_T_key, new_F_key]):
                    # set projection |b=0><b=0| to zero for T key
                    st_vec = self.cur_st_vec_dict[new_key]
                    # print("b, newkey=" + str(b) + "," + new_key)
                    if st_vec is not None:
                        slicex[tar_bit_pos] = b
                        st_vec[tuple(slicex)] = 0
                        tot_prob = StateVec.get_total_prob(st_vec)
                        # print('tot_prob=', tot_prob)
                        if tot_prob < 1e-8:
                            # this didn't work
                            # st_vec = None
                            self.cur_st_vec_dict[new_key] = None
                        slicex[tar_bit_pos] = slice(None)
                    # print(st_vec)

            # print(self.cur_st_vec_dict)
        else:
            assert False, 'unsupported measurement kind'

    def use_MP_Y(self, tar_bit_pos, trols, rad_angles):
        """
        This method should not be called. If called, it explains why it
        shouldn't be called.

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------

        """
        assert False, \
            "This class cannot simulate a circuit containing " \
            "raw multiplexors MP_Y. Work around: use first our" \
            "MultiplexorExpander class to expand " \
            "MP_Ys into simpler gates."

    def use_PHAS(self, angle_degs, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for PHAS.

        Parameters
        ----------
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneBitGates.phase_fac(angle_degs * np.pi/180)
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_P_PH(self, projection_bit,
                angle_degs, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for P_0 and P_1 phase factors.


        Parameters
        ----------
        projection_bit : int
            0 (resp. 1) for P_0 (resp. P_1) projection
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        fun = {
            0: OneBitGates.P_0_phase_fac,
            1: OneBitGates.P_1_phase_fac
        }
        gate = fun[projection_bit](angle_degs*np.pi/180)
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_PRINT(self, style, line_num):
        """
        Prints to screen the self.cur_st_vec_dict.

        Parameters
        ----------
        style : str
            style in wwhich to print
        line_num : int
            line number in eng & pic files in which PRINT command appears

        Returns
        -------
        None

        """
        print("\n*************************beginning PRINT output")
        print("PRINT line number=" + str(line_num))
        if style == "V1":
            st_vecs = self.cur_st_vec_dict
            StateVec.describe_st_vec_dict(st_vecs,
                                        # print_st_vec=True,
                                        show_probs=True)
        else:
            assert False, "unsupported PRINT style"
        print("****************************ending PRINT output")

    def use_ROT(self, axis,
                angle_degs, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for rot along axes x, y, or z.

        Parameters
        ----------
        axis : int
            1, 2, 3 for x, y, z
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneBitGates.rot_ax(angle_degs * np.pi / 180, axis)
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for rot along arbitrary axis.


        Parameters
        ----------
        angle_x_degs : float
        angle_y_degs : float
        angle_z_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneBitGates.rot(
                    angle_x_degs*np.pi/180,
                    angle_y_degs*np.pi/180,
                    angle_z_degs*np.pi/180)
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for sigx, sigy, sigz.

        Parameters
        ----------
        axis : int
            1, 2, 3 for x, y, z
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        fun = {
            1: OneBitGates.sigx,
            2: OneBitGates.sigy,
            3: OneBitGates.sigz
        }
        gate = fun[axis]()
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_SWAP(self, bit1, bit2, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_bit_swap().

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        self.evolve_by_controlled_bit_swap(bit1, bit2, controls)

if __name__ == "__main__":

    # use test = 0 if want to run all tests at once.
    test = 0
    # test = 3
    if test in [0, 1]:
        # test on circuit for a quantum fourier transform
        # (no loops, no internal measurements)
        sim = SEO_simulator('io_folder/sim_test1', 6, verbose=True)

    if test in [0, 2]:
        # test embedded loops
        sim = SEO_simulator('io_folder/sim_test2', 4, verbose=True)

    if test in [0, 3]:
        # test MEAS branching. Each kind 2 measurement doubles number of
        # branches
        sim = SEO_simulator('io_folder/sim_test3', 4, verbose=True)

