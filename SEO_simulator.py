import numpy as np
import copy as cp
import pprint as pp
from SEO_reader import *
from OneBitGates import *
import Utilities as ut


class SEO_simulator(SEO_reader):
    """
    This class simulates the evolution of a quantum state vector.

    This class has SEO_reader as a parent. Each line of an English file is
    read by the parent class and handed over to the use_ functions of this
    simulator class. The use functions multiply the current state vector by
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
    cur_st_vec_list : List[np.ndarray]
        current state vector list
    verbose : bool
        True if want to print a running commentary on console.

    num_ops : int
        number of operations. Lines inside a loop with 'reps' repetitions
        will count as 'reps' operations
    loop_to_cur_rep : dict[int, int]
        a dictionary mapping loop number TO current repetition
    just_jumped : bool
        flag used to alert when loop jumps from NEXT to LOOP
    line_count : int

    english_in : _io.TextIOWrapper
        file object for input text file that stores English description of
        circuit
    file_prefix : str
        beginning of the name of English file being scanned
    loop_to_start_line : dict[int, int]
        a dictionary mapping loop number TO loop line + 1
    loop_to_start_offset : dict[int, int]
        a dictionary mapping loop number TO offset of loop's start
    loop_to_reps : dict[int, int]
        a dictionary mapping loop number TO total number of repetitions of
        loop
    loop_queue : list[int]
        a queue of loops labelled by their id number
    num_bits : int
        number of qubits in whole circuit
    tot_num_lines : int
        number of lines in English file
    split_line : list[str]
        storage space for a list of strings obtained by splitting a line

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
            get this using the functions get_ground_st() or
            get_standard_basis_st()
        verbose : bool

        Returns
        -------

        """
        if init_st_vec is None:
            init_st_vec = SEO_simulator.get_ground_st(num_bits)
        self.cur_st_vec_list = [init_st_vec]
        self.verbose = verbose

        SEO_reader.__init__(self, file_prefix, num_bits)

    @staticmethod
    def get_ground_st(num_bits):
        """
        Returns ground state |0>|0>|0>...|0>, where |0> = [1,0]^t and |1> =
        [0,1]^t, t = transpose

        Parameters
        ----------
        num_bits : int

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([1 << num_bits], dtype=ty)
        mat[0] = 1
        mat = mat.reshape([2]*num_bits)
        return mat

    @staticmethod
    def get_standard_basis_st(spin_dir_list, zero_last=True):
        """
        If zero_last = True, Returns state ...|s2>|s1>|s0>, where
        spin_dir_list=[...,s2, s1, s0], s_j \in {0, 1} for all j, |0> = [1,
        0]^t and |1> = [0,1]^t, t = transpose

        Parameters
        ----------
        spin_dir_list : list[int]
        zero_last : bool
            True(False) if last(first) qubit is at position 0

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        num_bits = len(spin_dir_list)
        mat = np.zeros([1 << num_bits], dtype=ty)
        mat = mat.reshape([2]*num_bits)
        if zero_last:
            spin_dir_list = reversed(spin_dir_list)
        mat[tuple(spin_dir_list)] = 1
        return mat

    @staticmethod
    def get_total_prob(st_vec):
        """
        Returns total probability of state vector st_vec.

        Parameters
        ----------
        st_vec : np.ndarray
            state vector

        Returns
        -------
        float

        """
        return np.sum(np.real(st_vec*st_vec.conj()))

    @staticmethod
    def get_bit_probs(st_vec):
        """
        For a given state vector st_vec over num_qubits qubits, it returns a
        dictionary that maps each qubit to a pair (p, 1-p), where p is the
        probability that that particular qubit is 0, if the state of all
        other qubits is ignored.

        Parameters
        ----------
        st_vec : np.ndarray
            state vector

        Returns
        -------
        dict[int, (float, float)]

        """
        prob_dict = {}
        num_bits = st_vec.ndim
        # slicex is a portmanteau of slice index
        slicex = [slice(None)]*num_bits
        tot_prob = SEO_simulator.get_total_prob(st_vec)
        for k in range(num_bits):
            slicex[k] = 0
            vec = st_vec[tuple(slicex)]
            p = np.sum(np.real(vec*vec.conj()))/tot_prob
            prob_dict[k] = (p, 1-p)
            slicex[k] = slice(None)  # restore to all entries slice(None)
        return prob_dict

    def describe_fin_st(
            self, print_st_vec=False, do_pp=False, omit_zero_amps=False):
        """
        Prints a description of the final state vector

        Parameters
        ----------

        print_st_vec : bool
            if True, prints the final state vector (which may be huge. For n
            qubits, it has 2^n components.)

        do_pp : bool
            pp= pretty print. Only used if print_st_vec=True. For pp=False,
            it prints final state vector in usual numpy array print style.
            For pp=True, it prints final state vector as column of (index,
            array value) pairs.

        omit_zero_amps : bool
            If print_st_vec=True, pp=True and this parameter is True too,
            will omit states with zero amplitude

        Returns
        -------
        None

        """
        fin_st_vec = self.cur_st_vec_list[0]
        if print_st_vec:
            print('final state vector')
            if do_pp:
                print('(zero bit first in state tuple)')
                ut.pp_numpy_arr(fin_st_vec, omit_zero_amps)
            else:
                print(fin_st_vec)
        print('total probability of final state vector ' +
              '(=one if no measurements)=', self.get_total_prob(fin_st_vec))
        print('dictionary with key=qubit, value=final (P(0), P(1))')
        pp.pprint(self.get_bit_probs(fin_st_vec))

    def evolve_by_controlled_bit_swap(self, bit1, bit2, controls):
        """
        Evolve each branch of cur_st_vec_list by controlled bit swap.

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
        slicex = tuple(slicex)

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
        for br in range(len(self.cur_st_vec_list)):
            vec = self.cur_st_vec_list[br][slicex]
            self.cur_st_vec_list[br][slicex] = vec.transpose(perm)

    def evolve_by_controlled_one_bit_gate(self,
                tar_bit_pos, controls, one_bit_gate):
        """
        Evolve each branch of cur_st_vec_list by controlled one bit gate (
        from class OneBitGates). Note one_bit_gate is entered as np.ndarray.

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
        assert tar_bit_pos not in controls.bit_pos
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
        for br in range(len(self.cur_st_vec_list)):
            vec = self.cur_st_vec_list[br][vec_slicex]
            # Axes 1 of one_bit_gate and new_tar of vec are summed over. Axis
            #  0 of one_bit_gate goes to the front of all the axes of new vec.
            # Use transpose() to realign axes.
            vec = np.tensordot(one_bit_gate, vec, ([1], [new_tar]))
            self.cur_st_vec_list[br][vec_slicex] = np.transpose(vec, axes=perm)

    def finalize_next_line(self):
        """
        Prints running documentary as the end of the reading of each line.

        Returns
        -------

        """

        if self.verbose:
            print('\n')
            print(self.split_line)
            print('line number = ', self.line_count)
            print('operation = ', self.num_ops)
            for br in range(len(self.cur_st_vec_list)):
                print('#----------- BRANCH ' + str(br) + ':')
                print('tot_prob = ',
                      SEO_simulator.get_total_prob(self.cur_st_vec_list[br]))
                print('bit probs = ',
                      SEO_simulator.get_bit_probs(self.cur_st_vec_list[br]))
        
    def use_NOTA(self, bla_str):
        """
        Overrides the parent class use_ function. Does nothing.

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        pass

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

    def use_MEAS(self, kind, tar_bit_pos):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for MEAS.

        For kind 0 (resp., 1) measurements, it applies |0><0| (resp.,
        |1><1|) to each branch of cur_st_vec_list.

        For kind 3 measurements, it first creates a list concatenating x
        plus a deep copy of x, where x is the cur_st_vec_list. Next,
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
        list_len = len(self.cur_st_vec_list)
        # slicex = slice index
        slicex = [slice(None)]*self.num_bits
        # br = branch
        if kind == 0:
            for br in range(list_len):
                slicex[tar_bit_pos] = 1
                # set projection |1><1| to zero
                self.cur_st_vec_list[br][tuple(slicex)] = 0
                slicex[tar_bit_pos] = slice(None)
        elif kind == 1:
            for br in range(list_len):
                slicex[tar_bit_pos] = 0
                # set projection |0><0| to zero
                self.cur_st_vec_list[br][tuple(slicex)] = 0
                slicex[tar_bit_pos] = slice(None)
        elif kind == 2:
            self.cur_st_vec_list += cp.deepcopy(self.cur_st_vec_list)
            for br in range(list_len):
                slicex[tar_bit_pos] = 1
                # set projection |1><1| to zero
                self.cur_st_vec_list[br][tuple(slicex)] = 0
                slicex[tar_bit_pos] = slice(None)
            for br in range(list_len, 2*list_len):
                slicex[tar_bit_pos] = 0
                # set projection |0><0| to zero
                self.cur_st_vec_list[br][tuple(slicex)] = 0
                slicex[tar_bit_pos] = slice(None)
        else:
            assert False, 'unsupported measurement kind'

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

    def use_SIG(self, direction, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for sigx, sigy, sigz.

        Parameters
        ----------
        direction : int
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
        gate = fun[direction]()
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

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

    def use_ROT(self, direction,
                angle_degs, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for rot along axes x, y, or z.

        Parameters
        ----------
        direction : int
            1, 2, 3 for x, y, z
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneBitGates.rot_ax(angle_degs * np.pi/180, direction)
        self.evolve_by_controlled_one_bit_gate(tar_bit_pos, controls, gate)

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_bit_gate() for rot along arbitrary direction.


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
            "raw multiplexors. Work around: use first our" \
            "MultiplexorExpander application to expand " \
            "multiplexors into simpler gates "

if __name__ == "__main__":

    # use test = 0 if want to run all tests at once.
    test = 0
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

