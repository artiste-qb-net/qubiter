import copy as cp
# import pprint as pp
from qubiter.SEO_reader import *
from qubiter.OneQubitGate import *
from qubiter.StateVec import *
# import utilities_gen as ut

import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np


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
    measurement inserts a projector ``|0><0| = n = P_0`` at the target bit.
    A type 1 measurement inserts a projector ``|1><1| = nbar = P_1`` at the
    target bit. A type 2 measurement stores a copy of the state vector after
    ``|0><0|`` has been applied, and another copy after ``|1><1|`` has been
    applied.

    self.cur_st_vec_dict is a dictionary of strings (called branch keys) to
    state vectors StateVec on num_qbits qubits. We will refer to each state
    vec in the dict as a branch. Initially, this dict contains a single
    branch with branch key = "pure". A measurement MEAS of kinds 0 or 1 does
    not change the number of branches in the dict, but a measurement of kind
    2 doubles their number.

    Note that since projectors are not unitary matrices, the branches of
    cur_st_vec_dict are not expected have normalized state vectors as values
    except when there is only a single branch.

    If cur_st_vec_dict contains as values the states ``|br0>, |br1>, |br2>,
    ...``, then one can construct the density matrix of that state as ``rho
    = |br0><br0| + |br1><br1| + |br2><br2| + ...`` divided by a number so
    that trace(rho)=1. In other words, cur_st_vec_dict is just a particular
    way of storing the density matrix of a state. A state with a single
    branch is a pure state, but a state with more than one branch may not be.

    An item of cur_st_vec_dict may be key=some string, value=None. This
    means the state vector of that branch is zero.

    Attributes
    ----------
    cached_sts : dict[int, dict(str, StateVec)]
        A dictionary mapping an int to past values of self.cur_st_vec_dict.
        Used by use_PRINT() sometimes.
    cur_st_vec_dict : dict(str, StateVec)
        dictionary with key= branch_key string and value= StateVec|None. If
        there is a single item in dict, the dict represents a pure state and
        the key is "pure". If there is more than one item, the dict
        represents a mixed state and each branch key is a string that
        uniquely characterizes the measured controls. For example, if it has
        been measured previously (type 2 measurement only) that qubit 2 is
        True and qubit 4 is False, the branch key will be '4F2T'.
    lib : str
        tensor library. Either 'np' for numpy or 'tf' for tensorflow
    tensordot : function
    transpose : function
    use_tf : bool
        True iff using TensorFlow Eager. False by default

    """
    # class variables
    transpose = np.transpose
    tensordot = np.tensordot
    reshape = np.reshape

    # rrtucci: combines my java classes:
    # LineList, UnitaryMat, SEO_readerMu

    def __init__(self, file_prefix, num_qbits,
                 init_st_vec=None, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_qbits : int
        init_st_vec : StateVec
            Get this using the functions StateVec.get_ground_st_vec() or
            StateVec.get_standard_basis_st_vec().

        Returns
        -------

        """
        if StateVec.is_zero(init_st_vec):
            init_st_vec = StateVec.get_ground_st_vec(num_qbits)
        self.cur_st_vec_dict = {"pure": init_st_vec}
        self.cached_sts = {}
        self.use_tf = False
        self.lib = 'np'

        self.do_more_init_before_reading()

        SEO_reader.__init__(self, file_prefix, num_qbits, **kwargs)

    def do_more_init_before_reading(self):
        """
        Stub. Called in SEO_simulator.__init__ immediately before calling
        SEO_reader.__init__

        Returns
        -------
        None

        """
        pass

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
        a str key for self.cur_st_vec_dict.

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
        trols = Controls(self.num_qbits)
        trols.bit_pos_to_kind = dict(zip(bit_pos, kinds))
        trols.refresh_lists()
        return trols

    @staticmethod
    def get_br_key_with_new_link(br_key, new_bit_pos, new_kind):
        """
        Say new_bit_pos=2 and new_kind=True. This returns a new branch key
        with '2T' added to the end of the string br_key.

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

    def evolve_by_controlled_qbit_swap(self, bit1, bit2, controls):
        """
        Evolve each branch of cur_st_vec_dict by controlled bit swap iff the
        bit swap line is (1) outside of any IF_M block, or (2) it is inside
        such a block, and it satisfies self.mcase_trols.

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
            assert -1 < bit < self.num_qbits
            assert bit not in controls.bit_pos
 
        slicex = [slice(None)]*self.num_qbits
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
        for bit in range(self.num_qbits):
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
            if StateVec.is_zero(self.cur_st_vec_dict[br_key]):
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
                sub_arr = self.cur_st_vec_dict[br_key].arr[slicex]
                sub_arr = SEO_simulator.transpose(sub_arr, perm)

                # can't do array assignments with autograd or tensorflow
                # eager so achieve same result with other allowed tensor ops
                if self.use_tf or 'autograd.numpy' in sys.modules:
                    self.do_array_assignment_workaround(
                        br_key, slicex, sub_arr)
                    return

                self.cur_st_vec_dict[br_key].arr[slicex] = sub_arr

    def evolve_by_controlled_one_qbit_gate(self,
                tar_bit_pos, controls, one_qbit_gate):
        """
        Evolve each branch of cur_st_vec_dict by controlled one bit gate (
        from class OneQubitGate) iff the controlled one bit gate line is (1)
        outside of an IF_M block, or (2) it is inside such a block, and it
        satisfies self.mcase_trols. Note one_qbit_gate is entered as
        np.ndarray.

        Parameters
        ----------
        tar_bit_pos : int
            bit position of target of one bit gate.
        controls : Controls
        one_qbit_gate : np.ndarray

        Returns
        -------
        None

        """
        assert -1 < tar_bit_pos < self.num_qbits

        vec_slicex = [slice(None)]*self.num_qbits
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
        for bit in range(self.num_qbits):
            if bit == tar_bit_pos:
                new_tar = scout
            if bit not in controls.bit_pos:
                scout += 1

        perm_len = scout
        # example tar = 2
        # want to map [2, 0, 1, 3, 4] back to [0, 1, 2, 3, 4]

        # this didn't work
        # use perm 0=>2, 1=>0, 2=>1, 3=>3, 4=>4
        # perm = [new_tar] + list(range(new_tar))
        # perm += list(range(new_tar+1, perm_len))

        # use perm 2=>0, 0=>1, 1=>2, 3=>3, 4=>4
        perm = list(range(1, new_tar+1)) + [0]
        perm += list(range(new_tar+1, perm_len))

        # br = branch
        assert not(self.mcase_trols and not self.measured_bits)
        for br_key in self.cur_st_vec_dict.keys():
            if StateVec.is_zero(self.cur_st_vec_dict[br_key]):
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
                sub_arr = self.cur_st_vec_dict[br_key].arr[vec_slicex]
                # Axes 1 of one_qbit_gate and new_tar of vec are summed over.
                #  Axis 0 of one_qbit_gate goes to the front of all the axes
                # of new vec. Use transpose() to realign axes.
                sub_arr = SEO_simulator.tensordot(one_qbit_gate, sub_arr,
                                         ([1], [new_tar]))
                sub_arr = SEO_simulator.transpose(sub_arr, perm)

                # can't do array assignments with autograd or tensorflow
                # eager so achieve same result with other allowed tensor ops
                if self.use_tf or 'autograd.numpy' in sys.modules:
                    self.do_array_assignment_workaround(
                        br_key, vec_slicex, sub_arr)
                    return

                # original, if autograd is not being used
                self.cur_st_vec_dict[br_key].arr[vec_slicex] = sub_arr

    def do_array_assignment_workaround(self, br_key, slicex, sub_arr):
        """
        Internal function used in evolve_ methods iff autograd is on or
        use_tf is True. Should have same effect as

        self.cur_st_vec_dict[br_key].arr[slicex] = sub_arr

        Parameters
        ----------
        br_key : str
        slicex : tuple
        sub_arr : np.ndarray

        Returns
        -------
        None

        """
        test = False
        arr = self.cur_st_vec_dict[br_key].arr

        if test:
            # This creates a numpy copy of arr if arr is numpy.
            # If arr is tf, it converts to numpy then creates a numpy copy
            arr1_np = np.array(arr)
            # print('arr1 bef', arr1_np)
            arr1_np[slicex] = np.array(sub_arr)
            # print('arr1 aft', arr1_np)
        on_slicex = np.full(tuple(arr.shape), False)
        on_slicex[slicex] = True
        bigger_shape = [1]*self.num_qbits  # slicex is num_qbits long
        k = 0
        # print('wwwww', sub_arr.shape, slicex)
        for bit, kind in enumerate(slicex):
            if kind not in [0, 1]:
                bigger_shape[bit] = int(sub_arr.shape[k])
                k += 1
        one_on_slicex = on_slicex.astype(int)
        not_one_on_slicex = np.logical_not(on_slicex).astype(int)

        sub_arr = SEO_simulator.reshape(sub_arr, tuple(bigger_shape))

        arr = arr*not_one_on_slicex + sub_arr*one_on_slicex
        self.cur_st_vec_dict[br_key].arr = arr

        if test:
            # print('arr aft', arr)
            print('testing simulator ruse')
            assert np.linalg.norm(np.array(arr)-arr1_np) < 1e-6, \
                'sim ruse test failed'

    def convert_tensors_to_numpy(self):
        """
        This method is meant to be overridden to replace tensorflow tensors
        (or some other non-numpy tensor type) by numpy tensors.

        Returns
        -------

        """
        pass

    def convert_tensors_to_tf(self):
        """
        This method is meant to be overridden to replace numpy tensors by
        tensorflow tensors (or some other non-numpy tensor type).

        Returns
        -------

        """
        pass

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
            self.describe_st_vec_dict(  # print_st_vec=True,
                                    show_pp_probs=True)

    def describe_st_vec_dict(self, **kwargs):
        """
        Calls method with same name in class StateVec. It prints a
        description of the current state vector dictionary. Can call this at
        the end, after running through whole circuit.

        Parameters
        ----------
        kwargs : dict[]

        Returns
        -------
        None

        """
        self.convert_tensors_to_numpy()
        StateVec.describe_st_vec_dict(self.cur_st_vec_dict,
                                             **kwargs)
        self.convert_tensors_to_tf()

    def get_counts(self, num_shots, omit_zero_counts=True,
                   use_bin_labels=True, rand_seed=None):
        """
        This method calculates a probability distribution that we call pd
        from the current state vector if it is pure. (If the state vec is
        not pure, it calculates a density matrix from the
        self.cur_st_vec_dict. Then it extracts the diagonal of that density
        matrix. That diagonal must be a probability distribution that we
        call pd.) Then the method samples pd, num_shots times. The method
        returns the result of that sampling as an OrderedDict
        state_name_to_counts. Depending on the value of the flag
        use_bin_labels, the state names are a string '0', '1', '2', etc,
        or their binary representations followed by 'ZL', because the ZL
        convention is assumed.

        Parameters
        ----------
        num_shots : int
        omit_zero_counts : bool
        use_bin_labels : bool
        rand_seed : int

        Returns
        -------
        OrderedDict[str, int]

        """
        self.convert_tensors_to_numpy()

        if len(self.cur_st_vec_dict) == 1:
            # print('..,,mm', 'was here')
            pd = self.cur_st_vec_dict['pure'].get_pd()
        else:
            den_mat = StateVec.get_den_mat(self.num_qbits,
                                           self.cur_st_vec_dict)
            pd = StateVec.get_den_mat_pd(den_mat)
        # print('....,,,', pd.shape)
        obs_vec = StateVec.get_observations_vec(
            self.num_qbits, pd, num_shots, rand_seed)

        out = StateVec.get_counts_from_obs_vec(self.num_qbits, obs_vec,
                        use_bin_labels, omit_zero_counts)
        self.convert_tensors_to_tf()

        return out

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
        evolve_by_controlled_one_qbit_gate() for had2.


        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneQubitGate.had2(lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)

    def use_IF_M_beg(self, controls):
        """
        Do nothing.

        Parameters
        ----------
        controls : Controls

        Returns
        -------
        None

        """
        pass

    def use_IF_M_end(self):
        """
        Do nothing.

        Parameters
        ----------

        Returns
        -------
        None

        """
        pass

    def use_MEAS(self, tar_bit_pos, kind):
        """
        Overrides the parent class use_ function.

        For kind 0 (resp., 1) measurements, it applies ``P_0=|0><0|`` (
        resp., ``P_1=|1><1|``) to each branch of cur_st_vec_dict.

        For kind 2 measurements, it first doubles the number of branches in
        cur_st_vec_dict by adding a deep copy of each branch. Next,
        it applies P_0 to half of the branches of the dict and P_1 to the
        other half.

        Parameters
        ----------
        kind : int
        tar_bit_pos: int

        Returns
        -------
        None

        """
        self.convert_tensors_to_numpy()
        # slicex = slice index
        slicex = [slice(None)]*self.num_qbits
        # br = branch
        if kind in [0, 1]:
            b = 1 if kind == 0 else 0
            for br_key in self.cur_st_vec_dict:
                st_vec = self.cur_st_vec_dict[br_key]
                if not StateVec.is_zero(st_vec):
                    slicex[tar_bit_pos] = b
                    # set projection |b=0><b=0| to zero for kind=1
                    st_vec.arr[tuple(slicex)] = 0
                    tot_prob = st_vec.get_total_prob()
                    if tot_prob < 1e-8:
                        # this didn't work
                        # st_vec.arr = None
                        self.cur_st_vec_dict[br_key].arr = None
                    slicex[tar_bit_pos] = slice(None)
                # self.cur_st_vec_dict[br_key].arr = st_vec.arr
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
                    if not StateVec.is_zero(st_vec):
                        slicex[tar_bit_pos] = b
                        st_vec.arr[tuple(slicex)] = 0
                        tot_prob = st_vec.get_total_prob()
                        # print('tot_prob=', tot_prob)
                        if tot_prob < 1e-8:
                            # this didn't work
                            # st_vec.arr = None
                            self.cur_st_vec_dict[new_key].arr = None
                        slicex[tar_bit_pos] = slice(None)
                    # print(st_vec)

            # print(self.cur_st_vec_dict)
        else:
            assert False, 'unsupported measurement kind'
        self.convert_tensors_to_tf()

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

    def use_NOTA(self, bla_str):
        """
        Do nothing.

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        pass

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() for PHAS.

        Parameters
        ----------
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneQubitGate.phase_fac(angle_rads, lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)

    def use_P_PH(self, projection_bit,
                angle_rads, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() for P_0 and P_1 phase factors.


        Parameters
        ----------
        projection_bit : int
            0 (resp. 1) for P_0 (resp. P_1) projection
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        fun = {
            0: OneQubitGate.P_0_phase_fac,
            1: OneQubitGate.P_1_phase_fac
        }
        gate = fun[projection_bit](angle_rads, lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)

    def use_PRINT(self, style, line_num):
        """
        Prints to screen a description of self.cur_st_vec_dict.

        Parameters
        ----------
        style : str
            style in which to print
        line_num : int
            line number in eng & pic files in which PRINT command appears

        Returns
        -------
        None

        """
        self.convert_tensors_to_numpy()
        print("\n*************************beginning PRINT output")
        print("PRINT line number=" + str(line_num))
        st_vecs = self.cur_st_vec_dict
        StateVec.describe_st_vec_dict(st_vecs,
                                      **StateVec.get_style_dict(style))
        # must store copy or it will change
        self.cached_sts[line_num] = cp.deepcopy(st_vecs)
        print("****************************ending PRINT output")
        self.convert_tensors_to_tf()

    def use_ROTA(self, axis,
                 angle_rads, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() for rot along axes x, y, or z.

        Parameters
        ----------
        axis : int
            1, 2, 3 for x, y, z
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        # print('//////', angle_rads, axis)
        gate = OneQubitGate.rot_ax(angle_rads, axis, lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() for rot along arbitrary axis.


        Parameters
        ----------
        angle_x_rads : float
        angle_y_rads : float
        angle_z_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneQubitGate.rot(angle_x_rads, angle_y_rads, angle_z_rads,
                               lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() for sigx, sigy, sigz.

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
            1: OneQubitGate.sigx,
            2: OneQubitGate.sigy,
            3: OneQubitGate.sigz
        }
        gate = fun[axis](lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)

    def use_SWAP(self, bit1, bit2, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_qbit_swap().

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        self.evolve_by_controlled_qbit_swap(bit1, bit2, controls)

    def use_SWAY(self, bit1, bit2, controls, rads_list):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() 3 times.

        This relies on the fact that

        SWAY(0, 1) = SWAY(1, 0) =
        X---@
        @---U2
        X---@

        U2 = exp(j*(rads0 + rads1*sig_x))

        rads_list = [rads0, rads1]

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls
        rads_list : list[float]

        Returns
        -------
        None

        """
        controls1 = Controls.copy(controls)
        controls1.set_control(bit1, True, do_refresh=True)
        # print('m,m,', bit1, bit2)
        # print(",m,m", controls1.bit_pos_to_kind)

        controls2 = Controls.copy(controls)
        controls2.set_control(bit2, True, do_refresh=True)
        # print(",m,m", controls2.bit_pos_to_kind)

        assert len(rads_list) == 2
        rads0, rads1 = rads_list

        self.evolve_by_controlled_one_qbit_gate(
            bit2, controls1, OneQubitGate.sigx(lib=self.lib))
        self.evolve_by_controlled_one_qbit_gate(
            bit1, controls2, OneQubitGate.u2(rads0, rads1, 0.0, 0.0,
                                            lib=self.lib))
        self.evolve_by_controlled_one_qbit_gate(
            bit2, controls1, OneQubitGate.sigx(lib=self.lib))

    def use_U_2_(self, rads0, rads1, rads2, rads3,
                tar_bit_pos, controls):
        """
        Overrides the parent class use_ function. Calls
        evolve_by_controlled_one_qbit_gate() for arbitrary unitary 2-dim
        matrix.

        Parameters
        ----------
        rads0 : float
        rads1 : float
        rads2 : float
        rads3 : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        gate = OneQubitGate.u2(rads0, rads1, rads2, rads3, lib=self.lib)
        self.evolve_by_controlled_one_qbit_gate(tar_bit_pos, controls, gate)


if __name__ == "__main__":
    def main():
        # use test = 0 if want to run all tests at once.
        test = 0
        # test = 3
        if test in [0, 1]:
            # test on circuit for a quantum fourier transform
            # (no loops, no internal measurements)
            sim = SEO_simulator('sim_test1', 6, verbose=True)

        if test in [0, 2]:
            # test embedded loops
            sim = SEO_simulator('sim_test2', 4, verbose=True)

        if test in [0, 3]:
            # test MEAS branching. Each kind 2 measurement doubles number of
            # branches
            sim = SEO_simulator('sim_test3', 4, verbose=True)

    main()
