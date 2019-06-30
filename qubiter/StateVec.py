import pprint as pp
import scipy as sc
import qubiter.utilities_gen as ut
from collections import OrderedDict
import pandas as pan
from qubiter.Plotter import *

import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np
else:
    import autograd.numpy as np


class StateVec:
    """
    This class is a wrapper for its main attribute, a complex numpy array
    self.arr with shape [2]*num_bits. The class also provides functions for
    performing calculations on dictionaries of the objects of this class
    StateVec. The keys of these dictionaries of state vectors are strings
    that we call branch_keys, because they name "branches" in class
    SEO_simulation. This class also provides a function for constructing
    from such dictionaries of state vectors, a density matrix which is a 2
    dim square numpy array of dimension 2^num_bits.

    IMPORTANT: See docstring of method get_traditional_st_vec() for
    explanation of qubit ordering conventions and shape of self.arr
    
    Attributes
    ----------
    arr : np.ndarray
         a complex array of shape [2]*num_bits
    num_bits : int

    """
    def __init__(self, num_bits, arr=None):
        """
        Constructor

        Parameters
        ----------
        num_bits : int
        arr : np.ndarray

        Returns
        -------

        """
        self.num_bits = num_bits
        self.arr = arr
        if arr is not None:
            assert self.arr.shape == tuple([2]*self.num_bits)

    @staticmethod
    def is_zero(st_vec):
        """
        Returns True iff an object of this class is None or its parameter
        'arr' is None

        Parameters
        ----------
        st_vec : StateVec|None

        Returns
        -------
        bool

        """
        return st_vec is None or st_vec.arr is None

    def __str__(self):
        """
        Returns str(self.arr)

        Returns
        -------
        str

        """
        return str(self.arr)

    @staticmethod
    def get_ground_st_vec(num_bits):
        """
        Returns StateVec for the ground state |0>|0>|0>...|0>, where |0> = [
        1,0]^t and |1> = [0,1]^t, t = transpose

        Parameters
        ----------
        num_bits : int

        Returns
        -------
        StateVec

        """
        ty = np.complex128
        assert num_bits > 0
        arr = np.zeros([1 << num_bits], dtype=ty)
        arr[0] = 1
        arr = arr.reshape([2]*num_bits)
        return StateVec(num_bits, arr)

    @staticmethod
    def get_random_st_vec(num_bits, rand_seed=None):
        """
        Returns StateVec for random state \sum_b^n A(b^n)|b^n>, b^n \in {0,
        1}^n, where n=num_bits and \sum_b^n |A( b^n)|^2 = 1

        Parameters
        ----------
        num_bits : int
        rand_seed : int

        Returns
        -------
        StateVec

        """
        if rand_seed:
            np.random.seed(rand_seed)
        # returns array of random numbers in [0, 1] interval
        mat_phi = 2*np.pi*np.random.random(1 << num_bits)
        mat_r = np.random.random(1 << num_bits)
        arr = mat_r*(np.cos(mat_phi) + 1j*np.sin(mat_phi))
        magnitude = np.linalg.norm(arr)
        arr /= magnitude
        arr = arr.reshape([2]*num_bits)
        return StateVec(num_bits, arr)

    @staticmethod
    def get_standard_basis_st_vec(spin_dir_list, ZL=True):
        """
        If ZL = True, returns StateVec for state ...|s2>|s1>|s0>,
        where spin_dir_list=[...,s2, s1, s0], s_j \in {0, 1} for all j,
        |0> = [1, 0]^t and |1> = [0,1]^t, t = transpose. If ZL = False,
        same except spin_dir_list=reversed([...,s2, s1, s0]).

        Parameters
        ----------
        spin_dir_list : list[int]
        ZL : bool
            True(False) if last(first) entry of spin_dir_list refers to
            qubit 0

        Returns
        -------
        StateVec

        """
        num_bits = len(spin_dir_list)
        arr = np.zeros([1 << num_bits], dtype=np.complex128)
        arr = arr.reshape([2]*num_bits)
        if ZL:
            spin_dir_list = reversed(spin_dir_list)

        # print("spins", list(spin_dir_list))
        # print("arr", arr.shape)
        arr[tuple(spin_dir_list)] = 1
        return StateVec(num_bits, arr)

    def get_traditional_st_vec(self):
        """

        **IMPORTANT: Internally, self.arr in Qubiter has shape [2]*num_bits
        and assumes ZF convention because that way a numpy axis and a qubit
        number are the same thing. However, the traditional way of writing a
        state vector is as a column array of dimension 1<< num_bits in the
        ZL convention.**

        This function returns the traditional view. So it reshapes (
        flattens) the array and it reverses the axes (reversing axes takes
        it from ZF to ZL).

        The rows are always labelled 0, 1, 2, 3, ... or the binary
        representation thereof, regardless of whether ZL or ZF convention.
        One can go from digital to binary labels and vice versa
        using

        Examples
        --------
        >>> x = np.binary_repr(3, width=4)
        >>> x
        '0011'
        >>> int(x, 2)
        3

        Parameters
        ----------

        Returns
        -------
        np.array

        """
        perm = list(reversed(range(self.num_bits)))
        return np.transpose(self.arr, perm).flatten()

    @staticmethod
    def get_den_mat(num_bits, st_vec_dict):

        """
        Returns a density matrix (indexed in ZL convention) constructed from
        st_vec_dict which is a dict from strings to StateVec.

        The rows and columns are always labelled 0, 1, 2, .. or binary
        representation thereof, regardless of whether ZL or ZF convention.
        To switch between bin to dec representations of labels,
        see docstring of get_traditional_st_vec().

        Parameters
        ----------
        num_bits : int
        st_vec_dict : dict[str, StateVec]

        Returns
        -------
        np.ndarray

        """

        dim = 1 << num_bits
        den_mat = np.zeros((dim, dim), dtype=complex)
        # print(",,,", den_mat)
        for br_key in st_vec_dict:
            if StateVec.is_zero(st_vec_dict[br_key]):
                continue
            vec = st_vec_dict[br_key].get_traditional_st_vec()
            assert vec.shape == (dim,)
            den_mat += np.outer(vec, np.conj(vec))
            # print(',,..', br_key, vec)
        tr = np.trace(den_mat)
        assert abs(tr) > 1e-6
        return den_mat/tr

    @staticmethod
    def get_partial_tr(num_bits, den_mat, traced_bits_set):
        """
        Returns the partial trace of a density matrix den_mat. Traces over
        qubits in set traced_bits_set. To get full trace, just do np.trace(
        den_mat)

        Parameters
        ----------
        num_bits : int
        den_mat : np.ndarray
            if dim=2^num_bits, this function assumes that den_mat has shape
            (dim, dim) and that it's indexed in the ZL convention so qubit 0
            corresponds to axis num_bits-1.
        traced_bits_set : set[int]
             Set of qubits being traced over

        Returns
        -------
        np.ndarray

        """
        dim = 1 << num_bits
        assert den_mat.shape == (dim, dim)
        assert set(range(num_bits)) > traced_bits_set
        dm = den_mat.reshape([2]*(2*num_bits))
        # bit 0 corresponds to axis num_bits - 1
        traced_axes = [num_bits - 1 - k for k in traced_bits_set]
        num_traces = len(traced_axes)
        for k in range(num_traces):
            ax = traced_axes.pop(0)
            dm = np.trace(dm, axis1=ax, axis2=ax + num_bits - k)
            traced_axes = list(map(lambda x: (x if x <= ax else x-1),
                                           traced_axes))
        new_num_bits = num_bits - len(traced_bits_set)
        dim = 1 << new_num_bits
        dm = dm.reshape((dim, dim))
        return dm

    @staticmethod
    def get_impurity(den_mat):
        """
        Returns abs(trace(den_mat^2) -1). This is zero iff the density
        matrix den_mat represents a pure state. For example, for a pure
        state den_mat = |a><a|, den_mat^2 = den_mat = |a><a| so this
        quantity is indeed zero.

        Parameters
        ----------
        den_mat : np.nparray
            density matrix, shape=(dim, dim) where dim=2^num_bits

        Returns
        -------
        float

        """
        return abs(np.trace(np.dot(den_mat, den_mat)) - 1)

    @staticmethod
    def get_entropy(den_mat, method='eigen'):
        """
        Returns entropy of density matrix den_mat. Uses natural log for
        entropy.

        Parameters
        ----------
        den_mat : np.ndarray
            Density matrix. Eigenvalues must be non-negative and sum to 1
        method : str
            method used to calculate log of array. Either 'eigen' or 'pade'

        Returns
        -------
        float

        """
        ent = 0.0
        if method == 'eigen':
            evas = np.real(np.linalg.eigvalsh(den_mat))
            assert np.all(evas > -1e-6), evas
            assert abs(np.sum(evas) - 1) < 1e-6
            for val in evas:
                if val > 1e-6:
                    ent += - val*np.log(val)
        elif method == 'pade':
            ent = - np.trace(np.dot(den_mat, sc.linalg.logm(den_mat)))
        else:
            assert False, 'unsupported method for ' + \
                          'calculating entropy of a density matrix.'
        return ent

    @staticmethod
    def get_den_mat_pd(den_mat):
        """
        Returns the diagonal of den_mat (so indexed in ZL convention) .
        den_mat is expected to be a density matrix returned by get_den_mat()

        Parameters
        ----------
        den_mat : np.ndarray
            density matrix, shape=(dim, dim) where dim=2^num_bits, indexed
            in ZL convention.

        Returns
        -------
        np.ndarray

        """
        return np.real(np.diag(den_mat))

    def get_pd(self):
        """
        Returns copy of self.get_traditional_st_vec() with amplitudes
        replaced by probabilities. pd = probability distribution. So returns
        one column array indexed in ZL convention like the traditional state
        vec is. Doesn't check that the resulting array sums to 1.

        Parameters
        ----------

        Returns
        -------
        np.ndarray
            probability distribution of shape (2^num_bits,) IMP: will
            be indexed in ZL convention

        """
        x = self.get_traditional_st_vec()
        return np.real(x*np.conj(x))

    def get_mean_value_of_real_diag_mat(self, real_arr):
        """
        In Quantum Mechanics, one often needs to calculate the mean value of
        a Hermitian operator H, mean= <psi|H|psi>. Let H = U^\dag D U,
        where U is unitary and D is real diagonal matrix. If self = U|psi>,
        then this reduces to finding mean= <self|D|self>. So must decompose
        U into a SEO and evolve, using SEO_simulator, to the state U|psi>.

        Parameters
        ----------
        real_arr : np.ndarray
             a real array of shape=[2]^num_bits (same shape as self.arr). If
             flattened, real_arr contains the diagonal of the matrix D. If U
             is a Kronecker prod of 2-dim unitary matrices, the flattened
             real_arr can be obtained as Kronecker product of spinors, i.e.,
             shape=( 2, ) arrays.

        Returns
        -------
        float

        """
        return np.sum(np.real(np.conj(self.arr) * real_arr * self.arr))

    def get_total_prob(self):
        """
        Returns total probability of self.

        Parameters
        ----------

        Returns
        -------
        float

        """
        return np.sum(np.real(np.conj(self.arr)*self.arr))

    @staticmethod
    def get_observations_vec(num_bits, pd, num_shots, rand_seed=None):
        """
        vec = vector

        num_shots (number of shots) is often called number of trials or
        number of samples.

        For num_shots=1, this method returns an int (actually, a 1 X 1 array
        with an int in it) in range(1<<num_bits) chosen according to the
        probability distribution pd for num_bits qubits. If the output int
        were to be expressed in binary notation, its last, rightmost bit
        would be the measurement of the 0th qubit (because pd is assumed to
        be in ZL convention).

        For num_shots>1, the method returns an np.ndarray of shape (
        num_shots,) with the result of doing num_shots repetitions of what
        was done for num_shots=1.

        Does not assume that pd is normalized to 1.

        Parameters
        ----------
        num_bits : int
        pd : np.ndarray
            probability distribution of shape (2^num_bits,) IMP: assumed to
            be indexed in ZL convention
        num_shots : int
        rand_seed : int

        Returns
        -------
        np.ndarray
            shape (num_shots,)

        """
        if rand_seed:
            np.random.seed(rand_seed)
        len_pd = 1 << num_bits
        assert pd.shape == (len_pd,)
        tot_prob = np.sum(pd)
        p = pd
        if abs(tot_prob-1) > 1e-5:
            p = pd/tot_prob
        return np.random.choice(np.arange(0, len_pd), size=num_shots, p=p)

    @staticmethod
    def get_counts_from_obs_vec(num_bits, obs_vec,
                    use_bin_labels=True, omit_zero_counts=True):
        """
        This method takes as input an observations vector obs_vec such as
        returned by another method in this class, namely
        get_observations_vec(). This method returns an OrderedDict called
        state_name_to_count that maps the names of states to the number of
        times they occur in obs_vec. If use_bin_labels=True, state names are
        a string composed of a binary number that is num_bits long, followed
        by 'ZL' because ZL convention is assumed. If use_bin_labels=False,
        state names are '0', '1', '2', etc.

        Parameters
        ----------
        num_bits : int
        obs_vec : np.ndarray
        use_bin_labels : bool
        omit_zero_counts : bool

        Returns
        -------
        OrderedDict[str, int]

        """
        obs_list = list(obs_vec)
        num_states = 1 << num_bits
        state_name_to_count = OrderedDict()
        for s in range(num_states):
            s_count = obs_list.count(s)
            if use_bin_labels:
                # this returns string
                key = np.binary_repr(s, width=num_bits)
                key += 'ZL'
            else:
                key = str(s)
            if s_count > 0 or not omit_zero_counts:
                state_name_to_count[key] = s_count
        return state_name_to_count

    @staticmethod
    def get_empirical_pd_from_counts(num_bits, state_name_to_count):
        """
        This method takes as input "the counts dict" (i.e., an OrderedDict
        called state_name_to_count which is produced by another method in
        this class, namely get_counts_from_obs_vec()). This method returns
        an empirical probability distribution emp_pd calculated from the
        counts dict. emp_pd indices are ints referring to qubit states
        labelled in the ZL convention.

        Parameters
        ----------
        num_bits : int
        state_name_to_count : OrderedDict[str, int]

        Returns
        -------
        emp_pd : np.ndarray
            its shape is (1<<num_bits,)

        """
        emp_pd = np.zeros(shape=(1 << num_bits,), dtype=float)
        tot_counts = 0
        for st_name, count in state_name_to_count.items():
            # state name ends in ZL so trim last two chars
            pos = int(st_name[:-2], 2)
            emp_pd[pos] = count
            tot_counts += count
        return emp_pd/tot_counts

    @staticmethod
    def get_emp_state_vec_from_emp_pd(num_bits, emp_pd):
        """
        This method takes as input an empirical probability distribution
        emp_pd and it returns an empirical state vector calculated from
        emp_pd. This requires reshaping emp_pd to the shape [2]*num_bits,
        permuting its indices from the ZL to the ZF convention, and then
        taking the sqrt of the components to get an amplitude instead of a
        probability. All amplitudes of the output state vector are real
        though.

        Parameters
        ----------
        num_bits : int
        emp_pd : np.ndarray
            its shape is (1<<num_bits,)

        Returns
        -------
        StateVec

        """
        assert emp_pd.shape == (1 << num_bits,)
        arr = emp_pd.reshape(tuple([2]*num_bits))
        perm = list(reversed(range(num_bits)))
        arr = np.transpose(arr, perm)
        sqrt_probs = np.sqrt(arr)
        return StateVec(num_bits, sqrt_probs)

    @staticmethod
    def get_bit_probs(num_bits, pd):
        """
        Returns a list whose jth item is, for the jth qubit, the pair (p,
        1-p), where p is the probability that the jth qubit is 0, if the
        state of all other qubits is ignored.

        Does not assume that pd is normalized to 1.

        Parameters
        ----------
        num_bits : int
        pd : np.ndarray
            probability distribution of shape (2^num_bits,) IMP: assumed to
            be indexed in ZL convention

        Returns
        -------
        list[tuple[float, float]]

        """
        assert pd.shape == (1 << num_bits,)
        probs = []
        arr = pd.reshape([2] * num_bits)
        # tot_prob may not be one
        # if a measurement has been done
        tot_prob = np.sum(arr)
        if abs(tot_prob-1) > 1e-5:
            arr /= tot_prob
        # slicex is a portmanteau of slice index
        # print("state_vec_pd=\n", vec)
        slicex = [slice(None)]*num_bits
        # pd is assumed to be in ZL convention
        # so when reshape, bit k has axis= num_bits -1 - k
        for k in range(num_bits):
            k_axis = num_bits - 1 - k
            slicex[k_axis] = 0
            p = np.sum(arr[tuple(slicex)])
            probs.append((p, 1-p))
            slicex[k_axis] = slice(None)  # restore to all entries slice(None)
        # print(probs)
        return probs

    # this method is correct but superfluous, I think
    # @staticmethod
    # def get_bit_counts(bit_probs, num_shots, rand_seed=None):
    #     """
    #     num_shots (number of shots) is often called number of trials or
    #     number of samples.
    #
    #     Returns a list whose jth item is, for the jth qubit, the pair (
    #     count0, count1), where count0 + count1 = num_shots, and as num_shots
    #     ->infinity, (count0, count1)/num_shots tends to the probability pair
    #     bit_probs[j].
    #
    #     Parameters
    #     ----------
    #     bit_probs : list[tuple[float, float]]
    #     num_shots : int
    #     rand_seed : int
    #
    #     Returns
    #     -------
    #     list[tuple[int, int]]
    #
    #     """
    #     if rand_seed:
    #         np.random.seed(rand_seed)
    #     num_bits = len(bit_probs)
    #     counts = []
    #     for bit in range(num_bits):
    #         x1 = np.random.binomial(n=num_shots, p=bit_probs[bit][1])
    #         x0 = num_shots - x1
    #         counts.append((x0, x1))
    #     return counts

    def pp_arr_entries(self, omit_zero_amps=False,
                          show_pp_probs=False, ZL=True):
        """
        pp=pretty print. Prints for each entry of self.arr, a line of the
        form (i, j, k, ...) self.arr[i, j, k, ...], with zero bit last (
        resp., first) if ZL=True (resp., False).

        Parameters
        ----------
        omit_zero_amps : bool
            If True, will not list states with zero amplitude
        show_pp_probs : bool
            If True, will show probability of each amplitude

        ZL : bool
            If True, multi-index in ZL (Zero bit Last) convention. If False,
            multi-index in ZF (Zero bit First) convention.

        Returns
        -------
        None

        """
        for ind, x in np.ndenumerate(self.arr):
            index = ind
            label = 'ZF'
            if ZL:
                index = ind[::-1]  # this reverses order of tuple
                label = 'ZL'
            # turn index tuple into string and remove commas
            ind_str = str(index).replace(', ', '') + label
            mag = np.absolute(x)
            extra_str = ''
            if show_pp_probs:
                extra_str = ', prob=' + str(mag**2)
            if omit_zero_amps:
                if mag > 1E-6:
                    print(ind_str, x, extra_str)
            else:
                print(ind_str, x, extra_str)

    @staticmethod
    def get_style_dict(style):
        """
        Given a style string as input, this method returns a dict mapping
        various strings denoting parameters of the method
        StateVec.describe_self() to their bool values for the input style.

        Parameters
        ----------
        style : str

        Returns
        -------
        dict[str, bool]

        """
        vanilla = {
            'print_st_vec': False,
            'do_pp': False,
            'omit_zero_amps': False,
            'show_pp_probs': False,
            'ZL': True,
            'plot_st_vec_pd': False}
        if style == 'V1':
            out = vanilla
        elif style == 'ALL':
            out = {x: True for x in vanilla.keys()}
            out['plot_st_vec_pd'] = False
        elif style == 'ALL+':
            out = {x: True for x in vanilla.keys()}
        else:
            assert False, "unsupported style for StateVec.describe_self()"
        return out

    def describe_self(self, print_st_vec=False, do_pp=False,
                      omit_zero_amps=False, show_pp_probs=False, ZL=True,
                      plot_st_vec_pd=False):
        """
        Prints a description of self.

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
        show_pp_probs : bool
            If True, will show probability of each standard basis state
        ZL : bool
            If True, multi-index of ket in ZL (Zero bit Last) convention.
            If False, multi-index of ket in ZF (Zero bit First) convention.
        plot_st_vec_pd : bool
            If True, plots state vector's probability distribution

        Returns
        -------
        None

        """
        if self.arr is None:
            print("zero state vector")
            return
        if print_st_vec:
            print('state vector:')
            if do_pp:
                if not ZL:
                    print('ZF convention (Zero bit First in state tuple)')
                else:
                    print('ZL convention (Zero bit Last in state tuple)')
                self.pp_arr_entries(omit_zero_amps, show_pp_probs, ZL)
            else:
                print(self.arr)

        print('total probability of state vector ' +
              '(=one if no measurements)=', self.get_total_prob())

        print('dictionary with key=qubit, value=(Prob(0), Prob(1))')
        bit_probs = StateVec.get_bit_probs(self.num_bits, self.get_pd())
        pp.pprint(dict(enumerate(bit_probs)))
        if plot_st_vec_pd:
            st_vec_pd = self.get_pd()
            st_vec_pd_df = Plotter.get_pd_df(self.num_bits, st_vec_pd)
            Plotter.plot_probs_col(['st_vec_pd'], [st_vec_pd_df])

    @staticmethod
    def describe_st_vec_dict(st_vec_dict, **kwargs):
        """
        Calls describe_self() for each branch of st_vec_dict

        Parameters
        ----------
        st_vec_dict : dict[str, StateVec]
        kwargs : dict[]
            the keyword arguments of describe_self()

        Returns
        -------
        None

        """
        for br_key in st_vec_dict.keys():
            print("*********branch= " + br_key)
            if StateVec.is_zero(st_vec_dict[br_key]):
                print("zero state vector")
            else:
                st_vec_dict[br_key].describe_self(**kwargs)


if __name__ == "__main__":
    def main():
        num_bits = 3
        gs = StateVec(num_bits,
                      arr=StateVec.get_ground_st_vec(num_bits).arr)
        print('gs=\n', gs.arr)
        print("gs_trad=\n", gs.get_traditional_st_vec())

        S0100_ZL = StateVec(4,
            arr=StateVec.get_standard_basis_st_vec([0, 1, 0, 0], ZL=True).arr)
        print("S0100_ZL=\n", S0100_ZL.get_traditional_st_vec())

        S0100_ZF = StateVec(4,
            arr=StateVec.get_standard_basis_st_vec([0, 1, 0, 0], ZL=False).arr)
        print("S0100_ZF=\n", S0100_ZF.get_traditional_st_vec())

        st_vec0 = StateVec(num_bits,
            arr=StateVec.get_random_st_vec(num_bits).arr)
        st_vec1 = StateVec(num_bits,
            arr=StateVec.get_random_st_vec(num_bits).arr)
        st_vec_dict = {'br0': st_vec0,
                       'br1': st_vec1,
                       'br3': None}

        trad_st_vec = st_vec0.get_traditional_st_vec()

        den_mat = StateVec.get_den_mat(num_bits, st_vec_dict)
        print("den_mat\n", den_mat)
        print('trace_02 den_mat\n',
              StateVec.get_partial_tr(num_bits, den_mat, {0, 2}))
        print("impurity=", StateVec.get_impurity(den_mat))
        print("entropy=", StateVec.get_entropy(den_mat))
        den_mat_pd = StateVec.get_den_mat_pd(den_mat)
        print('den_mat_pd=', den_mat_pd)

        st_vec_pd = st_vec0.get_pd()
        bit_probs_vec = StateVec.get_bit_probs(num_bits, st_vec_pd)
        bit_probs_dm = StateVec.get_bit_probs(num_bits, den_mat_pd)

        # print("counts_dm=\n", StateVec.get_bit_counts(bit_probs_dm, 10))

        obs_vec = StateVec.get_observations_vec(
            num_bits, st_vec_pd, num_shots=20)
        print('observations vec\n', obs_vec)
        counts = StateVec.get_counts_from_obs_vec(num_bits, obs_vec)
        print('counts from obs_vec\n', counts)

        StateVec.describe_st_vec_dict(st_vec_dict,
                print_st_vec=True, do_pp=True,
                omit_zero_amps=False, show_pp_probs=True, ZL=True)
    main()
    import doctest
    doctest.testmod(verbose=True)
