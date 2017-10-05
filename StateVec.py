import numpy as np
import pprint as pp


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
        mat_phi = 2*np.pi*np.random.rand(1 << num_bits)
        mat_r = np.random.rand(1 << num_bits)
        arr = mat_r*(np.cos(mat_phi) + 1j*np.sin(mat_phi))
        magnitude = np.linalg.norm(arr)
        arr /= magnitude
        arr = arr.reshape([2]*num_bits)
        return StateVec(num_bits, arr)

    @staticmethod
    def get_standard_basis_st_vec(spin_dir_list, ZL=False):
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
        Internally, self.arr in Qubiter has shape [2]*num_bits and assumes
        ZF convention for axes indexes. However, the traditional way of
        writing a state vector is as a column array of dimension 1<<
        num_bits in the ZL convention. This function returns the traditional
        view. So it reshapes (flattens) the array and it reverses the axes (
        reversing axes takes it from ZF to ZL).

        The rows are always labelled 0, 1, 2, 3, ... or the binary
        representation thereof, regardless of whether ZL or ZF convention.
        One can go from digital to binary labels and vice versa using:

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

    def get_pd(self):
        """
        Returns copy of self.get_traditional_st_vec() with amplitudes
        replaced by probabilities. pd = probability distribution. So returns
        one column array indexed in ZL convention like the traditional state
        vec is.

        Parameters
        ----------

        Returns
        -------
        np.ndarray

        """
        x = self.get_traditional_st_vec()
        return np.real(x*np.conj(x))

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
    def get_partial_tr(num_bits, den_mat, tr_bits_set):
        """
        Returns the partial trace of a density matrix den_mat. Traces over
        qubits in set tr_bits_set. To get full trace, just do np.trace(
        den_mat)

        Parameters
        ----------
        num_bits : int
        den_mat : np.ndarray
            if dim=2^num_bits, this function assumes that den_mat has shape
            (dim, dim) and that it's indexed in the ZL convention so qubit 0
            corresponds to axis num_bits-1.
        tr_bits_set : set[int]
             Set of qubits being traced over

        Returns
        -------
        np.ndarray

        """
        dim = 1 << num_bits
        assert den_mat.shape == (dim, dim)
        assert set(range(num_bits)) > tr_bits_set
        dm = den_mat.reshape([2]*(2*num_bits))
        # largest axis must be traced first
        # bit 0 corresponds to axis num_bits - 1
        traced_axes = \
            sorted([num_bits - 1 - k for k in tr_bits_set], reverse=True)
        for k, ax in enumerate(traced_axes):
            dm = np.trace(dm, axis1=ax, axis2=ax + num_bits - k)
        return dm

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
    def get_entropy(den_mat):
        """
        Returns entropy of density matrix den_mat. Uses natural log for
        entropy.

        Parameters
        ----------
        den_mat : np.ndarray
            Density matrix. Real part of eigenvalues must be positive and
            sum to 1

        Returns
        -------
        float

        """
        evas = np.real(np.linalg.eigvals(den_mat))
        assert np.all(evas > 0)
        assert abs(np.sum(evas) - 1) < 1e-6
        if np.any(evas > 1 - 1e-6):
            ent = 0.0
        else:
            ent = - np.sum(evas * np.log(evas))
        return ent

    def get_total_prob(self):
        """
        Returns total probability of self.

        Parameters
        ----------

        Returns
        -------
        float

        """
        return np.sum(np.real(self.arr*self.arr.conj()))

    @staticmethod
    def get_bit_probs(num_bits, pd, normalize=False):
        """
        Returns a list whose jth item is, for the jth qubit, the pair (p,
        p_tot-p)/norm, where p is the probability that the jth qubit is 0,
        if the state of all other qubits is ignored. p_tot = np.sum( pd).
        norm=1 if normalize=False and norm=p_tot otherwise.

        Parameters
        ----------
        num_bits : int
        pd : np.ndarray
            probability distribution of shape (2^num_bits,) IMP: assumed to
            be indexed in ZL convention

        normalize : bool

        Returns
        -------
        list[tuple[float, float]]

        """
        assert 1 << num_bits == pd.shape[0]
        probs = []
        arr = pd.reshape([2] * num_bits)
        # p_total may not be one
        # if a measurement has been done
        p_total = np.sum(arr)
        # slicex is a portmanteau of slice index
        # print("state_vec_pd=\n", vec)
        slicex = [slice(None)]*num_bits
        # pd is assumed to be in ZL convention
        # so when reshape, bit k has axis= num_bits -1 - k
        for k in range(num_bits):
            k_axis = num_bits - 1 - k
            slicex[k_axis] = 0
            p = np.sum(arr[tuple(slicex)])
            norm = 1.0
            if normalize:
                norm = p_total
            probs.append((p/norm, (p_total-p)/norm))
            slicex[k_axis] = slice(None)  # restore to all entries slice(None)
        # print(probs)
        return probs

    @staticmethod
    def sample_bit_probs(bit_probs, num_samples, rand_seed=None):
        """
        Returns a list whose jth item is, for the jth qubit, the pair (
        count0, count1), where count0 + count1 = num_samples, and as
        num_samples ->infinity, (count0, count1)/num_samples tends to the
        probability pair bit_probs[j].

        Parameters
        ----------
        bit_probs : list[tuple[float, float]]
        num_samples : int
        rand_seed : int

        Returns
        -------
        list[tuple[int, int]]

        """
        if rand_seed:
            np.random.seed(rand_seed)
        num_bits = len(bit_probs)
        counts = []
        for bit in range(num_bits):
            x = np.random.choice(
                np.arange(0, 2), size=num_samples,
                p=bit_probs[bit])
            sum1 = np.sum(x)
            sum0 = num_samples - sum1
            counts.append((sum0, sum1))
        return counts

    def pp_arr_entries(self, omit_zero_amps=False,
                          show_probs=False, ZL=True):
        """
        pp=pretty print. Prints for each entry of self.arr, a line of the
        form (i, j, k, ...) self.arr[i, j, k, ...], with zero bit last (
        resp., first) if ZL=True (resp., False).

        Parameters
        ----------
        omit_zero_amps : bool
            If True, will not list states with zero amplitude
        show_probs : bool
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
            if show_probs:
                if omit_zero_amps:
                    if mag > 1E-6:
                        print(ind_str, x, ', prob=', mag**2)
                else:
                    print(ind_str, x, ', prob=', mag**2)
            else:
                if omit_zero_amps:
                    if mag > 1E-6:
                        print(ind_str, x)
                else:
                    print(ind_str, x)

    def describe_self(self, print_st_vec=False, do_pp=False,
            omit_zero_amps=False, show_probs=False, ZL=True):
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

        show_probs : bool
            If True, will show probability of each standard basis state

        ZL : bool
            If True, multi-index of ket in ZL (Zero bit Last) convention.
            If False, multi-index of ket in ZF (Zero bit First) convention.

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
                self.pp_arr_entries(omit_zero_amps, show_probs, ZL)
            else:
                print(self.arr)

        print('total probability of state vector ' +
              '(=one if no measurements)=', self.get_total_prob())

        print('dictionary with key=qubit, value=(Prob(0), Prob(1))')
        bit_probs = StateVec.get_bit_probs(self.num_bits, self.get_pd())
        pp.pprint(dict(enumerate(bit_probs)))

    @staticmethod
    def describe_st_vec_dict(st_vec_dict, **kwargs):
        """
        Calls describe_st_vec() for each branch of st_vec_dict

        Parameters
        ----------
        st_vec_dict : dict[str, StateVec]
        kwargs : dict
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
    # print("den_mat", den_mat)
    st_vec_pd = st_vec0.get_pd()
    den_mat_pd = StateVec.get_den_mat_pd(den_mat)
    bit_probs_vec = st_vec0.get_bit_probs(num_bits, st_vec_pd)
    bit_probs_dm = StateVec.get_bit_probs(num_bits, den_mat_pd)

    print("counts_dm=\n", StateVec.sample_bit_probs(bit_probs_dm, 10))
    print("impurity=", StateVec.get_impurity(den_mat))
    StateVec.describe_st_vec_dict(st_vec_dict,
            print_st_vec=True, do_pp=True,
            omit_zero_amps=False, show_probs=True, ZL=True)
