import numpy as np
import pprint as pp


class StateVec:
    """
    This class has no constructor. All its methods are static methods. It
    provides functions for constructing state vectors which are complex
    nummpy arrays with shape [2]*num_bits. The class also provides functions
    for performing calculations using as input: a single state vector,
    dictionaries of state vectors, and density matrices.

    The keys of these dictionaries of state vectors are strings that we call
    branch_keys. The class provides a function for constructing from such
    dictionaries of state vectors, a density matrix which is a 2 dim square
    numpy array of dimension 2^num_bits.

    """
    @staticmethod
    def get_ground_st_vec(num_bits):
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
    def get_random_st_vec(num_bits, rand_seed=None):
        """
        Returns random state \sum_b^n A(b^n)|b^n>, b^n \in {0,1}^n,
        where n=num_bits and \sum_b^n |A(b^n)|^2 = 1

        Parameters
        ----------
        num_bits : int
        rand_seed : int

        Returns
        -------
        np.ndarray

        """
        if rand_seed:
            np.random.seed(rand_seed)
        # returns array of random numbers in [0, 1] interval
        mat_phi = 2*np.pi*np.random.rand(1 << num_bits)
        mat_r = np.random.rand(1 << num_bits)
        mat = mat_r*(np.cos(mat_phi) + 1j*np.sin(mat_phi))
        magnitude = np.linalg.norm(mat)
        mat /= magnitude
        mat = mat.reshape([2]*num_bits)
        return mat

    @staticmethod
    def get_standard_basis_st_vec(spin_dir_list, ZL=False):
        """
        If ZL = True, Returns state ...|s2>|s1>|s0>, where
        spin_dir_list=[...,s2, s1, s0], s_j \in {0, 1} for all j, |0> = [1,
        0]^t and |1> = [0,1]^t, t = transpose

        Parameters
        ----------
        spin_dir_list : list[int]
        ZL : bool
            True(False) if last(first) qubit is at position 0

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        num_bits = len(spin_dir_list)
        mat = np.zeros([1 << num_bits], dtype=ty)
        mat = mat.reshape([2]*num_bits)
        if ZL:
            spin_dir_list = reversed(spin_dir_list)
        mat[tuple(spin_dir_list)] = 1
        return mat

    @staticmethod
    def get_traditional_st_vec(st_vec):
        """
        Internally, arrays in Qubiter assume ZF convention and state vectors
        have shape = [2]*num_bits. However, the traditional way of writing a
        state vector is as a column array of dimension 1<< num_bits in the
        ZL convention. This function returns the traditional view. So it
        reshapes (flattens) the array and it reverses the
        axes (reversing axes takes it from ZF to ZL).

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
        st_vec : np.array

        Returns
        -------
        np.array

        """
        num_bits = st_vec.ndim
        assert st_vec.shape == tuple([2]*num_bits)
        perm = list(reversed(range(num_bits)))
        return np.transpose(st_vec, perm).flatten()

    @staticmethod
    def get_den_mat(num_bits, st_vec_dict):
        """
        Returns a density matrix (in ZL convention) constructed from
        st_vec_dict which is a dict from strings to numpy arrays with shape
        [2]*num_bits.

        The rows and columns are always labelled 0, 1, 2, .. or binary
        representation thereof, regardless of whether ZL or ZF convention.
        For more info, see docstring of get_traditional_st_vec().

        Parameters
        ----------
        num_bits : int
        st_vec_dict : dict[str, np.ndarray]

        Returns
        -------
        np.ndarray

        """

        dim = 1 << num_bits
        den_mat = np.zeros((dim, dim), dtype=complex)
        # print(",,,", den_mat)
        for br_key in st_vec_dict:
            if st_vec_dict[br_key] is None:
                continue
            vec = StateVec.get_traditional_st_vec(st_vec_dict[br_key])
            assert vec.shape == (dim,)
            den_mat += np.outer(vec, np.conj(vec))
            # print(br_key, den_mat)
        tr = np.trace(den_mat)
        assert abs(tr) > 1e-6
        return den_mat/tr

    @staticmethod
    def get_st_vec_pd(st_vec):
        """
        Returns copy of st_vec.flatten() with amplitudes replaced by
        probabilities. pd = probability distribution

        Parameters
        ----------
        st_vec : np.ndarray

        Returns
        -------
        np.ndarray

        """
        # print(trad_st_vec)
        x = st_vec.flatten()
        return np.real(x*np.conj(x))

    @staticmethod
    def get_den_mat_pd(den_mat):
        """
        Returns the diagonal of den_mat. den_mat is expected to be a density
        matrix returned by get_den_mat()

        Parameters
        ----------
        den_mat : np.ndarray

        Returns
        -------
        np.ndarray

        """
        return np.real(np.diag(den_mat))

    @staticmethod
    def get_purity(den_mat):
        """
        Returns the 2-norm of get_den_mat^2 - get_den_mat. This
        is zero iff the density matrix represents a pure state

        Parameters
        ----------
        den_mat : np.nparray
            density matrix represented as square array

        Returns
        -------
        float

        """
        return np.linalg.norm(np.dot(den_mat, den_mat) - den_mat)

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
    def get_bit_probs(num_bits, st_vec_pd, normalize=False):
        """
        For a given state vector probability distribution st_vec_pd,
        which is a numpy array of shape (num_bits,), it returns a list whose
        jth item is, for the jth qubit, the pair (p, p_tot-p)/norm, where p
        is the probability that the jth qubit is 0, if the state of all
        other qubits is ignored. p_tot = np.sum(st_vec_pd). norm=1 if
        normalize=False and norm=p_tot otherwise.

        Parameters
        ----------
        num_bits : int
        st_vec_pd : np.ndarray
            state vector probability distribution
        normalize : True

        Returns
        -------
        list[tuple[float, float]]

        """
        assert 1 << num_bits == st_vec_pd.shape[0]
        probs = []
        vec = st_vec_pd.reshape([2]*num_bits)
        # p_total may not be one
        # if a measurement has been done
        p_total = np.sum(vec)
        # slicex is a portmanteau of slice index
        # print("state_vec_pd=\n", vec)
        slicex = [slice(None)]*num_bits
        for k in range(num_bits):
            slicex[k] = 0
            p = np.sum(vec[tuple(slicex)])
            norm = 1.0
            if normalize:
                norm = p_total
            probs.append((p/norm, (p_total-p)/norm))
            slicex[k] = slice(None)  # restore to all entries slice(None)
        # print(probs)
        return probs

    @staticmethod
    def get_bit_probs1(st_vec, normalize=False):
        """
        Returns same as get_bit_probs() but has st_vec instead of st_vec_pd
        as input.

        Parameters
        ----------
        st_vec : np.ndarray
        normalize : bool

        Returns
        -------
        list[tuple[float, float]]


        """
        num_bits = st_vec.shape[0]
        st_vec_pd = StateVec.get_st_vec_pd(st_vec)
        return StateVec.get_bit_probs(num_bits, st_vec_pd)

    @staticmethod
    def sample_bit_probs(bit_probs, num_samples, rand_seed=None):
        """
        Returns a list whose jth item is, for the jth qubit, the pair (
        count0, count1), where count0 + count1 = num_samples, and as
        num_samples ->infinity, count0/num_samples tends to the probability
        bit_probs[j].

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

    @staticmethod
    def pp_st_vec_entries(st_vec, omit_zero_amps=False,
                          show_probs=False, ZL=True):
        """
        pp=pretty print. Prints for each entry of the numpy array 'st_vec',
        a line of the form (i, j, k, ...) st_vec[i, j, k, ...], with zero
        bit last (resp., first) if ZL=True (resp., False).

        Parameters
        ----------
        st_vec : numpy.array
        omit_zero_amps : bool
            If True, will not list states with zero amplitude
        show_probs : bool
            If True, will show probability of each amplitude

        ZL : bool
            If False, multi-index in usual internal order, ZF (Zero bit
            First) convention. If True, multi-index in reverse of usual
            internal order, ZL (Zero bit Last) convention.

        Returns
        -------
        None

        """
        for ind, x in np.ndenumerate(st_vec):
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

    @staticmethod
    def describe_st_vec(
            st_vec, print_st_vec=False, do_pp=False,
            omit_zero_amps=False, show_probs=False, ZL=True):
        """
        Prints a description of the state vector st_vec

        Parameters
        ----------
        st_vec : np.ndarray|None

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
            If False, multi-index of ket in usual internal order, ZF (Zero
            bit First) convention. If True, multi-index of ket in reverse of
            usual internal order, ZL (Zero bit Last) convention.

        Returns
        -------
        None

        """
        if st_vec is None:
            print("zero state vector")
            return
        if print_st_vec:
            print('state vector:')
            if do_pp:
                if not ZL:
                    print('ZF convention (Zero bit First in state tuple)')
                else:
                    print('ZL convention (Zero bit Last in state tuple)')
                StateVec.pp_st_vec_entries(
                    st_vec, omit_zero_amps, show_probs, ZL)
            else:
                print(st_vec)

        print('total probability of state vector ' +
              '(=one if no measurements)=',
              StateVec.get_total_prob(st_vec))

        print('dictionary with key=qubit, value=(Prob(0), Prob(1))')
        st_vec_pd = StateVec.get_st_vec_pd(st_vec)
        bit_probs = StateVec.get_bit_probs(st_vec.ndim, st_vec_pd)
        pp.pprint(dict(enumerate(bit_probs)))

    @staticmethod
    def describe_st_vec_dict(st_vec_dict, **kwargs):
        """
        Calls describe_st_vec() for each branch of st_vec_dict

        Parameters
        ----------
        st_vec_dict : dict[str, np.ndarray]
        kwargs : dict
            the keyword arguments of describe_st_vec()

        Returns
        -------
        None

        """
        for br_key in st_vec_dict.keys():
            print("*********branch= " + br_key)
            StateVec.describe_st_vec(st_vec_dict[br_key], **kwargs)

if __name__ == "__main__":

    num_bits = 3
    gs = StateVec.get_ground_st_vec(num_bits)
    print('gs=\n', gs)
    gs_trad = StateVec.get_traditional_st_vec(gs)
    print("gs_trad=\n", gs_trad)

    S0100_ZL = StateVec.get_standard_basis_st_vec(
        [0, 1, 0, 0], ZL=True)
    print("S0100_ZL=\n", StateVec.get_traditional_st_vec(S0100_ZL))

    S0100_ZF = StateVec.get_standard_basis_st_vec(
        [0, 1, 0, 0], ZL=False)
    print("S0100_ZF=\n", StateVec.get_traditional_st_vec(S0100_ZF))

    st_vec0 = StateVec.get_random_st_vec(num_bits)
    st_vec1 = StateVec.get_random_st_vec(num_bits)
    st_vec_dict = {'br0': st_vec0,
                   'br1': st_vec1,
                   'br3': None}

    trad_st_vec = StateVec.get_traditional_st_vec(st_vec0)
    den_mat = StateVec.get_den_mat(num_bits, st_vec_dict)
    # print("den_mat", den_mat)
    st_vec_pd = StateVec.get_st_vec_pd(st_vec0)
    den_mat_pd = StateVec.get_den_mat_pd(den_mat)
    bit_probs_vec = StateVec.get_bit_probs(num_bits, st_vec_pd)
    bit_probs_dm = StateVec.get_bit_probs(num_bits, den_mat_pd)

    print("counts_dm=\n", StateVec.sample_bit_probs(bit_probs_dm, 10))
    print("purity=", StateVec.get_purity(den_mat))
    StateVec.describe_st_vec_dict(st_vec_dict,
            print_st_vec=True, do_pp=True,
            omit_zero_amps=False, show_probs=True, ZL=True)