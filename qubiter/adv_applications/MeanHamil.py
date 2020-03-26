from qubiter.StateVec import *
import qubiter.utilities_gen as utg
from qubiter.CodaSEO_writer import *
from qubiter.SEO_simulator import *
from qubiter.CktEmbedder import *
from openfermion.ops import QubitOperator


class MeanHamil:
    """
    This is an abstract class. The main purpose of this class is to evaluate
    the mean value of a Hamiltonian.

    The Hamiltonian hamil is stored as an object of QubitOperator (a class
    of the open-source lib OpenFermion). `terms` is an attribute of
    QubitOperator. hamil.terms is a dictionary that maps a term to a
    coefficient. A term represents a tensor product of Pauli matrices (a
    Pauli string) as a tuple of tuples of the form (bit_pos, action). An
    example of a term: ( (1, 'X'), (2, 'Y'))

    file_prefix identifies the location of an English file that specifies a
    quantum circuit. If init_st_vec=None, we assume that the initial state
    of that quantum circuit is the ground state (all qubits in state |0>).
    Let |psi> be the final state vector that evolves from that circuit. Let
    hamil be a Hamiltonian suitable for that circuit and stored as an object
    of QubitOperator (a class of the open-source lib OpenFermion). Then the
    mean value evaluated by this class is <psi|hamil|psi>.

    Subclasses of this class use different methods to evaluate this mean
    value. They might change the tensor lib (numpy, PyTorch, TensorFlow) or
    the device (native, Rigetti, etc.) or the simulator for a particular
    device. They might evaluate the mean value exactly or empirically.

    Attributes
    ----------
    all_var_nums : list[int]
        This is a list of all the non-functional placeholder variable numbers
    file_prefix : str
        Prefix to English file to be used in evaluating the mean hamil
    fun_name_to_fun : dict[str, function]
         This is a dict that maps function names to functions. Such
         functions are functional placeholders, meaning that their values
         are only decided at a later time. These functions do not vary
         during the minimization process.
    hamil : QubitOperator
        Hamiltonian
    init_st_vec : StateVec
        initial state vector
    num_qbits : int
        number of qubits
    num_samples : int
        number of samples (aka num_shots). If this is zero, the `|psi>` in
        `<psi|H|psi>` is calculated exactly from theory. If this is >0,
        the `|psi>` is calculated empirically from a number num_samples of
        "one-shot" experiments.
    simulator_name : str | None
        name of the simulator.

    """
    def __init__(self, file_prefix, num_qbits, hamil,
            all_var_nums, fun_name_to_fun, init_st_vec=None,
            simulator_name=None, num_samples=0):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_qbits : int
        hamil : QubitOperator
        all_var_nums : list[int]
        fun_name_to_fun : dict[str, function]
        init_st_vec : StateVec
        simulator_name : str
        num_samples : int

        Returns
        -------

        """
        self.file_prefix = file_prefix
        self.num_qbits = num_qbits
        self.hamil = hamil
        MeanHamil.check_hamil_is_herm(hamil)
        MeanHamil.check_hamil_is_in_range(hamil, num_qbits-1)
        self.all_var_nums = all_var_nums
        self.fun_name_to_fun = fun_name_to_fun
        self.init_st_vec = init_st_vec
        if self.init_st_vec is None:
            self.init_st_vec = StateVec.get_ground_st_vec(self.num_qbits)
        self.simulator_name = simulator_name
        self.num_samples = num_samples

    @staticmethod
    def check_hamil_is_herm(hamil):
        """
        Checks that the Hamiltonian hamil is a Hermitian operator. Emits
        warning and stops execution if it isn't.

        Parameters
        ----------
        hamil : QubitOperator

        Returns
        -------
        None

        """
        for term, coef in hamil.terms.items():
            coef = complex(coef)
            if abs(coef.imag) > 1e-8:
                assert False, 'The Hamiltonian should be Hermitian but it ' +\
                    "isn't. After being simplified by the " +\
                    'BosonOperator constructor, ' +\
                    'the coefficient of every term must be real.'

    @staticmethod
    def check_hamil_is_in_range(hamil, max_bit_pos):
        """
        Checks that the Hamiltonian hamil operates on range(max_bit_pos+1).

        Parameters
        ----------
        hamil : QubitOperator
        max_bit_pos : int

        Returns
        -------
        None

        """
        for term, coef in hamil.terms.items():
            for bit, action in term:
                assert bit in range(max_bit_pos+1)

    def get_real_vec(self, term):
        """
        Internal method that returns a numpy array, of shape [2]*num_qbits,
        that will be used as input to the method
        StateVec.get_mean_value_of_real_diag_mat()

        The input is a `term`. `terms` is an attribute of QubitOperator (a
        class in OpenFermion). terms is a dictionary that maps a term to a
        coefficient. A term represents a tensor product of Pauli matrices (a
        Pauli string) as a tuple of tuples of the form (bit_pos, action). An
        example of a term: ((1, 'X'), (2, 'Y'))

        Parameters
        ----------
        term : tuple

        Returns
        -------
        np.ndarray
            shape=[2]*num_qbits

        """
        arr_plus = np.array([1., 1.])
        arr_minus = np.array([1., -1.])
        arr_list = [arr_plus]*self.num_qbits
        for bit_pos, action in term:
            arr_list[bit_pos] = arr_minus
        real_arr = utg.kron_prod(arr_list)
        real_arr = np.reshape(real_arr, tuple([2]*self.num_qbits))
        return real_arr

    def get_mean_val(self, var_num_to_rads):
        """
        Abstract method. The main goal of subclasses of this class is to
        override this method.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------

        """
        assert False
