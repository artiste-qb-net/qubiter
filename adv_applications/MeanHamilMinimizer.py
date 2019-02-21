from adv_applications.CostMinimizer import *


class MeanHamilMinimizer(CostMinimizer):
    """
    This class is a child of class CostMinimizer. Like its parent,
    this class is also intended to be abstract and to be subclassed. For
    example, class MeanHamilMinimizer_naive is a child of this class.
    Whereas its parent embodies the essence of any cost minimizing object,
    this class embodies the essence of an object specifically designed to
    minimize a cost function which equals the mean value of a Hamiltonian.
     
    file_prefix identifies the location of an English file that specifies a 
    quantum circuit. Assuming that the initial state of that quantum circuit 
    is the ground state (all qubits in state |0>), let |psi> be the final 
    state vector that evolves from that circuit. Let hamil be a Hamiltonian 
    suitable for that circuit and stored as an object of QubitOperator, 
    which is a class of the open-source lib OpenFermion. Then the cost 
    function to be minimized is <psi|hamil|psi>. 
    
    Attributes
    ----------
    all_var_nums : list[int]
        this is a list of distinct ints that identify each continuous 
        variable (i.e., parameter, non-functional placeholder variable) on
        which the cost function depends.
    file_prefix : str
    fun_name_to_fun : dict[str, function]
         this is a dict that maps function names to functions. Such
         functions are functional placeholders, meaning that their values
         are only decided at a later time. These functions do not vary
         during the minimization process.
    hamil : QubitOperator
    init_var_num_to_rads : dict[int, float]
        this dictionary gives the initial values for the cost function being
        minimized. The dict maps variable numbers (int) to radians (float)
    num_bits : int
        number of qubits
    num_samples : int
        number of samples. If this is zero, the |psi> in <psi|H|psi> is
        calculated exactly. If this is >0, the |psi> is calculated
        empirically from a number num_samples of "one-shot" experiments.
    rand_seed : int
        random seed
      
    """

    def __init__(self, file_prefix, num_bits, hamil,
            init_var_num_to_rads, fun_name_to_fun,
            minimizer_fun, num_samples=0, rand_seed=None,
            print_hiatus=0, verbose=False, **mfun_kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        hamil : QubitOperator
        init_var_num_to_rads : dict[int, float]
        fun_name_to_fun : dict[str, function]
        minimizer_fun : function
        num_samples : int
        rand_seed : int
        print_hiatus : int
        verbose : bool
        mfun_kwargs : dict

        Returns
        -------

        """
        self.file_prefix = file_prefix
        self.num_bits = num_bits
        self.hamil = hamil
        MeanHamilMinimizer.check_hamil_is_herm(hamil)
        self.init_var_num_to_rads = init_var_num_to_rads
        self.all_var_nums, init_x_val = zip(*init_var_num_to_rads.items())
        self.fun_name_to_fun = fun_name_to_fun
        self.num_samples = num_samples
        self.rand_seed = rand_seed

        CostMinimizer.__init__(self, minimizer_fun, init_x_val,
                print_hiatus, verbose, **mfun_kwargs)

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
