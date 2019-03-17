from adv_applications.CostMinimizer import *
from StateVec import *


class MeanHamilMinimizer(CostMinimizer):
    """
    This class is a child of class CostMinimizer. Like its parent,
    this class is also intended to be abstract and to be subclassed. For
    example, class MeanHamilMinimizer_naive is a child of this class.
    Whereas its parent embodies the essence of any cost minimizing object,
    this class embodies the essence of an object specifically designed to
    minimize a cost function which equals the mean value of a Hamiltonian.
     
    file_prefix identifies the location of an English file that specifies a
    quantum circuit. If init_st_vec=None, we assume that the initial state
    of that quantum circuit is the ground state (all qubits in state |0>).
    Let |psi> be the final state vector that evolves from that circuit. Let
    hamil be a Hamiltonian suitable for that circuit and stored as an object
    of QubitOperator, which is a class of the open-source lib OpenFermion.
    Then the cost function to be minimized is <psi|hamil|psi>.
    
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
    init_st_vec : StateVec
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
            init_st_vec=None, num_samples=0, rand_seed=None,
            print_hiatus=0, verbose=False):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        hamil : QubitOperator
        init_var_num_to_rads : dict[int, float]
        fun_name_to_fun : dict[str, function]
        init_st_vec : StateVec
        num_samples : int
        rand_seed : int
        print_hiatus : int
        verbose : bool

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
        self.init_st_vec = init_st_vec
        if self.init_st_vec is None:
            self.init_st_vec = StateVec.get_ground_st_vec(self.num_bits)
        self.num_samples = num_samples
        self.rand_seed = rand_seed

        CostMinimizer.__init__(self, init_x_val, print_hiatus, verbose)

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

    def hamil_mean_val(self, var_num_to_rads):
        """
        This abstract method calculates the mean value of the Hamiltonian
        hamil.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        float

        """
        assert False

    def cost_fun(self, x_val):
        """
        This method wraps the static method hamil_mean_val() defined
        elsewhere in this class. This method will also print out whenever it
        is called, a report of the current values of x and cost.

        Parameters
        ----------
        x_val : tuple[float]

        Returns
        -------
        float

        """
        var_num_to_rads = dict(zip(self.all_var_nums, x_val))
        cost = self.hamil_mean_val(var_num_to_rads)

        self.cur_x_val = x_val
        self.cur_cost = cost
        self.broadcast_cost_fun_call()
        self.iter_count += 1

        return cost

    def find_min(self, interface, **kwargs):
        """
        This method finds min of cost function. It allows user to choose
        among several possible interfaces, namely, 'scipy', 'autograd',
        'pytorch', 'tflow'. Interface parameters can be passed in via kwargs.

        Parameters
        ----------
        interface : str
        kwargs : dict

        Returns
        -------
        OptimizeResult | None
            OptimizeResult is a class (basically an enum) defined in
            scipy.optimize to hold the output results of
            scipy.optimize.minimize

        """
        print('x_val~ (' +\
              ', '.join(['#' + str(k) for k in self.all_var_nums]) + ')')
        if interface == 'scipy':
            import scipy
            minimizer_fun = scipy.optimize.minimize
            opt_result = minimizer_fun(self.cost_fun,
                self.init_x_val, **kwargs)
            if self.verbose:
                print('*********final optimum result'
                      ' (final iter=' + str(self.iter_count) +\
                      '):\n', opt_result)
            return opt_result
        elif interface == 'autograd':
            assert False, 'not yet'
        elif interface == 'tflow':
            assert False, 'not yet'
        elif interface == 'pytorch':
            assert False, 'not yet'
        else:
            assert False, 'unsupported fin_min() interface'
