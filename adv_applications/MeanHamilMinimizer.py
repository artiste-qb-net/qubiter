from adv_applications.CostMinimizer import *
from StateVec import *
import utilities_gen as utg
from CodaSEO_writer import *
from SEO_simulator import *
from CktEmbedder import *
from openfermion.ops import QubitOperator
import sys


class MeanHamilMinimizer(CostMinimizer):
    """
    This class is a child of class CostMinimizer. Like its parent,
    this class is also intended to be abstract and to be subclassed. For
    example, classes MeanHamilMinimizer_native and
    MeanHamilMinimizer_rigetti are children of this class. Whereas its
    parent embodies the essence of any cost minimizing object, this class
    embodies the essence of an object specifically designed to minimize a
    cost function which equals the mean value of a Hamiltonian. In Qubiter's
    docstrings, we refer to this as the Mean Hamiltonian Minimization Problem.

    The qc history of this problem started with quantum chemists planning to
    use on a qc the phase estimation algo invented by Kitaev? (an algo that
    is also implemented in Qubiter) to estimate the energy levels (
    eigenvalues) of simple molecules, initially H2. Then a bunch of people
    realized, heck, rather than trying to estimate the eigenvalues of a
    Hamiltonian by estimating the phase changes it causes, we can estimate
    those eigenvalues more efficiently by estimating the mean value of that
    Hamiltonian as measured empirically on a qc. Basically, just the
    Rayleigh-Ritz method, one of the oldest tricks in the book. One of the
    first papers to propose this mean idea is
    https://arxiv.org/abs/1304.3061 Their algo is commonly referred to by
    the ungainly name VQE (Variational Quantum Eigensolver) VQE was
    originally applied to do quantum chemistry with a qc. But now Rigetti
    and others have renamed it hybrid quantum-classical quantum computing
    and pointed out that it's an algo that has wide applicability, not just
    to quantum chemistry.

    The idea behind hybrid quantum-classical is very simple. One has a
    classical box CBox and a quantum box QBox. The gates of QBox depend on N
    gate parameters. QBox sends info to CBox. CBox sends back to QBox N new
    gate parameters that will lower some cost function. This feedback
    process between CBox and QBox continues until the cost is minimized. The
    cost function is the mean value of a Hamiltonian which is estimated
    empirically from data obtained from the qc which resides inside the QBox.

    To minimize a function of N continuous parameters, one can use some
    methods like simulated annealing and Powell that do not require
    calculating derivatives, or one can use methods that do use derivatives.
    Another possible separation is between methods that don't care which
    local minimum they find, as long as they find one of them, and those
    methods that try to find the best local minimum of them all, the so
    called global minimum. Yet another separation is between methods that
    allow constraints and those that don't.

    Among the methods that do use derivatives, the so called gradient based
    methods only use the 1st derivative, whereas other methods use both
    first (Jacobian) and second (Hessian) derivatives. The performance of
    those that use both 1st and 2nd derivatives degrades quickly as N grows.
    Besides, calculating 2nd derivatives is very expensive. Hence, methods
    that use the 2nd derivatives are practically useless in the neural
    network field where N is usually very large. In that field, gradient
    based methods rule.

    A method that uses no derivatives is Powell. A gradient based method
    that is designed to have a fast convergence rate is the Conjugate
    Gradient (CG) method. Another gradient based method is back-propagation
    (BP). BP can be implemented as distributed computing much more easily
    than other gradient based methods so it is favored by the most popular
    computer programs for doing distributed AI, such as PyTorch and
    Tensorflow.

    The children of this class can perform minimization via various
    interfaces ('scipy', 'autograd', 'pytorch', 'tflow').

    Non-scipy interfaces implement backprop.

    The 'scipy' interface is a wrapper for the scipy function
    `scipy.optimize.minimize`. This scipy umbrella method implements many
    minimization methods, including Powell and CG.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    file_prefix identifies the location of an English file that specifies a
    quantum circuit. If init_st_vec=None, we assume that the initial state
    of that quantum circuit is the ground state (all qubits in state |0>).
    Let |psi> be the final state vector that evolves from that circuit. Let
    hamil be a Hamiltonian suitable for that circuit and stored as an object
    of QubitOperator (a class of the open-source lib OpenFermion). Then the
    cost function to be minimized is <psi|hamil|psi>.

    Attributes
    ----------
    all_var_nums : list[int]
        this is a list of distinct ints that identify each continuous
        variable (i.e., parameter, non-functional placeholder variable) on
        which the cost function depends. They are ordered in order of
        occurrence in the quantum circuit.
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
        number of samples (aka num_shots). If this is zero, the |psi> in
        <psi|H|psi> is calculated exactly from theory. If this is >0,
        the |psi> is calculated empirically from a number num_samples of
        "one-shot" experiments.
    rand_seed : int
        random seed
      
    """

    def __init__(self, file_prefix, num_bits, hamil,
            init_var_num_to_rads, fun_name_to_fun,
            init_st_vec=None, num_samples=0, rand_seed=None,
            print_hiatus=1, verbose=False):
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
        init_x_val = np.array(init_x_val)
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

    def emp_hamil_mean_val(self, var_num_to_rads):
        """
        This abstract method returns the empirically determined Hamiltonian
        mean value. Takes as input the values of placeholder variables.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        float

        """
        assert False

    def pred_hamil_mean_val(self, var_num_to_rads, num_fake_samples=0):
        """
        This method predicts the mean value of the Hamiltonian hamil using
        only Qubiter simulators and not data

        Parameters
        ----------
        var_num_to_rads : dict[int, float]
        num_fake_samples : int

        Returns
        -------
        float

        """
        # give it name unlikely to exist already
        fin_file_prefix = self.file_prefix + '99345125047'

        # hamil loop
        arr_plus = np.array([1., 1.])
        arr_minus = np.array([1., -1.])
        mean_val = 0
        for term, coef in self.hamil.terms.items():
            # we have checked before that coef is real
            coef = complex(coef).real

            # add measurement coda for this term of hamil
            # build real_arr from arr_list.
            # real_arr will be used at end of loop
            arr_list = [arr_plus]*self.num_bits
            bit_pos_to_xy_str = {}
            for bit_pos, action in term:
                arr_list[bit_pos] = arr_minus
                if action != 'Z':
                    bit_pos_to_xy_str[bit_pos] = action
            real_arr = utg.kron_prod(arr_list)
            real_arr = np.reshape(real_arr, tuple([2]*self.num_bits))
            wr = CodaSEO_writer(self.file_prefix,
                                fin_file_prefix, self.num_bits)
            wr.write_xy_measurements(bit_pos_to_xy_str)
            wr.close_files()

            # run simulation. get fin state vec
            vman = PlaceholderManager(
                var_num_to_rads=var_num_to_rads,
                fun_name_to_fun=self.fun_name_to_fun)
            # simulator will change init_st_vec so use
            # fresh copy of it each time
            init_st_vec = cp.deepcopy(self.init_st_vec)
            sim = SEO_simulator(fin_file_prefix, self.num_bits,
                                init_st_vec, vars_manager=vman)
            fin_st_vec = sim.cur_st_vec_dict['pure']
            # print('inside pred hamil in/out stvec',
            # self.init_st_vec, fin_st_vec)

            # get effective state vec
            if num_fake_samples:
                # if num_fake_samples !=0, then
                # sample qubiter-generated empirical prob dist
                pd = fin_st_vec.get_pd()
                obs_vec = StateVec.get_observations_vec(self.num_bits,
                        pd, num_fake_samples, rand_seed=self.rand_seed)
                counts_dict = StateVec.get_counts_from_obs_vec(self.num_bits,
                                                               obs_vec)
                emp_pd = StateVec.get_empirical_pd_from_counts(self.num_bits,
                                                               counts_dict)
                # print('mmmmmmmm,,,', np.linalg.norm(pd-emp_pd))
                emp_st_vec = StateVec.get_emp_state_vec_from_emp_pd(
                        self.num_bits, emp_pd)
                effective_st_vec = emp_st_vec
            else:  # num_fake_samples = 0
                effective_st_vec = fin_st_vec

            # add contribution to mean
            mean_val += coef*effective_st_vec.\
                    get_mean_value_of_real_diag_mat(real_arr)

        # create this coda writer in order to delete final files
        wr1 = CodaSEO_writer(self.file_prefix, fin_file_prefix, self.num_bits)
        wr1.delete_fin_files()

        return mean_val

    def cost_fun(self, x_val):
        """
        This method wraps the abstract method emp_hamil_mean_val() defined
        in a child class. This method will also print out, whenever it is
        called, a report of the current values of x and cost (and pred_cost
        if it is available).

        Parameters
        ----------
        x_val : np.ndarray

        Returns
        -------
        float

        """

        var_num_to_rads = dict(zip(self.all_var_nums, tuple(x_val)))
        cost = self.emp_hamil_mean_val(var_num_to_rads)

        self.cur_x_val = x_val
        self.cur_cost = cost
        self.broadcast_cost_fun_call()
        self.iter_count += 1

        return cost

    def pred_cost_fun(self, x_val):
        """
        Returns the cost, predicted from theory, rather than estimated from
        data as in cost_fun(). This method mimics the method cost_fun(),
        but that one wraps the abstract method emp_hamil_mean_val(). This
        one wraps the non-abstract method pred_hamil_mean_val() which is
        defined right in this class.

        Parameters
        ----------
        x_val : np.ndarray

        Returns
        -------
        float

        """
        # print('inside_pred_cost', x_val)
        var_num_to_rads = dict(zip(self.all_var_nums, tuple(x_val)))
        # print('mmmmmmmmmaaaaaa', var_num_to_rads)
        cost = self.pred_hamil_mean_val(var_num_to_rads)
        # print('bbvvv-cost', cost)
        return cost

    def find_min(self, interface, **kwargs):
        """
        This method finds minimum of cost function. It allows user to choose
        among several possible interfaces, namely, 'scipy', 'autograd',
        'pytorch', 'tflow'. Interface parameters can be passed in via kwargs.

        Non-scipy interfaces do backprop. It is known that the complexity of
        calculating forward  propagation and back propagation are about the
        same.

        kwargs (keyword arguments)
        interface = scipy
            the keyword args of scipy.optimize.minimize
        interface = autograd
            num_inter : float
                number of iterations (an iteration is every time call
                cost_fun)
            descent_rate : float
                positive float, constant that multiplies gradient of
                pred_cost_fun(). Often denoted as eta
            do_pred_cost : bool
                default=True, The gradient of the pred_cost_fun() is always
                calculated, but calculating the pred_cost_fun() itself is
                optional.

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
        print('x_val~ (' +
              ', '.join(['#' + str(k) for k in self.all_var_nums]) + ')')

        def pred_cost(xx):
            # self argument seems to confuse grad
            return self.pred_cost_fun(xx)
        if interface == 'scipy':
            import scipy
            minimizer_fun = scipy.optimize.minimize
            opt_result = minimizer_fun(self.cost_fun,
                self.init_x_val, **kwargs)
            if self.verbose:
                print('*********final optimum result'
                      ' (final step=' + str(self.iter_count) +
                      '):\n', opt_result)
            return opt_result
        elif interface == 'autograd':
            from autograd import grad

            assert 'num_iter' in kwargs, \
                "must pass-in keyword 'num_iter=' " \
                "if using autograd interface"
            num_iter = kwargs['num_iter']
            assert 'descent_rate' in kwargs, \
                "must pass-in keyword 'descent_rate=' " \
                "if using autograd interface"
            rate = kwargs['descent_rate']
            self.cur_x_val = self.init_x_val
            if 'do_pred_cost' in kwargs:
                do_pred_cost = kwargs['do_pred_cost']
            else:
                do_pred_cost = True
            for step in range(num_iter):
                xlist = list(self.cur_x_val)
                # print('mmbbb', self.cur_x_val, xlist)
                if do_pred_cost:
                    self.cur_pred_cost = pred_cost(self.cur_x_val)
                self.cur_cost = self.cost_fun(self.cur_x_val)
                # print('kkkhhh', grad(pred_cost)(self.cur_x_val))
                for dwrt in range(len(xlist)):
                    self.cur_x_val[dwrt] -= \
                        rate*grad(pred_cost)(self.cur_x_val)[dwrt]

        elif interface == 'pytorch':
            assert False, 'not yet'
        elif interface == 'tflow':
            assert False, 'not yet'
        else:
            assert False, 'unsupported find_min() interface'


if __name__ == "__main__":
    def main():
        print(5)
    main()
