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
    This class is a child of class CostMinimizer. It's purpose is to
    minimize a cost function which equals the mean value of a Hamiltonian.
    We refer to this task as the Mean Hamiltonian Minimization Problem.

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
    the ungainly name VQE (Variational Quantum Eigensolver). VQE was
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
    TensorFlow.

    Qubiter can perform minimization using various minlibs (minimization
    software libraries) such as 'scipy', 'autograd', 'pytorch', 'tflow'. It
    can also use various devices (aka simulators or backends), either
    virtual or real, to do the minimization.

    Non-scipy minlibs implement backprop.

    The 'scipy' minlib is a wrapper for the scipy function
    `scipy.optimize.minimize`. This scipy umbrella method implements many
    minimization methods, including Powell and CG.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    By a native device, we mean one that uses Qubiter native simulators like
    SEO_simulator.

    Attributes
    ----------
    all_var_nums : list[int]
        this is a list of distinct ints that identify each continuous
        variable (i.e., parameter, non-functional placeholder variable) on
        which the cost function depends. The ordering corresponds to the
        ordering of self.init_x_val
    emp_mhamil : MeanHamil
        Empirical mean Hamiltonian, used to evaluate cost
    init_x_val : nm.ndarray
        this array gives the initial values in radians for the cost function
        being minimized. The ordering corresponds to the ordering of
        self.all_var_nums
    pred_mhamil : MeanHamil
        Prediction mean Hamiltonian, used to evaluate pred_cost

    """

    def __init__(self, emp_mhamil, pred_mhamil,
                 all_var_nums, init_var_num_to_rads,
                 print_hiatus=1, verbose=False):
        """
        Constructor

        Parameters
        ----------
        emp_mhamil : MeanHamil
        pred_mhamil : MeanHamil
        all_var_nums : list[int]
        init_var_num_to_rads : dict[int, float]
        print_hiatus : int
        verbose : bool

        Returns
        -------

        """
        self.emp_mhamil = emp_mhamil
        self.pred_mhamil = pred_mhamil
        self.all_var_nums = all_var_nums
        assert emp_mhamil.all_var_nums == all_var_nums
        assert pred_mhamil.all_var_nums == all_var_nums
        init_x_val = [init_var_num_to_rads[k] for k in all_var_nums]
        self.init_x_val = np.array(init_x_val)

        CostMinimizer.__init__(self, init_x_val, print_hiatus, verbose)

    def cost_fun(self, x_val):
        """
        This method wraps self.emp_mhamil.get_mean_val(). This method will
        also print out, whenever it is called, a report of the current
        values of x and cost (and pred_cost if it is available).

        Parameters
        ----------
        x_val : np.ndarray

        Returns
        -------
        float

        """

        var_num_to_rads = dict(zip(self.all_var_nums, tuple(x_val)))
        cost = self.emp_mhamil.get_mean_val(var_num_to_rads)

        self.cur_x_val = x_val
        self.cur_cost = cost
        self.broadcast_cost_fun_call()
        self.iter_count += 1

        return cost

    def pred_cost_fun(self, x_val):
        """
        Returns the cost, predicted from theory, rather than estimated from
        data as in cost_fun(). This method mimics the method cost_fun(),
        but that one wraps self.emp_mhamil.get_mean_val(). This one wraps
        self.pred_mhamil.get_mean_val().

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
        assert self.pred_mhamil.num_samples == 0,\
            'predict cost with zero samples'
        cost = self.pred_mhamil.get_mean_val(var_num_to_rads)
        # print('bbvvv-cost', cost)
        return cost

    def find_min(self, minlib, **kwargs):
        """
        This method finds minimum of cost function. It allows user to choose
        among several possible minlibs, namely, 'scipy', 'autograd',
        'pytorch', 'tflow'. minlib parameters can be passed in via kwargs.

        Non-scipy minlibs do backprop. It is known that the complexity of
        calculating forward propagation and back propagation are about the
        same.

        kwargs (keyword arguments)
        minlib = scipy
            the keyword args of scipy.optimize.minimize
        minlib = autograd
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
        minlib : str
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
        if minlib == 'scipy':
            import scipy
            minimizer_fun = scipy.optimize.minimize
            opt_result = minimizer_fun(self.cost_fun,
                self.init_x_val, **kwargs)
            if self.verbose:
                print('*********final optimum result'
                      ' (final step=' + str(self.iter_count) +
                      '):\n', opt_result)
            return opt_result
        elif minlib == 'autograd':
            from autograd import grad

            assert 'num_iter' in kwargs, \
                "must pass-in keyword 'num_iter=' " \
                "if using autograd minlib"
            num_iter = kwargs['num_iter']
            assert 'descent_rate' in kwargs, \
                "must pass-in keyword 'descent_rate=' " \
                "if using autograd minlib"
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

        elif minlib == 'pytorch':
            assert False, 'not yet'
        elif minlib == 'tflow':
            assert False, 'not yet'
        else:
            assert False, 'unsupported find_min() minlib'


if __name__ == "__main__":
    def main():
        print(5)
    main()
