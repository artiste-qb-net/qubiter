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
    software libraries) such as 'scipy', 'autograd', 'tflow', 'pytorch'. It
    can also use various devices (aka simulators or backends), either
    virtual or real, to do the minimization. For example, tensorflow is a
    minlib and SEO_simulator_tf is a backend that is native to qubiter and
    uses tensorflow. By a native device, we mean one that uses Qubiter
    native simulators like SEO_simulator and SEO_simulator_tf.

    Non-scipy minlibs implement backprop.

    The 'scipy' minlib is a wrapper for the scipy function
    `scipy.optimize.minimize`. This scipy umbrella method implements many
    minimization methods, including Powell and CG.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

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
    targ_mhamil : MeanHamil
        Target mean Hamiltonian, used to evaluate targ_cost

    """

    def __init__(self, emp_mhamil, targ_mhamil,
                 all_var_nums, init_var_num_to_rads,
                 print_hiatus=1, verbose=False):
        """
        Constructor

        Parameters
        ----------
        emp_mhamil : MeanHamil
        targ_mhamil : MeanHamil
        all_var_nums : list[int]
        init_var_num_to_rads : dict[int, float]
        print_hiatus : int
        verbose : bool

        Returns
        -------

        """
        self.emp_mhamil = emp_mhamil
        self.targ_mhamil = targ_mhamil
        self.all_var_nums = all_var_nums
        assert emp_mhamil.all_var_nums == all_var_nums
        assert targ_mhamil.all_var_nums == all_var_nums
        init_x_val = [init_var_num_to_rads[k] for k in all_var_nums]
        self.init_x_val = np.array(init_x_val)

        CostMinimizer.__init__(self, init_x_val, print_hiatus, verbose)

    def cost_fun(self, x_val):
        """
        This method wraps self.emp_mhamil.get_mean_val(). This method will
        also print out, whenever it is called, a report of the current
        values of x and cost (and targ_cost if it is available).

        Parameters
        ----------
        x_val : np.ndarray

        Returns
        -------
        float

        """

        var_num_to_rads = dict(zip(self.all_var_nums, list(x_val)))
        cost = self.emp_mhamil.get_mean_val(var_num_to_rads)

        self.cur_x_val = x_val
        self.cur_cost = cost
        self.broadcast_cost_fun_call()
        self.iter_count += 1

        return cost

    def targ_cost_fun(self, x_val):
        """
        Returns the cost, predicted from theory, rather than estimated from
        data as in cost_fun(). This method mimics the method cost_fun(),
        but that one wraps self.emp_mhamil.get_mean_val(). This one wraps
        self.targ_mhamil.get_mean_val().

        Parameters
        ----------
        x_val : np.ndarray

        Returns
        -------
        float

        """
        # print('inside_targ_cost', x_val)
        var_num_to_rads = dict(zip(self.all_var_nums, list(x_val)))
        # print('mmmmmmmmmaaaaaa', var_num_to_rads)
        assert self.targ_mhamil.num_samples == 0,\
            'predict cost with zero samples'
        cost = self.targ_mhamil.get_mean_val(var_num_to_rads)
        # print('bbvvv-cost', cost)
        return cost

    def find_min(self, minlib, **kwargs):
        """
        This method finds minimum of cost function. It allows user to choose
        among several possible minlibs, namely, 'scipy', 'autograd',
        'tflow', 'pytorch'. minlib parameters can be passed in via kwargs.

        kwargs (keyword arguments)
        minlib = scipy
            the keyword args of scipy.optimize.minimize
        minlib = autograd, tflow
            num_inter : float
                number of iterations (an iteration is every time call
                cost_fun)
            descent_rate : float
                positive float, constant that multiplies gradient of
                cost function being minimized. Often denoted as eta

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

        def targ_cost(xx):
            # self argument seems to confuse grad
            return self.targ_cost_fun(xx)
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

            num_iter = utg.get_value(kwargs, 'num_iter')
            rate = utg.get_value(kwargs, 'descent_rate')

            self.cur_x_val = self.init_x_val

            for step in range(num_iter):
                xlist = list(self.cur_x_val)
                # print('mmbbb', self.cur_x_val, xlist)
                self.cur_targ_cost = targ_cost(self.cur_x_val)
                self.cur_cost = self.cost_fun(self.cur_x_val)
                # print('kkkhhh', grad(targ_cost)(self.cur_x_val))
                for dwrt in range(len(xlist)):
                    self.cur_x_val[dwrt] -= \
                        rate*grad(targ_cost)(self.cur_x_val)[dwrt]

        elif minlib == 'tflow':
            import tensorflow as tf
            assert tf.executing_eagerly()

            num_iter = utg.get_value(kwargs, 'num_iter')
            rate = utg.get_value(kwargs, 'descent_rate')

            self.init_x_val = tf.convert_to_tensor(self.init_x_val)
            self.cur_x_val = self.init_x_val

            # optimizer = tf.train.GradientDescentOptimizer(rate)

            for step in range(num_iter):
                with tf.GradientTape() as tape:
                    tape.watch(self.cur_x_val)
                    self.cur_targ_cost = targ_cost(self.cur_x_val)
                    # print('**********cccccccc', self.cur_x_val,
                    #       self.cur_targ_cost)
                self.cur_cost = self.cost_fun(self.cur_x_val)
                grads = tape.gradient(self.cur_targ_cost, self.cur_x_val)
                self.cur_x_val -= rate*grads
                # optimizer.apply_gradients(zip(grads, self.cur_x_val))

        elif minlib == 'pytorch':
            assert False, 'not yet'
        else:
            assert False, 'unsupported find_min() minlib'
    

if __name__ == "__main__":
    def main():
        print(5)
    main()
