from adv_applications.MeanHamilMinimizer import *


class MeanHamilMinimizer_native(MeanHamilMinimizer):
    """
    This class is a child of MeanHamilMinimizer. The purpose of this class
    is to implement an "all native' version of what is called in Qubiter the
    Mean Hamiltonian Minimization problem.

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

    Qubiter's class `MeanHamilMinimizer_native` can perform minimization via
    `scipy.optimize.minimize`. This scipy umbrella method implements many
    minimization methods, including Powell and CG.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    `MeanHamilMinimizer_native` does not call real physical hardware to do
    the simulation. Instead, it uses Qubiter's built-in simulator,
    `SEO_simulator`. That is why we call this class native.

    """

    def __init__(self, *args, **kwargs):
        """
        Constructor

        Parameters
        ----------
        args : list
            positional arguments of MeanHamilMinimizer
        kwargs : dict
            key-word arguments of MeanHamilMinimizer

        Returns
        -------

        """
        MeanHamilMinimizer.__init__(self, *args, **kwargs)

    def emp_hamil_mean_val(self, var_num_to_rads):
        """
        This method returns the empirically determined Hamiltonian mean
        value. Takes as input the values of placeholder variables. In this
        'native" case, we fake our data using SEO_simulator.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        float

        """
        return self.pred_hamil_mean_val(var_num_to_rads,
                                        num_fake_samples=self.num_samples)

if __name__ == "__main__":
    def main():
        num_bits = 4
        file_prefix = '../io_folder/mean_hamil_native_test'
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_Rx(2, rads=np.pi/7)
        wr.write_Rx(1, rads='#2*.5')
        wr.write_Rn(3, rads_list=['#1', '-#1*3', '#2'])
        wr.write_Rx(1, rads='-my_fun#2#1')
        wr.write_cnot(2, 3)
        wr.close_files()

        def my_fun(x, y):
            return x + .5*y
        fun_name_to_fun = {'my_fun': my_fun}

        init_var_num_to_rads = {1: 2.1, 2: 3.4}

        hamil = QubitOperator('X1 Y3 X1 Y1', .4) + QubitOperator('Y2 X1', .7)
        print('hamil=\n', hamil)

        minimizer_fun = scipy.optimize.minimize
        min_method = 'Powell'
        num_samples = 0
        rand_seed = 1234
        print_hiatus = 25
        verbose = False

        def case(**kwargs):
            return MeanHamilMinimizer_native(file_prefix, num_bits, hamil,
                                             init_var_num_to_rads, fun_name_to_fun,
                                             num_samples=num_samples, rand_seed=rand_seed,
                                             print_hiatus=print_hiatus,
                                             verbose=verbose).find_min(interface='scipy', **kwargs)
        print("*************************************")
        min_method = 'Powell'
        num_samples = 0
        case(method=min_method)
    main()
