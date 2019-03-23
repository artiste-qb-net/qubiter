from adv_applications.MeanHamilMinimizer import *


class MeanHamilMinimizer_native(MeanHamilMinimizer):
    """
    This class is a child of MeanHamilMinimizer.

    This class does not call real physical hardware, or someone else's
    simulator to calculate mean values. Instead, it uses Qubiter's built-in
    simulator, `SEO_simulator`. That is why we call this class native.

    Attributes
    ----------

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
        value. It takes as input the values of placeholder variables. In
        this "native" case, we use SEO_simulator to fake the data used to
        calculate the output mean value.

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
    import scipy
    from openfermion.ops import QubitOperator

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
                                             init_var_num_to_rads,
                                             fun_name_to_fun,
                                             num_samples=num_samples,
                                             rand_seed=rand_seed,
                                             print_hiatus=print_hiatus,
                                             verbose=verbose).find_min(
                                                interface='scipy', **kwargs)
        print("*************************************")
        min_method = 'Powell'
        num_samples = 0
        case(method=min_method)
    main()
