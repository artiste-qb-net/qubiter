import scipy
from openfermion.ops import QubitOperator
import utilities_gen as utg
from CodaSEO_writer import *
from SEO_simulator import *
from adv_applications.MeanHamilMinimizer import *


class MeanHamilMinimizer_naive(MeanHamilMinimizer):
    """
    This class is a child of MeanHamilMinimizer. The purpose of this class
    is to implement a naive but pedagogical version of what is called in
    Qubiter the Mean Hamiltonian Minimization problem.

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
    empirically on a qc which resides inside the QBox.

    To minimize a function of N continuous parameters, one can use some
    methods like simulated annealing that do not require calculating
    derivatives, or one can use methods that do use derivatives.

    The family of gradient descent methods (e.g., those that use auto
    differentiation and back-propagation) use only first derivatives.

    Other methods use both first (Jacobian) and second (Hessian)
    derivatives. Powell and Nelder-Mead are two methods of this kind. They
    are both implemented in the Numpy function `scipy.optimize.minimize`.

    Among the methods that do use derivatives, the performance of those that
    use both 1st and 2nd derivatives degrades quickly as N grows. Besides,
    calculating 2nd derivatives is very expensive. Hence, methods that use
    the 2nd derivatives are practically useless in the neural network field
    where N is usually very large. In that field, back-propagation rules.

    Qubiter's class `MeanHamilMinimizer_naive` is a wrapper class designed
    mainly around `scipy.optimize.minimize`. This scipy umbrella method
    implements many of those methods that use 1st and 2nd derivatives. That
    is why I call this class naive. So why write and study this class at
    all? Because it is pedagogical and allows one to set up some of the
    structures used to implement more advanced methods that do use
    back-propagation.

    """

    def __init__(self, *args, **kwargs):
        """

        Parameters
        ----------
        args : list
        kwargs : dict

        Returns
        -------

        """
        MeanHamilMinimizer.__init__(self, *args, **kwargs)

    def find_min(self):
        """
        This method wraps the method scipy.optimize.minimize

        Returns
        -------
        OptimizeResult
            OptimizeResult is a class (basically an enum) defined in
            scipy.optimize to hold the output results of
            scipy.optimize.minimize

        """
        opt_result = self.minimizer_fun(self.cost_fun,
            self.init_x_val, **self.mfun_kwargs)
        if self.verbose:
            print('*********final optimum result'
                  ' (final iter=' + str(self.iter_count) + '):\n', opt_result)
        return opt_result

    @staticmethod
    def hamil_mean_val(var_num_to_rads, fun_name_to_fun, 
                       file_prefix, num_bits, hamil,
                       num_samples=0, rand_seed=None):
        """
        This method calculates the mean value of the Hamiltonian hamil. See
        docstrings of class MeanHamilMinimizer for explanations of the
        various arguments of this method

        Parameters
        ----------
        var_num_to_rads : dict[int, float]
        fun_name_to_fun : dict[str, function]
        file_prefix : str
        num_bits : int
        hamil : QubitOperator
        num_samples : int
        rand_seed : int

        Returns
        -------
        float

        """
        # give it name unlikely to exist already
        fin_file_prefix = file_prefix + '99345125047'
        init_st_vec = StateVec.get_ground_st_vec(num_bits)
        arr_1 = np.array([1., 1.])
        arr_z = np.array([1., -1.])
        mean_val = 0
        for term, coef in hamil.terms.items():
            # we have checked before that coef is real
            coef = complex(coef).real
            bit_pos_to_xy_str = {}
            arr_list = [arr_1]*num_bits
            for bit_pos, action in term:
                arr_list[bit_pos] = arr_z
                if action != 'Z':
                    bit_pos_to_xy_str[bit_pos] = action
            wr = CodaSEO_writer(file_prefix, fin_file_prefix, num_bits)
            wr.write_xy_measurements(bit_pos_to_xy_str)
            wr.close_files()
            
            vman = PlaceholderManager(
                var_num_to_rads=var_num_to_rads,
                fun_name_to_fun=fun_name_to_fun)
            sim = SEO_simulator(fin_file_prefix, num_bits, init_st_vec,
                                vars_manager=vman)
            fin_st_vec = sim.cur_st_vec_dict['pure']
            real_vec = utg.kron_prod(arr_list)
            real_vec = np.reshape(real_vec, tuple([2]*num_bits))
            if not num_samples:
                mean_val += coef*fin_st_vec.\
                        get_mean_value_of_real_diag_mat(real_vec).real
            else:  # if num_samples !=0, then
                # sample qubiter-generated empirical prob dist
                pd = fin_st_vec.get_pd()
                obs_vec = StateVec.get_observations_vec(num_bits,
                        pd, num_samples, rand_seed=rand_seed)
                counts_dict = StateVec.get_counts_from_obs_vec(num_bits,
                                                               obs_vec)
                emp_pd = StateVec.get_empirical_pd_from_counts(num_bits,
                                                               counts_dict)
                # print('mmmmmmmm,,,', np.linalg.norm(pd-emp_pd))
                emp_st_vec = StateVec.get_emp_state_vec_from_emp_pd(
                        num_bits, emp_pd)
                mean_val += coef*emp_st_vec.\
                        get_mean_value_of_real_diag_mat(real_vec).real

        # create this coda writer in order to delete final files
        wr1 = CodaSEO_writer(file_prefix, fin_file_prefix, num_bits)
        wr1.close_and_del_fin_files()

        return mean_val

    def cost_fun(self, x_val):
        """
        This method wraps the static method hamil_mean_val() defined
        elsewhere in this class. This method will also print out whenever it
        is asked a report of the current values in x and cost.

        Parameters
        ----------
        x_val : tuple[float]

        Returns
        -------
        float

        """
        var_num_to_rads = dict(zip(self.all_var_nums, x_val))
        cost = MeanHamilMinimizer_naive.hamil_mean_val(
            var_num_to_rads, self.fun_name_to_fun,
            self.file_prefix, self.num_bits, self.hamil,
            self.num_samples, self.rand_seed)
        
        self.cur_x_val = x_val
        self.cur_cost = cost
        self.iter_count += 1
        self.broadcast_cost_fun_call()
            
        return cost

if __name__ == "__main__":
    def main():
        num_bits = 4
        file_prefix = 'io_folder/mean_hamil_naive_test'
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

        init_var_num_to_rads = {1: np.pi/6, 2: np.pi/3}
        fun_name_to_fun = {'my_fun': my_fun}

        hamil = QubitOperator('X1 Y3 X1 Y1', .4) + QubitOperator('Y2 X1', .7)
        print('hamil=\n', hamil)

        minimizer_fun = scipy.optimize.minimize
        min_method = 'Powell'
        num_samples = 0
        rand_seed = 1234
        print_hiatus = 25
        verbose = False

        def case():
            return MeanHamilMinimizer_naive(file_prefix, num_bits, hamil,
                init_var_num_to_rads, fun_name_to_fun,
                minimizer_fun, num_samples, rand_seed,
                print_hiatus, verbose, method=min_method).find_min()
        num_samples = 0
        case()
        print("*************************************")
        num_samples = 1000
        case()
    main()
