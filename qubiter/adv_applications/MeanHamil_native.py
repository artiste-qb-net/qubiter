from qubiter.adv_applications.MeanHamil import *
from qubiter.SEO_simulator_tf import *


class MeanHamil_native(MeanHamil):
    """
    This class is a child of MeanHamil.

    This class does not call real physical hardware, or someone else's
    simulator to calculate mean values. Instead, it uses Qubiter's built-in
    simulators, such as `SEO_simulator` and `SEO_simulator_tf`. That is why
    we call this class native.

    Attributes
    ----------
    list_of_supported_sims : list[str]
        list of the names of simulators supported by this class.
        self.simulator_name must be in this list.
    use_tf : bool

    """
    # class variable
    list_of_supported_sims = ['SEO_simulator',
                              'SEO_simulator_tf']

    def __init__(self, *args, **kwargs):
        """
        Constructor

        Parameters
        ----------
        args : list
            positional arguments of MeanHamil
        kwargs : dict
            key-word arguments of MeanHamil

        Returns
        -------

        """
        MeanHamil.__init__(self, *args, **kwargs)
        assert self.simulator_name in MeanHamil_native.\
            list_of_supported_sims
        self.use_tf = (self.simulator_name == 'SEO_simulator_tf')
        if self.use_tf:
            assert tf.executing_eagerly()

    def get_mean_val(self, var_num_to_rads):
        """
        This method predicts the mean value of the Hamiltonian hamil using
        only Qubiter simulators.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        float

        """
        # give it name unlikely to exist already
        fin_file_prefix = self.file_prefix + '99345125047'

        # hamil loop
        mean_val = 0
        for term, coef in self.hamil.terms.items():
            # we have checked before that coef is real
            coef = complex(coef).real

            # add measurement coda for this term of hamil
            wr = CodaSEO_writer(self.file_prefix,
                                fin_file_prefix, self.num_qbits)
            bit_pos_to_xy_str =\
                {bit: action for bit, action in term if action != 'Z'}
            wr.write_xy_measurements(bit_pos_to_xy_str)
            wr.close_files()

            # run simulation. get fin state vec
            vman = PlaceholderManager(
                var_num_to_rads=var_num_to_rads,
                fun_name_to_fun=self.fun_name_to_fun)
            # simulator will change init_st_vec so use
            # fresh copy of it each time
            init_st_vec = cp.deepcopy(self.init_st_vec)
            if self.simulator_name == 'SEO_simulator':
                sim = SEO_simulator(fin_file_prefix, self.num_qbits,
                                    init_st_vec, vars_manager=vman)
            elif self.simulator_name == 'SEO_simulator_tf':
                init_st_vec.arr = tf.convert_to_tensor(init_st_vec.arr)
                sim = SEO_simulator_tf(fin_file_prefix, self.num_qbits,
                                    init_st_vec, vars_manager=vman)
            else:
                assert False, 'unsupported native simulator'

            fin_st_vec = sim.cur_st_vec_dict['pure']

            # print('********bbbvvvvvv',
            # self.init_st_vec.arr, fin_st_vec.arr)

            # get effective state vec
            if self.num_samples:
                # if num_samples !=0, then
                # sample qubiter-generated empirical prob dist
                pd = fin_st_vec.get_pd()
                obs_vec = StateVec.get_observations_vec(self.num_qbits,
                        pd, self.num_samples)
                counts_dict = StateVec.get_counts_from_obs_vec(self.num_qbits,
                                                               obs_vec)
                emp_pd = StateVec.get_empirical_pd_from_counts(self.num_qbits,
                                                               counts_dict)
                # print('mmmmmmmm,,,', np.linalg.norm(pd-emp_pd))
                emp_st_vec = StateVec.get_emp_state_vec_from_emp_pd(
                        self.num_qbits, emp_pd)
                effective_st_vec = emp_st_vec
            else:  # num_samples = 0
                effective_st_vec = fin_st_vec

            # add contribution to mean
            real_arr = self.get_real_vec(term)
            if not self.use_tf:
                mean_val += coef*effective_st_vec.\
                        get_mean_value_of_real_diag_mat(real_arr)
            else:
                real_arr = tf.convert_to_tensor(real_arr, dtype=tf.complex128)
                arr = effective_st_vec.arr
                mean_val += coef*tf.reduce_sum(
                    tf.math.real(tf.math.conj(arr) * real_arr * arr))


        # create this writer in order to delete final files
        wr1 = SEO_writer(fin_file_prefix,
                CktEmbedder(self.num_qbits, self.num_qbits))
        wr1.delete_files()

        return mean_val


if __name__ == "__main__":
    def main():
        print(5)
    main()
