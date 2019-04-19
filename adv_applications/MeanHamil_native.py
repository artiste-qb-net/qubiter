from adv_applications.MeanHamil import *


class MeanHamil_native(MeanHamil):
    """
    This class is a child of MeanHamil.

    This class does not call real physical hardware, or someone else's 
    simulator to calculate mean values. Instead, it uses Qubiter's built-in 
    simulators, such as `SEO_simulator`. That is why we call this class 
    native. 

    Attributes
    ----------
    list_of_supported_sims : list[str]
        list of the names of simulators supported by this class.
        self.simulator_name must be in this list.

    """

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
        # can add to list of supported simulators in future
        self.list_of_supported_sims = ['SEO_simulator']
        assert self.simulator_name in self.list_of_supported_sims

    def get_mean_val(self, var_num_to_rads):
        """
        This method predicts the mean value of the Hamiltonian hamil using
        only Qubiter simulators and no data

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
                                fin_file_prefix, self.num_bits)
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
                sim = SEO_simulator(fin_file_prefix, self.num_bits,
                                    init_st_vec, vars_manager=vman)
            else:
                assert False, 'unsupported native simulator'
            fin_st_vec = sim.cur_st_vec_dict['pure']
            # print('inside pred hamil in/out st_vec',
            # self.init_st_vec, fin_st_vec)

            # get effective state vec
            if self.num_samples:
                # if num_samples !=0, then
                # sample qubiter-generated empirical prob dist
                pd = fin_st_vec.get_pd()
                obs_vec = StateVec.get_observations_vec(self.num_bits,
                        pd, self.num_samples)
                counts_dict = StateVec.get_counts_from_obs_vec(self.num_bits,
                                                               obs_vec)
                emp_pd = StateVec.get_empirical_pd_from_counts(self.num_bits,
                                                               counts_dict)
                # print('mmmmmmmm,,,', np.linalg.norm(pd-emp_pd))
                emp_st_vec = StateVec.get_emp_state_vec_from_emp_pd(
                        self.num_bits, emp_pd)
                effective_st_vec = emp_st_vec
            else:  # num_samples = 0
                effective_st_vec = fin_st_vec

            # add contribution to mean
            real_arr = self.get_real_vec(term)
            mean_val += coef*effective_st_vec.\
                    get_mean_value_of_real_diag_mat(real_arr)

        # create this writer in order to delete final files
        wr1 = SEO_writer(fin_file_prefix, self.num_bits)
        wr1.delete_files()

        return mean_val

if __name__ == "__main__":
    def main():
        print(5)
    main()
