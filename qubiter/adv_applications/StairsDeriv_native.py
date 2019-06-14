from qubiter.adv_applications.StairsDeriv import *
from qubiter.SEO_Lista import *
import itertools as it


class StairsDeriv_native(StairsDeriv):
    """
    This class is a child of StairsDeriv. Its main purpose is to override
    the method get_mean_val() of its abstract parent class StairsDeriv. In
    this class, the simulation necessary to evaluate the output of
    get_mean_val() is done by native, Qubiter simulators.

    Attributes
    ----------

    """
    def __init__(self, deriv_gate_str, gate_str_to_rads_list,
                 file_prefix, parent_num_bits, hamil, **kwargs):
        """
        Constructor

        Parameters
        ----------
        deriv_gate_str : str
        gate_str_to_rads_list : dict[str, list[float|str]]
        file_prefix : str
        parent_num_bits : int
        hamil : QubitOperator
        kwargs : dict
            key-word arguments of MeanHamil

        Returns
        -------

        """
        StairsDeriv.__init__(self, deriv_gate_str, gate_str_to_rads_list,
                 file_prefix, parent_num_bits, hamil, **kwargs)

    def get_mean_val(self, var_num_to_rads):
        """
        This method returns a list partials_list consisting of 4 floats
        which are the partial derivatives wrt the 4 possible derivative
        directions ( deriv_direc), of the multi-controlled gate U specified
        by self.deriv_gate_str.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        list[float]

        """
        partials_list = [0., 0., 0., 0.]
        # number of bits with (i.e., including) ancilla
        num_bits_w_anc = self.num_bits
        for has_neg_polarity, deriv_direc in it.product(
                *[[False, True], range(4)]):
            if self.deriv_gate_str == 'prior':
                if has_neg_polarity:
                    has_neg_polarity = None
                else:
                    continue  # this skips iteration in loop
            for dpart_name in StairsDeriv.dpart_dict[deriv_direc]:
                emb = CktEmbedder(num_bits_w_anc, num_bits_w_anc)
                wr = StairsDerivCkt_writer(self.deriv_gate_str,
                    has_neg_polarity, deriv_direc, dpart_name,
                        self.gate_str_to_rads_list, self.file_prefix, emb)
                t_list = self.gate_str_to_rads_list[self.deriv_gate_str]
                coef_of_dpart = StairsDerivCkt_writer.\
                    get_coef_of_dpart(t_list, deriv_direc,
                                      dpart_name, var_num_to_rads)
                fun_name_to_fun = StairsDerivCkt_writer.\
                    get_fun_name_to_fun(t_list, deriv_direc, dpart_name)
                lili = SEO_Lista.eng_file_to_line_list(wr.file_prefix,
                                                       num_bits_w_anc)
                lista = SEO_Lista(num_bits_w_anc, line_list=lili)
                len_lista_in = len(lista.line_list)
                for term, coef in self.hamil.terms.items():
                    # we have checked before that coef is real
                    coef = complex(coef).real
                    # print('nnnnnbbbbb', term)
                    new_term = tuple(list(term) + [(num_bits_w_anc-1, 'X')])
                    # print('jjjjjjj', new_term)

                    # throw out previous coda
                    lista = lista[:len_lista_in]

                    # add measurement coda for this term of hamil
                    # and for X at ancilla
                    bit_pos_to_xy_str =\
                        {bit: action for bit, action in new_term
                         if action != 'Z'}
                    lista.add_xy_meas_coda(bit_pos_to_xy_str)

                    # run simulation. get fin state vec
                    vman = PlaceholderManager(
                        var_num_to_rads=var_num_to_rads,
                        fun_name_to_fun=fun_name_to_fun)
                    # simulator will change init_st_vec so use
                    # fresh copy of it each time
                    init_st_vec = cp.deepcopy(self.init_st_vec)
                    sim = lista.simulate(init_st_vec=init_st_vec,
                                         vars_manager=vman)
                    fin_st_vec = sim.cur_st_vec_dict['pure']

                    # get effective state vec
                    if self.num_samples:
                        # if num_samples !=0, then
                        # sample qubiter-generated empirical prob dist
                        pd = fin_st_vec.get_pd()
                        obs_vec = StateVec.get_observations_vec(
                            num_bits_w_anc, pd, self.num_samples)
                        counts_dict = StateVec.get_counts_from_obs_vec(
                            num_bits_w_anc, obs_vec)
                        emp_pd = StateVec.get_empirical_pd_from_counts(
                            num_bits_w_anc, counts_dict)
                        # print('mmmmmmmm,,,', np.linalg.norm(pd-emp_pd))
                        emp_st_vec = StateVec.\
                            get_emp_state_vec_from_emp_pd(num_bits_w_anc,
                                                          emp_pd)
                        effective_st_vec = emp_st_vec
                    else:  # num_samples = 0
                        effective_st_vec = fin_st_vec

                    # add contribution to mean
                    real_arr = self.get_real_vec(new_term)
                    mean_val_change = coef*effective_st_vec.\
                            get_mean_value_of_real_diag_mat(real_arr)
                    mean_val_change *= coef_of_dpart
                    if has_neg_polarity:
                        mean_val_change *= -1
                    partials_list[deriv_direc] += mean_val_change
        return partials_list


if __name__ == "__main__":
    def main():
        num_bits = 4
        parent_num_bits = num_bits - 1  # one bit for ancilla

        # u2_bit_to_higher_bits = None
        u2_bit_to_higher_bits = {0: [2], 1: [2], 2: []}
        gate_str_to_rads_list = StairsCkt_writer.\
            get_gate_str_to_rads_list(parent_num_bits,
                '#int', rads_const=np.pi/2,
                u2_bit_to_higher_bits=u2_bit_to_higher_bits)
        pp.pprint(gate_str_to_rads_list)

        deriv_gate_str = list(gate_str_to_rads_list.keys())[2]

        file_prefix = '../io_folder/stairs_deriv_native_test'

        hamil = QubitOperator('X1 Y0 X1 Y1', .4) +\
            QubitOperator('Y2 X1', .7)

        der = StairsDeriv_native(deriv_gate_str,
                                 gate_str_to_rads_list, file_prefix,
                                 parent_num_bits, hamil)

        var_num_to_rads = StairsCkt_writer.get_var_num_to_rads(
            gate_str_to_rads_list, 'const', rads_const=np.pi/2)

        partials_list = der.get_mean_val(var_num_to_rads)
        print('partials_list=', partials_list)
    main()
