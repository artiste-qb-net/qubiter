import copy as cp

from adv_applications.MeanHamilMinimizer import *
from device_specific.Qubiter_to_RigettiPyQuil import *
from device_specific.Cloud_rigetti import *
import utilities_gen as utg

from openfermion.ops import QubitOperator

from pyquil import Program, Pragma
from pyquil.gates import *


class MeanHamilMinimizer_rigetti(MeanHamilMinimizer):
    """
    
    Attributes
    ----------
    device_name : str
    do_resets : bool
    len_pg_in : int
    pg : Program
    translator : Qubiter_to_RigettiPyQuil
    
    """

    def __init__(self, device_name, do_resets=True,
                 *args, **kwargs):
        """
        Constructor

        Do in constructor as much hamil indep stuff as possible so don't
        have to redo it with every call to cost fun

        Parameters
        ----------
        device_name : str
        do_resets : bool
        args : list
        kwargs : dict

        Returns
        -------

        """
        MeanHamilMinimizer.__init__(self, *args, **kwargs)
        self.device_name = device_name
        self.do_resets = do_resets

        # this creates a file with all Pyquil gates that
        # are independent of hamil, with free parameters
        self.translator = Qubiter_to_RigettiPyQuil(
            self.file_prefix, self.num_bits, prelude_str='', ending_str='')

        # pg prelude
        self.pg = Program()
        self.pg += Pragma('INITIAL_REWIRING', ['"PARTIAL"'])
        if self.do_resets:
            self.pg += RESET()
        qubits = list(self.pg.get_qubits())
        assert len(qubits) == self.num_bits
        ro = self.pg.declare('ro', 'BIT', self.num_bits)
        s = ''
        for var_num in self.all_var_nums:
            vname = self.translator.vprefix + str(var_num)
            s += vname
            s += ' = pg.declare("'
            s += vname
            s += '", memory_type="REAL")\n'
        eval(s)

        # pg operations independent of hamil
        with open(self.translator.aqasm_path, 'r') as fi:
            lines = fi.readlines()
        for line in lines:
            line = line.strip('\n')
            if line:
                eval(line)

        self.len_pg_in = len(self.pg)

    @staticmethod
    def add_xy_meas_coda_to_program(prog, bit_pos_to_xy_str):
        """
        
        Parameters
        ----------
        prog : Program
        bit_pos_to_xy_str : dict[int, str]

        Returns
        -------
        None

        """
        for bit_pos, xy_str in bit_pos_to_xy_str.items():
            if xy_str == 'X':
                # exp(-i*sigy*pi/4)*sigz*exp(i*sigy*pi/4) = sigx
                prog += RY(-np.pi/2, bit_pos)
            elif xy_str == 'Y':
                # exp(i*sigx*pi/4)*sigz*exp(-i*sigx*pi/4) = sigy
                prog += RX(np.pi/2, bit_pos)
            else:
                assert False, "Unsupported qbit measurement. '" + \
                            xy_str + "' Should be either 'X' or 'Y'"

    def hamil_mean_val(self, var_num_to_rads):
        """
        This method calculates the mean value of the Hamiltonian hamil.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        float

        """
        # hamil loop
        qubits = list(self.pg.get_qubits())
        arr_1 = np.array([1., 1.])
        arr_z = np.array([1., -1.])
        mean_val = 0
        for term, coef in self.hamil.terms.items():
            # we have checked before that coef is real
            coef = complex(coef).real
            
            # add measurement coda for this term of hamil
            # build real_vec from arr_list.
            # real_vec will be used at end of loop
            arr_list = [arr_1]*self.num_bits
            bit_pos_to_xy_str = {}
            for bit_pos, action in term:
                arr_list[bit_pos] = arr_z
                if action != 'Z':
                    bit_pos_to_xy_str[bit_pos] = action
            real_vec = utg.kron_prod(arr_list)
            real_vec = np.reshape(real_vec, tuple([2]*self.num_bits))
            # reset pg to initial length
            self.pg = self.pg[:self.len_pg_in]
            MeanHamilMinimizer_rigetti.add_xy_meas_coda_to_program(
                self.pg, bit_pos_to_xy_str)
            # request measurements
            for i, q in enumerate(qubits):
                self.pg += MEASURE(q, ro[i])

            # send and receive from cloud, get obs_vec
            qc = Cloud_rigetti.get_qc(self.device_name)
            vprefix = self.translator.vprefix
            var_name_to_rads = {vprefix + str(vnum): rads
                for vnum, rads in var_num_to_rads.items()}
            bitstrings = Cloud_rigetti.get_bitstrings(pg, qc,
                        var_name_to_rads, self.num_samples)
            obs_vec = Cloud_rigetti.obs_vec_from_bitstrings(
                    bitstrings, self.num_bits)
            
            # go from obs_vec to effective state vec
            counts_dict = StateVec.get_counts_from_obs_vec(self.num_bits,
                                                           obs_vec)
            emp_pd = StateVec.get_empirical_pd_from_counts(self.num_bits,
                                                           counts_dict)
            emp_st_vec = StateVec.get_emp_state_vec_from_emp_pd(
                    self.num_bits, emp_pd)
            effective_st_vec = emp_st_vec
            
            # add contribution to mean
            mean_val += coef*effective_st_vec.\
                    get_mean_value_of_real_diag_mat(real_vec).real

        return mean_val

if __name__ == "__main__":
    def main():
        pass
