from adv_applications.MeanHamil import *
from adv_applications.StairsDerivCkt_writer import *


class StairsDeriv(MeanHamil):
    """
    This abstract class is a child of abstract class MeanHamil. The main
    purpose of its children classes is to override the method get_mean_val()
    of its parent class MeanHamil. The override returns the partial
    derivatives of a quantum cost function defined from a stairs circuit (an
    object of of class StairsCkt.writer) and a hamil (an object of the
    Openfermion class QubitOperator).


    Attributes
    ----------
    deriv_direc_to_dpart_range : dict[int, list[str]]
    deriv_gate_str : str
    gate_str_to_rads_list : dict[str, list[float|str]]

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
        all_var_nums = StairsCkt_writer.\
            get_all_var_nums(gate_str_to_rads_list)
        fun_name_to_fun = {}  # filled by get_mean_val()
        MeanHamil.__init__(self, file_prefix, parent_num_bits+1, hamil,
            all_var_nums, fun_name_to_fun, **kwargs)
        # ancilla at bit pos parent_num_bits
        MeanHamil.check_hamil_is_in_range(hamil, parent_num_bits-1)
        self.deriv_gate_str = deriv_gate_str
        self.gate_str_to_rads_list = gate_str_to_rads_list

        self.deriv_direc_to_dpart_range = {
            0: ['single'],
            1: ['1', 's', '1s'],
            2: ['1', 's', '1s'],
            3: ['1', 's', '1s']
        }

    def get_mean_val(self, var_num_to_rads):
        """
        Abstract method. The main goal of subclasses of this class is to
        override this method.

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------

        """
        assert False

if __name__ == "__main__":
    def main():
        print(5)
    main()
