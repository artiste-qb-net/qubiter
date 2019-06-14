from qubiter.adv_applications.StairsDeriv_rigetti import *
import pprint as pp


class StairsAllDeriv_rigetti(StairsDeriv_rigetti):
    """
    This class is a child of StairsDeriv_rigetti. For the parent class,
    the get_mean_val() method returns a list of 4 partial derivatives
    belonging to a particular gate string (a gate_str is a key in
    gate_str_to_rads_list). For this class, get_mean_val() returns an
    ordered dictionary mapping each gate_str to its 4 partials.

    Attributes
    ----------
    deriv_gate_str : str

    """
    def __init__(self, qc, gate_str_to_rads_list,
                 file_prefix, parent_num_bits, hamil, **kwargs):
        """
        Constructor

        Parameters
        ----------
        qc : QuantumComputer
        gate_str_to_rads_list : dict[str, list[float|str]]
        file_prefix : str
        parent_num_bits : int
        hamil : QubitOperator
        kwargs : dict
            key-word arguments of MeanHamil

        Returns
        -------

        """
        deriv_gate_str = 'dummy'
        StairsDeriv_rigetti.__init__(self, qc, deriv_gate_str,
                                    gate_str_to_rads_list, file_prefix,
                                    parent_num_bits, hamil, **kwargs)

    def get_mean_val(self, var_num_to_rads):
        """
        This method returns an ordered dictionary gate_str_to_partials_list
        mapping each gate_str to its 4 partial derivatives. This method
        calls the get_mean_val() of the parent class for all possible
        gate_str.

        The output dictionary of this method can be converted to a numpy
        array using
        StairDerivCkt_writer.make_array_from_gate_str_to_rads_list()

        Parameters
        ----------
        var_num_to_rads : dict[int, float]

        Returns
        -------
        dict[str, list[float]]

        """
        gate_str_to_partials_list = OrderedDict()
        for gate_str in self.gate_str_to_rads_list.keys():
            self.deriv_gate_str = gate_str
            gate_str_to_partials_list[gate_str] =\
                StairsDeriv_rigetti.get_mean_val(self, var_num_to_rads)
        return gate_str_to_partials_list

if __name__ == "__main__":
    def main():
        print(5)
    main()
