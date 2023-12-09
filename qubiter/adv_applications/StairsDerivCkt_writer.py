from qubiter.SEO_writer import *
from qubiter.adv_applications.StairsCkt_writer import *


class StairsDerivCkt_writer(SEO_writer):
    """
    This class is a subclass of `SEO_writer`. It writes several intermediary
    stairs derivative circuits that will be used in class
    `StairsDeriv_native` for calculating the gradients of a quantum cost
    function (mean hamiltonian).

    Suppose U = exp[i*(t_0 + t_1*sigx + t_2*sigy + t_3*sigz)], where sigx,
    sigy, sigz are the Pauli matrices and t_r for r in range(4) are 4 real
    parameters. To take the derivative wrt t_r of a given multi-controlled
    gate U in a stairs circuit, we need to evaluate several circuits (we
    call them dparts, which stands for derivative parts). Say, for instance,
    that ``GATE= @---O---+---U``. To calculate d/dt_r GATE(t_0, t_1, t_2, t_3),
    for r=0,1, 2, 3, we need to calculate a new circuit wherein the GATE in
    the parent circuit is replaced by::

        sum_k  c_k  @---@---O---+---U_k

    (which is said to have `has_neg_polarity`=False) and::

        sum_k  c_k  @---@---O---+---U_k
                    Z---@---O   |   |

    (which is said to have `has_neg_polarity`=True)

    Also, some extra stuff (a coda) must be appended to the end of the
    parent stairs circuit.

    Note that an extra "ancilla" qbit has been added (as the new last qubit)
    to the parent stairs circuit being differentiated. So if the parent
    stairs circuit has a number `parent_num_qbits` of qubits, then the one
    written by this class has that many qubits plus one.

    The index r which is in range(4) is called the derivative direction (
    `deriv_direc`)

    `gate_str_to_rads_list` is the same as for the parent stairs circuit.

    `deriv_gate_str` is a well formed gate_str that specifies which U is
    being differentiated

    The index k is given as a string called `dpart_name` ("dpart" stands for
    derivative part).

    The coefficients c_k can be obtained via the method get_coef_of_dpart()

    Each U_k is a U(2) matrix itself, and its 4 parameters are defined in
    terms of the parameters tlist=[t_0, t_1, t_2, t_3] of the U(tlist) being
    differentiated, via 4 functions of tlist. These functions can be
    obtained via the method get_fun_name_to_fun().

    Attributes
    ----------
    deriv_direc : int
        in range(4)
    deriv_gate_str : str
    dpart_name : str
    gate_str_to_rads_list : dict[str, list[float]]
    has_neg_polarity : bool


    """

    def __init__(self, deriv_gate_str, has_neg_polarity,
                 deriv_direc, dpart_name,
                 gate_str_to_rads_list,
                 file_prefix, emb, **kwargs):
        """
        Constructor

        This constructor writes English and Picture files but it doesn't
        close those files after writing them. You must do that yourself
        using close_files().


        Parameters
        ----------
        deriv_gate_str : str
        has_neg_polarity : bool
        deriv_direc : int
            in range(4)
        dpart_name : str
        gate_str_to_rads_list : dict[int, list[float]]
        file_prefix : str
        emb : CktEmbedder
        kwargs : dict
            key-word arguments of SEO_writer

        Returns
        -------

        """
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        self.deriv_gate_str = deriv_gate_str
        self.has_neg_polarity = has_neg_polarity
        self.deriv_direc = deriv_direc
        self.dpart_name = dpart_name
        self.gate_str_to_rads_list = gate_str_to_rads_list

        assert deriv_gate_str in gate_str_to_rads_list.keys()
        assert deriv_direc in range(4)
        if deriv_gate_str == 'prior':
            assert has_neg_polarity is None

        self.write()

    def write(self):
        """
        This method writes English and Picture files for a quantum circuit
        Der. Der is used in calculating the derivative of a parent stairs
        circuit Pa with respect to one of 4 parameters for each of the
        multi-controlled U's of Pa.

        Returns
        -------

        """
        num_qbits = self.emb.num_qbits_bef
        anc_bit_pos = num_qbits-1

        for gate_str, rads_list in self.gate_str_to_rads_list.items():
            u2_pos = self.get_u2_pos(gate_str)
            trols = StairsCkt_writer.get_controls_from_gate_str(
                num_qbits, gate_str)
            if gate_str == self.deriv_gate_str:
                # add control on ancilla qubit
                trols.bit_pos_to_kind[anc_bit_pos] = True
                trols.refresh_lists()
                deriv_rads_list = \
                    self.gate_str_to_rads_list[self.deriv_gate_str]
                if not isinstance(deriv_rads_list[0], str):
                    t_list = self.gate_str_to_rads_list[self.deriv_gate_str]
                    rads_list = StairsDerivCkt_writer.\
                        get_rads_list_of_dpart(t_list,
                                               self.deriv_direc,
                                               self.dpart_name)
                else:
                    rads_list = ['rads' + deriv_rads_list[k][1:] +
                                 ''.join(deriv_rads_list) for k in range(4)]
                self.write_H(anc_bit_pos)
            else:
                rads_list = self.gate_str_to_rads_list[gate_str]
            self.write_controlled_one_qbit_gate(u2_pos, trols,
                    OneQubitGate.u2, rads_list)
            if gate_str == self.deriv_gate_str and gate_str != 'prior' and \
                    self.has_neg_polarity:
                # remove ancilla qubit control that was added previously
                del trols.bit_pos_to_kind[anc_bit_pos]
                trols.refresh_lists()
                # add controlled sigz for neg X term
                self.write_controlled_one_qbit_gate(anc_bit_pos, trols,
                    OneQubitGate.sigz)

    def get_u2_pos(self, gate_str):
        """
        Given a well formed gate_str (one of the keys of
        gate_str_to_rads_list), this method returns the bit position of the
        U(2) matrix.

        Class StairsCkt_writer has a method of same name but which returns a
        different value for the same gate_str input. This is due to the fact
        that the circuit generated by this class has an extra ancilla qubit
        compared with its parent circuit.

        Parameters
        ----------
        gate_str : str

        Returns
        -------
        int

        """
        # last qubit position, num_qbits-1, is reserved for ancilla qubit
        # "prior" gate_str has U2 at bit pos num_qbits-2
        num_qbits = self.emb.num_qbits_bef
        if gate_str == 'prior':
            return num_qbits - 2
        else:
            return num_qbits - 2 - len(gate_str) // 2

    @staticmethod
    def float_t_list(t_list, var_num_to_rads):
        """
        Internal method that checks whether t_list is of form list[str] or
        list[float]. In first case, it uses var_num_to_rads to replace
        t_list by a list[float].

        Parameters
        ----------
        t_list : list[float|str]
        var_num_to_rads : list[int, float]

        Returns
        -------
        list[float]

        """

        if any([isinstance(t, str) for t in t_list]):
            assert var_num_to_rads is not None,\
                't_list has #int strings, items not all floats,' \
                ' so must provide a var_num_to_rads'
            t_list = [var_num_to_rads[int(t[1:])] for t in t_list]
        return t_list

    @staticmethod
    def get_rads_list_of_dpart(t_list, deriv_direc, dpart_name,
                               var_num_to_rads=None):
        """
        This method returns the rads list (list of 4 floats) for a given
        given: 1. the rads list `t_list` of the U that is being
        differentiated, 2. the dpart (derivative part), and 3. which of 4
        directions deriv_direc the derivative is wrt.

        var_num_to_rads is used if self wrote the English file with #int
        string symbols for placeholder variables. var_num_to_rads is used to
        float those strings.

        Parameters
        ----------
        t_list : list[float | str]
        deriv_direc : int
            int in range(4)
        dpart_name : str
        var_num_to_rads : dict[int, float] | None

        Returns
        -------
        list[float]

        """
        t_list = StairsDerivCkt_writer.\
            float_t_list(t_list, var_num_to_rads)
        t_vec = t_list[1:]
        t = np.linalg.norm(np.array(t_vec))

        if deriv_direc == 0:
            assert dpart_name == 'single'
            rads_array = np.array(t_list)
            rads_array[0] += np.pi/2
        else:  # deriv_direc in [1, 2, 3]
            rads_array = np.zeros((4,), dtype=float)
            rads_array[0] = t_list[0]
            if dpart_name == '1':
                if t > 1e-6:
                    for k in range(1, 4):
                        rads_array[k] += (np.pi/2 + t)*t_list[k]/t
            elif dpart_name == 's':
                rads_array[deriv_direc] += np.pi/2
            elif dpart_name == '1s':
                if t > 1e-6:
                    for k in range(1, 4):
                        rads_array[k] += (np.pi/2)*t_list[k]/t
            else:
                assert False, 'unsupported deriv part name'
        return list(rads_array)

    @staticmethod
    def get_coef_of_dpart(t_list,
                          deriv_direc,
                          dpart_name,
                          var_num_to_rads=None):
        """
        This method returns coefficient of dpart (derivative part), either
        p1, ps or -p1*ps

        var_num_to_rads is used if self wrote the English file with #int
        string symbols for placeholder variables. var_num_to_rads is used to
        float those strings.

        Parameters
        ----------
        t_list : list[float | str]
        deriv_direc : int
            int in range(4)
        dpart_name : str
        var_num_to_rads : dict[int, float] | None

        Returns
        -------
        float

        """
        # t stands for theta not time
        t_list = StairsDerivCkt_writer.\
            float_t_list(t_list, var_num_to_rads)

        # space components
        t_vec = t_list[1:]
        t = np.linalg.norm(np.array(t_vec))
        if t < 1e-6:
            ps = 1.
            p1 = 0.
        else:
            ps = np.sin(t)/t
            p1 = t_list[deriv_direc]/t

        if deriv_direc == 0:
            assert dpart_name == 'single'
            coef = 1.
        else:  # deriv_direc in [1, 2, 3]
            if dpart_name == '1':
                coef = p1
            elif dpart_name == 's':
                coef = ps
            elif dpart_name == '1s':
                coef = -p1*ps
            else:
                assert False, 'unsupported deriv part name'
        return coef

    @staticmethod
    def get_fun_name_to_fun(t_list, deriv_direc, dpart_name):
        """
        This method returns a dictionary fun_name_to_fun mapping the
        function name to function, for all functions defined by this class.

        Parameters
        ----------
        t_list : list[float | str]
        deriv_direc : int
        dpart_name : str

        Returns
        -------
        dict[str, function]

        """
        fun_name_to_fun = {}
        if isinstance(t_list[0], str):
            for k in range(4):
                fun_name = 'rads' + t_list[k][1:]
                fun_name_to_fun[fun_name] =\
                    (lambda t0, t1, t2, t3, k1=k,
                     x2=deriv_direc, x3=dpart_name:
                     StairsDerivCkt_writer.
                     get_rads_list_of_dpart([t0, t1, t2, t3], x2, x3)[k1])
        return fun_name_to_fun


if __name__ == "__main__":
    def main():
        num_qbits = 4
        parent_num_qbits = num_qbits - 1  # one bit for ancilla
        gate_str_to_rads_list = StairsCkt_writer.\
            get_gate_str_to_rads_list(
                parent_num_qbits, '#int', rads_const=np.pi/2)

        file_prefix = 'stairs_deriv_writer_test'
        emb = CktEmbedder(num_qbits, num_qbits)
        deriv_gate_str = list(gate_str_to_rads_list.keys())[2]
        for deriv_direc, dpart_name, has_neg_polarity in \
                [(0, 'single', None), (3, 's', True)]:
            wr = StairsDerivCkt_writer(deriv_gate_str,
                                       has_neg_polarity,
                                       deriv_direc,
                                       dpart_name,
                                       gate_str_to_rads_list,
                                       file_prefix, emb)
            wr.close_files()
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%")
            wr.print_eng_file()
            wr.print_pic_file()
    main()
