from Controls import *
from SEO_reader import *
from SEO_writer import *
from quantum_CSD_compiler.UnitaryMat import *
import utilities_gen as ut


class Qubiter_to_AnyQasm(SEO_reader):
    """
    This abstract class is a child of SEO_reader. It reads an input English
    file and writes an AnyQasm file that is a translation of the input
    English file into the AnyQasm language. If the option
    write_qubiter_files is set to True, this class will also write new
    English and Picture files that are in 1-1 onto line correspondence with
    the output AnyQasm file.

    The input English file that is read can only have lines of the following
    types or else the program will abort with an error message:

    1. single qubit rotations (HAD2, SIGX, SIGY, SIGZ, ROTX, ROTY, ROTZ or
    ROTN with no controls)

    2. simple CNOTs (SIGX with a single True control). Call them c->t=(c,
    t) if c is the control and t the target. (c, t) must be allowed by
    'c_to_tars'.

    3. NOTA or PRINT lines. PRINT lines are commented out.

    If you have an English file that contains lines that are more
    complicated than this (because, for example, they contain rotations with
    one or more controls attached, or because a CNOT is not allowed
    according to 'c_to_tars'), you can use the expander classes
    CGateExpander, DiagUnitaryExpander, MultiplexorExpander,
    and ForbiddenCNotExpander to expand the circuit to an equivalent albeit
    longer circuit that satisfies constraints 1, 2, 3.

    This class can handle a chip with any number of qubits.

    This class halts execution if it encounters a CNOT that is disallowed
    according to the input 'c_to_tars'. 'c_to_tars' varies with chip. Some
    'c_to_tars's are listed in the files 'chip_couplings_...' found in same
    folder as this file. If c_to_tars = None, the class assumes any CNOT is
    possible.

    Footnote: Some AnyQasm's distinguish between quantum registers qreg and
    classical registers creg. Qubiter does not use cregs because it uses the
    classical memory of your Linux PC instead. AnyQasm has an intricate set
    of commands for measurements. Qubiter has a complete set of measurement
    commands too (see MEAS in Rosetta stone). The AnyQasm and Qubiter
    measurement commands can obviously be translated into each other. We
    leave that part of the translation to a future version of this class.

    Attributes
    ----------
    c_to_tars : dict[int, list[int]]
        a dictionary mapping j in range(num_bits) to a list, possibly empty,
        of the physically allowed targets of qubit j, when j is the control
        of a CNOT. If c_to_tars = None, the class assumes any CNOT is
        possible.
    file_prefix : str
    num_bits : int
    qasm_name : str
        the name of the qasm language, for example, IBMqasm. Used as ending
        of file name, between '_' and '.txt'
    qasm_out : _io.TextIOWrapper
        This output stream is used to write a qasm file based on the input
        English file.
    qbtr_wr : SEO_writer
        A SEO_writer object created iff write_qubiter_files is True.
    var_nums_list : list[int]
        list of all the distinct variable numbers encountered
    vname : str
        all variables in qasm file will be called vname + an int
    write_qubiter_files : bool
        The class always writes an AnyQasm text file based on the input
        English file that is read. Iff this is True, the class also writes
        English and Picture files in 1-1 line correspondence with the output
        AnyQasm file


    """
    def __init__(self, file_prefix, num_bits, qasm_name='',
                 c_to_tars=None, write_qubiter_files=False, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        qasm_name : str
        c_to_tars : dict[int, list[int]]|None
        write_qubiter_files : bool

        Returns
        -------

        """
        self.file_prefix = file_prefix
        self.num_bits = num_bits

        vman = PlaceholderManager(eval_all_vars=False)
        re = SEO_reader(file_prefix, num_bits, vars_manager=vman,
                        write_log=True)
        self.var_nums_list = re.vars_manager.var_nums_list

        self.qasm_name = qasm_name
        self.vname = 'rads_'
        self.c_to_tars = c_to_tars
        self.write_qubiter_files = write_qubiter_files

        self.qasm_out = open(file_prefix + '_' + qasm_name + '.txt', 'wt')

        self.qbtr_wr = None
        if write_qubiter_files:
            emb = CktEmbedder(num_bits, num_bits)
            out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
            self.qbtr_wr = SEO_writer(out_file_prefix, emb)

        self.write_prelude()

        SEO_reader.__init__(self, file_prefix, num_bits, **kwargs)

        self.write_ending()

        self.qasm_out.close()
        if write_qubiter_files:
            self.qbtr_wr.close_files()

    def write_prelude(self):
        """
        abstract function, writes AnyQasm's opening statements before calls
        to use_ methods for gates.

        Returns
        -------
        None

        """

        assert False

    def write_ending(self):
        """
        abstract function, writes AnyQasm's ending statements after calls to
        use_ methods for gates.

        Returns
        -------
        None

        """
        assert False

if __name__ == "__main__":
    def main():
        print(5)
    main()
