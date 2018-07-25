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

    3. NOTA lines

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
    verbose : bool
    write_qubiter_files : bool
        The class always writes an AnyQasm text file based on the input
        English file that is read. Iff this is True, the class also writes
        English and Picture files in 1-1 line correspondence with the output
        AnyQasm file


    """
    def __init__(self, file_prefix, qasm_name,
                 num_bits, c_to_tars=None, verbose=False,
                 write_qubiter_files=False, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        qasm_name : str
        num_bits : int
        c_to_tars : dict[int, list[int]]|None
        verbose : bool
        write_qubiter_files : bool
        kwargs : dict[]

        Returns
        -------
        None

        """
        self.file_prefix = file_prefix
        self.qasm_name = qasm_name
        self.num_bits = num_bits
        self.c_to_tars = c_to_tars
        self.verbose = verbose
        self.write_qubiter_files = write_qubiter_files

        self.qasm_out = open(file_prefix + '_' + qasm_name + '.txt', 'wt')

        self.qbtr_wr = None
        if write_qubiter_files:
            emb = CktEmbedder(num_bits, num_bits)
            out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
            self.qbtr_wr = SEO_writer(out_file_prefix, emb, **kwargs)

        self.write_prelude()

        SEO_reader.__init__(self, file_prefix, num_bits, verbose)

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

    def use_DIAG(self, trols, rad_angles):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        assert False, "No DIAG lines allowed"

    def use_IF_M_beg(self, controls):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        controls : Controls

        Returns
        -------
        None

        """
        assert False, "No IF_M{ lines allowed"

    def use_IF_M_end(self):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------

        Returns
        -------
        None

        """
        assert False, "No }IF_M lines allowed"

    def use_LOOP(self, loop_num, reps):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        loop_num : int
        reps : int

        Returns
        -------
        None

        """
        assert False, "No LOOP lines allowed"

    def use_MEAS(self, tar_bit_pos, kind):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        kind : int
        tar_bit_pos : int

        Returns
        -------
        None

        """
        assert False, "No MEAS lines allowed"

    def use_MP_Y(self, tar_bit_pos, trols, rad_angles):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        assert False, "No MP_Y lines allowed"

    def use_NEXT(self, loop_num):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        assert False, "No NEXT lines allowed"

    def use_P_PH(self, projection_bit, angle_degs, tar_bit_pos, controls):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        projection_bit : int
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False, "No P0PH or P1PH lines allowed"

    def use_PRINT(self, style, line_num):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        style : str
        line_num : int

        Returns
        -------
        None

        """
        assert False, "No PRINT lines allowed"

    def use_SWAP(self, bit1, bit2, controls):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False, "No SWAP lines allowed"
