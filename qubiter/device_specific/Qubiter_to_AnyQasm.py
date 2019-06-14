from qubiter.Controls import *
from qubiter.SEO_reader import *
from qubiter.SEO_writer import *
from qubiter.UnitaryMat import *
import qubiter.utilities_gen as ut


class Qubiter_to_AnyQasm(SEO_reader):
    """
    This abstract class is a child of SEO_reader. It reads an input English
    file and writes an AnyQasm file that is a translation of the input
    English file into the AnyQasm language. If the flag write_qubiter_files
    is set to True, this class will also write new English and Picture files
    that are in 1-1 onto line correspondence with the output AnyQasm file.

    Footnote: Some AnyQasm's distinguish between quantum registers qreg and
    classical registers creg. Qubiter does not use cregs because it uses the
    classical memory of your Linux PC instead. AnyQasm has an intricate set
    of commands for measurements. Qubiter has a complete set of measurement
    commands too (see MEAS in Rosetta stone). The AnyQasm and Qubiter
    measurement commands can obviously be translated into each other. We
    leave that part of the translation to a future version of this class.

    This class can run in either a strict or a non-strict mode depending on
    the flag `strict_mode`, which equals False in default mode. In the
    strict mode, the set of gates allowed is constrained to a small but
    universal set that is specified below, and that is allowed in any target
    qasm. In the non-strict mode, more gates are allowed that depend on
    specific target qasm. In the strict mode, the program will end if you
    try to use gates that are not allowed. In the non-strict mode,
    the program will end if you try to use gates for a target language that
    have not been implemented yet in the Qubiter class targeting that
    language, often because the target language doesn't support those gates.
    
    Will refer to target qasm as AnyQasm or aqasm

    Next we give a description of the strict_mode:

    In the strict mode, the input English file that is read can only have
    lines of the following types or else the program will abort with an
    error message:

    1. single qubit rotations (HAD2, SIGX, SIGY, SIGZ, ROTX, ROTY,
    ROTZ or ROTN with no controls)

    2. simple CNOTs (SIGX with a single True control). Call them c->t=(
    c, t) if c is the control and t the target. (c, t) must be allowed
    by 'c_to_tars'.

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


    Attributes
    ----------
    all_fun_names : list[str]
        a list of all the distinct function names encountered in circuit
    all_var_nums : list[int]
        a list of all distinct numbers of the variables encountered in circuit
    aqasm_name : str
        the name of the aqasm language, for example, IBMqasm. Used as ending
        of file name, between '_' and '.txt'
    aqasm_path : str
        path to aqasm file
    aqasm_out : _io.TextIOWrapper
        This output stream is used to write an aqasm file based on the input
        English file.
    c_to_tars : dict[int, list[int]]
        a dictionary mapping j in range(num_bits) to a list, possibly empty,
        of the physically allowed targets of qubit j, when j is the control
        of a CNOT. If c_to_tars = None, the class assumes any CNOT is
        possible.
    file_prefix : str
    num_bits : int
    qbtr_wr : SEO_writer
        A SEO_writer object created iff write_qubiter_files is True.
    strict_mode : bool
    vprefix : str
        all variables in aqasm file will be called vprefix + an int
    write_qubiter_files : bool
        The class always writes an AnyQasm text file based on the input
        English file that is read. Iff this is True, the class also writes
        English and Picture files in 1-1 line correspondence with the output
        AnyQasm file


    """
    def __init__(self, file_prefix, num_bits, aqasm_name='',
            strict_mode=False, c_to_tars=None, write_qubiter_files=False,
                 vars_manager=None, aqasm_ftype='txt',
                 prelude_str=None, ending_str=None, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        aqasm_name : str
        strict_mode : bool
        c_to_tars : dict[int, list[int]]|None
        write_qubiter_files : bool
        vars_manager : PlaceholderManager
        aqasm_ftype : str
            file type of output aqasm file. If this equals 'txt', name of
            aqasm file will end in '.txt'
        prelude_str : str | None
            string to write as prelude to aqasm file. If None, then the
            override method of self.write_prelude() is called
        ending_str : str | None
            string to write as ending to aqasm file. If None, then the
            override method of self.write_ending() is called

        Returns
        -------

        """
        self.file_prefix = file_prefix
        self.num_bits = num_bits

        vman = PlaceholderManager(eval_all_vars=False)
        rdr = SEO_reader(file_prefix, num_bits, vars_manager=vman,
                        write_log=True)
        self.all_var_nums = rdr.vars_manager.all_var_nums
        self.all_fun_names = rdr.vars_manager.all_fun_names

        self.aqasm_name = aqasm_name
        self.strict_mode = strict_mode
        self.vprefix = 'rads'
        self.c_to_tars = c_to_tars
        self.write_qubiter_files = write_qubiter_files

        self.aqasm_path = file_prefix +\
                          '_' + aqasm_name + '.' + aqasm_ftype
        self.aqasm_out = open(self.aqasm_path, 'wt')

        self.qbtr_wr = None
        if write_qubiter_files:
            emb = CktEmbedder(num_bits, num_bits)
            out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
            self.qbtr_wr = SEO_writer(out_file_prefix, emb)

        if prelude_str is not None:
            self.write(prelude_str)
        else:
            self.write_prelude()

        vman1 = PlaceholderManager(eval_all_vars=False)
        SEO_reader.__init__(self, file_prefix, num_bits,
                            vars_manager=vman1, **kwargs)

        if ending_str is not None:
            self.write(ending_str)
        else:
            self.write_ending()

        self.aqasm_out.close()
        if write_qubiter_files:
            self.qbtr_wr.close_files()

    def write(self, s):
        """
        Writes string s to aqasm and qubiter out files

        Parameters
        ----------
        s : str

        Returns
        -------
        None

        """
        self.aqasm_out.write(s + '\n')

        if self.write_qubiter_files:
            lines = s.split('\n')
            for line in lines:
                self.qbtr_wr.write_NOTA(line)

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

    def new_var_name(self, var_name, coda='', strict=False):
        """
        Starts by asserting that var_name is a legal variable name.

        If var_name is not functional, this method replaces # in var_name by
        self.vprefix and adds coda to end of string.  For example,
        if self.vprefix='rads' and var_name='-#2*.5", then output is
        '-rads2*.5' + coda

        If var_name is functional, this method replaces each # in var_name
        by self.vprefix, adds commas and parenthesis, and adds coda to end
        of string. For example, if self.vprefix='rads' and
        var_name='-fun#1#2', then output is '-fun(rads1, rads2)' + coda

        The above applies only if strict=False. In the strict mode, only an
        empty coda is allowed for functional placeholders. For
        non-functional placeholders, if var_name contains an *, then the str
        after the * and the coda are merged using eval().

        Parameters
        ----------
        var_name : str
        coda : str
        strict : bool

        Returns
        -------
        str

        """
        assert PlaceholderManager.is_legal_var_name(var_name)
        if not PlaceholderManager.is_functional_var(var_name):
            new_coda = coda
            if coda:
                assert len(coda) > 1, "illegal coda: " + coda
            star_pos = var_name.find("*")
            if strict and star_pos != -1 and coda:
                assert coda[0] == '*', "A coda must start with * " +\
                    "in strict mode. Got coda: " + coda
                fac1 = var_name[star_pos+1:]
                fac2 = coda[1:]
                fac12 = fac1 + '*' + fac2
                try:
                    new_coda = '*' + str(eval(fac12))
                except:
                    assert False, 'cannot eval "' +\
                                  fac12 + '" to merge "' +\
                            var_name + '" and "' + coda + '"'
            end_pos = len(var_name)
            if strict and star_pos != -1:
                end_pos = star_pos

            if var_name[0] == "#":
                return self.vprefix + var_name[1:end_pos] + new_coda
            else:  # starts with -#
                return "-" + self.vprefix + var_name[2:end_pos] + new_coda

        else:
            if strict:
                assert not coda, "functional placeholders cannot " +\
                    'have scaling factors in strict mode'
            first_hash_pos = var_name.find('#')
            nums_strings = var_name[first_hash_pos+1:].split('#')
            arg_str = '('
            for num_str in nums_strings:
                arg_str += self.vprefix + num_str + ', '
            arg_str = arg_str[:-2] + ')'
            return var_name[:first_hash_pos] + arg_str + coda

    def print_aqasm_file(self):
        """
        Prints aqasm file created by constructor.

        Returns
        -------

        """
        with open(self.aqasm_path) as f:
            print(f.read())

if __name__ == "__main__":
    def main():
        print(5)
    main()
