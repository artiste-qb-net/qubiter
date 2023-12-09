from qubiter.Controls import *
from qubiter.SEO_pre_reader import *
from qubiter.PlaceholderManager import *
from qubiter.LoopyPlaceholderManager import *
import qubiter.utilities_gen as utg
import os

import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np


class SEO_reader(SEO_pre_reader):
    """
    This class inherits from the class SEO_pre_reader. It's an abstract
    class because it has a bunch of use_ methods that must be overridden by
    a child class. This class reads each line of an English file, parses it,
    and sends the info obtained to a use_ method for further processing. One
    very important child of this class is SEO_simulator which uses each line
    of the English file to evolve by one further step a quantum state vector.

    See the docstring for the class SEO_writer for more info about English
    files.

    The main use of this "abstract" class is as an intermediate class that
    is a parent class to a child class. However, this class is not 100%
    abstract. An object of it can be created without getting an error
    message iff the write_log input parameter, which is False by default,
    is set to True instead. In that case, creating an object of the class
    will produce a log file about the English file that it reads. That log
    file will contain useful information about the English file, like its
    number of lines, its number of elementary ops, its number of CNOT
    operations (SIGX with one control), etc.

    Attributes
    ----------
    english_in : _io.TextIOWrapper
        file object for input text file that stores English description of
        circuit
    line_count : int
    loop_to_cur_rep : dict[int, int]
        a dictionary mapping loop number TO current repetition
    mcase_trols if outside if_m block, else a control specifying the current
        if_m case
    measured_bits : list(int)
        list of bits that have been measured with type 2 measurement and
        haven't been reset to ``|0>`` or ``|1>``
    num_cnots : int
    num_ops : int
    split_line : list[str]
    vars_manager : PlaceholderManager
        handles variables indicated by #int in the English file being read
    verbose : bool
    write_log : bool
    xfile_num : int
        You can ignore this if the English file has no loops. This number is
        -1 by default. If you change it to something non-negative, the class
        looks in the io_folder for a "Loop xfile" with a very specific name
        that mentions the xfile_num: file_prefix + '_' + str( num_qbits) +
        '_loop' + str(xfile_num) + ".py". If the class can't find that xfile
        , it will abort. Otherwise, it tries to exec() that xfile. A Loop
        xfile is a file that you write yourself from a template file called
        "Loop File" that is generated by the classes LoopFileGenerator and
        LoopyPlaceholder. The Loop File has the same name as the Loop xfile,
        except that the Loop xfile has that additional xfile_num before the
        ".py". Read docstrings and main() of classes LoopFileGenerator and
        LoopyPlaceholder for more info and examples illustrating how to use
        Loop Files and Loop xfiles.

    """

    def __init__(self, file_prefix, num_qbits, vars_manager=None,
                 verbose=False, write_log=False, xfile_num=-1):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_qbits : int
        vars_manager : PlaceholderManager
        verbose : bool
        write_log : bool
        xfile_num : int

        Returns
        -------

        """
        SEO_pre_reader.__init__(self, file_prefix, num_qbits)
        self.split_line = None
        self.vars_manager = vars_manager
        if vars_manager is None:
            self.vars_manager = PlaceholderManager()

        self.verbose = verbose
        self.write_log = write_log
        if write_log:
            self.vars_manager.eval_all_vars = False
        self.xfile_num = xfile_num
        if xfile_num >= 0:
            assert not self.vars_manager.var_num_to_rads, "we don't "\
                "allow a var_num_to_rads and a loop xfile simultaneously"
            assert not self.vars_manager.fun_name_to_fun, "we don't "\
                "allow a fun_name_to_fun and a loop xfile simultaneously"
            self.fill_history_lists_by_executing_loop_xfile()
        self.measured_bits = []
        self.mcase_trols = None

        self.english_in = open(utg.preface(
            file_prefix + '_' + str(num_qbits) + '_eng.txt'), 'rt')

        self.loop_to_cur_rep = {loop_num: 0 for
                                loop_num in self.loop_to_nreps.keys()}

        self.num_ops = 0
        self.num_cnots = 0
        self.line_count = 0

        while not self.english_in.closed:
            self.next_line()

        if write_log:
            self.do_log()

    @staticmethod
    def xed_file_prefix(file_prefix):
        """
        Xed file_prefix. Returns file_prefix + '_X1', assuming that '_X' +
        str(k) for some integer k is not already the ending of file_prefix.
        If it is, then the ending is changed to '_X' + str( k+1). Classes
        that use this are called "expanders" throughout Qubiter. The X
        stands for expanded.

        Parameters
        ----------
        file_prefix : str

        Returns
        -------
        str

        """
        k = file_prefix.rfind('_X')
        if k == -1:
            out_file_prefix = file_prefix + '_X1'
        else:
            out_file_prefix = file_prefix[:k+2] + \
                              str(int(file_prefix[k+2:])+1)
        return out_file_prefix

    def do_log(self):
        """
        Write a log file and print info on console too.

        Returns
        -------
        None

        """
        log = open(utg.preface(
            self.file_prefix + '_' + str(self.num_qbits) + '_log.txt'), 'wt')
        s = ''
        s += "Number of lines in file = " + str(self.tot_num_lines) + '\n'
        s += "Number of Elem. Ops = " + str(self.num_ops) + '\n'
        s += "Number of CNOTS (SIGX with single control) = " + \
            str(self.num_cnots) + '\n'

        s += "List of distinct variable numbers encountered "
        s += "(length=" + str(len(self.vars_manager.all_var_nums)) + ')=\n'
        s += str(self.vars_manager.all_var_nums) + "\n"

        s += "List of distinct function names encountered "
        s += "(length=" + str(len(self.vars_manager.all_fun_names)) + ')=\n'
        s += str(self.vars_manager.all_fun_names) + "\n"

        log.write(s)
        if self.verbose:
            print(s)

        log.close()

    def get_log_file_path(self, rel=False):
        """
        Returns path (relative if rel is True, absolute if rel is False) of
        log file.

        Parameters
        ----------
        rel : bool

        Returns
        -------
        str

        """
        rel_path = self.file_prefix + '_' + str(self.num_qbits) + '_log.txt'
        return rel_path if rel else utg.preface(rel_path)

    def print_log_file(self):
        """
        Prints log file.

        Returns
        -------
        None

        """
        path = self.get_log_file_path(rel=True)
        with open(utg.preface(path)) as f:
            print(f.read())

    def degs_str_to_rads(self, degs_str):
        """
        Wrapper for function of same name in PlaceholderManager.

        Parameters
        ----------
        degs_str : str

        Returns
        -------
        float | str

        """
        return self.vars_manager.degs_str_to_rads(degs_str, self.line_count)

    def fill_history_lists_by_executing_loop_xfile(self):
        """
        This method is called by the constructor of this class, iff the user
        enters a valid (non-negative) xfile number. Before using this
        method, user is expected to have generated a Loop File from an
        English file via classes LoopFileGenerator and
        LoopyPlaceholderManager, and created a Loop xfile by editing that
        Loop File. This class executes the Loop xfile to fill its history
        dictionaries (the ones that end in _hist).

        Returns
        -------
        None

        """
        assert self.xfile_num >= 0, \
            "user entered xfile number must be a non-negative int"
        all_var_nums = []
        all_fun_names = []
        var_num_to_hist = defaultdict(list)
        fun_name_to_hist = defaultdict(list)
        xfile_name = self.file_prefix + '_' + str(self.num_qbits) +\
            '_loop' + str(self.xfile_num) + '.py'
        try:
            loopx_in = open(utg.preface(xfile_name), 'rt')
        except IOError:
            print("Expected to find but didn't find a file named\n" +
                  xfile_name)
            exit()
        # var_dict = \
        #     {
        #         'all_var_nums': all_var_nums,
        #         'all_fun_names': all_fun_names,
        #         'var_num_to_hist': var_num_to_hist,
        #         'fun_name_to_hist': fun_name_to_hist
        #     }
        exec(loopx_in.read())
        self.vars_manager.all_var_nums = all_var_nums
        self.vars_manager.all_fun_names = all_fun_names
        self.vars_manager.var_num_to_hist = var_num_to_hist
        self.vars_manager.fun_name_to_hist = fun_name_to_hist

        # print('--------before resolving')
        # print(all_var_nums, all_fun_names)
        # print(var_num_to_hist)
        # print(fun_name_to_hist)
        self.vars_manager.resolve_all_histories()
        # print('--------after resolving')
        # print(var_num_to_hist)
        # print(fun_name_to_hist)

    def next_line(self):
        """
        Analyze the inputted line. Send info to use_ methods labelled by
        first four letters of line) for further use.

        Parameters
        ----------

        Returns
        -------
        None

        """
        line = self.english_in.readline()
        if not line or not line.strip():
            self.english_in.close()
            return

        self.split_line = line.split()
        line_name = self.split_line[0]
        self.num_ops += 1
        self.line_count += 1

        if line_name == "DIAG":
            # example:
            # DIAG IF 2:1 1:0  0T BY 30.0 10.5 11.0 83.1

            BY_pos = -1
            for k in range(len(self.split_line)):
                if self.split_line[k] == 'BY':
                    BY_pos = k
                    break
            trol_tokens = self.split_line[2: BY_pos]
            ang_tokens = self.split_line[BY_pos + 1: len(self.split_line)]
            trols = self.read_multi_controls(trol_tokens)
            rad_angles = [self.degs_str_to_rads(ang_tokens[k])
                          for k in range(len(ang_tokens))]
            self.use_DIAG(trols, rad_angles)

        elif line_name == "HAD2":
            # example:
            # HAD2 AT 1 IF 3F 2T

            tar_bit_pos = int(self.split_line[2])
            controls = self.read_TF_controls(self.split_line[4:])
            self.use_HAD2(tar_bit_pos, controls)

        elif line_name == "IF_M(":
            # don't count IF_M(<controls>){ as operation
            self.num_ops -= 1

            # example:
            # IF_M( 3F 2T ){
            self.mcase_trols = self.read_TF_controls(
                self.split_line[1:-1])
            for bit in self.mcase_trols.bit_pos:
                assert bit in self.measured_bits, \
                    "IF_M() argument mentions a qubit that" \
                    " hasn't been measured yet"
            self.use_IF_M_beg(self.mcase_trols)

        elif line_name == "}IF_M":
            # don't count }IF_M as operation
            self.num_ops -= 1
            self.mcase_trols = None
            self.use_IF_M_end()

        elif line_name == "LOOP":
            # don't count LOOP as operation
            self.num_ops -= 1

            # example:
            # LOOP 5 NREPS= 2

            loop_num = int(self.split_line[1])
            nreps = int(self.split_line[3])
            self.use_LOOP(loop_num, nreps)

        elif line_name == "MEAS":
            # example:
            # MEAS  0  AT  5
            # MEAS  1  AT  5
            # MEAS  2  AT  5

            kind = int(self.split_line[1])
            tar_bit_pos = int(self.split_line[3])
            if kind == 2:
                # don't measure same bit twice
                assert tar_bit_pos not in self.measured_bits,\
                    "attempting to measure (kind=2) same qubit twice"
                self.measured_bits.append(tar_bit_pos)
            self.use_MEAS(tar_bit_pos, kind)

        elif line_name == "MP_Y":
            # example:
            # MP_Y AT 3 IF 2:1 1:0  0T BY 30.0 10.5 11.0 83.1

            tar_bit_pos = int(self.split_line[2])
            BY_pos = -1
            for k in range(len(self.split_line)):
                if self.split_line[k] == 'BY':
                    BY_pos = k
                    break
            trol_tokens = self.split_line[4: BY_pos]
            ang_tokens = self.split_line[BY_pos + 1: len(self.split_line)]
            trols = self.read_multi_controls(trol_tokens)
            rad_angles = [self.degs_str_to_rads(ang_tokens[k])
                for k in range(len(ang_tokens))]
            self.use_MP_Y(tar_bit_pos, trols, rad_angles)

        elif line_name == "NEXT":
            # don't count NEXT as operation
            self.num_ops -= 1

            # example:
            # NEXT 5

            loop_num = int(self.split_line[1])
            self.use_NEXT(loop_num)

        elif line_name == 'NOTA':
            # don't count NOTA as operation
            self.num_ops -= 1

            # example:
            # NOTA  "I love you Mary."

            self.use_NOTA(line[4:].strip())

        elif line_name == "PHAS":
            # example:
            # PHAS 42.7 AT 1 IF 3F 2T

            angle_rads = self.degs_str_to_rads(self.split_line[1])
            tar_bit_pos = int(self.split_line[3])
            controls = self.read_TF_controls(self.split_line[5:])
            self.use_PHAS(angle_rads, tar_bit_pos, controls)

        elif line_name == "P0PH":
            self.read_P_phase_factor(0)
        elif line_name == "P1PH":
            self.read_P_phase_factor(1)

        elif line_name == "PRINT":
            # don't count PRINT as operation
            self.num_ops -= 1

            # example:
            # PRINT V1
            assert len(self.split_line) == 2, \
                "PRINT line must contain style str"
            self.use_PRINT(self.split_line[1], self.line_count)

        elif line_name == "ROTX":
            self.read_ROT(1)
        elif line_name == "ROTY":
            self.read_ROT(2)
        elif line_name == "ROTZ":
            self.read_ROT(3)
        elif line_name == "ROTN":
            # example:
            # ROTN 42.7 30.2 78.5 AT 1 IF 3F 2T

            angle_x_rads = self.degs_str_to_rads(self.split_line[1])
            angle_y_rads = self.degs_str_to_rads(self.split_line[2])
            angle_z_rads = self.degs_str_to_rads(self.split_line[3])
            tar_bit_pos = int(self.split_line[5])
            controls = self.read_TF_controls(self.split_line[7:])
            self.use_ROTN(angle_x_rads, angle_y_rads, angle_z_rads,
                             tar_bit_pos, controls)
        elif line_name == "SIGX":
            self.read_SIG(1)
        elif line_name == "SIGY":
            self.read_SIG(2)
        elif line_name == "SIGZ":
            self.read_SIG(3)

        elif line_name == "SWAP":
            # example:
            # SWAP 1 0 IF 3F 2T

            bit1 = int(self.split_line[1])
            bit2 = int(self.split_line[2])
            controls = self.read_TF_controls(self.split_line[4:])
            self.use_SWAP(bit1, bit2, controls)
        elif line_name == "SWAY":
            # example:
            # SWAY 0 1 BY 25.1 42.7 IF 3F 2T

            bit1 = int(self.split_line[1])
            bit2 = int(self.split_line[2])
            rads_list = [self.degs_str_to_rads(self.split_line[k])
                                               for k in range(4, 6)]
            controls = self.read_TF_controls(self.split_line[7:])
            self.use_SWAY(bit1, bit2, controls, rads_list)
        elif line_name == "U_2_":
            # example:
            # U_2_  25.1 42.7 30.2 78.5 AT 1 IF 3F 2T

            rads0 = self.degs_str_to_rads(self.split_line[1])
            rads1 = self.degs_str_to_rads(self.split_line[2])
            rads2 = self.degs_str_to_rads(self.split_line[3])
            rads3 = self.degs_str_to_rads(self.split_line[4])
            tar_bit_pos = int(self.split_line[6])
            controls = self.read_TF_controls(self.split_line[8:])
            self.use_U_2_(rads0, rads1, rads2, rads3,
                             tar_bit_pos, controls)
        else:
            assert False, \
                "reading an unsupported line kind: " + line_name

        self.finalize_next_line()

    def finalize_next_line(self):
        """
        Useful for intercepting the end of each call to next_line().

        Returns
        -------

        """
        if self.verbose:
            print('line_num, operation =', self.line_count, self.num_ops)

    def read_multi_controls(self, tokens, allow_only_TF=False):
        """
        Given a list of tokens of the form:
        * an int followed by either T or F,
        * int, colon, int,

        construct a control out of it.

        Parameters
        ----------
        tokens : list[str]
        allow_only_TF : bool

        Returns
        -------
        Controls

        """
        # safe to use when no "IF"
        # when no "IF", will return controls with _numControls=0
        controls = Controls(self.num_qbits)
        if tokens:
            for t in tokens:
                t_end = t[-1]
                if allow_only_TF:
                    assert t_end in ['T', 'F']
                if t_end == 'T':
                    controls.set_control(int(t[:-1]), True)
                elif t_end == 'F':
                    controls.set_control(int(t[:-1]), False)
                else:
                    k1, k2 = t.split(':')
                    controls.set_control(int(k1), int(k2))
            controls.refresh_lists()
        return controls

    def read_TF_controls(self, tokens):
        """
        Same as read_multi_controls() but only allows T/F kind controls.

        Parameters
        ----------
        tokens : list[str]

        Returns
        -------
        Controls

        """
        return self.read_multi_controls(tokens, allow_only_TF=True)

    def read_P_phase_factor(self, projection_bit):
        """
        Collect useful info from P0PH or P1PH split_line and forward it to
        use_ method.

        Parameters
        ----------
        projection_bit : int
        Returns
        -------
        None
        """
        # example:
        # P0PH 42.7 AT 1 IF 3F 2T
        # P1PH 42.7 AT 1 IF 3F 2T

        angle_rads = self.degs_str_to_rads(self.split_line[1])
        tar_bit_pos = int(self.split_line[3])
        controls = self.read_TF_controls(self.split_line[5:])
        assert projection_bit in [0, 1]
        self.use_P_PH(projection_bit,
                          angle_rads, tar_bit_pos, controls)

    def read_ROT(self, axis):
        """
        Collect useful info from ROTX, ROTY, or ROTZ split_line and forward
        it to use_ method.

        Parameters
        ----------
        axis : int

        Returns
        -------
        None

        """
        # example:
        # ROTX 42.7 AT 1 IF 3F 2T
        # ROTY 42.7 AT 1 IF 3F 2T
        # ROTZ 42.7 AT 1 IF 3F 2T

        angle_rads = self.degs_str_to_rads(self.split_line[1])
        tar_bit_pos = int(self.split_line[3])
        controls = self.read_TF_controls(self.split_line[5:])
        self.use_ROTA(axis, angle_rads, tar_bit_pos, controls)

    def read_SIG(self, axis):
        """
        Collect useful info from SIGX, SIGY, or SIGZ split_line and forward
        it to use_ method.

        Parameters
        ----------
        axis : int

        Returns
        -------
        None

        """
        # example:
        # SIGX AT 1 IF 3F 2T
        # SIGY AT 1 IF 3F 2T
        # SIGZ AT 1 IF 3F 2T

        tar_bit_pos = int(self.split_line[2])
        controls = self.read_TF_controls(self.split_line[4:])
        assert axis in [1, 2, 3]
        if axis == 1 and len(controls.bit_pos) == 1:
            self.num_cnots += 1

        self.use_SIG(axis, tar_bit_pos, controls)

    def use_DIAG(self, trols, rad_angles):
        """
        Abstract use_ method that must be overridden by child class. 

        Parameters
        ----------
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'DIAG not used'

    def use_HAD2(self, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'HAD2 not used'

    def use_IF_M_beg(self, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'IF_M(){ not used'

    def use_IF_M_end(self):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, '}IF_M not used'

    def use_LOOP(self, loop_num, nreps):
        """
        Don't override this unless you know what you are doing and have very
        good reasons. It has been carefully set up to deal properly with
        embedded loops.

        Parameters
        ----------
        loop_num : int
        nreps : int

        Returns
        -------
        None

        """
        pass

    def use_MEAS(self, tar_bit_pos, kind):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        kind : int
        tar_bit_pos : int

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'MEAS not used'

    def use_MP_Y(self, tar_bit_pos, trols, rad_angles):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'MP_Y not used'

    def use_NEXT(self, loop_num):
        """
        Don't override this unless you know what you are doing and have very
        good reasons. It has been carefully set up to deal properly with
        embedded loops.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        cur_rep = self.loop_to_cur_rep[loop_num]
        if cur_rep < self.loop_to_nreps[loop_num]-1:
            self.english_in.seek(self.loop_to_start_offset[loop_num])
            self.line_count = self.loop_to_start_line[loop_num] - 1
            self.loop_to_cur_rep[loop_num] += 1

        else:  # cur_rep = self.loop_to_nreps[loop_num]-1
            self.loop_to_cur_rep[loop_num] = 0

    def use_NOTA(self, bla_str):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'NOTA not used'

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'PHAS not used'

    def use_P_PH(self, projection_bit, angle_rads, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        projection_bit : int
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'P0PH or P1PH not used'

    def use_PRINT(self, style, line_num):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        style : str
        line_num : int

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'PRINT not used'

    def use_ROTA(self, axis, angle_rads, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        axis : int
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'ROTX, ROTY or ROTZ not used'

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        angle_x_rads : float
        angle_y_rads : float
        angle_z_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'ROTN not used'

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        axis : int
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'SIGX, SIGY or SIGZ not used'

    def use_SWAP(self, bit1, bit2, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'SWAP not used'

    def use_SWAY(self, bit1, bit2, controls, rads_list):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls
        rads_list: list[float | str]

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'SWAY not used'

    def use_U_2_(self, rads0, rads1, rads2, rads3,
                tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        rads0 : float
        rads1 : float
        rads2 : float
        rads3 : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if self.write_log:
            return
        assert False, 'U_2_ not used'


if __name__ == "__main__":
    def main():
        file_prefix = 'expansions_examples_X1'
        num_qbits = 3
        SEO_reader(file_prefix, num_qbits, write_log=True)
    main()
