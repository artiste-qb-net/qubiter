from Controls import *
from SEO_pre_reader import *
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


    Attributes
    ----------
    num_ops : int
        number of operations. Lines inside a loop with 'reps' repetitions
        will count as 'reps' operations
    loop_to_cur_rep : dict[int, int]
        a dictionary mapping loop number TO current repetition
    just_jumped : bool
        flag used to alert when loop jumps from NEXT to LOOP
    line_count : int

    english_in : _io.TextIOWrapper
        file object for input text file that stores English description of
        circuit
    file_prefix : str
        beginning of the name of English file being scanned
    loop_to_start_line : dict[int, int]
        a dictionary mapping loop number TO loop line + 1
    loop_to_start_offset : dict[int, int]
        a dictionary mapping loop number TO offset of loop's start
    loop_to_reps : dict[int, int]
        a dictionary mapping loop number TO total number of repetitions of
        loop
    loop_queue : list[int]
        a queue of loops labelled by their id number
    num_bits : int
        number of qubits in whole circuit
    tot_num_lines : int
        number of lines in English file
    split_line : list[str]
        storage space for a list of strings obtained by splitting a line


    """
    
    def __init__(self, file_prefix, num_bits):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int

        Returns
        -------

        """
        SEO_pre_reader.__init__(self, file_prefix, num_bits)

        self.english_in = open(
            file_prefix + '_' + str(num_bits) + '_eng.txt', 'rt')

        self.loop_to_cur_rep = {loop_num: 0 for
                loop_num in self.loop_to_reps.keys()}

        self.num_ops = 0
        self.line_count = 0
        self.just_jumped = False

        while not self.english_in.closed:
            self.next_line()

        self.write_log()

    def write_log(self):
        """
        Write a log file and print info on console too.

        Returns
        -------
        None

        """
        log_out = open(
            self.file_prefix + '_' + str(self.num_bits) + '_log.txt', 'wt')

        s = "Number of lines in file = " + str(self.tot_num_lines) + '\n'
        log_out.write(s)
        print(s)

        s = "Number of Elem. Ops = " + str(self.num_ops) + '\n'
        log_out.write(s)
        print(s)

        log_out.close()

    def next_line(self):
        """
        Analyze the inputted line. Send info to abstract use_ methods (
        labelled by first four letters of line) for further use.


        Parameters
        ----------

        Returns
        -------
        None

        """
        line = self.english_in.readline()
        if not line:
            self.english_in.close()
            return

        self.split_line = line.split()
        line_name = self.split_line[0]
        self.num_ops += 1
        self.line_count += 1
        self.just_jumped = False

        if line_name == "LOOP":
            # don't count LOOP as operation
            self.num_ops -= 1

            # example:
            # LOOP 5 REPS: 2

            loop_num = int(self.split_line[1])
            reps = int(self.split_line[3])
            self.use_LOOP(loop_num, reps)

        elif line_name == "NEXT":
            # don't count NEXT as operation
            self.num_ops -= 1

            # example:
            # NEXT 5

            loop_num = int(self.split_line[1])
            self.use_NEXT(loop_num)

        elif line_name == 'NOTA':
            # don't count NOTA as operation

            # example:
            # NOTA  "I love you Mary."
            self.num_ops -= 1
            self.use_NOTA(line[4:])

        elif line_name == "SWAP":
            # example:
            # SWAP 1 0 IF 3F 2T

            bit1 = int(self.split_line[1])
            bit2 = int(self.split_line[2])
            controls = self.read_TF_controls(self.split_line[4:])
            self.use_SWAP(bit1, bit2, controls)

        elif line_name == "MEAS":
            # example:
            # MEAS  0  AT  5
            # MEAS  1  AT  5
            # MEAS  2  AT  5

            kind = int(self.split_line[1])
            tar_bit_pos = int(self.split_line[3])
            self.use_MEAS(kind, tar_bit_pos)

        elif line_name == "PHAS":
            # example:
            # PHAS 42.7 AT 3 IF 3F 2T

            angle_degs = float(self.split_line[1])
            tar_bit_pos = int(self.split_line[3])
            controls = self.read_TF_controls(self.split_line[5:])
            self.use_PHAS(angle_degs, tar_bit_pos, controls)

        elif line_name == "P0PH":
            self.read_P_phase_factor(0)
        elif line_name == "P1PH":
            self.read_P_phase_factor(1)
        elif line_name == "SIGX":
            self.read_SIG(1)
        elif line_name == "SIGY":
            self.read_SIG(2)
        elif line_name == "SIGZ":
            self.read_SIG(3)
        elif line_name == "HAD2":
            # example:
            # HAD2 AT 1 IF 3F 2T

            tar_bit_pos = int(self.split_line[2])
            controls = self.read_TF_controls(self.split_line[4:])
            self.use_HAD2(tar_bit_pos, controls)

        elif line_name == "ROTX":
            self.read_ROT(1)
        elif line_name == "ROTY":
            self.read_ROT(2)
        elif line_name == "ROTZ":
            self.read_ROT(3)
        elif line_name == "ROTN":
            # example:
            # ROTN 42.7 30.2 78.5 AT 3 IF 3F 2T

            angle_x_degs = float(self.split_line[1])
            angle_y_degs = float(self.split_line[2])
            angle_z_degs = float(self.split_line[3])
            tar_bit_pos = int(self.split_line[5])
            controls = self.read_TF_controls(self.split_line[7:])
            self.use_ROTN(angle_x_degs, angle_y_degs, angle_z_degs,
                             tar_bit_pos, controls)

        elif line_name == "MP_Y":
            # example:
            # MP_Y AT 3 IF 2(1 1(0  0T BY 30.0 10.5 11.0 83.1

            tar_bit_pos = int(self.split_line[2])
            BY_pos = -1
            for k in range(len(self.split_line)):
                if self.split_line[k] == 'BY':
                    BY_pos = k
                    break
            trol_tokens = self.split_line[4: BY_pos]
            ang_tokens = self.split_line[BY_pos + 1: len(self.split_line)]
            trols = self.read_multi_controls(trol_tokens)
            rad_angles = [float(ang_tokens[k])*np.pi/180
                for k in range(len(ang_tokens))]
            self.use_MP_Y(tar_bit_pos, trols, rad_angles)
        else:
            assert False, \
                "reading an unsupported line kind: " + line_name

        self.finalize_next_line()
        if self.just_jumped:
            self.line_count = self.loop_to_start_line[loop_num] - 1

    def finalize_next_line(self):
        """
        Useful for intercepting the end of each call to read_line()

        Returns
        -------

        """
        print('line_num, operation =', self.line_count, self.num_ops)

    def read_multi_controls(self, tokens, allow_only_TF=False):
        """
        Given a list of tokens of the form:
        * an int followed by either T or F,
        * int, opening parenthesis, int,

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
        controls = Controls(self.num_bits)
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
                    k1, k2 = t.split('(')
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
        abstract use_ method.

        Parameters
        ----------
        projection_bit : int

        Returns
        -------
        None

        """
        # example:
        # P0PH 42.7 AT 3 IF 2T
        # P1PH 42.7 AT 3 IF 2T

        angle_degs = float(self.split_line[1])
        tar_bit_pos = int(self.split_line[3])
        controls = self.read_TF_controls(self.split_line[5:])
        assert projection_bit in [0, 1]
        self.use_P_PH(projection_bit,
                          angle_degs, tar_bit_pos, controls)

    def read_SIG(self, direction):
        """
        Collect useful info from SIGX, SIGY, or SIGZ split_line and forward
        it to abstract use_ method.

        Parameters
        ----------
        direction : int

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
        assert direction in [1, 2, 3]
        self.use_SIG(direction, tar_bit_pos, controls)

    def read_ROT(self, direction):
        """
        Collect useful info from ROTX, ROTY, or ROTZ split_line and forward
        it to abstract use_ method.

        Parameters
        ----------
        direction : int

        Returns
        -------
        None

        """
        # example:
        # ROTX 42.7 AT 3 IF 3F 2T
        # ROTY 42.7 AT 3 IF 3F 2T
        # ROTZ 42.7 AT 3 IF 3F 2T

        angle_degs = float(self.split_line[1])
        tar_bit_pos = int(self.split_line[3])
        controls = self.read_TF_controls(self.split_line[5:])
        self.use_ROT(direction,
                         angle_degs, tar_bit_pos, controls)

    def use_LOOP(self, loop_num, reps):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        loop_num : int
        reps : int

        Returns
        -------
        None

        """
        pass

    def use_NEXT(self, loop_num):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        if self.loop_to_cur_rep[loop_num] < self.loop_to_reps[loop_num]-1:
            self.english_in.seek(self.loop_to_start_offset[loop_num])
            self.just_jumped = True
            self.loop_to_cur_rep[loop_num] += 1
        else:
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
        assert False

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
        assert False

    def use_MEAS(self, kind, tar_bit_pos):
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
        assert False

    def use_PHAS(self, angle_degs, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False

    def use_P_PH(self, projection_bit, angle_degs, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

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
        assert False

    def use_SIG(self, direction, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        direction : int
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False

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
        assert False

    def use_ROT(self, direction, angle_degs, tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        direction : int
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        Abstract use_ method that must be overridden by child class.

        Parameters
        ----------
        angle_x_degs : float
        angle_y_degs : float
        angle_z_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False

    def use_MP_Y(self, tar_bit_pos, trols, rad_angles):
        """

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------

        """
        assert False

if __name__ == "__main__":
    print(5)
