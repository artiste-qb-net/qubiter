from CktEmbedder import *
from Controls import *
from OneBitGates import *
import re


class SEO_writer:
    """
    The constructor of this class opens an English and a Picture file. Every
    other method of the class writes, each time it is called, a single line
    in each of those 2 files.

    Note SEO stands for Sequence of Elementary Operations.

    So what are English and Picture files?

    We use 3 types of files to characterize a single quantum circuit (in the
    gate model): (1) an English file (2) a Picture file (3) a Log file.

    Log files are written by class SEO_reader, whereas English and Picture
    files are written by this class, SEO_writer.

    A log file just contains useful information like the number of lines of
    the English file (same as that of Picture file) and their number of
    operations.

    The lines of an English and Picture file are in 1-1 correspondence,
    each line representing a single operation (e.g., a multi-controlled one
    qubit gate or a multi-controlled 2 qubit swap).

    The English file (resp., Picture file) contains complete (resp.,
    partial, in the form of an ASCII picture) info about the operation
    specified by each line.

    In English and Picture files, time flows downward.

    The class SEO_writer gives the bool option zero_bit_first. When this
    argument is set to True (resp., False), the Picture file shows the zero
    qubit first (resp., last), and the remaining qubits in consecutive
    order. Picture files written with zero bit first (resp., last) are
    labelled prefix + '_ZFpict.text' (reps., prefix + '_ZLpict.txt'). The
    zero_bit_first choice does not affect the English file.

    See the following and earlier arXiv papers by R.R.Tucci for more info on
    English and Picture files.

        http://arxiv.org/abs/1004.2205
        "Quibbs, a Code Generator for Quantum Gibbs Sampling"

    The following images from the above paper are stored in the same folder
    as this python file.

        ZLpic_file_examples.png
        eng_file_examples.png

    These 2 images give examples of lines in analytic/Picture/English
    formats. The Picture file png follows the ZL convention.

    3 kinds (called 0, 1, 2) of measurements MEAS are allowed. A type 0
    measurement inserts a projector |0><0| = n = P_0 at the target bit. A
    type 1 measurement inserts a projector |1><1| = nbar = P_1 at the target
    bit. A type 2 measurement stores a copy of the state vector after |0><0|
    has been applied, and another copy after |1><1| has been applied.

    Attributes
    ----------
    emb : CktEmbedder
    english_out : _io.TextIOWrapper
        file object for output text file that stores English description of
        circuit
    picture_out : _io.TextIOWrapper
        file object for output text file that stores ASCII Picture
        description of circuit
    file_prefix : str
        beginning of the name of both English and Picture files
    line_counter : int
    zero_bit_first : bool

    """

    def __init__(self, file_prefix, emb, zero_bit_first=False,
                english_out=None, picture_out=None):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        zero_bit_first : bool
        english_out : _io.TextIOWrapper
        picture_out : _io.TextIOWrapper
        emb : CktEmbedder

        Returns
        -------

        """

        self.line_counter = 0
        self.file_prefix = file_prefix
        self.emb = emb
        self.zero_bit_first = zero_bit_first

        if english_out is None:
            self.english_out = open(
                file_prefix + '_' +
                str(self.emb.num_bits_aft) + '_eng.txt', 'wt')
        else:
            self.english_out = english_out

        if picture_out is None:
            self.picture_out = open(
                file_prefix + '_' + str(self.emb.num_bits_aft) +
                ('_ZF' if zero_bit_first else '_ZL') + 'pic.txt', 'wt')
        else:
            self.picture_out = picture_out

    def close_files(self):
        """
        Closes English and Picture files that were opened by the constructor.

        Returns
        -------
        None

        """
        self.english_out.close()
        self.picture_out.close()

    def write_pic_line(self, pic_line):
        """
        Writes a line in the Picture file using either the ZF or ZL
        conventions. pic_line is originally written in ZL format, so this
        method does nothing to pic_line if ZL option is chosen but reverses
        order of gates if ZF option chosen.

        Parameters
        ----------
        pic_line : str

        Returns
        -------
        None

        """
        # example:
        # Ry--R---<--->---Rz
        if self.zero_bit_first:
            pic_line = pic_line.strip()
            # delimiter is -- or --- or 2 spaces or 3 spaces
            nodes = re.split('-{2,3}| {2,3}', pic_line)
            dash_or_space = [pic_line[k*4 + 3] for k in range(len(nodes)-1)]
            nodes = list(reversed(nodes))
            dash_or_space = list(reversed(dash_or_space))
            for kk in range(len(nodes)):
                if nodes[kk] == '<':
                    nodes[kk] = '>'
                elif nodes[kk] == '>':
                    nodes[kk] = '<'
            new_line = ''
            k = 0
            for nd in nodes:
                if k < len(nodes) - 1:
                    if len(nd) == 1:
                        new_line += (nd + dash_or_space[k]*3)
                    elif len(nd) == 2:
                        new_line += (nd + dash_or_space[k]*2)
                else:  # k = len(nodes) - 1
                    new_line += nd
                k += 1
        else:
            new_line = pic_line

        self.picture_out.write(new_line)

    def write_NOTA(self, bla_str):
        """
        Writes a 'NOTA' line in eng & pic files. As the name implies,
        a NOTA is just a note or comment such as "I love you Mary". It is
        not a gate.

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        s = "NOTA\t" + bla_str + '\n'
        self.english_out.write(s)
        self.picture_out.write(s)

    def write_LOOP(self, loop_num, reps):
        """
        Writes a 'LOOP' line in eng & pic files. The gates between a LOOP
        line and its partner NEXT line are to be repeated a number of times
        called reps.

        Parameters
        ----------
        loop_num : int
        reps : int

        Returns
        -------
        None

        """
        s = "LOOP\t" + str(loop_num) + "\tREPS:\t" + str(reps) + '\n'
        self.english_out.write(s)
        self.picture_out.write(s)

    def write_NEXT(self, loop_num):
        """
        Writes a 'NEXT' line in eng & pic files.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        s = "NEXT\t" + str(loop_num) + '\n'
        self.english_out.write(s)
        self.picture_out.write(s)

    def write_controlled_bit_swap(self, bit1, bit2, trols):
        """
        Writes a line in eng & pic files for a 'SWAP' with >= 0 controls. 

        Parameters
        ----------
        bit1 : int
        bit2 : int
            bit1 and bit2 are the positions of the swapped bits.
        trols : Controls

        Returns
        -------
        None

        """

        # preamble, same for all 3 controlled gate methods
        self.line_counter += 1
        assert not self.english_out.closed
        assert not self.picture_out.closed

        num_bits_bef = self.emb.num_bits_bef
        num_bits_aft = self.emb.num_bits_aft
        # aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        aft_trols = trols.new_embedded_self(self.emb)
        # add extra controls if there are any
        extra_dict = self.emb.extra_controls.bit_pos_to_kind
        if extra_dict:
            aft_trols.bit_pos_to_kind.update(extra_dict)
            aft_trols.refresh_lists()

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)
        # end of preamble

        assert bit1 != bit2, "swapped bits must be different"
        assert -1 < bit1 < num_bits_bef
        assert -1 < bit2 < num_bits_bef
        x = [self.emb.aft(bit1), self.emb.aft(bit2)]
        big = max(x)
        small = min(x)

        # english file
        self.english_out.write("SWAP\t" + str(big) + "\t" + str(small))
        self.english_out.write("\tIF\t" if num_controls != 0 else "\n")

        # list bit-positions in decreasing order
        for c in range(num_controls):
            self.english_out.write(
                str(aft_trols.bit_pos[c]) +
                ("T" if aft_trols.kinds[c] == True else "F") +
                ("\n" if (c == num_controls - 1) else "\t"))

        # picture file
        pic_line = ""
        biggest = big
        smallest = small
        if num_controls != 0:
            biggest = max(aft_trols.bit_pos[0], big)
            smallest = min(aft_trols.bit_pos[num_controls-1], small)

        # k a bit position
        for k in range(num_bits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_big = (k == big)
            is_small = (k == small)
            is_control = False
            control_kind = False
            tres = "   " if (k == smallest) else "---"
            # dos = "  " if (k == smallest) else "--"

            for c in range(c_int, num_controls, +1):
                if k == aft_trols.bit_pos[c]:
                    is_control = True
                    control_kind = aft_trols.kinds[c]
                    c_int += 1
                    break

            if is_control:
                pic_line += ('@' if control_kind else 'O') + tres
            else:  # control not found
                if is_big:
                    pic_line += "<" + tres
                elif is_small:
                    pic_line += ">" + tres
                else:
                    pic_line += "+" + tres

        for k in range(smallest-1, -1, -1):
            pic_line += "|   "
            
        self.write_pic_line(pic_line)
        self.picture_out.write("\n")

    def write_controlled_one_bit_gate(
            self, tar_bit_pos, trols, one_bit_gate_fun, fun_arg_list=None):
        """
        Writes a line in eng & pic files for a one bit gate (from class 
        OneBitGates) with >= 0 controls. 

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        one_bit_gate_fun : function mapping Any->np.ndarray
        fun_arg_list : list[]

        Returns
        -------
        None

        """

        # preamble, same for all 3 controlled gate methods
        self.line_counter += 1
        assert not self.english_out.closed
        assert not self.picture_out.closed

        # num_bits_bef = self.emb.num_bits_bef
        num_bits_aft = self.emb.num_bits_aft
        aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        aft_trols = trols.new_embedded_self(self.emb)

        # add extra controls if there are any
        extra_dict = self.emb.extra_controls.bit_pos_to_kind
        if extra_dict:
            aft_trols.bit_pos_to_kind.update(extra_dict)
            aft_trols.refresh_lists()

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)
        # end of preamble

        assert aft_tar_bit_pos not in aft_trols.bit_pos,\
            "target bit cannot be a control bit"

        if fun_arg_list is None:
            fun_arg_list = []

        # english file
        if one_bit_gate_fun == OneBitGates.meas:
            kind = fun_arg_list[0]
            assert kind in [0, 1, 2], "unsupported measurement kind"
            self.english_out.write('MEAS\t' + str(kind) + '\t')
        elif one_bit_gate_fun == OneBitGates.phase_fac:
            self.english_out.write("PHAS\t" +
            str(fun_arg_list[0]*180/np.pi))
        elif one_bit_gate_fun == OneBitGates.P_0_phase_fac:
            self.english_out.write("P0PH\t" +
                str(fun_arg_list[0]*180/np.pi))
        elif one_bit_gate_fun == OneBitGates.P_1_phase_fac:
            self.english_out.write("P1PH\t" +
                str(fun_arg_list[0]*180/np.pi))
        elif one_bit_gate_fun == OneBitGates.sigx:
            self.english_out.write("SIGX")
        elif one_bit_gate_fun == OneBitGates.sigy:
            self.english_out.write("SIGY")
        elif one_bit_gate_fun == OneBitGates.sigz:
            self.english_out.write("SIGZ")
        elif one_bit_gate_fun == OneBitGates.had2:
            self.english_out.write("HAD2")
        elif one_bit_gate_fun == OneBitGates.rot_ax:
            ang_rads = fun_arg_list[0]
            axis = fun_arg_list[1]
            if axis == 1:
                self.english_out.write("ROTX\t" + str(ang_rads*180/np.pi))
            elif axis == 2:
                self.english_out.write("ROTY\t" + str(ang_rads*180/np.pi))
            elif axis == 3:
                self.english_out.write("ROTZ\t" + str(ang_rads*180/np.pi))
            else:
                assert False
        elif one_bit_gate_fun == OneBitGates.rot:
            x_degs = fun_arg_list[0]*180/np.pi
            y_degs = fun_arg_list[1]*180/np.pi
            z_degs = fun_arg_list[2]*180/np.pi
            self.english_out.write("ROTN\t" +
                str(x_degs) + "\t" + str(y_degs) + "\t" + str(z_degs))
        else:
            assert False, "writing an unsupported controlled gate"

        self.english_out.write("\tAT\t" + str(aft_tar_bit_pos) +
                               ("\tIF\t" if num_controls != 0 else "\n"))

        # list bit-positions in decreasing order
        for c in range(num_controls):
            self.english_out.write(
                str(aft_trols.bit_pos[c]) +
                ("T" if aft_trols.kinds[c] == True else "F") +
                ("\n" if c == num_controls - 1 else "\t"))

        # picture file
        pic_line = ""
        biggest = aft_tar_bit_pos
        smallest = aft_tar_bit_pos
        if num_controls != 0:
            biggest = max(aft_trols.bit_pos[0], aft_tar_bit_pos)
            smallest = min(aft_trols.bit_pos[num_controls-1], aft_tar_bit_pos)

        # k a bit position
        for k in range(num_bits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_target = (k == aft_tar_bit_pos)
            is_control = False
            control_kind = False
            tres = "   " if (k == smallest) else "---"
            dos = "  " if (k == smallest) else "--"

            # c_int starts at last value
            for c in range(c_int, num_controls, +1):
                if k == aft_trols.bit_pos[c]:
                    is_control = True
                    control_kind = aft_trols.kinds[c]
                    c_int += 1
                    break

            if is_control:
                pic_line += ("@" if control_kind else "O") + tres
            else:  # is not control
                if not is_target:  # is not control or target
                    pic_line += "+" + tres
                else:  # is target
                    if one_bit_gate_fun == OneBitGates.meas:
                        kind = fun_arg_list[0]
                        if kind == 0:
                            pic_line += "M0" + dos
                        elif kind == 1:
                            pic_line += 'M1' + dos
                        else:
                            pic_line += 'M' + tres
                    elif one_bit_gate_fun == OneBitGates.phase_fac:
                        pic_line += "Ph" + dos
                    elif one_bit_gate_fun == OneBitGates.P_0_phase_fac:
                        pic_line += "OP" + dos
                    elif one_bit_gate_fun == OneBitGates.P_1_phase_fac:
                        pic_line += "@P" + dos
                    elif one_bit_gate_fun == OneBitGates.sigx:
                        pic_line += "X" + tres
                    elif one_bit_gate_fun == OneBitGates.sigy:
                        pic_line += "Y" + tres
                    elif one_bit_gate_fun == OneBitGates.sigz:
                        pic_line += "Z" + tres
                    elif one_bit_gate_fun == OneBitGates.had2:
                        pic_line += "H" + tres
                    elif one_bit_gate_fun == OneBitGates.rot_ax:
                        axis = fun_arg_list[1]
                        if axis == 1:
                            pic_line += "Rx" + dos
                        elif axis == 2:
                            pic_line += "Ry" + dos
                        elif axis == 3:
                            pic_line += "Rz" + dos
                        else:
                            assert False
                    elif one_bit_gate_fun == OneBitGates.rot:
                        pic_line += "R" + tres
                    else:
                        assert False, "writing an unsupported controlled gate"

        for k in range(smallest-1, -1, -1):
            pic_line += "|   "
            
        self.write_pic_line(pic_line)
        self.picture_out.write("\n")

    def write_controlled_multiplexor_gate(self, tar_bit_pos,
            trols, with_minus, rad_angles):
        """
        Writes a line in eng & pic files for a multiplexor 'MP_Y' with
        >= 0 controls.

        The definition of multiplexors and how they are specified in both
        English and Picture files is described in the references given in
        the docstring of this class.

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        with_minus : bool
            When this is True, angles are multiplied by -1.
            Useful if need complex conjugate of multiplexor.
        rad_angles : list[float]

        Returns
        -------
        None

        """

        # preamble, same for all 3 controlled gate methods
        self.line_counter += 1
        assert not self.english_out.closed
        assert not self.picture_out.closed

        # num_bits_bef = self.emb.num_bits_bef
        num_bits_aft = self.emb.num_bits_aft
        aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        aft_trols = trols.new_embedded_self(self.emb)
        # add extra controls if there are any
        extra_dict = self.emb.extra_controls.bit_pos_to_kind
        if extra_dict:
            aft_trols.bit_pos_to_kind.update(extra_dict)
            aft_trols.refresh_lists()

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)
        # end of preamble
        
        assert aft_tar_bit_pos not in aft_trols.bit_pos,\
            "target bit cannot be a control bit"

        num_int_controls = aft_trols.get_num_int_controls()
        assert num_int_controls != 0, \
            "multiplexor with no half-moon controls"
        num_angles = len(rad_angles)
        assert num_angles == (1 << num_int_controls),\
            "wrong number of multiplexor angles"
        
        # english file
        self.english_out.write(
            "MP_Y\tAT\t" + str(aft_tar_bit_pos) + "\tIF\t")
        
        # list bit-positions in decreasing order
        for c in range(num_controls):
            x = aft_trols.kinds[c]
            kind_str = ""
            # int is subclass of bool
            # so isinstance(x, int) will be true for x bool too!
            if not isinstance(x, bool):
                kind_str = "(" + str(x)
            elif not x:
                kind_str = "F"
            elif x:
                kind_str = "T"
            self.english_out.write(
                str(aft_trols.bit_pos[c]) + kind_str + "\t")
        # use BY to indicate end of controls
        self.english_out.write("\tBY\t")
        for k in range(num_angles):
            self.english_out.write(
                str((-1 if with_minus else 1)*rad_angles[k]*180/np.pi) +
                ("\n" if k == (num_angles-1) else "\t"))
        
        # picture file
        pic_line = ""
        biggest = aft_tar_bit_pos
        smallest = aft_tar_bit_pos
        if num_controls != 0:
            biggest = max(aft_trols.bit_pos[0], aft_tar_bit_pos)
            smallest = min(aft_trols.bit_pos[num_controls-1], aft_tar_bit_pos)

        # k a bit position
        for k in range(num_bits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_target = (k == tar_bit_pos)
            is_control = False
            control_kind = False
            tres = "   " if (k == smallest) else "---"
            dos = "  " if (k == smallest) else "--"

            # c_int starts at last value
            for c in range(c_int, num_controls, +1):
                if k == aft_trols.bit_pos[c]:
                    is_control = True
                    control_kind = aft_trols.kinds[c]
                    c_int += 1
                    break

            if is_control:
                if not isinstance(control_kind, bool):
                    pic_line += "@O" + dos
                elif not control_kind:
                    pic_line += "O" + tres
                elif control_kind:
                    pic_line += "@" + tres
            else:  # is not control
                if not is_target:  # is not control nor target
                    pic_line += "+" + tres
                else:  # is target
                    pic_line += "Ry" + dos
        
        for k in range(smallest-1, -1, -1):
            pic_line += "|   "
            
        self.write_pic_line(pic_line)
        self.picture_out.write("\n")
        
    def write_bit_swap(self, bit1, bit2):
        """
        Write a line in eng & pic files for a 'SWAP' with no controls.

        Parameters
        ----------
        bit1 : int
        bit2 : int

        Returns
        -------
        None

        """
        trols = Controls(2)  # dummy with zero controls
        self.write_controlled_bit_swap(bit1, bit2, trols)

    def write_one_bit_gate(
            self, tar_bit_pos, one_bit_gate_fun, fun_arg_list=None):
        """
        Write a line in eng & pic files for a one qubit gate (from class
        OneBitGates) with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        one_bit_gate_fun : function
        fun_arg_list : list[]

        Returns
        -------
        None

        """
        trols = Controls(2)  # dummy with zero controls
        self.write_controlled_one_bit_gate(
            tar_bit_pos, trols, one_bit_gate_fun, fun_arg_list)

    def write_measurement(self, tar_bit_pos, kind):
        """
        Writes a 'MEAS' line in eng & pic files. This denotes a measurement
        step. We allow 3 kinds of measurements (0, 1, 2) at a target bit.

        Parameters
        ----------
        tar_bit_pos : int
        kind : int
            either 0, 1 or 2

        Returns
        -------
        None

        """
        self.write_one_bit_gate(tar_bit_pos, OneBitGates.meas, [kind])

    def write_global_phase_fac(self, ang_rads):
        """
        Write a line in eng & pic files for a global phase factor 'PHAS'
        with no controls.

        Parameters
        ----------
        ang_rads : float

        Returns
        -------
        None

        """
        tar_bit_pos = 0  # anyone will do
        trols = Controls(2)  # dummy with zero controls
        gate_fun = OneBitGates.phase_fac
        self.write_controlled_one_bit_gate(
            tar_bit_pos, trols, gate_fun, [ang_rads])
        
    def write_multiplexor_gate(
            self, tar_bit_pos, controls, with_minus, rad_angles):
        """
        Write a line in eng & pic files for a multiplexor 'MP_Y' with
        no T/F controls.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls
        with_minus : bool
        rad_angles : list[float]

        Returns
        -------
        None

        """
        num_controls = len(controls.bit_pos)
        assert num_controls == controls.get_num_int_controls(),\
            "some of the controls of this multiplexor are not half-moons"
        self.write_controlled_multiplexor_gate(
            tar_bit_pos, controls, with_minus, rad_angles)

if __name__ == "__main__":
    num_bits = 5
    emb = CktEmbedder(num_bits, num_bits)
    trols = Controls(num_bits)
    trols.bit_pos_to_kind = {3: True, 4: False}
    trols.refresh_lists()
    ang_rads = 30*np.pi/180

    for zf in [False, True]:
        wr = SEO_writer('io_folder//wr_test', emb, zero_bit_first=zf)

        wr.write_NOTA('zero bit first = ' + str(zf))

        wr.write_LOOP(10, 15)

        wr.write_NEXT(10)

        wr.write_controlled_bit_swap(0, 2, trols)

        wr.write_bit_swap(1, 2)

        tar_bit_pos = 2
        for kind in [0, 1, 2]:
            wr.write_measurement(tar_bit_pos, kind)

        gate = OneBitGates.phase_fac
        wr.write_controlled_one_bit_gate(2, trols, gate, [ang_rads])

        wr.write_global_phase_fac(30*np.pi/180)

        gate = OneBitGates.P_0_phase_fac
        wr.write_controlled_one_bit_gate(2, trols, gate, [ang_rads])

        gate = OneBitGates.P_1_phase_fac
        wr.write_controlled_one_bit_gate(2, trols, gate, [ang_rads])

        gate = OneBitGates.sigx
        wr.write_controlled_one_bit_gate(2, trols, gate)

        gate = OneBitGates.sigy
        wr.write_controlled_one_bit_gate(2, trols, gate)

        gate = OneBitGates.sigz
        wr.write_controlled_one_bit_gate(2, trols, gate)

        gate = OneBitGates.had2
        wr.write_controlled_one_bit_gate(2, trols, gate)

        gate = OneBitGates.rot_ax
        wr.write_controlled_one_bit_gate(2, trols, gate, [ang_rads, 1])

        gate = OneBitGates.rot_ax
        wr.write_controlled_one_bit_gate(2, trols, gate, [ang_rads, 2])

        gate = OneBitGates.rot_ax
        wr.write_controlled_one_bit_gate(2, trols, gate, [ang_rads, 3])

        gate = OneBitGates.rot
        wr.write_controlled_one_bit_gate(2, trols, gate,
            [ang_rads/3, ang_rads*2/3, ang_rads])

        gate = OneBitGates.sigx
        wr.write_one_bit_gate(2, gate)

        tar_bit_pos = 0
        trols1 = Controls(num_bits)
        trols1.bit_pos_to_kind = {1: 0, 2: 1, 3: True, 4: False}
        trols1.refresh_lists()
        wr.write_controlled_multiplexor_gate(tar_bit_pos, trols1, False,
            [ang_rads/3, ang_rads*2/3, ang_rads, ang_rads*4/3])

        trols2 = Controls(num_bits)
        trols2.bit_pos_to_kind = {1: 0, 2: 1}
        trols2.refresh_lists()
        wr.write_multiplexor_gate(tar_bit_pos, trols2, False,
            [ang_rads/3, ang_rads*2/3, ang_rads, ang_rads*4/3])

        wr.close_files()
