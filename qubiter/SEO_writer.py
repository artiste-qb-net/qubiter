from qubiter.CktEmbedder import *
from qubiter.Controls import *
from qubiter.OneQubitGate import *
import re
import qubiter.utilities_gen as utg
from qubiter.PlaceholderManager import *
import sys
import os
from IPython.display import HTML, display
if 'autograd.numpy' not in sys.modules:
    import numpy as np


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

    The class SEO_writer gives the bool option ZL. When this argument is set
    to True (resp., False), the Picture file shows the zero qubit last (
    resp., first), and the remaining qubits in consecutive order. Picture
    files written with zero bit last (resp., first) are labelled prefix +
    '_ZLpict.text' (resp., prefix + '_ZFpict.txt'). The ZL choice does not
    affect the English file.

    See the following and earlier arXiv papers by R.R.Tucci for more info on
    English and Picture files.

        http://arxiv.org/abs/1004.2205
        "Quibbs, a Code Generator for Quantum Gibbs Sampling"

    The following pdf is stored in the same folder as this python file.

        qubiter_rosetta_stone.pdf

    This pdf gives examples of lines in analytic/Picture/English formats.
    The Picture file examples follow the ZL convention.

    3 kinds (called 0, 1, 2) of measurements MEAS are allowed. A type 0
    measurement inserts a projector ``|0><0| = n = P_0`` at the target bit.
    A type 1 measurement inserts a projector ``|1><1| = nbar = P_1`` at the
    target bit. A type 2 measurement stores a copy of the state vector after
    ``|0><0|`` has been applied, and another copy after ``|1><1|`` has been
    applied.

    If a vertical wire hasn't been measured as type 2 measurement,
    it is drawn in pic file as "|";  otherwise, it is drawn as ":".

    Attributes
    ----------
    emb : CktEmbedder
    english_out : _io.TextIOWrapper
        file object for output text file that stores English description of
        circuit
    file_prefix : str
        beginning of the name of both English and Picture files
    gate_line_counter : int
    indentation : int
        Starts at 0, Grows by 4 at end of each write_LOOP and decrease by 4
        at beginning of each write_NEXT
    measured_bits : list(int)
        list of bits that have been measured with type 2 measurement and
        haven't been reset to ``|0>`` or ``|1>``
    picture_out : _io.TextIOWrapper
        file object for output text file that stores ASCII Picture
        description of circuit
    ZL : bool

    """

    def __init__(self, file_prefix, emb, ZL=True,
                english_out=None, picture_out=None):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        emb : CktEmbedder
        ZL : bool
        english_out : _io.TextIOWrapper
        picture_out : _io.TextIOWrapper


        Returns
        -------

        """
        self.gate_line_counter = 0
        self.file_prefix = file_prefix
        self.emb = emb
        self.ZL = ZL
        self.measured_bits = []

        if english_out is None and file_prefix:
            self.english_out = open(utg.preface(self.get_eng_file_path(
                rel=True)), 'wt')
        else:
            self.english_out = english_out

        if picture_out is None and file_prefix:
            self.picture_out = open(utg.preface(self.get_pic_file_path(
                rel=True)), 'wt')
        else:
            self.picture_out = picture_out

        self.indentation = 0

    def get_eng_file_path(self, rel=False):
        """
        Returns path (relative if rel is True, absolute if rel is False) of
        English file.

        Attributes
        ----------
        rel : bool

        Returns
        -------
        str

        """
        rel_path = utg.get_eng_file_rel_path(self.file_prefix,
                                             self.emb.num_qbits_aft)
        # print("..,,mmmm", rel_path)
        return rel_path if rel else utg.preface(rel_path)

    def get_pic_file_path(self, rel=False):
        """
        Returns path (relative if rel is True, absolute if rel is False) of
        Picture file.

        Attributes
        ----------
        rel : bool

        Returns
        -------
        str

        """
        rel_path = utg.get_pic_file_rel_path(self.file_prefix,
                                             self.emb.num_qbits_aft,
                                             ZL=self.ZL)
        return rel_path if rel else utg.preface(rel_path)

    def close_files(self):
        """
        Closes English and Picture files that were opened by the constructor.

        Returns
        -------
        None

        """
        # print("bbnnnvvv", 'closing files')
        self.english_out.close()
        self.picture_out.close()

    def delete_files(self):
        """
        This method will delete the English and Picture files. The files are
        closed before being deleted in case that hasn't been done yet.
        Closing a file a second time does nothing.

        Returns
        -------
        None

        """
        self.close_files()
        import os
        os.remove(self.get_eng_file_path())
        os.remove(self.get_pic_file_path())

    @staticmethod
    def gen_html_from_eng_or_pic_file(f):
        """
        This function returns a string for an html table generated from an
        English or Picture file. The table has two columns.
        The first column gives the line numbers starting from 1.

        Parameters
        ----------
        f : file
            f is file object returned by open()

        Returns
        -------
        str

        """
        k = 1
        line = f.readline()
        all_lines = ''
        td0 = "<td style='border-right:1px solid red;'>"
        td1 = "<td style='text-align:left;'>"
        while line:
            all_lines += td0 + str(k) + "</td>" + td1 +\
                         "<pre>" + line.strip() + "</pre>" +\
                         '</td></tr>'
            k += 1
            line = f.readline()
        table = "<table style='font-family:monospace'><tr>" +\
                all_lines + '</table>'
        # print(table)
        return table

    def print_eng_file(self, jup=False):
        """
        Prints English file.

        Parameters
        ----------
        jup : bool

            If jup=False, it prints text. Otherwise, it draws in a jupyter
            notebook a table with line numbers starting at 1

        Returns
        -------
        None

        """
        with open(utg.preface(self.get_eng_file_path(rel=True)), 'r') as f:
            if not jup:
                print(f.read())
            else:
                dis_obj = HTML(SEO_writer.gen_html_from_eng_or_pic_file(f))
                display(dis_obj)

    def print_pic_file(self, jup=False):
        """
        Prints Picture file.

        Parameters
        ----------
        jup : bool
            If jup=False, it prints text. Otherwise, it draws in a jupyter
            notebook a table with line numbers starting at 1

        Returns
        -------
        None

        """
        with open(utg.preface(self.get_pic_file_path(rel=True)), 'r') as f:
            if not jup:
                print(f.read())
            else:
                dis_obj = HTML(SEO_writer.gen_html_from_eng_or_pic_file(f))
                display(dis_obj)

    def colonize(self, pic_line):
        """
        This function returns new version of pic_line. Every "|" wire is
        replaced by ":" colon wire and vice versa iff the wire is at a bit
        position that has been measured (type 2 measurement) in the past and
        not reset. This function assumes pic_line in ZL convention so must
        call this function before calling write_ZF_or_ZL_pic_line()

        Parameters
        ----------
        pic_line : str

        Returns
        -------
        str

        """
        num_qbits = self.emb.num_qbits_aft
        li = list(pic_line)
        for bit in range(0, num_qbits):
            m = 4*(num_qbits - 1 - bit)
            if bit in self.measured_bits:
                if li[m] == "|":
                    li[m] = ":"
            else:
                if li[m] == ":":
                    li[m] = "|"

        return "".join(li)

    def write_ZF_or_ZL_pic_line(self, pic_line):
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
        if not self.ZL:
            pic_line = pic_line.strip()
            # delimiter is -- or --- or 2 spaces or 3 spaces
            nodes = re.split('-{2,3}| {2,3}', pic_line)
            dash_or_space = [pic_line[k*4 + 3] for k in range(len(nodes)-1)]
            nodes = list(reversed(nodes))
            dash_or_space = list(reversed(dash_or_space))
            for kk in range(len(nodes)):
                if nodes[kk] == '<':
                    nodes[kk] = '>'
                elif nodes[kk] == '<<':
                    nodes[kk] = '>>'
                elif nodes[kk] == '>':
                    nodes[kk] = '<'
                elif nodes[kk] == '>>':
                    nodes[kk] = '<<'
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

    def write_IF_M_beg(self, trols):
        """
        Writes a 'IF_M( <controls> ){' line in eng & pic files.

        Parameters
        ----------
        trols : Controls

        Returns
        -------
        None

        """
        aft_trols = trols.new_embedded_self(self.emb)
        s = "IF_M(\t"
        num_controls = len(aft_trols.bit_pos)
        for c in range(num_controls):
            s += str(aft_trols.bit_pos[c]) + \
                 ("T" if aft_trols.kinds[c] == True else "F") + \
                 ("\t){\n" if (c == num_controls - 1) else "\t")

        self.english_out.write(' '*self.indentation + s)
        self.picture_out.write(' '*self.indentation + s)

    def write_IF_M_end(self):
        """
        Writes an '}IF_M' line in eng and pic files.

        Parameters
        ----------

        Returns
        -------
        None

        """
        s = "}IF_M\n"
        self.english_out.write(' '*self.indentation + s)
        self.picture_out.write(' '*self.indentation + s)

    def write_LOOP(self, loop_num, nreps):
        """
        Writes a 'LOOP' line in eng & pic files. The gates between a LOOP
        line and its partner NEXT line are to be repeated a number of times
        called nreps.

        Parameters
        ----------
        loop_num : int
        nreps : int

        Returns
        -------
        None

        """
        s = "LOOP\t" + str(loop_num) + "\tNREPS=\t" + str(nreps) + '\n'
        self.english_out.write(' '*self.indentation + s)
        self.picture_out.write(' '*self.indentation + s)
        self.indentation += 4

    def write_MEAS(self, tar_bit_pos, kind):
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

        # num_qbits_bef = self.emb.num_qbits_bef
        num_qbits_aft = self.emb.num_qbits_aft
        aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        assert kind in [0, 1, 2], "unsupported measurement kind"
        if kind == 2:
            self.measured_bits.append(aft_tar_bit_pos)

        # english file
        s = 'MEAS\t' + str(kind) + '\tAT\t' + str(aft_tar_bit_pos) + '\n'
        self.english_out.write(' '*self.indentation + s)

        # picture file
        pic_line = ""
        biggest = aft_tar_bit_pos
        smallest = aft_tar_bit_pos
        # k a bit position
        for k in range(num_qbits_aft-1, biggest, -1):
            pic_line += "|   "
        if kind == 0:
            pic_line += "M0  "
        elif kind == 1:
            pic_line += 'M1  '
        else:
            pic_line += 'M   '

        for k in range(smallest-1, -1, -1):
            pic_line += "|   "

        pic_line = self.colonize(pic_line)
        self.write_ZF_or_ZL_pic_line(pic_line)
        self.picture_out.write("\n")

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
        self.indentation -= 4
        s = "NEXT\t" + str(loop_num) + '\n'
        self.english_out.write(' '*self.indentation + s)
        self.picture_out.write(' '*self.indentation + s)

    def write_NOTA(self, bla_str, permission=True):
        """
        Writes a 'NOTA' line in eng & pic files. As the name implies,
        a NOTA is just a note or comment such as "I love you Mary". It is
        not a gate.

        Parameters
        ----------
        bla_str : str
        permission : bool
            General permission, useful for turning off whole batches of NOTAs

        Returns
        -------
        None

        """
        if permission:
            s = "NOTA\t" + bla_str.rstrip() + '\n'
            self.english_out.write(' '*self.indentation + s)
            self.picture_out.write(' '*self.indentation + s)

    def write_PRINT(self, style):
        """
        Writes a 'PRINT' line in eng & pic files.

        Parameters
        ----------
        style : str

        Returns
        -------
        None

        """
        s = "PRINT\t" + style + '\n'
        self.english_out.write(' '*self.indentation + s)
        self.picture_out.write(' '*self.indentation + s)

    def rads_to_degs_str(self, rads):
        """
        This method returns

        str(rads*180/pi) if isinstance(rads, float)

        rads if is_legal_var_name(rads) (this implies rads is str)

        aborts otherwise.

        The method is only used inside this class so I am making it
        non-static even though it doesn't use self.

        Parameters
        ----------
        rads : float | str 

        Returns
        -------
        str

        """
        # print("--nnn", type(rads))
        # np.float types are different from float!!!
        if isinstance(rads, (int, float, np.floating)):
            # print("--nnn", str(rads*180/np.pi))
            return '{0:0.6f}'.format(rads*180/np.pi)
        else:
            assert PlaceholderManager.is_legal_var_name(rads), \
                "attempting to write an illegal variable name: '" +\
                str(rads) + "'"
            return rads

    def write_controlled_preamble(self, trols):
        """
        This is an internal function, used as preamble to all methods that
        are named write_controlled_...()

        Parameters
        ----------
        trols : Controls

        Returns
        -------
        Controls

        """
        self.gate_line_counter += 1
        assert not self.english_out.closed
        assert not self.picture_out.closed

        aft_trols = trols.new_embedded_self(self.emb)
        # add extra controls if there are any
        extra_dict = self.emb.extra_controls.bit_pos_to_kind
        if extra_dict:
            aft_trols.bit_pos_to_kind.update(extra_dict)
            aft_trols.refresh_lists()

        return aft_trols

    def write_controlled_one_qbit_gate(
            self, tar_bit_pos, trols, one_qbit_gate_fun, fun_arg_list=None):
        """
        Writes a line in eng & pic files for a one bit gate (from class 
        OneQubitGate) with >= 0 controls. 

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        one_qbit_gate_fun : function
            maps Any->np.ndarray
        fun_arg_list : list

        Returns
        -------
        None

        """
        aft_trols = self.write_controlled_preamble(trols)

        # num_qbits_bef = self.emb.num_qbits_bef
        num_qbits_aft = self.emb.num_qbits_aft
        aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)

        assert aft_tar_bit_pos not in aft_trols.bit_pos,\
            "target bit cannot be a control bit"

        if fun_arg_list is None:
            fun_arg_list = []

        # english file
        self.english_out.write(' '*self.indentation)
        if one_qbit_gate_fun == OneQubitGate.had2:
            self.english_out.write("HAD2")
        elif one_qbit_gate_fun == OneQubitGate.phase_fac:
            self.english_out.write("PHAS\t" +
            self.rads_to_degs_str(fun_arg_list[0]))
        elif one_qbit_gate_fun == OneQubitGate.P_0_phase_fac:
            self.english_out.write("P0PH\t" +
                self.rads_to_degs_str(fun_arg_list[0]))
        elif one_qbit_gate_fun == OneQubitGate.P_1_phase_fac:
            self.english_out.write("P1PH\t" +
                self.rads_to_degs_str(fun_arg_list[0]))
        elif one_qbit_gate_fun == OneQubitGate.rot_ax:
            ang_rads = fun_arg_list[0]
            axis = fun_arg_list[1]
            degs_str = self.rads_to_degs_str(ang_rads)
            if axis == 1:
                self.english_out.write("ROTX\t" + degs_str)
            elif axis == 2:
                self.english_out.write("ROTY\t" + degs_str)
            elif axis == 3:
                self.english_out.write("ROTZ\t" + degs_str)
            else:
                assert False
        elif one_qbit_gate_fun == OneQubitGate.rot:
            x_degs = self.rads_to_degs_str(fun_arg_list[0])
            y_degs = self.rads_to_degs_str(fun_arg_list[1])
            z_degs = self.rads_to_degs_str(fun_arg_list[2])
            self.english_out.write("ROTN\t" +
                x_degs + "\t" + y_degs + "\t" + z_degs)
        elif one_qbit_gate_fun == OneQubitGate.sigx:
            self.english_out.write("SIGX")
        elif one_qbit_gate_fun == OneQubitGate.sigy:
            self.english_out.write("SIGY")
        elif one_qbit_gate_fun == OneQubitGate.sigz:
            self.english_out.write("SIGZ")
        elif one_qbit_gate_fun == OneQubitGate.u2:
            ph_degs = self.rads_to_degs_str(fun_arg_list[0])
            x_degs = self.rads_to_degs_str(fun_arg_list[1])
            y_degs = self.rads_to_degs_str(fun_arg_list[2])
            z_degs = self.rads_to_degs_str(fun_arg_list[3])
            self.english_out.write("U_2_\t" + ph_degs + "\t" +
                x_degs + "\t" + y_degs + "\t" + z_degs)
        else:
            assert False, "writing an unsupported controlled gate\n" +\
                            one_qbit_gate_fun.__name__

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
        for k in range(num_qbits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_target = (k == aft_tar_bit_pos)
            is_control = False
            control_kind = False
            tres = ' '*3 if (k == smallest) else "---"
            dos = ' '*2 if (k == smallest) else "--"

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
                    if one_qbit_gate_fun == OneQubitGate.had2:
                        pic_line += "H" + tres
                    elif one_qbit_gate_fun == OneQubitGate.phase_fac:
                        pic_line += "Ph" + dos
                    elif one_qbit_gate_fun == OneQubitGate.P_0_phase_fac:
                        pic_line += "OP" + dos
                    elif one_qbit_gate_fun == OneQubitGate.P_1_phase_fac:
                        pic_line += "@P" + dos
                    elif one_qbit_gate_fun == OneQubitGate.rot_ax:
                        axis = fun_arg_list[1]
                        if axis == 1:
                            pic_line += "Rx" + dos
                        elif axis == 2:
                            pic_line += "Ry" + dos
                        elif axis == 3:
                            pic_line += "Rz" + dos
                        else:
                            assert False
                    elif one_qbit_gate_fun == OneQubitGate.rot:
                        pic_line += "R" + tres
                    elif one_qbit_gate_fun == OneQubitGate.sigx:
                        pic_line += "X" + tres
                    elif one_qbit_gate_fun == OneQubitGate.sigy:
                        pic_line += "Y" + tres
                    elif one_qbit_gate_fun == OneQubitGate.sigz:
                        pic_line += "Z" + tres
                    elif one_qbit_gate_fun == OneQubitGate.u2:
                        rads1 = fun_arg_list[1]
                        rads2 = fun_arg_list[2]
                        rads3 = fun_arg_list[3]

                        def both_are_small(r1, r2):
                            if isinstance(r1, str) or isinstance(r2, str):
                                return False
                            if abs(r1) > 1e-6 or abs(r2) > 1e-6:
                                return False
                            return True

                        small_12 = both_are_small(rads1, rads2)
                        small_23 = both_are_small(rads2, rads3)
                        small_13 = both_are_small(rads1, rads3)

                        if small_12 and small_23:
                            pic_line += "Ph" + dos
                        elif small_23:
                            pic_line += "Ux" + dos
                        elif small_13:
                            pic_line += "Uy" + dos
                        elif small_12:
                            pic_line += "Uz" + dos
                        else:
                            pic_line += "U" + tres
                    else:
                        assert False, \
                            "writing an unsupported controlled gate\n" +\
                            one_qbit_gate_fun.__name__

        for k in range(smallest-1, -1, -1):
            pic_line += "|   "

        pic_line = self.colonize(pic_line)
        self.write_ZF_or_ZL_pic_line(pic_line)
        self.picture_out.write("\n")

    def write_controlled_qbit_swap(self, bit1, bit2, trols, rads_list=None):
        """
        If rads_list=None, this method writes a line in eng & pic files for
        a 'SWAP' with >= 0 controls

        NOTE: SWAP is qbit symmetric: SWAP(0,1) = SWAP(1,0)

        If rads_list is not None and equals a list of 2 angles, rads_list=[
        rads0, rads1], this method writes a generalization of SWAP that I
        call SWAY (just to have a verb that is close to swap)

        SWAY  =
        [1  0  0]
        [0  U2 0]
        [0  0  1]

        where U2 is the most general 2-dim unitary matrix satisfying

        sigx U2 sigx = U2.

        If U2 is parametrized as

        U2 = exp(j*(rads0 + rads1*sigx + rads2*sigy + rads3*sigz))

        then SWAY is qbit symmetric (SWAY(0,1)=SWAY(1,0)) iff sigx U2 sigx
        = U2 iff rads2=rads3=0.

        SWAY includes SWAP, sqrt(SWAP), iSWAP, sqrt(iSWAP), PWAP,
        sqrt(PSWAP) etc.


        Parameters
        ----------
        bit1 : int
        bit2 : int
            bit1 and bit2 are the positions of the swapped bits.
        trols : Controls
        rads_list : list[float | str ] | None

        Returns
        -------
        None

        """
        aft_trols = self.write_controlled_preamble(trols)

        num_qbits_bef = self.emb.num_qbits_bef
        num_qbits_aft = self.emb.num_qbits_aft
        # aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)

        assert bit1 != bit2, "swapped bits must be different"
        assert -1 < bit1 < num_qbits_bef
        assert -1 < bit2 < num_qbits_bef
        x = [self.emb.aft(bit1), self.emb.aft(bit2)]
        big = max(x)
        small = min(x)

        use_sway = False
        if rads_list is not None:
            assert len(rads_list) == 2
            use_sway = True

        # english file
        gate_name = 'SWAP'
        if use_sway:
            gate_name = 'SWAY'
        self.english_out.write(' '*self.indentation + gate_name +
                               '\t' + str(big) + "\t" + str(small))
        if use_sway:
            self.english_out.write('\tBY')
            for k in range(2):
                deg_str = self.rads_to_degs_str(rads_list[k])
                self.english_out.write('\t' + deg_str)
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
        for k in range(num_qbits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_big = (k == big)
            is_small = (k == small)
            is_control = False
            control_kind = False
            tres = ' '*3 if (k == smallest) else "---"
            dos = ' '*2 if (k == smallest) else "--"

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
                    if use_sway:
                        pic_line += "<<" + dos
                    else:
                        pic_line += "<" + tres
                elif is_small:
                    if use_sway:
                        pic_line += ">>" + dos
                    else:
                        pic_line += ">" + tres
                else:
                    pic_line += "+" + tres

        for k in range(smallest-1, -1, -1):
            pic_line += "|   "

        pic_line = self.colonize(pic_line)
        self.write_ZF_or_ZL_pic_line(pic_line)
        self.picture_out.write("\n")

    def write_controlled_multiplexor_gate(self,
                tar_bit_pos, trols, rad_angles):
        """
        Writes a line in eng & pic files for a multiplexor 'MP_Y' with >= 0
        controls of either int (integer, intrinsic) or True/False kind.


        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        aft_trols = self.write_controlled_preamble(trols)

        # num_qbits_bef = self.emb.num_qbits_bef
        num_qbits_aft = self.emb.num_qbits_aft
        aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)
        
        assert aft_tar_bit_pos not in aft_trols.bit_pos,\
            "target bit cannot be a control bit"

        num_int_controls = aft_trols.get_num_int_controls()
        assert num_int_controls != 0, \
            "multiplexor with no half-moon controls"
        num_angles = len(rad_angles)
        assert num_angles == (1 << num_int_controls),\
            "wrong number of multiplexor angles"
        
        # english file
        self.english_out.write(' '*self.indentation +
            "MP_Y\tAT\t" + str(aft_tar_bit_pos) + "\tIF\t")
        
        # list bit-positions in decreasing order
        for c in range(num_controls):
            x = aft_trols.kinds[c]
            kind_str = ""
            # bool is subclass of int
            # so isinstance(x, int) will be true if x is bool!
            if not isinstance(x, bool):
                kind_str = ":" + str(x)
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
                self.rads_to_degs_str(rad_angles[k]) +
                ("\n" if k == (num_angles-1) else "\t"))
        
        # picture file
        pic_line = ""
        biggest = aft_tar_bit_pos
        smallest = aft_tar_bit_pos
        if num_controls != 0:
            biggest = max(aft_trols.bit_pos[0], aft_tar_bit_pos)
            smallest = min(aft_trols.bit_pos[num_controls-1], aft_tar_bit_pos)

        # k a bit position
        for k in range(num_qbits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_target = (k == tar_bit_pos)
            is_control = False
            control_kind = False
            tres = ' '*3 if (k == smallest) else "---"
            dos = ' '*2 if (k == smallest) else "--"

            # c_int starts at last value
            for c in range(c_int, num_controls, +1):
                if k == aft_trols.bit_pos[c]:
                    is_control = True
                    control_kind = aft_trols.kinds[c]
                    c_int += 1
                    break

            if is_control:
                if not isinstance(control_kind, bool):
                    pic_line += "%" + tres
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

        pic_line = self.colonize(pic_line)
        self.write_ZF_or_ZL_pic_line(pic_line)
        self.picture_out.write("\n")

    def write_controlled_diag_unitary_gate(self, trols, rad_angles):
        """
        Writes a line in eng & pic files for a diagonal unitary gate 'DIAG'
        with >= 0 controls of either int (integer, intrinsic) or True/False
        kind.

        Parameters
        ----------
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        aft_trols = self.write_controlled_preamble(trols)

        # num_qbits_bef = self.emb.num_qbits_bef
        num_qbits_aft = self.emb.num_qbits_aft
        # aft_tar_bit_pos = self.emb.aft(tar_bit_pos)

        # number of controls may be zero
        num_controls = len(aft_trols.bit_pos)

        num_int_controls = aft_trols.get_num_int_controls()
        assert num_int_controls != 0, \
            "d-unitary with no half-moon controls"
        num_angles = len(rad_angles)
        assert num_angles == (1 << num_int_controls),\
            "wrong number of d-unitary angles"

        # english file
        self.english_out.write(' '*self.indentation +
                               "DIAG\tIF\t")

        # list bit-positions in decreasing order
        for c in range(num_controls):
            x = aft_trols.kinds[c]
            kind_str = ""
            # bool is subclass of int
            # so isinstance(x, int) will be true if x is bool!
            if not isinstance(x, bool):
                kind_str = ":" + str(x)
            elif not x:
                kind_str = "F"
            elif x:
                kind_str = "T"
            self.english_out.write(
                str(aft_trols.bit_pos[c]) + kind_str + "\t")
        # use BY to indicate end of controls
        self.english_out.write("\tBY\t")
        for k in range(num_angles):
            # print('nmnmkkk', k, num_angles, rad_angles[k])
            self.english_out.write(
                self.rads_to_degs_str(rad_angles[k]) +
                ("\n" if k == (num_angles-1) else "\t"))

        # picture file
        pic_line = ""
        biggest = aft_trols.bit_pos[0]
        smallest = aft_trols.bit_pos[num_controls-1]

        # k a bit position
        for k in range(num_qbits_aft-1, biggest, -1):
            pic_line += "|   "

        c_int = 0
        for k in range(biggest, smallest-1, -1):
            is_control = False
            control_kind = False
            tres = ' '*3 if (k == smallest) else "---"

            # c_int starts at last value
            for c in range(c_int, num_controls, +1):
                if k == aft_trols.bit_pos[c]:
                    is_control = True
                    control_kind = aft_trols.kinds[c]
                    c_int += 1
                    break

            if is_control:
                if not isinstance(control_kind, bool):
                    pic_line += "%" + tres
                elif not control_kind:
                    pic_line += "O" + tres
                elif control_kind:
                    pic_line += "@" + tres
            else:  # is not control
                pic_line += "+" + tres

        for k in range(smallest-1, -1, -1):
            pic_line += "|   "

        pic_line = self.colonize(pic_line)
        self.write_ZF_or_ZL_pic_line(pic_line)
        self.picture_out.write("\n")

    def write_one_qbit_gate(
            self, tar_bit_pos, one_qbit_gate_fun, fun_arg_list=None):
        """
        Write a line in eng & pic files for a one qubit gate (from class
        OneQubitGate) with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        one_qbit_gate_fun : function
        fun_arg_list : list

        Returns
        -------
        None

        """
        trols = Controls(2)  # dummy with zero controls
        self.write_controlled_one_qbit_gate(
            tar_bit_pos, trols, one_qbit_gate_fun, fun_arg_list)

    def write_qbit_swap(self, bit1, bit2, rads_list=None):
        """
        Write a line in eng & pic files for a 'SWAP' if rads_list=None or
        'SWAY' if rads_list!=None, with no controls.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        rads_list : list[float | str] | None

        Returns
        -------
        None

        """
        trols = Controls(2)  # dummy with zero controls
        self.write_controlled_qbit_swap(bit1, bit2, trols, rads_list)

    def write_H(self, tar_bit_pos):
        """
        Writes HAD2 with no controls.

        Parameters
        ----------
        tar_bit_pos : int

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.had2)

    def write_Rx(self, tar_bit_pos, rads):
        """
        writes ROTX = exp(1j*rads*sig_x) with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        rads : float

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.rot_ax, [rads, 1])

    def write_Ry(self, tar_bit_pos, rads):
        """
        Writes ROTY = exp(1j*rads*sig_y) with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        rads : float

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.rot_ax, [rads, 2])

    def write_Rz(self, tar_bit_pos, rads):
        """
        Writes ROTZ = exp(1j*rads*sig_z) with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        rads : float

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.rot_ax, [rads, 3])

    def write_Rn(self, tar_bit_pos, rads_list):
        """
        Writes

        ROTN = exp(1j*(rads_x*sig_x + rads_y*sig_y + rads_z*sig_z))

        with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        rads_list : list[float]
            [rads_x, rads_y, rads_z]

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.rot, rads_list)

    def write_S(self, tar_bit_pos, herm=False):
        """
        Writes P1PH = exp(1j*P_1*pi/2) or its Hermitian with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        herm : bool

        Returns
        -------
        None

        """
        sign = +1
        if herm:
            sign = -1
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.P_1_phase_fac,
                                [sign*np.pi/2])

    def write_T(self, tar_bit_pos, herm=False):
        """
        Writes P1PH = exp(1j*P_1*pi/4) or its Hermitian with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        herm : bool

        Returns
        -------
        None

        """
        sign = +1
        if herm:
            sign = -1
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.P_1_phase_fac,
                                [sign*np.pi/4])

    def write_U2(self, tar_bit_pos, rads_list):
        """
        Writes

        UN_2= exp(1j*(rads0 + rads1*sig_x + rads2*sig_y + rads3*sig_z))

        with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        rads_list : list[float]
            [rads0, rads1, rads2, rads3]

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.u2, rads_list)

    def write_X(self, tar_bit_pos):
        """
        Writes SIGX with no controls.

        Parameters
        ----------
        tar_bit_pos : int

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.sigx)

    def write_Y(self, tar_bit_pos):
        """
        Writes SIGY with no controls.

        Parameters
        ----------
        tar_bit_pos : int

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.sigy)

    def write_Z(self, tar_bit_pos):
        """
        Writes SIGZ with no controls.

        Parameters
        ----------
        tar_bit_pos : int

        Returns
        -------
        None

        """
        self.write_one_qbit_gate(tar_bit_pos, OneQubitGate.sigz)

    def write_cnot(self, control_bit, target_bit, kind=True):
        """
        Writes a simple singly controlled not. If kind=True (resp. False),
        cnot fires when control is ``|1>`` (resp. ``|0>``).

        Parameters
        ----------
        control_bit : int
        target_bit : int
        kind : bool

        Returns
        -------
        None

        """
        num_qbits = self.emb.num_qbits_aft
        trols = Controls.new_single_trol(num_qbits, control_bit, kind)
        self.write_controlled_one_qbit_gate(target_bit, trols,
            OneQubitGate.sigx)

    def write_cz(self, control_bit, target_bit, kind=True):
        """
        Writes a simple singly controlled Z. If kind=True (resp. False),
        cz fires when control is ``|1>`` (resp. ``|0>``).

        Parameters
        ----------
        control_bit : int
        target_bit : int
        kind : bool

        Returns
        -------
        None

        """
        num_qbits = self.emb.num_qbits_aft
        trols = Controls.new_single_trol(num_qbits, control_bit, kind)
        self.write_controlled_one_qbit_gate(target_bit, trols,
            OneQubitGate.sigz)

    def write_c_P1PH(self, control_bit, target_bit, rads=np.pi, kind=True):
        """
        Writes a simple singly controlled P1PH. If kind=True (resp. False),
        c_P1PH fires when control is ``|1>`` (resp. ``|0>``). When kind=
        True and rads=p1, c_P1PH equals ``(-1)^{n(t)n(c)} = sigz(t)^{n(c)}``
        where c is the control and t is the target. This is often called a
        controlled Z, and denoted by Cz.

        Parameters
        ----------
        control_bit : int
        target_bit : int
        rads : float
        kind : bool

        Returns
        -------
        None

        """
        num_qbits = self.emb.num_qbits_aft
        trols = Controls.new_single_trol(num_qbits, control_bit, kind)
        self.write_controlled_one_qbit_gate(target_bit, trols,
            OneQubitGate.P_1_phase_fac, [rads])

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
        gate_fun = OneQubitGate.phase_fac
        self.write_controlled_one_qbit_gate(
            tar_bit_pos, trols, gate_fun, [ang_rads])
        
    def write_multiplexor_gate(self, tar_bit_pos, controls, rad_angles):
        """
        Write a line in eng & pic files for a multiplexor 'MP_Y' with
        no T/F controls.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        num_controls = len(controls.bit_pos)
        assert num_controls == controls.get_num_int_controls(),\
            "some of the controls of this multiplexor are not half-moons"
        self.write_controlled_multiplexor_gate(
            tar_bit_pos, controls, rad_angles)

    def write_diag_unitary_gate(self, controls, rad_angles):
        """
        Write a line in eng & pic files for a diagonal unitary 'DIAG' with
        no T/F controls.

        Parameters
        ----------
        controls : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        num_controls = len(controls.bit_pos)
        assert num_controls == controls.get_num_int_controls(),\
            "some of the controls of this diagonal unitary are not half-moons"
        self.write_controlled_diag_unitary_gate(controls, rad_angles)


if __name__ == "__main__":
    def main():
        num_qbits = 5
        emb = CktEmbedder(num_qbits, num_qbits)
        trols = Controls(num_qbits)
        trols.bit_pos_to_kind = {3: True, 4: False}
        trols.refresh_lists()
        ang_rads = 30*np.pi/180

        for ZL in [False, True]:
            wr = SEO_writer('wr_test', emb, ZL=ZL)

            wr.write_NOTA('zero bit last = ' + str(ZL))

            wr.write_IF_M_beg(trols)
            wr.write_IF_M_end()

            wr.write_LOOP(10, 15)
            wr.write_NEXT(10)

            tar_bit_pos = 1
            for kind in [0, 1, 2]:
                wr.write_MEAS(tar_bit_pos, kind)

            wr.write_PRINT('F2')

            wr.write_controlled_qbit_swap(0, 2, trols)

            wr.write_qbit_swap(1, 2)

            gate = OneQubitGate.phase_fac
            wr.write_controlled_one_qbit_gate(2, trols, gate, [ang_rads])

            wr.write_global_phase_fac(30*np.pi/180)

            gate = OneQubitGate.P_0_phase_fac
            wr.write_controlled_one_qbit_gate(2, trols, gate, [ang_rads])

            gate = OneQubitGate.P_1_phase_fac
            wr.write_controlled_one_qbit_gate(2, trols, gate, [ang_rads])

            gate = OneQubitGate.sigx
            wr.write_controlled_one_qbit_gate(2, trols, gate)

            gate = OneQubitGate.sigy
            wr.write_controlled_one_qbit_gate(2, trols, gate)

            gate = OneQubitGate.sigz
            wr.write_controlled_one_qbit_gate(2, trols, gate)

            gate = OneQubitGate.had2
            wr.write_controlled_one_qbit_gate(2, trols, gate)

            gate = OneQubitGate.rot_ax
            wr.write_controlled_one_qbit_gate(2, trols, gate, [ang_rads, 1])

            gate = OneQubitGate.rot_ax
            wr.write_controlled_one_qbit_gate(2, trols, gate, [ang_rads, 2])

            gate = OneQubitGate.rot_ax
            wr.write_controlled_one_qbit_gate(2, trols, gate, [ang_rads, 3])

            gate = OneQubitGate.rot
            wr.write_controlled_one_qbit_gate(2, trols, gate,
                [ang_rads/3, ang_rads*2/3, ang_rads])

            gate = OneQubitGate.sigx
            wr.write_one_qbit_gate(2, gate)

            wr.write_cnot(2, 1)

            tar_bit_pos = 0
            trols1 = Controls(num_qbits)
            trols1.bit_pos_to_kind = {1: 0, 2: 1, 3: True, 4: False}
            trols1.refresh_lists()
            wr.write_controlled_multiplexor_gate(tar_bit_pos, trols1,
                [ang_rads/3, ang_rads*2/3, ang_rads, ang_rads*4/3])

            trols2 = Controls(num_qbits)
            trols2.bit_pos_to_kind = {1: 0, 2: 1}
            trols2.refresh_lists()
            wr.write_multiplexor_gate(tar_bit_pos, trols2,
                [ang_rads/3, ang_rads*2/3, ang_rads, ang_rads*4/3])

            wr.close_files()


    main()
