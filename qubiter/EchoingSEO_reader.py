from qubiter.SEO_reader import *
from qubiter.OneQubitGate import *
from qubiter.SEO_writer import *


class EchoingSEO_reader(SEO_reader):
    """
    This class is a child of SEO_reader. The class reads any previously
    created Qubiter English file and it writes new English & Picture files
    wherein every line is echoed faithfully, except maybe the qubits are
    permuted.

    The constructor of this class takes as input wr which should be a
    SEO_writer or child thereof. If wr is a plain SEO_writer writing to a
    file_prefix different to the file_prefix being read, this class will
    generate a Picture File & English file starting from only an English
    file.

    This class has many uses. Here are some:

    1. If given as input a pictureless English file, it can be used to draw
    an ASCII picture of the input English file.

    2. If wr permutes the qubits, it can be used to write new Picture File
    & English files that have permuted qubits with respect to the input
    English file.

    3. It is the parent class to the Expander classes MultiplexorExpander
    and DiagUnitaryExpander, both of which echo every line except those
    starting with MP_Y and DIAG, respectively, and therefore need only
    override the functions use_MP_Y() and use_DIAG, respectively, of this
    class.

    Attributes
    ----------
    wr : SEO_writer
        This object of SEO_writer is called inside every use_  function to
        do some writing in the output files.

    """

    def __init__(self, file_prefix, num_qbits, wr,
                 vars_manager=None, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_qbits : int
        wr : SEO_writer
        vars_manager : PlaceholderManager

        Returns
        -------

        """
        self.wr = wr

        SEO_reader.__init__(self, file_prefix, num_qbits,
                            vars_manager=vars_manager, **kwargs)

        self.wr.close_files()

    def use_DIAG(self, controls, rad_angles):
        """
        This function echoes a DIAG line; i.e., it transcribes the line from
        the input English file to the output English & Picture files.

        Parameters
        ----------
        controls : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """

        self.wr.write_controlled_diag_unitary_gate(controls, rad_angles)

    def use_HAD2(self, tar_bit_pos, controls):
        """
        This function echoes a HAD2 line.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                           OneQubitGate.had2)

    def use_IF_M_beg(self, controls):
        """
        This function echoes IF_M_beg line.

        Parameters
        ----------
        controls : Controls

        Returns
        -------
        None

        """
        self.wr.write_IF_M_beg(controls)

    def use_IF_M_end(self):
        """
        This function echoes IF_M_end line

        Parameters
        ----------

        Returns
        -------
        None

        """
        self.wr.write_IF_M_end()

    def use_LOOP(self, loop_num, nreps):
        """
        This function echoes a LOOP line.

        Parameters
        ----------
        loop_num : int
        nreps : int

        Returns
        -------
        None

        """
        self.wr.write_LOOP(loop_num, nreps)

    def use_MEAS(self, tar_bit_pos, kind):
        """
        This function echoes a MEAS line.

        Parameters
        ----------
        kind : int
        tar_bit_pos : int

        Returns
        -------
        None

        """
        self.wr.write_MEAS(tar_bit_pos, kind)

    def use_MP_Y(self, tar_bit_pos, controls, rad_angles):
        """
        This function echoes an MP_Y line.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """

        self.wr.write_controlled_multiplexor_gate(
                    tar_bit_pos, controls, rad_angles)

    def use_NEXT(self, loop_num):
        """
        This function echoes a NEXT line.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        self.wr.write_NEXT(loop_num)

    def use_NOTA(self, bla_str):
        """
        This function echoes a NOTA line.

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        self.wr.write_NOTA(bla_str)

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        This function echoes a PHAS line.

        Parameters
        ----------
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
            OneQubitGate.phase_fac, [angle_rads])

    def use_PRINT(self, style, line_num):
        """
        This function echoes PRINT line

        Parameters
        ----------
        style : str
        line_num : int

        Returns
        -------
        None

        """
        self.wr.write_PRINT(style)

    def use_P_PH(self, projection_bit, angle_rads, tar_bit_pos, controls):
        """
        This function echoes a P0PH or P1PH line.

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
        if projection_bit == 0:
            u2_fun = OneQubitGate.P_0_phase_fac
        elif projection_bit == 1:
            u2_fun = OneQubitGate.P_1_phase_fac
        else:
            assert False

        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                                              u2_fun, [angle_rads])

    def use_ROTA(self, axis, angle_rads, tar_bit_pos, controls):
        """
        This function echoes a ROTX, ROTY or ROTZ line.

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
        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                           OneQubitGate.rot_ax, [angle_rads, axis])

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        This function echoes a ROTN line.

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
        rad_ang_list = [angle_x_rads, angle_y_rads, angle_z_rads]
        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                           OneQubitGate.rot, rad_ang_list)

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        This function echo a SIGX, SIGY or SIGZ line.

        Parameters
        ----------
        axis : int
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if axis == 1:
            u2_fun = OneQubitGate.sigx
        elif axis == 2:
            u2_fun = OneQubitGate.sigy
        elif axis == 3:
            u2_fun = OneQubitGate.sigz
        else:
            assert False

        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                                              u2_fun)

    def use_SWAP(self, bit1, bit2, controls):
        """
        This function echoes a SWAP line.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        self.wr.write_controlled_qbit_swap(bit1, bit2, controls)

    def use_SWAY(self, bit1, bit2, controls, rads_list):
        """
        This function echoes a SWAY line.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls
        rads_list : list[float | str]

        Returns
        -------
        None

        """
        self.wr.write_controlled_qbit_swap(bit1, bit2, controls, rads_list)

    def use_U_2_(self, rads0, rads1, rads2, rads3,
                tar_bit_pos, controls):
        """
        This function echoes a U_2_ line.

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
        rad_ang_list = [rads0, rads1, rads2, rads3]
        self.wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                           OneQubitGate.u2, rad_ang_list)

    def do_log(self):
        """
        This class does a "flat" reading of the input file; i.e.,
        the reading does not respect loop structure. Hence, we won't let it
        write a log file, for if we did, it would be incorrect. A correct
        log file can always be obtained by creating a SEO_reader object.

        Returns
        -------
        None

        """
        pass

    @staticmethod
    def pic_file_from_eng_file(file_prefix, num_qbits, ZL=True):
        """
        This function reads an English file with file prefix = file_prefix
        and it writes a Picture file for it with the same file prefix.

        Parameters
        ----------
        file_prefix : str
        num_qbits : int
        ZL : bool

        Returns
        -------
        None

        """
        end_str = '_' + str(num_qbits) + '_eng.txt'
        file_prefix_tempo = file_prefix + '_tempo'
        from shutil import copyfile
        copyfile(utg.preface(file_prefix + end_str),
                 utg.preface(file_prefix_tempo + end_str))

        emb = CktEmbedder(num_qbits, num_qbits)
        # English out file must different from English in file because one
        # can't read a file at the same time one is writing to it
        wr = SEO_writer(file_prefix, emb, ZL=ZL)
        vman = PlaceholderManager(eval_all_vars=False)
        EchoingSEO_reader(file_prefix_tempo, num_qbits, wr, vars_manager=vman)

        import os
        os.remove(utg.preface(file_prefix_tempo + end_str))


if __name__ == "__main__":
    def main():
        file_prefix_in = 'echo_test'
        file_prefix_out = 'echo_test_perm'
        num_qbits = 6

        # permute qubits by advancing their positions by 1
        bit_map = [1, 2, 3, 4, 5, 0]
        emb = CktEmbedder(num_qbits, num_qbits, bit_map)
        wr = SEO_writer(file_prefix_out, emb)
        EchoingSEO_reader(file_prefix_in, num_qbits, wr)

        EchoingSEO_reader.pic_file_from_eng_file(file_prefix_in, num_qbits)
    main()
