from SEO_reader import *
from OneBitGates import *
from SEO_writer import *


class EchoingSEO_reader(SEO_reader):
    """
    This class is a child of SEO_reader. The class reads any previously 
    created Qubiter English file and it writes new English & Picture files 
    wherein every line is echoed faithfully. The purpose of this class is to 
    act as a parent class to the Expander classes MultiplexorExpander and 
    DiagUnitaryExpander, both of which echo every line except those starting 
    with MP_Y and DIAG, respectively, and therefore need only override the 
    functions use_MP_Y() and use_DIAG, respectively, of this class. 

    The constructor of this class takes as input wr which should be a
    SEO_writer or child thereof. If wr is a plain SEO_writer writing to a
    file_prefix different to the file_prefix being read, this class will
    generate a Picture File & English file starting from only an English
    file. So it is an ASCII drawing class with an English file as input.

    Attributes
    ----------
    gbit_list : list(int)
        Only needed for some expansion styles, this is a list of grounded bits
    wr : SEO_writer
        This object of SEO_writer is called inside every use_  function to
        do some writing in the output files.

    """

    def __init__(self, file_prefix, num_bits, wr, style, gbit_list=None):
        """

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        wr : SEO_writer
        style : str
        gbit_list : list(int)

        Returns
        -------
        None

        """
        self.gbit_list = gbit_list

        self.wr = wr

        SEO_reader.__init__(self, file_prefix, num_bits)

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
        self.wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                           OneBitGates.had2)

    def use_LOOP(self, loop_num, reps):
        """
        This function echoes a LOOP line.

        Parameters
        ----------
        loop_num : int
        reps : int

        Returns
        -------
        None

        """
        self.wr.write_LOOP(loop_num, reps)

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

    def use_PHAS(self, angle_degs, tar_bit_pos, controls):
        """
        This function echoes a PHAS line.

        Parameters
        ----------
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        ang_rads = angle_degs*np.pi/180
        self.wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
            OneBitGates.phase_fac, [ang_rads])

    def use_P_PH(self, projection_bit, angle_degs, tar_bit_pos, controls):
        """
        This function echoes a P0PH or P1PH line.

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
        if projection_bit == 0:
            u2_fun = OneBitGates.P_0_phase_fac
        elif projection_bit == 1:
            u2_fun = OneBitGates.P_1_phase_fac
        else:
            assert False

        ang_rads = angle_degs*np.pi/180
        self.wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                                              u2_fun, [ang_rads])

    def use_ROT(self, axis, angle_degs, tar_bit_pos, controls):
        """
        This function echoes a ROTX, ROTY or ROTZ line.

        Parameters
        ----------
        axis : int
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        rad_ang = angle_degs*np.pi/180
        self.wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                           OneBitGates.rot_ax, [rad_ang, axis])

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        This function echoes a ROTN line.

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
        rad_ang_list = list(np.array([angle_x_degs,
                                     angle_y_degs,
                                     angle_z_degs])*np.pi/180)
        self.wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                           OneBitGates.rot, rad_ang_list)

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
            u2_fun = OneBitGates.sigx
        elif axis == 2:
            u2_fun = OneBitGates.sigy
        elif axis == 3:
            u2_fun = OneBitGates.sigz
        else:
            assert False

        self.wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
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
        self.wr.write_controlled_bit_swap(bit1, bit2, controls)

    def write_log(self):
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

if __name__ == "__main__":
    print(5)