from CGateSEO_writer import *
from SEO_reader import *


class CktExpander(SEO_reader):
    """
    Qubiter English and Picture files allow single lines that represent U(2)
    matrices of numerous types with 0, 1 or multiple controls of either the
    n or nbar kind. If we say a gate is controlled, it may have 1 or more
    controls (it might be singly or multiply controlled). This class is a
    child of SEO_reader. The class reads any previously created Qubiter
    English file and it writes new English & Picture files wherein every
    line of the original English file is expanded into a sequence of (1)
    single qubit rotations and (2) simple CNOTs. Such a class is useful
    because many quantum computers (for example, IBM Quantum Experience) can
    only do (1) and (2).

    If the input English file has file_prefix as file prefix, then the
    output English & Picture files have as file prefix file_prefix + '_X1',
    assuming that '_X' + str(k) for some integer k is not already the ending
    of file_prefix. If it is, then the ending is changed to '_X' + str(
    k+1).

    Global phase factors are ignored, so expansions equal the original line
    up to a phase factor.

    You can get a count of the number of CNOTs in the expanded file by
    creating an object of the class SEO_reader.

    Attributes
    ----------
    wr : CGateSEO_writer
        This object of CGateSEO_writer, created in the class constructor,
        is called inside every use_  function to do some writing in the output
        files.

    """

    def __init__(self, file_prefix, num_bits, verbose=False):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        verbose : bool

        Returns
        -------
        None

        """
        self.verbose = verbose

        # temporary embedder
        emb = CktEmbedder(num_bits, num_bits)
        out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
        self.wr = CGateSEO_writer(out_file_prefix, emb,
            one_line=False, expand_1c_u2=True)

        SEO_reader.__init__(self, file_prefix, num_bits, verbose)

        self.wr.close_files()

    def two_embs_for_c_bit_swap(self, bit1, bit2, controls):
        """
        This internal function returns two CktEmbedder objects called emb0,
        emb1 that are used to write an expansion for a controlled swap of
        bits bit1, bit2.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        (CktEmbedder, CktEmbedder)

        """
        num_trols = len(controls.kinds)
        num_bits_bef = num_trols + 2
        num_bits_aft = self.num_bits

        bit_map0 = [0]*num_bits_bef
        for k in range(num_trols):
            bit_map0[k] = controls.bit_pos[num_trols-k-1]
        bit_map0[-2] = bit2
        bit_map0[-1] = bit1

        bit_map1 = bit_map0.copy()
        bit_map1[-2] = bit1
        bit_map1[-1] = bit2

        emb0 = CktEmbedder(num_bits_bef, num_bits_aft, bit_map0)
        emb1 = CktEmbedder(num_bits_bef, num_bits_aft, bit_map1)

        return emb0, emb1

    def emb_for_c_u2(self, tar_bit_pos, controls):
        """
        This internal function returns a CktEmbedder object called emb that
        is used to write an expansion for a controlled U(2) matrix with
        target at tar_bit_pos.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        CktEmbedder

        """
        num_trols = len(controls.kinds)
        num_bits_bef = num_trols + 1
        num_bits_aft = self.num_bits
        bit_map = [0]*num_bits_bef
        for k in range(num_trols):
            bit_map[k] = controls.bit_pos[num_trols-k-1]
        bit_map[-1] = tar_bit_pos
        emb = CktEmbedder(num_bits_bef, num_bits_aft, bit_map)
        return emb

    def write_gate_name(self, name, num_trols):
        """
        This function is used solely in the verbose mode. It announces with
        a NOTA comment line the beginning of a gate expansion.

        Parameters
        ----------
        name : str
        num_trols : int

        Returns
        -------
        None

        """
        skip0 = num_trols == 0 and not name == "SWAP"
        skip1 = num_trols == 1 and name == "SIGX"
        if skip0 or skip1:
            pass
        else:
            if self.verbose:
                self.wr.write_NOTA('--------expansion of ' + name)

    def use_DIAG(self, trols, rad_angles):
        """
        Returns error message if input circuit contains diagonal unitaries
        DIAG.

        Parameters
        ----------
        trols : Controls
        rad_angles : list[float]

        Returns
        -------

        """
        assert False, "This circuit contains diagonal unitaries DIAG." \
            "You should expand them first using class DiagUnitaryExpander"

    def use_HAD2(self, tar_bit_pos, controls):
        """
        This function expands a HAD2 line; i.e., it reads the line from the
        input English file and writes an expansion of it in the output
        English & Picture files.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        self.write_gate_name("HAD2", len(controls.kinds))

        self.wr.emb = self.emb_for_c_u2(tar_bit_pos, controls)

        self.wr.write(controls.kinds,
                      OneBitGates.had2)

    def use_LOOP(self, loop_num, reps):
        """
        This function echoes a LOOP line; i.e., it transcribes the line from
        the input English file to the output English & Picture files.

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
        This function echoes a MEAS line; i.e., it transcribes the line from
        the input English file to the output English & Picture files.

        Parameters
        ----------
        kind : int
        tar_bit_pos : int

        Returns
        -------
        None

        """
        self.wr.emb = CktEmbedder(self.num_bits, self.num_bits)
        # print("----", tar_bit_pos, self.wr.emb.bit_map)
        self.wr.write_MEAS(tar_bit_pos, kind)

    def use_MP_Y(self, tar_bit_pos, trols, rad_angles):
        """
        Returns error message if input circuit contains multiplexors MP_Y.

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------

        """
        assert False, "This circuit contains multiplexors." \
            "You should expand them first using class MultiplexorExpander"

    def use_NEXT(self, loop_num):
        """
        This function echoes a NEXT line; i.e., it transcribes the line from
        the input English file to the output English & Picture files.

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
        This function echoes a NOTA line; i.e., it transcribes the line from
        the input English file to the output English & Picture files.

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
        This function expands a PHAS line; i.e., it reads the line from the
        input English file and writes an expansion of it in the output
        English & Picture files.

        Parameters
        ----------
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        self.write_gate_name("PHAS", len(controls.kinds))

        self.wr.emb = self.emb_for_c_u2(tar_bit_pos, controls)
        ang_rads = angle_degs*np.pi/180
        self.wr.write(controls.kinds,
                      OneBitGates.phase_fac, [ang_rads])

    def use_P_PH(self, projection_bit, angle_degs, tar_bit_pos, controls):
        """
        This function expands a P0PH or P1PH line; i.e., it reads the line
        from the input English file and writes an expansion of it in the
        output English & Picture files.

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
        name = ''
        if projection_bit == 0:
            u2_fun = OneBitGates.P_0_phase_fac
            name = 'P0PH'
        elif projection_bit == 1:
            u2_fun = OneBitGates.P_1_phase_fac
            name = 'P1PH'
        else:
            assert False
        self.write_gate_name(name, len(controls.kinds))

        self.wr.emb = self.emb_for_c_u2(tar_bit_pos, controls)
        ang_rads = angle_degs*np.pi/180
        self.wr.write(controls.kinds, u2_fun, [ang_rads])

    def use_ROT(self, axis, angle_degs, tar_bit_pos, controls):
        """
        This function expands a ROTX, ROTY or ROTZ line; i.e., it reads the
        line from the input English file and writes an expansion of it in
        the output English & Picture files.

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
        name = ''
        if axis == 1:
            name = 'ROTX'
        elif axis == 2:
            name = 'ROTY'
        elif axis == 3:
            name = 'ROTZ'
        else:
            assert False
        self.write_gate_name(name, len(controls.kinds))

        self.wr.emb = self.emb_for_c_u2(tar_bit_pos, controls)

        rad_ang = angle_degs*np.pi/180
        self.wr.write(controls.kinds,
                      OneBitGates.rot_ax, [rad_ang, axis])

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        This function expands a ROTN line; i.e., it reads the line from the
        input English file and writes an expansion of it in the output
        English & Picture files.

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
        self.write_gate_name("ROTN", len(controls.kinds))

        self.wr.emb = self.emb_for_c_u2(tar_bit_pos, controls)

        rad_ang_list = list(np.array([angle_x_degs,
                                     angle_y_degs,
                                     angle_z_degs])*np.pi/180)
        self.wr.write(controls.kinds,
                      OneBitGates.rot,
                      rad_ang_list)

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        This function expands a SIGX, SIGY or SIGZ line; i.e., it reads the
        line from the input English file and writes an expansion of it in
        the output English & Picture files.

        Parameters
        ----------
        axis : int
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        name = ''
        if axis == 1:
            u2_fun = OneBitGates.sigx
            name = 'SIGX'
        elif axis == 2:
            u2_fun = OneBitGates.sigy
            name = 'SIGY'
        elif axis == 3:
            u2_fun = OneBitGates.sigz
            name = 'SIGZ'
        else:
            assert False
        self.write_gate_name(name, len(controls.kinds))

        self.wr.emb = self.emb_for_c_u2(tar_bit_pos, controls)
        self.wr.write(controls.kinds, u2_fun)

    def use_SWAP(self, bit1, bit2, controls):
        """
        This function expands a SWAP line; i.e., it reads the line from the
        input English file and writes an expansion of it in the output
        English & Picture files.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        self.write_gate_name("SWAP", len(controls.kinds))

        emb0, emb1 = self.two_embs_for_c_bit_swap(bit1, bit2, controls)
        num_trols = len(controls.kinds)

        self.wr.emb = emb0  # intialize emb

        # insert opening Hadamards for controls equal to n_bar = |0><0|
        self.wr.write_hads(controls.kinds)

        self.wr.write([True] * (num_trols + 1), OneBitGates.sigx)

        self.wr.emb = emb1  # change emb
        self.wr.write([True] * (num_trols + 1), OneBitGates.sigx)
        self.wr.emb = emb0  # restore emb

        self.wr.write([True] * (num_trols + 1), OneBitGates.sigx)

        # insert closing Hadamards for controls equal to n_bar = |0><0|
        self.wr.write_hads(controls.kinds, herm_conj=True)

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
    xer = CktExpander('io_folder/fou_test', 6, verbose=True)
    xer = CktExpander('io_folder/fou_test_X1', 6, verbose=True)
    xer = CktExpander('io_folder/ph_est_test', 8, verbose=True)
    xer = CktExpander('io_folder/sim_test2', 4, verbose=True)
    # write log file for sim_test2
    SEO_reader('io_folder/sim_test2', 4)
    xer = CktExpander('io_folder/sim_test3', 4, verbose=True)
