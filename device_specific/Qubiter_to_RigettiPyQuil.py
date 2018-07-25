from device_specific.Qubiter_to_AnyQasm import *
import device_specific.chip_couplings_rigetti as cc


class Qubiter_to_RigettiPyQuil(Qubiter_to_AnyQasm):
    """
    See docstring of parent class Qubiter_to_AnyQasm

    If input c_to_tars = None, all CNOTs allowed.

    References
    ----------
    1. https://github.com/rigetticomputing/pyquil

    Attributes
    ----------

    """
    def __init__(self, *args, **kwargs):
        """
        Constructor

        Parameters
        ----------

        Returns
        -------
        None

        """
        Qubiter_to_AnyQasm.__init__(self, *args, **kwargs)

    def write_prelude(self):
        """
        Writes PyQuil opening statements before calls to use_ methods for
        gates.

        Returns
        -------
        None

        """

        str_list = []
        str_list.append('from from pyquil.quil import Program\n')
        str_list.append('from pyquil.gates import import X, Y, Z, H, CNOT\n')
        str_list.append('from pyquil.gates import import RX, RY, RZ\n')
        str_list.append('\n')
        str_list.append('\n')
        str_list.append('pg = Program()\n')

        for st in str_list:
            self.qasm_out.write(st)
        if self.write_qubiter_files:
            for st in str_list:
                self.qbtr_wr.write_NOTA(st[:-1])

    def write_ending(self):
        """
        Writes PyQuil ending statements after calls to use_ methods for gates.

        Returns
        -------
        None

        """

        pass

    def use_HAD2(self, tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: HAD2 with no controls.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert len(controls.bit_pos) == 0
        self.qasm_out.write("pg.inst(H(" + str(tar_bit_pos) + "))\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                   OneBitGates.had2)

    def use_NOTA(self, bla_str):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: NOTA

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        self.qasm_out.write("# " + bla_str + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(bla_str)

    def use_PHAS(self, angle_degs, tar_bit_pos, controls):
        """
        If called for a controlled phase, this function will halt execution
        of program. If it's just a global phase with no controls,
        the function will comment the phase out in the output files (PyQuil
        and output Qubiter English and Picture files.) and move on to the
        next line.

        Parameters
        ----------
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        if controls.bit_pos_to_kind:
            assert False, "No PHAS lines with controls allowed"
        else:
            bla_str = 'PHAS\t' + str(angle_degs) + '\tAT\t' + str(tar_bit_pos)
            self.qasm_out.write("# " + bla_str + "\n")
            if self.write_qubiter_files:
                self.qbtr_wr.write_NOTA(bla_str)

    def use_PRINT(self, style, line_num):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: PRINT
        
        Parameters
        ----------
        style : str
        line_num : int

        Returns
        -------
        None

        """
        str1 = 'PRINT\t' + style
        self.qasm_out.write("# " + str1 + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(str1)
            
    def use_ROT(self, axis, angle_degs, tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: ROTX, ROTY or ROTZ with no controls.

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
        assert len(controls.bit_pos) == 0
        rad_ang = angle_degs*np.pi/180

        line_str = "pg.inst("
        if axis == 1:
            line_str += "RX("
        elif axis == 2:
            line_str += "RY("
        elif axis == 3:
            line_str += "RZ("
        else:
            assert False
        line_str += str(2*rad_ang) + ', '
        line_str += str(tar_bit_pos) + "))\n"
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [rad_ang, axis])

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: ROTN with no controls.

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
        assert len(controls.bit_pos) == 0
        rad_ang_list = list(np.array([angle_x_degs,
                                     angle_y_degs,
                                     angle_z_degs])*np.pi/180)

        arr = OneBitGates.rot(*rad_ang_list)
        delta, left_rads, center_rads, right_rads = \
            UnitaryMat.u2_zyz_decomp(arr)
        phi_le = ut.centered_rads(-2*left_rads)
        phi_ce = ut.centered_rads(-2*center_rads)
        phi_ri = ut.centered_rads(-2*right_rads)
        end_str = ', ' + str(tar_bit_pos) + '))\n'

        line_str = 'pg.inst(RZ(' + str(phi_ri) + end_str
        self.qasm_out.write(line_str)

        line_str = 'pg.inst(RY(' + str(phi_ce) + end_str
        self.qasm_out.write(line_str)

        line_str = 'pg.inst(RZ(' + str(phi_le) + end_str
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [right_rads, 3])
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [center_rads, 2])
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [left_rads, 3])

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: SIGX, SIGY or SIGZ with no controls, or else SIGX with one
        True control (i.e., simple CNOT).

        Parameters
        ----------
        axis : int
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        num_trols = len(controls.bit_pos)
        assert num_trols in [0, 1]
        if num_trols == 1:
            assert axis == 1
            assert controls.bit_pos_to_kind[controls.bit_pos[0]] == True
        if num_trols == 0:
            line_str = "pg.inst("
            if axis == 1:
                line_str += "X("
            elif axis == 2:
                line_str += "Y("
            elif axis == 3:
                line_str += "Z("
            else:
                assert False
            line_str += str(tar_bit_pos) + "))\n"
            self.qasm_out.write(line_str)
            if self.write_qubiter_files:
                if axis == 1:
                    u2_fun = OneBitGates.sigx
                elif axis == 2:
                    u2_fun = OneBitGates.sigy
                elif axis == 3:
                    u2_fun = OneBitGates.sigz
                else:
                    assert False

                self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos,
                                    controls, u2_fun)
        else:  # num_trols == 1
            tar_pos = tar_bit_pos
            trol_pos = controls.bit_pos[0]
            if not self.c_to_tars or tar_pos in self.c_to_tars[trol_pos]:
                line_str = 'pg.inst(CNOT('
                line_str += str(trol_pos) + ', '
                line_str += str(tar_pos) + '))\n'
                self.qasm_out.write(line_str)
                if self.write_qubiter_files:
                    self.qbtr_wr.write_cnot(trol_pos, tar_pos)
            else:
                assert False, "Forbidden CNOT detected: " \
                    + str(trol_pos) + "->" + str(tar_pos) \
                    + " in line " + str(self.line_count) \
                    + ". Use class ForbiddenCNotExpander " \
                    + "before attempting translation to PyQuil."

if __name__ == "__main__":
    import device_specific.chip_couplings_rigetti as rig

    def main():
        file_prefix = "../io_folder/qbtr2rigetti_test"
        qasm_name = 'RigPyQuil'
        num_bits = 6
        c_to_tars = rig.rigetti20_c_to_tars
        Qubiter_to_RigettiPyQuil(file_prefix, qasm_name,
                num_bits, c_to_tars, write_qubiter_files=True)

    main()
