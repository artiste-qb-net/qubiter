from device_specific.Qubiter_to_AnyQasm import *
from device_specific.QbitPlanarLattice import *
import device_specific.chip_couplings_google as cc


class Qubiter_to_GoogleCirq(Qubiter_to_AnyQasm):
    """
    See docstring of parent class Qubiter_to_AnyQasm

    If input c_to_tars = None, all CNOTs allowed. If c_to_tars = 'do_fill',
    class fills c_to_tars.

    References
    ----------
    1. https://github.com/quantumlib/Cirq

    Attributes
    ----------
    lattice : QbitPlanarLattice

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
        self.lattice = QbitPlanarLattice(cc.BRISTLECONE_GRID)
        Qubiter_to_AnyQasm.__init__(self, *args, **kwargs)

    def write_prelude(self):
        """
        Writes Google's cirq opening statements before calls to use_ methods
        for gates.

        Returns
        -------
        None

        """
        if self.c_to_tars == 'do_fill':
            self.c_to_tars = self.lattice.get_c_to_tars()

        str_list = []
        str_list.append('import cirq\n')
        str_list.append('from cirq.devices import GridQubit\n')
        str_list.append('from cirq.ops import X, Y, Z, H, CNOT\n')
        str_list.append('from cirq.ops import RotXGate, RotYGate, RotZGate\n')
        str_list.append('\n')
        str_list.append('\n')
        str_list.append('ckt = cirq.Circuit()\n')

        for st in str_list:
            self.qasm_out.write(st)
        if self.write_qubiter_files:
            for st in str_list:
                self.qbtr_wr.write_NOTA(st[:-1])

    def write_ending(self):
        """
        Writes Google's cirq ending statements after calls to use_ methods
        for gates.

        Returns
        -------
        None

        """

        pass

    def use_HAD2(self, tar_bit_pos, controls):
        """
        Writes line in Cirq file corresponding to an English file line
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
        row, col = self.lattice.one2two(tar_bit_pos)
        self.qasm_out.write("ckt.append(H(GridQubit(" +
                            str(row) + ', ' +
                            str(col) + ")))\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                   OneBitGates.had2)

    def use_NOTA(self, bla_str):
        """
        Writes line in Cirq file corresponding to an English file line
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
        the function will comment the phase out in the output files (Cirq 
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
        Writes line in Cirq file corresponding to an English file line
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
        Writes line in Cirq file corresponding to an English file line
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

        row, col = self.lattice.one2two(tar_bit_pos)
        line_str = "ckt.append("
        if axis == 1:
            line_str += "RotXGate("
        elif axis == 2:
            line_str += "RotYGate("
        elif axis == 3:
            line_str += "RotZGate("
        else:
            assert False
        line_str += 'rads=' + str(2*rad_ang) + ').on('
        line_str += 'GridQubit(' + str(row) + ', ' + str(col) + ")))\n"
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [rad_ang, axis])

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        Writes line in Cirq file corresponding to an English file line
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
        row, col = self.lattice.one2two(tar_bit_pos)
        end_str = ').on(GridQubit(' + str(row) + ', ' + str(col) + ')))\n'

        line_str = 'ckt.append(RotZGate(rads=' + str(phi_ri) + end_str
        self.qasm_out.write(line_str)

        line_str = 'ckt.append(RotYGate(rads=' + str(phi_ce) + end_str
        self.qasm_out.write(line_str)

        line_str = 'ckt.append(RotZGate(rads=' + str(phi_le) + end_str
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
        Writes line in Cirq file corresponding to an English file line
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
            row, col = self.lattice.one2two(tar_bit_pos)
            line_str = "ckt.append("
            if axis == 1:
                line_str += "X("
            elif axis == 2:
                line_str += "Y("
            elif axis == 3:
                line_str += "Z("
            else:
                assert False
            line_str += 'GridQubit(' + str(row) + ', ' + str(col) + ")))\n"
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
                row_tar, col_tar = self.lattice.one2two(tar_pos)
                row_trol, col_trol = self.lattice.one2two(trol_pos)
                line_str = 'ckt.append(CNOT('
                line_str += 'GridQubit(' +\
                            str(row_trol) + ', ' + str(col_trol) + '), '
                line_str += 'GridQubit(' +\
                            str(row_tar) + ', ' + str(col_tar) + ')))\n'
                self.qasm_out.write(line_str)
                if self.write_qubiter_files:
                    self.qbtr_wr.write_cnot(trol_pos, tar_pos)
            else:
                assert False, "Forbidden CNOT detected: " \
                    + str(trol_pos) + "->" + str(tar_pos) \
                    + " in line " + str(self.line_count) \
                    + ". Use class ForbiddenCNotExpander " \
                    + "before attempting translation to Cirq."

if __name__ == "__main__":

    def main():
        file_prefix = "../io_folder/qbtr2google_test"
        qasm_name = 'GooCirq'
        num_bits = 5
        c_to_tars = 'do_fill'  # filled by constructor
        Qubiter_to_GoogleCirq(file_prefix, qasm_name,
                num_bits, c_to_tars, write_qubiter_files=True)

    main()
