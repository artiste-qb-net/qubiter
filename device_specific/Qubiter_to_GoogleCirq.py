from device_specific.Qubiter_to_AnyQasm import *
from device_specific.QbitPlanarLattice import *
import device_specific.chip_couplings_google as cc
import utilities_gen as ug


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
    def __init__(self, file_prefix, num_bits, **kwargs):
        """
        Constructor

        Parameters
        ----------

        Returns
        -------

        """
        self.lattice = QbitPlanarLattice(cc.BRISTLECONE_GRID)
        Qubiter_to_AnyQasm.__init__(self, file_prefix, num_bits, **kwargs)

    def write_prelude(self):
        """
        Writes Cirq opening statements before calls to use_ methods for gates.

        Returns
        -------
        None

        """
        if self.c_to_tars == 'do_fill':
            self.c_to_tars = self.lattice.get_c_to_tars()

        s = 'import cirq\n'
        s += 'from cirq.devices import GridQubit\n'
        s += 'from cirq.ops import X, Y, Z, H, Rx, Ry, Rz\n'
        s += 'from cirq.ops import CNOT, CZ, SWAP\n\n\n'
        s += 'ckt = cirq.Circuit()\n'
        for var_num in self.all_var_nums:
            vname = self.vprefix + str(var_num)
            s += vname
            s += ' = cirq.Symbol("'
            s += vname
            s += '")\n'
        s = s.strip()
        self.qasm_out.write(s)
        self.qasm_out.write('\n\n')

        if self.write_qubiter_files:
            lines = s.split('\n')
            for line in lines:
                self.qbtr_wr.write_NOTA(line)
        self.qbtr_wr.write_NOTA('')

    def write_ending(self):
        """
        Writes Cirq ending statements after calls to use_ methods for gates.

        Returns
        -------
        None

        """

        pass

    def bit2str(self, bit_pos):
        """
        Returns a string of form 'GridQubit(' ... ')'

        Parameters
        ----------
        bit_pos : int

        Returns
        -------
        str

        """

        row, col = self.lattice.one2two(bit_pos)
        return 'GridQubit(' + str(row) + ', ' + str(col) + ')'

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
        self.qasm_out.write("ckt.append(H(" +
                            self.bit2str(tar_bit_pos) + "))\n")
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

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        If called for a controlled phase, this function will halt execution 
        of program. If it's just a global phase with no controls, 
        the function will comment the phase out in the output files (Cirq 
        and output Qubiter English and Picture files.) and move on to the 
        next line. 

        Parameters
        ----------
        angle_rads : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        def degs_str(x):
            return x if isinstance(x, str) else str(x*180/np.pi)

        if controls.bit_pos_to_kind:
            assert False, "No PHAS lines with controls allowed"
        else:
            bla_str = 'PHAS\t' + degs_str(angle_rads) +\
                      '\tAT\t' + str(tar_bit_pos)
            self.qasm_out.write("# " + bla_str + "\n")
            if self.write_qubiter_files:
                self.qbtr_wr.write_NOTA(bla_str)


    def use_P_PH(self, projection_bit, angle_rads, tar_bit_pos, controls):
        """
        0Writes line in Cirq file corresponding to an English file line of
        type: P0PH or P1PH with 0 or 1 controls.

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
        assert not self.strict_mode
        num_trols = len(controls.bit_pos)
        assert num_trols in [0, 1]

        line_str = "ckt.append("
        if num_trols == 0:
            assert projection_bit == 1, \
                'exp(j*P_0*alp) not implemented in Cirq. ' +\
                'You  can use exp(j*P_0*alp)=sig_x*exp(j*P_1*alp)*sig_x'
            line_str += 'Z'

        else:  # num_trols == 1
            trol_bit_pos = controls.bit_pos[0]
            trol_type = controls.bit_pos_to_kind[trol_bit_pos]
            second_bit = 1 if trol_type else 0
            if projection_bit == 0:
                if second_bit == 0:
                    line_str += 'CPHASE00('
                    assert False, 'this gate not implemented in Cirq'
                else:
                    line_str += 'CPHASE01('
                    assert False, 'this gate not implemented in Cirq'
            elif projection_bit == 1:
                if second_bit == 0:
                    line_str += 'CPHASE10('
                    assert False, 'this gate not implemented in Cirq'
                else:
                    line_str += 'CZ'
            else:
                assert False
        if isinstance(angle_rads, float):
            cirq_turns = angle_rads/np.pi
        elif isinstance(angle_rads, str):
            cirq_turns = self.new_var_name(angle_rads, "/np.pi")
        else:
            assert False
        line_str += "**" + str(cirq_turns)
        if num_trols == 0:
            line_str += '.on(' + self.bit2str(tar_bit_pos)
        else:  # num_trols == 1:
            line_str += '.on(' + self.bit2str(controls.bit_pos[0])
            line_str += ', ' + self.bit2str(tar_bit_pos)
        line_str += "))\n"
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            if projection_bit == 0:
                u2_fun = OneBitGates.P_0_phase_fac
            elif projection_bit == 1:
                u2_fun = OneBitGates.P_1_phase_fac
            else:
                assert False

            self.qbtr_wr.write_controlled_one_bit_gate(
                tar_bit_pos, controls, u2_fun, [angle_rads])

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

    def use_ROT(self, axis, angle_rads, tar_bit_pos, controls):
        """
        Writes line in Cirq file corresponding to an English file line
        of type: ROTX, ROTY or ROTZ with no controls.

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
        assert len(controls.bit_pos) == 0

        line_str = "ckt.append("
        if axis == 1:
            line_str += "Rx("
        elif axis == 2:
            line_str += "Ry("
        elif axis == 3:
            line_str += "Rz("
        else:
            assert False

        if isinstance(angle_rads, float):
            cirq_rads = angle_rads*(-2)
        elif isinstance(angle_rads, str):
            cirq_rads = self.new_var_name(angle_rads, "*(-2)")
        else:
            assert False

        line_str += 'rads=' + str(cirq_rads) + ').on('
        line_str += self.bit2str(tar_bit_pos) + "))\n"
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [angle_rads, axis])

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        Writes line in Cirq file corresponding to an English file line
        of type: ROTN with no controls.

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
        assert len(controls.bit_pos) == 0

        rad_ang_list = [angle_x_rads, angle_y_rads, angle_z_rads]
        assert ug.all_floats(rad_ang_list), \
            "With Cirq, ROTN with any of its 3 angles variable is " +\
            "not allowed. Workaround: can use 3 rotations of type " +\
            "Rx, Ry or Rz with variable angles."
        arr = OneBitGates.rot(*rad_ang_list)
        delta, left_rads, center_rads, right_rads = \
            UnitaryMat.u2_zyz_decomp(arr)

        end_str = ').on(' + self.bit2str(tar_bit_pos) + '))\n'

        line_str = 'ckt.append(Rz(rads=' + str(-2*right_rads) + end_str
        self.qasm_out.write(line_str)

        line_str = 'ckt.append(Ry(rads=' + str(-2*center_rads) + end_str
        self.qasm_out.write(line_str)

        line_str = 'ckt.append(Rz(rads=' + str(-2*left_rads) + end_str
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
            line_str = "ckt.append("
            if axis == 1:
                line_str += "X("
            elif axis == 2:
                line_str += "Y("
            elif axis == 3:
                line_str += "Z("
            else:
                assert False
            line_str += self.bit2str(tar_bit_pos) + "))\n"
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
                line_str = 'ckt.append(CNOT('
                line_str += self.bit2str(trol_pos) + ', '
                line_str += self.bit2str(tar_pos) + '))\n'
                self.qasm_out.write(line_str)
                if self.write_qubiter_files:
                    self.qbtr_wr.write_cnot(trol_pos, tar_pos)
            else:
                assert False, "Forbidden CNOT detected: " \
                    + str(trol_pos) + "->" + str(tar_pos) \
                    + " in line " + str(self.line_count) \
                    + ". Use class ForbiddenCNotExpander " \
                    + "before attempting translation to Cirq."

    def use_SWAP(self, bit1, bit2, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: SWAP with no controls.


        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        assert not self.strict_mode
        assert len(controls.bit_pos) == 0

        line_str = 'ckt.append(SWAP('
        line_str += self.bit2str(bit1)
        line_str += ", " + self.bit2str(bit2) + "))\n"
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_bit_swap(bit1, bit2, controls)

if __name__ == "__main__":

    def main():
        file_prefix = "../io_folder/qbtr2google_test"
        qasm_name = 'GooCirq'
        num_bits = 5
        c_to_tars = 'do_fill'  # filled by constructor
        Qubiter_to_GoogleCirq(file_prefix, num_bits, qasm_name=qasm_name,
                              c_to_tars=c_to_tars, write_qubiter_files=True)

    main()
