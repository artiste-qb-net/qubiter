from qubiter.device_specific.Qubiter_to_AnyQasm import *
import qubiter.utilities_gen as utg


class Qubiter_to_IBMqasm(Qubiter_to_AnyQasm):
    """
    See docstring of parent class Qubiter_to_AnyQasm.

    If input c_to_tars = None, all CNOTs allowed.

    References
    ----------
    1. https://github.com/Qiskit

    Attributes
    ----------

    """
    def __init__(self, file_prefix, num_qbits, **kwargs):
        """
        Constructor

        Parameters
        ----------

        Returns
        -------

        """
        Qubiter_to_AnyQasm.__init__(self, file_prefix, num_qbits, **kwargs)

    def write_prelude(self):
        """
        Writes IBM qasm opening statements before calls to use_ methods for
        gates.

        Returns
        -------
        None

        """
        s = "import numpy as np\n"
        s += 'from qiskit import QuantumCircuit\n'
        s += 'from qiskit import ClassicalRegister, QuantumRegister\n'
        s += 'from qiskit import execute\n\n\n'
        s += "q = QuantumRegister(" + str(self.num_qbits) + ", 'q')\n"
        s += 'ckt = QuantumCircuit(q)'
        self.write(s + '\n')

    def write_ending(self):
        """
        Writes IBM qasm ending statements after calls to use_ methods for
        gates.

        Returns
        -------
        None

        """
        pass

    @staticmethod
    def qasm_line_for_rot(arr, tar_bit_pos):
        """
        This function returns a string for an IBM qasm file line for a one
        qubit rotation.

        Parameters
        ----------
        arr : np.array
            the matrix of the rotation
        tar_bit_pos : int
            target bit position at which rotation occurs.

        Returns
        -------
        str

        """
        # U(\theta,\phi,\lambda) := R_z(\phi)R_y(\theta)R_z(\lambda)
        # R_y(\theta)=\mathrm{exp}(-i\theta Y/2)
        # R_z(\phi)=\mathrm{exp}(-i\phi Z/2)

        # excerpt from qelib1.inc

        # // 3-parameter 2-pulse single qubit gate
        # gate u3(theta,phi,lambda) q { U(theta,phi,lambda) q; }
        # // 2-parameter 1-pulse single qubit gate
        # gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }
        # // 1-parameter 0-pulse single qubit gate
        # gate u1(lambda) q { U(0,0,lambda) q; }
        # // controlled-NOT
        # gate cx c,t { CX c,t; }
        # // idle gate (identity)
        # gate id a { U(0,0,0) a; }

        delta, left_rads, center_rads, right_rads = \
            UnitaryMat.u2_zyz_decomp(arr)
        phi = -2*left_rads
        theta = -2*center_rads
        lam = -2*right_rads
        if abs(phi) < 1e-6 and abs(theta) < 1e-6:
            line_str = "ckt.u1(" + str(lam)
        elif abs(theta - np.pi/2) < 1e-6:
            line_str = "ckt.u2(" + str(phi) + ", " + str(lam)
        else:
            line_str = "ckt.u3(" + str(theta) + ", " + \
                       str(phi) + ", " + str(lam)
        line_str += "), q[" + str(tar_bit_pos) + "])\n"
        return line_str

    def use_HAD2(self, tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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
        self.aqasm_out.write("ckt.h(q[" + str(tar_bit_pos) + "])\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                   OneQubitGate.had2)

    def use_NOTA(self, bla_str):
        """
        Writes line in IBM qasm file corresponding to an English file line
        of type: NOTA

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        self.aqasm_out.write("# " + bla_str + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(bla_str)

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        If called for a controlled phase, this function will halt execution
        of program. If it's just a global phase with no controls,
        the function will comment the phase out in the output files (IBM
        qasm and output Qubiter English and Picture files.) and move on to
        the next line.

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
            self.aqasm_out.write("# " + bla_str + "\n")
            if self.write_qubiter_files:
                self.qbtr_wr.write_NOTA(bla_str)

    def use_P_PH(self, projection_bit, angle_rads, tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
        of type: P0PH or P1PH with 0 or 1 controls.

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
        assert isinstance(angle_rads, float),\
            'At present, IBM qasm does not support variable angles'
        assert not self.strict_mode
        num_trols = len(controls.bit_pos)
        assert num_trols in [0, 1]

        line_str = "ckt."
        if num_trols == 0:
            assert projection_bit == 1, \
                'exp(j*P_0*alp) not implemented in IBM qasm. ' +\
                'You  can use exp(j*P_0*alp)=sig_x*exp(j*P_1*alp)*sig_x'
            line_str += 'u1('

        else:  # num_trols == 1
            trol_bit_pos = controls.bit_pos[0]
            trol_type = controls.bit_pos_to_kind[trol_bit_pos]
            second_bit = 1 if trol_type else 0
            if projection_bit == 0:
                if second_bit == 0:
                    line_str += 'CPHASE00('
                    assert False, 'this gate not implemented in IBM qasm'
                else:
                    line_str += 'CPHASE01('
                    assert False, 'this gate not implemented in IBM qasm'
            elif projection_bit == 1:
                if second_bit == 0:
                    line_str += 'CPHASE10('
                    assert False, 'this gate not implemented in IBM qasm'
                else:
                    line_str += 'cu1('
            else:
                assert False
        if isinstance(angle_rads, float):
            ibm_rads = angle_rads
        elif isinstance(angle_rads, str):
            ibm_rads = self.new_var_name(angle_rads)
        else:
            assert False
        line_str += str(ibm_rads)
        if num_trols == 0:
            line_str += ', q[' + str(tar_bit_pos)
        else:  # num_trols == 1:
            line_str += ', q[' + str(controls.bit_pos[0])
            line_str += '], q[' + str(tar_bit_pos)
        line_str += "])\n"
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            if projection_bit == 0:
                u2_fun = OneQubitGate.P_0_phase_fac
            elif projection_bit == 1:
                u2_fun = OneQubitGate.P_1_phase_fac
            else:
                assert False

            self.qbtr_wr.write_controlled_one_qbit_gate(
                tar_bit_pos, controls, u2_fun, [angle_rads])

    def use_PRINT(self, style, line_num):
        """
        Writes line in IBM qasm file corresponding to an English file line
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
        self.aqasm_out.write("# " + str1 + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(str1)

    def use_ROTA(self, axis, angle_rads, tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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
        assert isinstance(angle_rads, float),\
            'At present, IBM qasm does not support variable angles'
        assert len(controls.bit_pos) == 0

        arr = OneQubitGate.rot_ax(angle_rads, axis)
        line_str = Qubiter_to_IBMqasm.qasm_line_for_rot(arr, tar_bit_pos)
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                               OneQubitGate.rot_ax, [angle_rads, axis])

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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
        assert utg.all_floats(rad_ang_list), \
            'At present, IBM qasm does not support variable angles'
        arr = OneQubitGate.rot(*rad_ang_list)
        line_str = Qubiter_to_IBMqasm.qasm_line_for_rot(arr, tar_bit_pos)
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                               OneQubitGate.rot, rad_ang_list)

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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
            if axis == 1:
                prefix = "ckt.x(q["
            elif axis == 2:
                prefix = "ckt.y(q["
            elif axis == 3:
                prefix = "ckt.z(q["
            else:
                assert False
            self.aqasm_out.write(prefix + str(tar_bit_pos) + "])\n")
            if self.write_qubiter_files:
                if axis == 1:
                    u2_fun = OneQubitGate.sigx
                elif axis == 2:
                    u2_fun = OneQubitGate.sigy
                elif axis == 3:
                    u2_fun = OneQubitGate.sigz
                else:
                    assert False

                self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos,
                                    controls, u2_fun)
        else:  # num_trols == 1
            tar_pos = tar_bit_pos
            trol_pos = controls.bit_pos[0]
            if not self.c_to_tars or tar_pos in self.c_to_tars[trol_pos]:
                self.aqasm_out.write("ckt.cx(q[" + str(trol_pos) + "], "
                                    "q[" + str(tar_pos) + "])\n")
                if self.write_qubiter_files:
                    self.qbtr_wr.write_cnot(trol_pos, tar_pos)
            else:
                assert False, "Forbidden CNOT detected: " \
                    + str(trol_pos) + "->" + str(tar_pos) \
                    + " in line " + str(self.line_count) \
                    + ". Use class ForbiddenCNotExpander " \
                    + "before attempting translation to IBM qasm."

    def use_SWAP(self, bit1, bit2, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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

        line_str = 'ckt.swap(q[' + \
                   str(bit1) + "], q[" + str(bit2) + "])\n"
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_qbit_swap(bit1, bit2, controls)


if __name__ == "__main__":
    import qubiter.device_specific.chip_couplings_ibm as ibm

    def main():
        file_prefix = "qbtr2ibm_test"
        aqasm_name = 'IBMqasm'
        num_qbits = 5
        c_to_tars = ibm.ibmq5YorktownTenerife_c_to_tars
        Qubiter_to_IBMqasm(file_prefix, num_qbits, aqasm_name=aqasm_name,
                           c_to_tars=c_to_tars, write_qubiter_files=True)

    main()
