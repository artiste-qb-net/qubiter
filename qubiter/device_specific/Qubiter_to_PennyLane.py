from qubiter.device_specific.Qubiter_to_AnyQasm import *
import qubiter.utilities_gen as utg


class Qubiter_to_PennyLane(Qubiter_to_AnyQasm):
    """
    See docstring of parent class Qubiter_to_AnyQasm

    If input c_to_tars = None, all CNOTs and CZs are allowed.

    References
    ----------
    1. https://github.com/XanaduAI/pennylane

    Attributes
    ----------
    fun_defs_path : str
        path to a py file that defines all the distinct functions used in
        functional placeholders in the English file being translated to
        PennyLane.
    indentation : int
        internal int that keeps track of indentation. Starts at 0
    qnode_name : str
        name to be given to qnode in output PennyLane file. The whole
        English File will be in the body of a SINGLE qnode.
    rotn_has_been_defined : bool
        If the English file uses ROTN **with** placeholder variables,
        then the first time, and only the 1st time, that ROTN is used,
        the function rot() defining ROTN as an array is included in the body
        of the qnode. This internal boolean flag helps to insure that the
        def of rot( ) is included only once in the body of the qnode.


    """
    def __init__(self, file_prefix, num_qbits, qnode_name='qnode',
                 fun_defs_path=None, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_qbits : int
        qnode_name : str
        fun_defs_path : str
        rotn_has_been_defined : bool

        Returns
        -------

        """
        self.qnode_name = qnode_name
        self.fun_defs_path = fun_defs_path
        self.rotn_has_been_defined = False
        self.indentation = 0
        Qubiter_to_AnyQasm.__init__(self, file_prefix, num_qbits,
                                    aqasm_ftype='py', **kwargs)

    def write_prelude(self):
        """
        Writes PennyL opening statements before calls to use_ methods for
        gates.

        Returns
        -------
        None

        """
        vars_str = ""
        for num in self.all_var_nums:
            vars_str += self.vprefix + str(num) + ', '
        vars_str = vars_str[:-2]

        s = 'import pennylane as qml\n\n\n'
        s += 'def ' + self.qnode_name + '(' + vars_str + '):\n'
        self.indentation += 4
        s += ' '*self.indentation
        s += '# distinct fun names in functional placeholders=\n'
        s += ' '*self.indentation + '# '
        s += str(self.all_fun_names) + '\n'
        if self.fun_defs_path:
            with open(utg.preface(self.fun_defs_path), 'r') as fi:
                fi_lines = fi.readlines()
            for line in fi_lines:
                s += ' '*self.indentation + line
        s = s.rstrip()
        self.write(s)

    def write_ending(self):
        """
        Writes PennyL ending statements after calls to use_ methods for gates.

        Returns
        -------
        None

        """
        s = ' '*self.indentation + 'return qml.expval.Hermitian(hamil)'
        self.write(s)

    def use_HAD2(self, tar_bit_pos, controls):
        """
        Writes line in PennyL file corresponding to an English file line
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
        self.aqasm_out.write(' '*self.indentation +
            "qml.Hadamard(" + str(tar_bit_pos) + ")\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_H(tar_bit_pos)

    def use_NOTA(self, bla_str):
        """
        Writes line in PennyL file corresponding to an English file line
        of type: NOTA

        Parameters
        ----------
        bla_str : str

        Returns
        -------
        None

        """
        self.aqasm_out.write(' '*self.indentation +
                            "# " + bla_str + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(bla_str)

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        If called for a controlled phase, this function will halt execution
        of program. If it's just a global phase with no controls,
        the function will comment the phase out in the output files (PennyL
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
            self.aqasm_out.write("# " + bla_str + "\n")
            if self.write_qubiter_files:
                self.qbtr_wr.write_NOTA(bla_str)

    def use_P_PH(self, projection_bit, angle_rads, tar_bit_pos, controls):
        """
        Writes line in PennyL file corresponding to an English file line of
        type: P1PH with 0 controls.

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
        assert num_trols == 0

        line_str = ' '*self.indentation
        assert projection_bit == 1, \
            'exp(j*P_0*alp) not implemented in PennyL. ' +\
            'You  can use exp(j*P_0*alp)=sig_x*exp(j*P_1*alp)*sig_x'
        line_str += 'qml.PhaseShift('
        if isinstance(angle_rads, float):
            penny_rads = angle_rads
        elif isinstance(angle_rads, str):
            penny_rads = self.new_var_name(angle_rads)
        else:
            assert False
        line_str += str(penny_rads)
        line_str += ', wires=' + str(tar_bit_pos)
        line_str += ")\n"
        self.aqasm_out.write(line_str)

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
        Writes line in PennyL file corresponding to an English file line
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
        self.aqasm_out.write(' '*self.indentation +
                            "# " + str1 + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(str1)

    def use_ROTA(self, axis, angle_rads, tar_bit_pos, controls):
        """
        Writes line in PennyL file corresponding to an English file line
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

        line_str = ' '*self.indentation + "qml."
        if axis == 1:
            line_str += "RX("
        elif axis == 2:
            line_str += "RY("
        elif axis == 3:
            line_str += "RZ("
        else:
            assert False

        if isinstance(angle_rads, float):
            penny_rads = angle_rads*(-2)
        elif isinstance(angle_rads, str):
            penny_rads = self.new_var_name(angle_rads, "*(-2)")
        else:
            assert False

        line_str += str(penny_rads) + ', wires='
        line_str += str(tar_bit_pos) + ")\n"
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [angle_rads, axis])

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        Writes line in PennyL file corresponding to an English file line
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
        if not self.rotn_has_been_defined:
            import inspect
            # gives tuple of list so need 0th component
            # first line is @staticmethod, omit it
            lines = inspect.getsourcelines(OneBitGates.rot)[0][1:]
            # print(",999999999999999,0", lines)
            s = ''
            for line in lines:
                # print('9kkkkkk', line)
                line = line[4:] if line != '\n' else line
                s += ' '*self.indentation + str(line)
            self.rotn_has_been_defined = True
            self.aqasm_out.write(s)
            s = s.rstrip()
            if self.write_qubiter_files:
                lines = s.split('\n')
                for line in lines:
                    self.qbtr_wr.write_NOTA(line)

        penny_rad_ang_list = []
        for rads in rad_ang_list:
            if isinstance(rads, float):
                penny_rad_ang_list.append(rads*(-2))
            elif isinstance(rads, str):
                penny_rad_ang_list.append(self.new_var_name(rads, "*(-2)"))
            else:
                assert False

        line_str = ' '*self.indentation + 'qml.QubitUnitary(rot('
        for rads in penny_rad_ang_list:
            line_str += str(rads) + ', '
        line_str = line_str[:-2] + '), '
        line_str += 'wires=' + str(tar_bit_pos) + ')\n'
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos,
                    controls, OneBitGates.rot, rad_ang_list)

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        Writes line in PennyL file corresponding to an English file line of
        type: SIGX, SIGY or SIGZ with no controls, or else SIGX with one
        True control (i.e., simple CNOT), or else SIGZ with one True control.

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
            assert axis in [1, 3]
            assert controls.bit_pos_to_kind[controls.bit_pos[0]] == True
        if axis == 1:
            u2_fun = OneBitGates.sigx
        elif axis == 2:
            u2_fun = OneBitGates.sigy
        elif axis == 3:
            u2_fun = OneBitGates.sigz
        else:
            assert False
        if num_trols == 0:
            line_str = ' '*self.indentation + "qml."
            if axis == 1:
                line_str += "PauliX("
            elif axis == 2:
                line_str += "PauliY("
            elif axis == 3:
                line_str += "PauliZ("
            else:
                assert False
            line_str += str(tar_bit_pos) + ")\n"
            self.aqasm_out.write(line_str)
            if self.write_qubiter_files:
                self.qbtr_wr.write_controlled_one_bit_gate(
                    tar_bit_pos, controls, u2_fun)
        else:  # num_trols == 1
            tar_pos = tar_bit_pos
            trol_pos = controls.bit_pos[0]
            if not self.c_to_tars or tar_pos in self.c_to_tars[trol_pos]:
                line_str = ' '*self.indentation
                if axis == 1:
                    line_str += 'qml.CNOT(wires=['
                elif axis == 3:
                    line_str += 'qml.CZ(wires=['
                else:
                    assert False, 'unsupported axis'
                line_str += str(trol_pos) + ', '
                line_str += str(tar_pos) + '])\n'
                self.aqasm_out.write(line_str)
                if self.write_qubiter_files:
                    self.qbtr_wr.write_controlled_one_bit_gate(
                        tar_bit_pos, controls, u2_fun)
            else:
                assert False, "Forbidden CNOT detected: " \
                    + str(trol_pos) + "->" + str(tar_pos) \
                    + " in line " + str(self.line_count) \
                    + ". Use class ForbiddenCNotExpander " \
                    + "before attempting translation to PennyL."

    def use_SWAP(self, bit1, bit2, controls):
        """
        Writes line in PennyL file corresponding to an English file line
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

        line_str = ' '*self.indentation + \
                   'qml.SWAP(wires=[' + \
                   str(bit1) + ", " + str(bit2) + "])\n"
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_bit_swap(bit1, bit2, controls)


if __name__ == "__main__":

    def main1():
        file_prefix = "qbtr2penny_test1"
        num_qbits = 3
        emb = CktEmbedder(num_qbits, num_qbits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_H(0)
        wr.write_X(1)
        wr.write_Y(1)
        wr.write_Z(1)
        wr.write_cnot(0, 1)
        wr.write_cz(0, 1)
        wr.write_bit_swap(1, 0)
        wr.write_Rx(2, rads=np.pi)
        wr.write_Ry(2, rads=np.pi)
        wr.write_Rz(2, rads=np.pi)
        wr.write_one_bit_gate(1, OneBitGates.P_1_phase_fac, [np.pi])
        wr.write_Rn(0, rads_list=[np.pi, np.pi, np.pi])
        wr.close_files()

        qasm_name = 'PennyL'
        qnode_name = 'Turing'
        Qubiter_to_PennyLane(file_prefix, num_qbits,
                qnode_name,
                aqasm_name=qasm_name,
                write_qubiter_files=True)

    def main2():
        file_prefix = "qbtr2penny_test2"
        num_qbits = 4
        emb = CktEmbedder(num_qbits, num_qbits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_Rx(2, rads=np.pi/7)
        wr.write_Rx(1, rads='#2*.5')
        wr.write_Rn(3, rads_list=['#1', '-#1*3', '#2'])
        wr.write_Rx(1, rads='-my_fun#2#1')
        wr.write_cnot(2, 3)
        wr.close_files()

        aqasm_name = 'PennyL'
        fun_defs_path = 'qbtr2penny_test2_fun_defs.py'
        qnode_name = 'Feynman'
        Qubiter_to_PennyLane(file_prefix, num_qbits,
                qnode_name,
                fun_defs_path,
                aqasm_name=aqasm_name,
                write_qubiter_files=True)
    main1()
    main2()
