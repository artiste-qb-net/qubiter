from qubiter.device_specific.Qubiter_to_AnyQasm import *
import qubiter.device_specific.chip_couplings_rigetti as cc
import qubiter.utilities_gen as utg


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
    def __init__(self, file_prefix, num_qbits, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_qbits : int

        Returns
        -------

        """
        Qubiter_to_AnyQasm.__init__(self, file_prefix, num_qbits, **kwargs)

    def write_prelude(self):
        """
        Writes PyQuil opening statements before calls to use_ methods for
        gates.

        Returns
        -------
        None

        """
        s = 'from pyquil import Program, get_qc\n'
        s += 'from pyquil.gates import *\n\n\n'
        s += 'pg = Program()\n'
        for var_num in self.all_var_nums:
            vname = self.vprefix + str(var_num)
            s += vname
            s += ' = pg.declare("'
            s += vname
            s += '", memory_type="REAL")\n'
        s += 'ro = pg.declare("ro", memory_type="BIT", memory_size='
        s += str(self.num_qbits)
        s += ')'
        self.write(s + '\n')

    def write_ending(self):
        """
        Writes PyQuil ending statements after calls to use_ methods for gates.

        Returns
        -------
        None

        """

        self.write('\n')

        s = ''
        for k in range(self.num_qbits):
            s += 'pg.MEASURE('
            s += str(k)
            s += ', ro['
            s += str(k)
            s += '])\n'
        s = s.strip()
        self.write(s)

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
        self.aqasm_out.write("pg += H(" + str(tar_bit_pos) + ")\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                   OneQubitGate.had2)

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
        self.aqasm_out.write("# " + bla_str + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(bla_str)

    def use_PHAS(self, angle_rads, tar_bit_pos, controls):
        """
        If called for a controlled phase, this function will halt execution
        of program. If it's just a global phase with no controls,
        the function will comment the phase out in the output files (PyQuil
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
        Writes line in PyQuil file corresponding to an English file line of
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

        line_str = "pg += "
        if num_trols == 0:
            assert projection_bit == 1, \
                'exp(j*P_0*alp) not implemented in pyquil. ' +\
                'You  can use exp(j*P_0*alp)=sig_x*exp(j*P_1*alp)*sig_x'
            line_str += 'PHASE('

        else:  # num_trols == 1
            trol_bit_pos = controls.bit_pos[0]
            trol_type = controls.bit_pos_to_kind[trol_bit_pos]
            second_bit = 1 if trol_type else 0
            if projection_bit == 0:
                if second_bit == 0:
                    line_str += 'CPHASE00('
                else:
                    line_str += 'CPHASE01('
            elif projection_bit == 1:
                if second_bit == 0:
                    line_str += 'CPHASE10('
                else:
                    line_str += 'CPHASE('
            else:
                assert False
        if isinstance(angle_rads, float):
            quil_rads = angle_rads
        elif isinstance(angle_rads, str):
            quil_rads = self.new_var_name(angle_rads)
        else:
            assert False
        line_str += str(quil_rads)
        if num_trols == 0:
            line_str += ', ' + str(tar_bit_pos)
        else:  # num_trols == 1:
            line_str += ', ' + str(controls.bit_pos[0])
            line_str += ', ' + str(tar_bit_pos)
        line_str += ")\n"
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
        self.aqasm_out.write("# " + str1 + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(str1)

    def use_ROTA(self, axis, angle_rads, tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
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

        line_str = "pg += "
        if axis == 1:
            line_str += "RX("
        elif axis == 2:
            line_str += "RY("
        elif axis == 3:
            line_str += "RZ("
        else:
            assert False
      
        if isinstance(angle_rads, float):
            quil_rads = angle_rads*(-2)
        elif isinstance(angle_rads, str):
            quil_rads = self.new_var_name(angle_rads, "*(-2)")
        else:
            assert False

        line_str += str(quil_rads) + ', '
        line_str += str(tar_bit_pos) + ")\n"
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos, controls,
                               OneQubitGate.rot_ax, [angle_rads, axis])

    def use_ROTN(self, angle_x_rads, angle_y_rads, angle_z_rads,
                tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
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
            "With Pyquil, ROTN with any of its 3 angles variable is " +\
            "not allowed. Workaround: can use 3 rotations of type " +\
            "Rx, Ry or Rz with variable angles."
        arr = OneQubitGate.rot(*rad_ang_list)
        delta, left_rads, center_rads, right_rads = \
            UnitaryMat.u2_zyz_decomp(arr)
                  
        end_str = ', ' + str(tar_bit_pos) + ')\n'

        line_str = 'pg += RZ(' + str(-2*right_rads) + end_str
        self.aqasm_out.write(line_str)

        line_str = 'pg += RY(' + str(-2*center_rads) + end_str
        self.aqasm_out.write(line_str)

        line_str = 'pg += RZ(' + str(-2*left_rads) + end_str
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos,
                    controls, OneQubitGate.rot_ax, [right_rads, 3])
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos,
                    controls, OneQubitGate.rot_ax, [center_rads, 2])
            self.qbtr_wr.write_controlled_one_qbit_gate(tar_bit_pos,
                    controls, OneQubitGate.rot_ax, [left_rads, 3])

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
            line_str = "pg += "
            if axis == 1:
                line_str += "X("
            elif axis == 2:
                line_str += "Y("
            elif axis == 3:
                line_str += "Z("
            else:
                assert False
            line_str += str(tar_bit_pos) + ")\n"
            self.aqasm_out.write(line_str)
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
                line_str = 'pg += CNOT('
                line_str += str(trol_pos) + ', '
                line_str += str(tar_pos) + ')\n'
                self.aqasm_out.write(line_str)
                if self.write_qubiter_files:
                    self.qbtr_wr.write_cnot(trol_pos, tar_pos)
            else:
                assert False, "Forbidden CNOT detected: " \
                    + str(trol_pos) + "->" + str(tar_pos) \
                    + " in line " + str(self.line_count) \
                    + ". Use class ForbiddenCNotExpander " \
                    + "before attempting translation to PyQuil."

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

        line_str = 'pg += SWAP(' + \
                   str(bit1) + ", " + str(bit2) + ")\n"
        self.aqasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_qbit_swap(bit1, bit2, controls)

    def use_U_2_(self, rads0, rads1, rads2, rads3,
                tar_bit_pos, controls):
        """
        Writes line in PyQuil file corresponding to an English file line
        of type: U_2_ with no controls.

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
        assert len(controls.bit_pos) == 0
        # we drop global phase
        self.use_ROTN(rads1, rads2, rads3, tar_bit_pos, controls)


if __name__ == "__main__":
    import qubiter.device_specific.chip_couplings_rigetti as rig

    def main():
        file_prefix = "qbtr2rigetti_test"
        aqasm_name = 'RigPyQuil'
        num_qbits = 6
        c_to_tars = rig.rigetti20_c_to_tars
        Qubiter_to_RigettiPyQuil(file_prefix, num_qbits, aqasm_name=aqasm_name,
                c_to_tars=c_to_tars, write_qubiter_files=True)

    main()
