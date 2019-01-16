from device_specific.Qubiter_to_AnyQasm import *
import device_specific.chip_couplings_rigetti as cc
import utilities_gen as ug


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
    def __init__(self, file_prefix, num_bits, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int

        Returns
        -------

        """
        Qubiter_to_AnyQasm.__init__(self, file_prefix, num_bits, **kwargs)

    def write_prelude(self):
        """
        Writes PyQuil opening statements before calls to use_ methods for
        gates.

        Returns
        -------
        None

        """

        s = ''
        s += 'from from pyquil import Program, get_qc\n'
        s += 'from pyquil.gates import *\n\n\n'
        s += 'pg = Program()\n'
        for var_num in self.var_nums_list:
            s += 'pg.declare("'
            s += self.vname
            s += str(var_num)
            s += '", memory_type="REAL")\n'
        s += 'pg.declare("ro", memory_type="BIT", memory_size='
        s += str(self.num_bits)
        s += ')'
        self.qasm_out.write(s)
        self.qasm_out.write('\n')

        if self.write_qubiter_files:
            lines = s.split('\n')
            for line in lines:
                self.qbtr_wr.write_NOTA(line)
        self.qbtr_wr.write_NOTA('')

    def write_ending(self):
        """
        Writes PyQuil ending statements after calls to use_ methods for gates.

        Returns
        -------
        None

        """

        s = ''
        for k in range(self.num_bits):
            s += 'pg.MEASURE('
            s += str(k)
            s += 'ro['
            s += str(k)
            s += '])\n'
        self.qasm_out.write(s)

        if self.write_qubiter_files:
            lines = s.split('\n')
            for line in lines:
                self.qbtr_wr.write_NOTA(line)

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

    @staticmethod
    def use_ROT_rads_map(rads):
        """
        Asserts that rads is a float. This method maps a qubiter rads to a 
        pyquil rads 
        
        Parameters
        ----------
        rads : float

        Returns
        -------
        float

        """
        assert isinstance(rads, float)
        return 2*rads
            
    def use_ROT(self, axis, angle_rads, tar_bit_pos, controls):
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

        line_str = "pg.inst("
        if axis == 1:
            line_str += "RX("
        elif axis == 2:
            line_str += "RY("
        elif axis == 3:
            line_str += "RZ("
        else:
            assert False
      
        if isinstance(angle_rads, float):
            quil_rads = \
                Qubiter_to_RigettiPyQuil.use_ROT_rads_map(angle_rads)
        elif isinstance(angle_rads, str):
            if angle_rads[0] == '#':
                quil_rads = self.vname + angle_rads[1:]
            else:  # starts with -#
                quil_rads = '-' + self.vname + angle_rads[2:]
        else:
            assert False

        line_str += str(quil_rads) + ', '
        line_str += str(tar_bit_pos) + "))\n"
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [angle_rads, axis])
            
    @staticmethod
    def use_ROTN_rads_map(rad_ang_list):
        """
        Asserts that rad_ang_list is list of floats. This method maps a 
        qubiter rad_ang_list to a pyquil rad_ang_list 
        
        Parameters
        ----------
        rad_ang_list : list[float]

        Returns
        -------
        list[float]

        """

        assert ug.all_floats(rad_ang_list)

        arr = OneBitGates.rot(*rad_ang_list)
        delta, left_rads, center_rads, right_rads = \
            UnitaryMat.u2_zyz_decomp(arr)
        quil_rads_L = -2*left_rads
        quil_rads_C = -2*center_rads
        quil_rads_R = -2*right_rads
        return quil_rads_L, quil_rads_C, quil_rads_R


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
        if ug.all_floats(rad_ang_list):
            quil_rads_L, quil_rads_C, quil_rads_R = \
                Qubiter_to_RigettiPyQuil.use_ROTN_rads_map(rad_ang_list)
        elif ug.all_strings(rad_ang_list):
            quil_rads_L, quil_rads_C, quil_rads_R = \
                    ['rads_' + x[1:] for x in rad_ang_list]
        else:
            assert False, "For PyQuil, angles of ROTN must " \
                          "all be either floats or variables"
                  
        end_str = ', ' + str(tar_bit_pos) + '))\n'

        line_str = 'pg.inst(RZ(' + str(quil_rads_R) + end_str
        self.qasm_out.write(line_str)

        line_str = 'pg.inst(RY(' + str(quil_rads_C) + end_str
        self.qasm_out.write(line_str)

        line_str = 'pg.inst(RZ(' + str(quil_rads_L) + end_str
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [-quil_rads_R/2, 3])
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [-quil_rads_C/2, 2])
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [-quil_rads_L/2, 3])

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
        Qubiter_to_RigettiPyQuil(file_prefix, num_bits, qasm_name=qasm_name,
                c_to_tars=c_to_tars, write_qubiter_files=True)

    main()
