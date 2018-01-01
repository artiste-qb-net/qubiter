from Controls import *
from SEO_reader import *
from SEO_writer import *
from quantum_CSD_compiler.UnitaryMat import *
from ForbiddenCNotExpander import *
from ChipCouplingsFitter import *
import utilities_gen as ut


class Qubiter_to_IBMqasm2(SEO_reader):
    """
    This class is a child of SEO_reader. It reads an input English file and
    writes an IBM qasm2 file that is a translation of the input English file
    into the IBM qasm2 language. If the option write_qubiter_files is set to
    True, this class will also write new English and Picture files that are
    in 1-1 onto line correspondence with the output IBM qasm file.

    The input English file that is read can only have lines of the following
    types or else the program will abort with an error message:

    1. single qubit rotations (HAD2, SIGX, SIGY, SIGZ, ROTX, ROTY, ROTZ or
    ROTN with no controls)

    2. simple CNOTs (SIGX with a single True control). Call them c->t=(c,
    t) if c is the control and t the target. (c, t) must be an item of the
    input 'c_to_t'.

    3. NOTA lines

    If you have an English file that contains lines that are more
    complicated than this (because, for example, they contain rotations with
    one or more controls attached, or because a CNOT is not allowed
    according to `c_to_t`), you can use the expander classes CGateExpander,
    DiagUnitaryExpander, MultiplexorExpander, and ForbiddenCNotExpander to
    expand the circuit to an equivalent albeit longer circuit that satisfies
    constraints 1, 2, 3.

    This class can handle an IBM chip with any number of qubits.

    This class halts execution if it encounters a CNOT that is disallowed
    according to the input `c_to_t`. `c_to_t` varies with IBM chip. Some
    `c_to_t`s are listed in the file `ibm_chip_couplings.py` found in same
    folder as this file.

    Footnote: QASM distinguishes between quantum registers qreg and
    classical registers creg. Qubiter does not use cregs because it uses the
    classical memory of your Linux PC instead. QASM has an intricate set of
    commands for measurements. Qubiter has a complete set of measurement
    commands too (see MEAS in Rosetta stone). The QASM and Qubiter
    measurement commands can obviously be translated into each other. We
    leave that part of the translation to a future version of this class.

    References
    ----------
    1. https://github.com/IBMResearch/python-sdk-quantum-experience
    2. https://github.com/IBMQuantum/QASM

    Attributes
    ----------
    c_to_t : tuple[tuple(int,int)]
        Pairs of qubits that are physically connected so they can be the two
        ends of an elementary CNOT. Order of qubits matters: first entry of
        tuple is control and second is target of a possible CNOT.
    qasm_out : _io.TextIOWrapper
        This output stream is used to write a qasm file based on the input
        English file.
    qbtr_wr : SEO_writer
        A SEO_writer object created iff write_qubiter_files is True.
    write_qubiter_files : bool
        The class always writes a qasm text file based on the input English
        file that is read. Iff this is True, the class also writes English
        and Picture files


    """
    def __init__(self, file_prefix, num_bits, c_to_t, verbose=False,
                 write_qubiter_files=False, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        c_to_t : tuple(tuple(int, int))
        verbose : bool
        write_qubiter_files : bool
        kwargs : dict[]

        Returns
        -------
        None

        """
        self.c_to_t = c_to_t
        self.targets = ForbiddenCNotExpander.get_targets_from_c_to_t(num_bits, c_to_t)

        self.write_qubiter_files = write_qubiter_files

        self.qasm_out = open(file_prefix + '_qasm2.txt', 'wt')
        self.qasm_out.write("IBMQASM 2.0;\n")
        self.qasm_out.write('include "qelib1.inc";\n')
        self.qasm_out.write("qreg q[5];\n")
        self.qasm_out.write("creg c[5];\n")
        self.qbtr_wr = None
        if write_qubiter_files:
            emb = CktEmbedder(num_bits, num_bits)
            out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
            self.qbtr_wr = SEO_writer(out_file_prefix, emb, **kwargs)
            self.qbtr_wr.write_NOTA("IBMQASM 2.0;")
            self.qbtr_wr.write_NOTA('include "qelib1.inc";')
            self.qbtr_wr.write_NOTA("qreg q[5];")
            self.qbtr_wr.write_NOTA("creg c[5];")

        SEO_reader.__init__(self, file_prefix, num_bits, verbose)

        self.qasm_out.write("measure q -> c;\n")
        if write_qubiter_files:
            self.qbtr_wr.write_NOTA("measure q -> c;")

        self.qasm_out.close()
        if write_qubiter_files:
            self.qbtr_wr.close_files()

    @staticmethod
    def qasm_line_for_rot(arr, tar_bit_pos):
        """
        This function returns a string for a qasm file line for a one qubit
        rotation.

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
        phi = ut.centered_rads(-2*left_rads)
        theta = ut.centered_rads(-2*center_rads)
        lam = ut.centered_rads(-2*right_rads)
        if abs(phi) < 1e-6 and abs(theta) < 1e-6:
            line_str = "u1(" + str(lam)
        elif abs(theta - np.pi/2) < 1e-6:
            line_str = "u2(" + str(phi) + ", " + str(lam)
        else:
            line_str = "u3(" + str(theta) + ", " + str(phi) + ", " + str(lam)
        line_str += ")  q[" + str(tar_bit_pos) + "];\n"
        return line_str

    def use_DIAG(self, trols, rad_angles):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        assert False, "No DIAG lines allowed"

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
        self.qasm_out.write("h   q[" + str(tar_bit_pos) + "];\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                   OneBitGates.had2)

    def use_LOOP(self, loop_num, reps):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        loop_num : int
        reps : int

        Returns
        -------
        None

        """
        assert False, "No LOOP lines allowed"

    def use_MEAS(self, tar_bit_pos, kind):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        kind : int
        tar_bit_pos : int

        Returns
        -------
        None

        """
        assert False, "No MEAS lines allowed"

    def use_MP_Y(self, tar_bit_pos, trols, rad_angles):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        tar_bit_pos : int
        trols : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """
        assert False, "No MP_Y lines allowed"

    def use_NEXT(self, loop_num):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        assert False, "No NEXT lines allowed"

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
        self.qasm_out.write("// " + bla_str + "\n")
        if self.write_qubiter_files:
            self.qbtr_wr.write_NOTA(bla_str)

    def use_PHAS(self, angle_degs, tar_bit_pos, controls):
        """
        If called, this function will halt execution of program.

        Parameters
        ----------
        angle_degs : float
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False, "No PHAS lines allowed"

    def use_P_PH(self, projection_bit, angle_degs, tar_bit_pos, controls):
        """
        If called, this function will halt execution of program.

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
        assert False, "No P0PH or P1PH lines allowed"

    def use_ROT(self, axis, angle_degs, tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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

        arr = OneBitGates.rot_ax(rad_ang, axis)
        line_str = Qubiter_to_IBMqasm2.qasm_line_for_rot(arr, tar_bit_pos)
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot_ax, [rad_ang, axis])

    def use_ROTN(self, angle_x_degs, angle_y_degs, angle_z_degs,
                tar_bit_pos, controls):
        """
        Writes line in IBM qasm file corresponding to an English file line
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
        line_str = Qubiter_to_IBMqasm2.qasm_line_for_rot(arr, tar_bit_pos)
        self.qasm_out.write(line_str)

        if self.write_qubiter_files:
            self.qbtr_wr.write_controlled_one_bit_gate(tar_bit_pos, controls,
                               OneBitGates.rot, rad_ang_list)

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
                prefix = "x   q["
            elif axis == 2:
                prefix = "y   q["
            elif axis == 3:
                prefix = "z   q["
            else:
                assert False
            self.qasm_out.write(prefix + str(tar_bit_pos) + "];\n")
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
            if tar_pos in self.targets[trol_pos]:
                self.qasm_out.write("cx  q[" + str(trol_pos) + "], "
                                    "q[" + str(tar_pos) + "];\n")
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
        If called, this function will halt execution of program.

        Parameters
        ----------
        bit1 : int
        bit2 : int
        controls : Controls

        Returns
        -------
        None

        """
        assert False, "No SWAP lines allowed"

if __name__ == "__main__":
    file_prefix = "../io_folder/qbtr2ibm_test"
    num_bits = 5
    import for_IBM_devices.ibm_chip_couplings as ibm
    c_to_t = ibm.ibmqx2_edges
    q2i = Qubiter_to_IBMqasm2(file_prefix,
            num_bits, c_to_t, write_qubiter_files=True)
