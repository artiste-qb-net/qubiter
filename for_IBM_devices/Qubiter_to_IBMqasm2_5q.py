from Controls import *
from SEO_reader import *
from SEO_writer import *
from quantum_CSD_compiler.UnitaryMat import *
import Utilities as ut


class Qubiter_to_IBMqasm2_5q(SEO_reader):
    """
    This class is a child of SEO_reader. It reads an input English file and
    writes an IBM qasm2 for 5 qubits file which is a line by line
    translation of the input English file into the IBM qasm2 language. If
    the option write_qubiter_files is set to True, this class will also
    write new English and Picture files that are in 1-1 line correspondence
    with the output qasm file.

    The input English file that is read can only have lines of the following
    types or else the program will abort with an error message:

    1. single qubit rotations (HAD2, SIGX, SIGY, SIGZ, ROTX, ROTY, ROTZ or 
    ROTN with no controls) 

    2. simple CNOTs (SIGX with a single True control)

    3. NOTA lines

    If you have an English file that contains lines that are more
    complicated than this (because, for example, they contain rotations with
    one or more controls attached), you can use the expander classes
    CktExpander, DiagUnitaryExpander, MultiplexorExpander, to expand the
    circuit to an equivalent albeit longer circuit that satisfies
    constraints 1, 2, 3.

    This class expects exactly 5 qubits, call them 0, 1, .., 4. The input
    English file circuit can contain CNOTs between ANY pair of qubits and
    with any qubit as target.

    This class conforms with the most recent IBM qc (IBM Quantum Experience
    QASM2.0). The current IBM qc is not fully connected (not all pairs of
    qubits are physically connected). Furthermore, not all qubits can be
    targets of an elementary CNOT: only qubits 1, 2 and 4 can be targets.

    If an elementary CNOT is not allowed because its ends are disconnected
    or its target is forbidden, this class will replace that elementary CNOT
    by a compound CNOT; i.e., a sequence of 1 or 4 allowed elementary CNOTs
    (and a bunch of Hadamards) that is equivalent to the original elementary
    CNOT. Next we discuss how these compound CNOTs are defined.

    Note that the positions of the target X and control @ of a CNOT can be
    swapped with Hadamard matrices H. For example,

    X---@

    can be replaced by:

    H   H
    @---X
    H   H

    Hence, a CNOT between 2 connected ends but having a disallowed target
    can be replaced by 4 Hadamards and an elementary CNOT with the same ends
    and with an allowed target.

    Suppose qubits a and b are disconnected and we want to implement this:

    a   2   b
    X---+---@

    This class replaces this problematic elementary CNOT by:

    a   2   b
    X---@   |
    |   X---@
    X---@   |
    |   X---@

    =

    H   |   |
    |   H   |
    @---X
    |   H   |
    |   X---@
    |   H   |
    @---X
    |   H   |
    |   X---@
    H   |   |

    (This is allowed because qubit 2 is connected to all other qubits,
    including a and b).

    Footnote: QASM has "measure" operations and distinguishes between
    quantum registers qreg and classical registers creg. Qubiter does not
    use cregs because it uses the classical memory of your Linux PC instead.
    QASM has an intricate set of commands for measurements. Qubiter has a
    complete set of measurement commands too (see MEAS in Rosetta stone).
    The QASM and Qubiter measurement commands can obviously be translated
    into each other. We leave that part of the translation to a future
    version of this class.

    References
    ----------
    1. https://github.com/IBMResearch/python-sdk-quantum-experience
    2. https://github.com/IBMQuantum/QASM

    Attributes
    ----------
    allowed_tars : list[int]
        Allowed targets. Qubits that are equipped in hardware to be targets
        of an elementary CNOT. IBM qc currently has 5 qubits 0, 1, ...,
        4. Out of those, only 1, 2 and 4 can be targets.
    connections : list[tuple(int,int)]
        Pairs of qubits that are physically connected so they can be the two
        ends of an elementary CNOT. Order of qubits in pairs is irrelevant.
        This picture indicates which qubits of the current IBM qc are
        connected:

            4     0
            | \ / |
            |  2  |
            | / \ |
            3     1

    qasm_out : _io.TextIOWrapper
        This output stream is used to write a qasm file based on the input
        English file.
    qbtr_wr : SEO_writer
        A SEO_writer object created iff write_qubiter_files is True.
    write_qubiter_files : bool
        The class always writes a qasm text file based on the input English
        file that is read. Iff this is True, the class also writes English
        and Picture files.

    """

    def __init__(self, file_prefix, num_bits, verbose=False,
                 write_qubiter_files=False, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        verbose : bool
        write_qubiter_files : bool
        kwargs : dict[]

        Returns
        -------
        None

        """
        assert num_bits == 5
        self.connections = [(4, 3),
                           (4, 2),
                           (3, 2),
                           (0, 2),
                           (2, 1),
                           (0, 1)]
        self.allowed_tars = [1, 2, 4]

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
            
    def are_connected(self, x, y):
        """
        This function returns true iff qubits x and y are connected in the
        hardware.

        Parameters
        ----------
        x : int
        y : int

        Returns
        -------
        bool

        """
        for a, b in self.connections:
            if (x, y) == (a, b) or (x, y) == (b, a):
                return True
        return False

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
        line_str = Qubiter_to_IBMqasm2_5q.qasm_line_for_rot(arr, tar_bit_pos)
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
        line_str = Qubiter_to_IBMqasm2_5q.qasm_line_for_rot(arr, tar_bit_pos)
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
            
        def cx_fun(c, t):
            self.qasm_out.write("cx  q[" + str(c) + "], q[" + str(t) + "];\n")

        def h_fun(t):
            self.qasm_out.write("h   q[" + str(t) + "];\n")
        
        def cx_fun1(c, t):
            trols = Controls.new_knob(5, c, True)
            self.qbtr_wr.write_controlled_one_bit_gate(t,
                trols, OneBitGates.sigx)

        def h_fun1(t):
            self.qbtr_wr.write_one_bit_gate(t, OneBitGates.had2)
            
        if num_trols == 0:
            prefix = ""
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
            if self.are_connected(trol_pos, tar_pos):
                if tar_pos in self.allowed_tars:
                    cx_fun(trol_pos, tar_pos)
                else:
                    h_fun(tar_pos)
                    h_fun(trol_pos)
                    cx_fun(tar_pos, trol_pos)
                    h_fun(trol_pos)
                    h_fun(tar_pos)                   
                if self.write_qubiter_files:
                    if tar_pos in self.allowed_tars:
                        cx_fun1(trol_pos, tar_pos)
                    else:
                        h_fun1(tar_pos)
                        h_fun1(trol_pos)
                        cx_fun1(tar_pos, trol_pos)
                        h_fun1(trol_pos)
                        h_fun1(tar_pos)
            else:
                # a   2   b
                # X---+---@
                # replaced by

                # H   |   |
                # |   H   |
                # @---X
                # |   H   |
                # |   X---@
                # |   H   |
                # @---X
                # |   H   |
                # |   X---@
                # H   |   |
                a = tar_pos
                b = trol_pos

                h_fun(a)
                h_fun(2)
                cx_fun(a, 2)
                h_fun(2)
                cx_fun(b, 2)
                h_fun(2)
                cx_fun(a, 2)
                h_fun(2)
                cx_fun(b, 2)
                h_fun(a)
                
                if self.write_qubiter_files:
                    h_fun1(a)
                    h_fun1(2)
                    cx_fun1(a, 2)
                    h_fun1(2)
                    cx_fun1(b, 2)
                    h_fun1(2)
                    cx_fun1(a, 2)
                    h_fun1(2)
                    cx_fun1(b, 2)
                    h_fun1(a)

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
    file_prefix = "io_folder/qbtr2ibm_test"
    q2i = Qubiter_to_IBMqasm2_5q(file_prefix, 5, write_qubiter_files=True)
