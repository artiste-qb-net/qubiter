# from chemistry.SuzukiSEO_writer import *
import Utilities as ut
from SEO_writer import *
# import itertools as it
from chemistry.CombisInterlacer import *


class MolEvolOpSEO_writer(SEO_writer):
    """
    Molecular Evolution Operator SEO writer. (1) This class and (2) a list
    of parameters constructed with the class MolEvolOpData or a child
    thereof, are passed as arguments to the constructor of
    PhaseEstSEO_writer. The role of this class is to write the "atom"
    circuit used by PhaseEstSEO_writer.

    So what is this atom circuit? It's a Trotter approximation of the
    evolution operator exp(-i t H) where the Hamiltonian H depends on the
    molecule being considered and is discussed further in the docstring for
    the class MolEvolOpData.

    A combi or combination is a list of ints. For example, a 4 bit combi is
    a list of 4 ints. A bunch is a list of combis.

    The Trotter approximation of an evolution operator exp(-i t H) is given
    by

    exp(-i t H) =
        [ exp(-i dt H_1) exp(-i dt H_2)...exp(-i dt H_n)]^(t/dt)

    We refer to t/dt as the number of Trotter cycles. We refer to each H_k
    as a spoke (of the wheel that is being cycled). A spoke has a core and a
    sheath. The core is a bunch of interlaced combis. The sheath is a
    unitary transformation V (aka as the left endcap) on left side of the
    core and the Hermitian conjugate of V (aka the right endcap) on the
    right side of the core. We choose the combis in the bunch of each spoke
    to be interlaced, because then all the combis of the bunch can be
    performed in parallel. Unlike the combis in the core bunch,
    the operations in the sheath of the spoke cannot be parallelized in any
    obvious way. However, some gates in the sheaths of adjacent spokes will
    cancel out. The class JWDummbBellNeutralizer tries do such cancellations.

    We will include one ancilla as the last qubit. Hence the total number of
    qubits of the circuit writen by this class will equal the number of
    orbitals plus 1.

    Parameters
    ----------
    global_theta : float
        global phase factor in radians
    diag_bit_to_theta : dict[int: float]
        for diagonal part of evolution operator, dictionary of int to theta
        angle in radians
    diag_2bits_to_theta : dict[list(int): float]
        for diagonal part of evolution operator, dictionary of list of 2
        ints to theta angle in radians
    _2bits_to_theta : dict[list(int): float]
        for non-diagonal part of evolution operator, dictionary of list of 2
        ints to theta angle in radians
    _3bits_to_theta : dict[list(int): float]
        for non-diagonal part of evolution operator, dictionary of list of 3
        ints to theta angle in radians
    _4bits_to_3thetas : dict[list(int): list(float)]
        for non-diagonal part of evolution operator, dictionary of list of 4
        ints to list of 3 theta angles in radians
    num_orbitals : int
        number of qubits = number of orbitals + 1. One ancilla to implement
        JW tails
    num_trot_cycles : int
        number of Trotter cycles
    total_time : float
        total time of evolution, t in exp( - i t H)
    do_notas : bool
        setting this to False turns off many NOTAS at once
    loop_count : int
    """

    def __init__(self, arg_list, do_write, file_prefix, emb, **kwargs):
        """
        Constructor

        Parameters
        ----------
        arg_list : list
        do_write : bool
        file_prefix : str
        emb : CktEmbedder
        kwargs : dict[]

        Returns
        -------

        """
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        [
            self.do_notas,
            self.num_orbitals,
            self.num_trot_cycles,
            self.total_time,
            self.global_theta,
            self.diag_bit_to_theta,
            self.diag_2bits_to_theta,
            self._2bits_to_theta,
            self._3bits_to_theta,
            self._4bits_to_3thetas
        ] = arg_list

        self.loop_count = 0

        if do_write:
            self.write()

    def write_stair(self, tar_bits, trol_bit, dir):
        """
        This function writes a stair (sequence) of controlled Pauli matrix (
        sigx, sigy or sigz, decided by the value of dir=1,2, 3) gates. All
        gates in the sequence have their control at trol_bit. The list
        tar_bits gives the locations of the Pauli matrices.

        example of controlled sigx stair. Assuming ZL, tar_bits = [0, 2]
        and trol_bit = 3

        |   |   @---+---+---X
        |   |   @---X   |   |

        Parameters
        ----------
        tar_bits : list[int]
            target bits
        trol_bit : int
            control bit
        dir : int
            direction, either 1,2 or 3

        Returns
        -------
        None

        """
        trol = Controls.new_knob(self.num_orbitals+1, trol_bit, True)
        if dir == 1:
            fun = OneBitGates.sigx
        elif dir == 2:
            fun = OneBitGates.sigy
        elif dir == 3:
            fun = OneBitGates.sigz
        else:
            assert False
        for k in range(len(tar_bits)):
            self.write_controlled_one_bit_gate(
                tar_bits[k], trol, fun, [])

    def write_stair_herm(self, tar_bits, trol_bit, dir):
        """
        This function writes the Hermitian of write_stair()

        Parameters
        ----------
        tar_bits : list[int]
            target bits
        trol_bit : int
            control bit
        dir : int
            direction, either 1,2 or 3

        Returns
        -------
        None

        """
        trol = Controls.new_knob(self.num_orbitals+1, trol_bit, True)
        if dir == 1:
            fun = OneBitGates.sigx
        elif dir == 2:
            fun = OneBitGates.sigy
        elif dir == 3:
            fun = OneBitGates.sigz
        else:
            assert False

        for k in reversed(range(len(tar_bits))):
            self.write_controlled_one_bit_gate(
                tar_bits[k], trol, fun, [])

    def write_jw_dumbbell_stair(self, z_bits, trol_bit):
        """
        This function writes a stair (sequence) of dumbbells. A dumbbell is
        a gate (-1)^{n(k1)n(k2)} for qubits k1 and k2, where n( k) = |1><1|
        at k = P_1(k). In this situation, the dumbbells arise from using the
        Jordan Wigner (JW) transformation to model fermionic degrees of
        freedom. All gates in the sequence have their control at trol_bit.
        The list z_bits gives the locations of the sigz's.

        example of JW dumbbell stair. Assuming ZL, z_bits = [0, 2]
        and trol_bit = 3

        |   |   @---+---+---Z
        |   |   @---Z   |   |

        Parameters
        ----------
        z_bits : list[int]
        trol_bit : int

        Returns
        -------
        None

        """
        self.write_NOTA("Begin JW dumbbell stair", self.do_notas)
        self.write_stair(z_bits, trol_bit, 3)
        self.write_NOTA("End JW dumbbell stair", self.do_notas)

    def write_jw_dumbbell_stair_herm(self, z_bits, trol_bit):
        """
        This function writes the Hermitian of write_jw_dumbbell_stair()

        Parameters
        ----------
        z_bits : list[int]
        trol_bit : int

        Returns
        -------
        None

        """
        self.write_NOTA("Begin JW dumbbell stair", self.do_notas)
        self.write_stair_herm(z_bits, trol_bit, 3)
        self.write_NOTA("End JW dumbbell stair", self.do_notas)

    def write_controlled_sigx_stair(self, x_bits, trol_bit):
        """
        This function writes a stair (sequence) of controlled sigx gates.
        All gates in the sequence have their control at trol_bit. The list
        x_bits gives the locations of the sigx's.

        example of controlled sigx stair. Assuming ZL, x_bits = [0, 2]
        and trol_bit = 3

        |   |   @---+---+---X
        |   |   @---X   |   |

        Parameters
        ----------
        x_bits : list[int]
        trol_bit : int

        Returns
        -------
        None

        """
        self.write_NOTA("Begin controlled sigx_stair", self.do_notas)
        self.write_stair(x_bits, trol_bit, 1)
        self.write_NOTA("End controlled sigx_stair", self.do_notas)

    def write_controlled_sigx_stair_herm(self, x_bits, trol_bit):
        """
        This function writes the Hermitian of write_controlled_sigx_stair()

        Parameters
        ----------
        x_bits : list[int]
        trol_bit : int

        Returns
        -------
        None

        """
        self.write_NOTA("Begin controlled sigx stair", self.do_notas)
        self.write_stair_herm(x_bits, trol_bit, 1)
        self.write_NOTA("End controlled sigx stair", self.do_notas)

    def write_2bit_combi(self, bits, theta):
        """
        This function writes, without the begin and end caps, a gate
        representing

        exp(i sum_{k1 < k2} theta(k1, k2) n(k1)n(k2))

        Caps are written by other functions.

        Parameters
        ----------
        bits : list[int]
        theta : float

        Returns
        -------
        None

        """

        assert len(bits) == 2
        assert bits[0] < bits[1]

        trol0 = Controls.new_knob(self.num_orbitals+1, bits[0], True)
        self.write_controlled_one_bit_gate(
            bits[1], trol0, OneBitGates.rot_ax, [theta, 1])

    def write_3bit_combi(self, bits, theta):
        """
        This function writes, without the begin and end caps, a gate
        representing

        exp(i sum_{k1 < k2 , k3 != k1, k2}
            theta(k1, k2, k3) a^\dag(k1) a(k2) n(k3) + h.c.)

        Caps are written by other functions.

        Parameters
        ----------
        bits : list[int]
        theta : float

        Returns
        -------
        None

        """
        assert len(bits) == 3
        assert bits[0] < bits[1]

        trols02 = Controls(self.num_orbitals+1)
        trols02.bit_pos_to_kind = {bits[0]: True, bits[2]: True}
        trols02.refresh_lists()
        self.write_controlled_one_bit_gate(
            bits[1], trols02, OneBitGates.rot_ax, [theta, 1])

    def write_4bit_combi(self, bits, thetas):
        """
        This function writes, without the begin and end caps, a gate 
        representing 

        exp(i sum_{k1 < k2 < k3 < k4}
            theta_0(k1, k2, k3, k4) a^\dag(k1) a^\dag(k2) a(k3) a(k4) + h.c.
            theta_1(k1, k2, k3, k4) a^\dag(k1) a(k2) a^\dag(k3) a(k4) + h.c.
            theta_2(k1, k2, k3, k4) a^\dag(k1) a(k2) a(k3) a^\dag(k4) + h.c.

        Parameters
        ----------
        bits : list[int]
        thetas : list[float]

        Returns
        -------
        None

        """
        assert len(bits) == 4
        assert len(thetas) == 3
        assert bits[0] < bits[1] < bits[2] < bits[3]

        for theta_index in range(0, 3):
            theta = ut.centered_rads(thetas[theta_index])
            if abs(theta) > ut.TOL:
                trols = Controls(self.num_orbitals+1)
                for r in range(0, 3):
                    if r == 2 - theta_index:
                        trols.set_control(bits[r], False)
                    else:
                        trols.set_control(bits[r], True)
                trols.refresh_lists()
                self.write_controlled_one_bit_gate(
                    bits[3], trols, OneBitGates.rot_ax, [theta, 1])

    def get_2bit_bunch_coverage(self, bunch):
        """
        the coverage of a 2 bit combi [a, b] is defined as the set of ints k
        such that a <= k <= b. The coverage of a bunch of 2 bit combis is
        defined as the union of the coverages of all the combis in the bunch.

        This function returns a list of ints that represents the coverage of
        the input bunch of 2 bit combis.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        list[int]

        """

        is_covered = [False]*self.num_orbitals

        for combi in bunch:
            assert combi[0] < combi[1]
            for k in range(combi[0], combi[1]+1):
                is_covered[k] = True
        coverage = [k for k in range(self.num_orbitals) if is_covered[k]]
        return coverage

    def get_3bit_bunch_coverage(self, bunch):
        """
        Let 'bunch' be a bunch of 3 bit combis and 'bunch2' be a bunch of 2
        bit combis obtained from bunch by dropping the third int of each 3
        bit combi. This function returns a list of ints that represents the
        coverage of bunch2.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        list[int]

        """
        return self.get_2bit_bunch_coverage(bunch)

    def get_4bit_bunch_coverage(self, bunch):
        """
        Let 'bunch' be a bunch of 4 bit combis and 'bunch2' be a bunch of 2
        bit combis obtained from bunch by breaking each 4 bit combi [a, b,
        c, d] into two combis [a, b] and [c, d]. This function returns a
        list of ints that represents the coverage of bunch2.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        list[int]

        """

        is_covered = [False]*self.num_orbitals

        for combi in bunch:
            assert combi[0] < combi[1] < combi[2] < combi[3]
            for k in range(combi[0], combi[1]+1):
                is_covered[k] = True
            for k in range(combi[2], combi[3]+1):
                is_covered[k] = True
        coverage = [k for k in range(self.num_orbitals) if is_covered[k]]
        return coverage

    def write_diag_1bit_spoke(self):
        """
        This function writes a spoke consisting of a sigz rotation at each
        qubit. This spoke is a diagonal matrix.

        Returns
        -------
        None

        """

        self.write_NOTA("****Begin diag 1 bit spoke", self.do_notas)

        self.write_global_phase_fac(self.global_theta)

        for k, theta in self.diag_bit_to_theta.items():
            self.write_one_bit_gate(
                k, OneBitGates.rot_ax, [theta, 3])

        self.write_NOTA("****End diag 1 bit spoke", self.do_notas)

    def write_diag_2bit_spoke(self, bunch):
        """
        This function writes a spoke consisting of a bunch of interlaced
        2 qubit gates. This spoke is a diagonal matrix.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        None

        """

        self.write_NOTA("****Begin diag 2bit spoke", self.do_notas)

        for combi in bunch:
            if combi in self.diag_2bits_to_theta:
                trols = Controls.new_knob(self.num_orbitals+1, combi[1], True)
                theta = self.diag_2bits_to_theta[combi]
                self.write_controlled_one_bit_gate(
                    combi[0], trols, OneBitGates.rot_ax, [theta, 3])

        self.write_NOTA("****End diag 2bit spoke", self.do_notas)

    def write_2bit_spoke(self, bunch):
        """
        This function writes a spoke whose core consists of a bunch of
        interlaced 2 qubit gates.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        None

        """

        self.write_NOTA("****Begin 2 bit spoke", self.do_notas)

        coverage = self.get_2bit_bunch_coverage(bunch)

        self.write_jw_dumbbell_stair(coverage, self.num_orbitals)
        for combi in bunch:
            self.write_controlled_sigx_stair(combi[0:1], combi[1])

        self.write_one_bit_gate(self.num_orbitals, OneBitGates.sigx)
        for combi in bunch:
            self.write_2bit_combi(
                combi, self._2bits_to_theta[combi])

        for combi in reversed(bunch):
            self.write_controlled_sigx_stair_herm(combi[0:1], combi[1])
        self.write_jw_dumbbell_stair_herm(coverage, self.num_orbitals)

        self.write_NOTA("****End 2 bit spoke", self.do_notas)

    def write_3bit_spoke(self, bunch):
        """
        This function writes a spoke whose core consists of a bunch of
        interlaced 3 qubit gates.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        None

        """
        self.write_NOTA("****Begin 3 bit spoke", self.do_notas)

        coverage = self.get_3bit_bunch_coverage(bunch)

        self.write_jw_dumbbell_stair(coverage, self.num_orbitals)
        for combi in bunch:
            self.write_controlled_sigx_stair(combi[0:1], combi[1])

        self.write_one_bit_gate(self.num_orbitals, OneBitGates.sigx)
        for combi in bunch:
            self.write_3bit_combi(combi, self._3bits_to_theta[combi])

        for combi in reversed(bunch):
            self.write_controlled_sigx_stair_herm(combi[0:1], combi[1])
        self.write_jw_dumbbell_stair_herm(coverage, self.num_orbitals)

        self.write_NOTA("****End 3 bit spoke", self.do_notas)

    def write_4bit_spoke(self, bunch):
        """
        This function writes a spoke whose core consists of a bunch of
        interlaced 4 qubit gates.

        Parameters
        ----------
        bunch : tuple[tuple[int]]

        Returns
        -------
        None

        """

        self.write_NOTA("****Begin 4 bit spoke", self.do_notas)

        coverage = self.get_4bit_bunch_coverage(bunch)

        self.write_jw_dumbbell_stair(coverage, self.num_orbitals)
        for combi in bunch:
            self.write_controlled_sigx_stair(combi[0:3], combi[3])

        self.write_one_bit_gate(self.num_orbitals, OneBitGates.sigx)
        for combi in bunch:                    
            self.write_4bit_combi(combi, self._4bits_to_3thetas[combi])

        for combi in bunch:
            self.write_controlled_sigx_stair_herm(combi[0:3], combi[3])
        self.write_jw_dumbbell_stair_herm(coverage, self.num_orbitals)

        self.write_NOTA("****End 4 bit spoke", self.do_notas)

    def write_trotterized_evol_op(self):
        """
        This function writes a Trotter approximation of an evolution
        operator U = exp(-i t H)

        Returns
        -------
        None

        """
        assert self.num_orbitals+1 == self.emb.num_bits_bef
        self.loop_count -= 1
        self.write_LOOP(self.loop_count, self.num_trot_cycles)
        
        self.write_diag_1bit_spoke()

        ci = CombisInterlacer(n=self.num_orbitals,
            combi_list=self.diag_2bits_to_theta.keys())
        for bunch in ci.bunch_gen():
            self.write_diag_2bit_spoke(bunch)

        ci = CombisInterlacer(n=self.num_orbitals,
            combi_list=self._2bits_to_theta.keys())
        for bunch in ci.bunch_gen():
            self.write_2bit_spoke(bunch)

        ci = CombisInterlacer(n=self.num_orbitals,
            combi_list=self._3bits_to_theta.keys())
        for bunch in ci.bunch_gen():
            self.write_3bit_spoke(bunch)

        ci = CombisInterlacer(n=self.num_orbitals,
            combi_list=self._4bits_to_3thetas.keys())
        for bunch in ci.bunch_gen():
            self.write_4bit_spoke(bunch)
            
        self.write_NEXT(self.loop_count)

    def write(self):
        """
        This function writes a Trotter approximation of an evolution
        operator U = exp(-i t H)


        Parameters
        ----------

        Returns
        -------
        None

        """
        self.write_trotterized_evol_op()

    def write_pow(self, power):
        """
        This function writes the nth power of an evolution operator U = exp(
        -i t H). Raising to the nth power is achieved by a LOOP with n
        repetitions.


        Parameters
        ----------
        power : int

        Returns
        -------
        None

        """

        self.write_LOOP(power, power)
        self.write()
        self.write_NEXT(power)

from chemistry.MolEvolOpData import *
if __name__ == "__main__":
    data = MolEvolOpData()
    data.set_test_arg_list()

    arg_list = data.get_arg_list()
    file_prefix = "chem_io_folder//trotter_evol"
    bit_map = list(range(5))
    emb = CktEmbedder(5, 5, bit_map)
    wr = MolEvolOpSEO_writer(arg_list,
        False,
        file_prefix,
        emb)
    wr.write_pow(12)
