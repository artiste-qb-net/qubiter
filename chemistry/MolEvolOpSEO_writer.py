import Utilities as ut
from PhaseEstSEO_writer import AtomWriter
from Controls import *
from CktEmbedder import *
from OneBitGates import *
# import itertools as it
from chemistry.CombisInterlacer import *
from chemistry.MolEvolOpData import *


class MolEvolOpSEO_writer(AtomWriter):
    """
    Molecular Evolution Operator SEO writer. An object of this class is
    passed as an argument to the constructor of PhaseEstSEO_writer. The role
    of this class is to write the "atom" circuit used by PhaseEstSEO_writer.

    So what is this atom circuit? It's a Trotter or Trotter-Suzuki
    approximation of the evolution operator exp(-i t H) where the
    Hamiltonian H depends on the molecule being considered and is discussed
    further in the docstring for the class MolEvolOpData.

    The Trotter approximation of an evolution operator exp(-i t H) is given
    by

    exp(-i t H) = [E_1 E_2...E_{N-1}E_N]^(t/dt)

    where E_k = exp(-i dt H_k). The more accurate and symmetrical 2cnd order
    Trotter-Suzuki approximation is defined by

    exp(-i t H) = [A_1 A_2...A_{N-1}A_N   A_N k_{N-1}...A_2, A_1]^(t/dt)

    where A_k = sqrt(E_k)

    We refer to t/dt as the number of Trotter cycles. We refer to each E_k
    or A_k as a spoke (of the wheel that is being cycled).

    A combi or combination is a list of ints. For example, a 4 bit combi is
    a list of 4 ints. A bunch is a list of combis. A bunch of combis will be
    said to be interlaced if the combis are disjoint as sets.

    A spoke has a core and a sheath. The core is a bunch of interlaced
    combis. The sheath is a unitary transformation V (aka as the left
    endcap) on left side of the core and the Hermitian conjugate of V (aka
    the right endcap) on the right side of the core.

    We choose the combis in the bunch of each spoke to be interlaced,
    because then all the combis of the bunch can be performed in parallel.
    Unlike the combis in the core bunch, the operations in the sheath of the
    spoke cannot be parallelized in any obvious way. However, some gates in
    the sheaths of adjacent spokes will cancel out. The class
    DummbBellNeutralizer tries do such cancellations.

    We will use the last qubit as an ancilla to implement JW tails. Hence
    the total number of qubits of the circuit written by this class will
    equal the number of orbitals plus 1.

    For more information about the technical details underpinning this class
    and its dependents, see the file `ph-est-chem.pdf' in the same directory
    as this file. That pdf is a paper entitled ``Quantum Circuit for
    Estimating the Ground State Energy of Molecules", by R.R.Tucci

    Attributes
    ----------
    data : MolEvolOpData
    loop_count : int

    """

    def __init__(self, data, do_write, **kwargs):
        """
        Constructor

        Parameters
        ----------
        data : MolEvolOpData
        do_write : bool
        kwargs : dict[]

        Returns
        -------

        """
        AtomWriter.__init__(self, do_write, **kwargs)
        self.data = data

        self.loop_count = 0

        if do_write:
            self.write()

    def write_stair(self, tar_bits, trol_bit, axis):
        """
        This function writes a stair (sequence) of controlled Pauli matrix (
        sigx, sigy or sigz, decided by the value of axis=1, 2, 3) gates. All
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
        axis : int
            axis, either 1,2 or 3

        Returns
        -------
        None

        """
        trol = Controls.new_knob(self.data.num_orbitals+1, trol_bit, True)
        if axis == 1:
            fun = OneBitGates.sigx
        elif axis == 2:
            fun = OneBitGates.sigy
        elif axis == 3:
            fun = OneBitGates.sigz
        else:
            assert False
        for k in range(len(tar_bits)):
            self.write_controlled_one_bit_gate(
                tar_bits[k], trol, fun, [])

    def write_stair_herm(self, tar_bits, trol_bit, axis):
        """
        This function writes the Hermitian of write_stair()

        Parameters
        ----------
        tar_bits : list[int]
            target bits
        trol_bit : int
            control bit
        axis : int
            axis, either 1,2 or 3

        Returns
        -------
        None

        """
        trol = Controls.new_knob(self.data.num_orbitals+1, trol_bit, True)
        if axis == 1:
            fun = OneBitGates.sigx
        elif axis == 2:
            fun = OneBitGates.sigy
        elif axis == 3:
            fun = OneBitGates.sigz
        else:
            assert False

        for k in reversed(range(len(tar_bits))):
            self.write_controlled_one_bit_gate(
                tar_bits[k], trol, fun, [])

    def write_jw_dumbbell_stair(self, z_bits, trol_bit):
        """
        This function writes a stair (sequence) of dumbbells. A dumbbell is
        a gate (-1)^{n(k1)n(k2)} for qubits k1 and k2, where n(k) = |1><1|
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
        self.write_NOTA("Begin JW dumbbell stair", self.data.do_notas)
        self.write_stair(z_bits, trol_bit, 3)
        self.write_NOTA("End JW dumbbell stair", self.data.do_notas)

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
        self.write_NOTA("Begin JW dumbbell stair", self.data.do_notas)
        self.write_stair_herm(z_bits, trol_bit, 3)
        self.write_NOTA("End JW dumbbell stair", self.data.do_notas)

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
        self.write_NOTA("Begin controlled sigx_stair", self.data.do_notas)
        self.write_stair(x_bits, trol_bit, 1)
        self.write_NOTA("End controlled sigx_stair", self.data.do_notas)

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
        self.write_NOTA("Begin controlled sigx stair", self.data.do_notas)
        self.write_stair_herm(x_bits, trol_bit, 1)
        self.write_NOTA("End controlled sigx stair", self.data.do_notas)

    def write_2bit_ckt(self, bits, theta):
        """
        This function writes, without the left and right endcaps, a gate
        representing

        exp(i sum_{k1 < k2} theta(k1, k2) a^\dag(k1)a(k2))

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

        trol0 = Controls.new_knob(self.data.num_orbitals+1, bits[0], True)
        self.write_controlled_one_bit_gate(
            bits[1], trol0, OneBitGates.rot_ax, [theta, 1])

    def write_3bit_ckt(self, bits, theta):
        """
        This function writes, without the left and right endcaps, a gate
        representing

        exp(i \sum_{k1 < k2} \sum{k3 != k1, k2}
            theta(k1, k2, k3) a^\dag(k1) a(k2) n(k3) + h.c.)

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

        trols02 = Controls(self.data.num_orbitals+1)
        trols02.bit_pos_to_kind = {bits[0]: True, bits[2]: True}
        trols02.refresh_lists()
        self.write_controlled_one_bit_gate(
            bits[1], trols02, OneBitGates.rot_ax, [theta, 1])

    def write_4bit_ckt(self, bits, thetas):
        """
        This function writes, without the left and right endcaps, a gate
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
            if theta_index == 1:
                sign = +1
            else:
                sign = -1
            if abs(theta) > 1e-6:
                trols = Controls(self.data.num_orbitals+1)
                for r in range(0, 3):
                    if r == 2 - theta_index:
                        trols.set_control(bits[r], False)
                    else:
                        trols.set_control(bits[r], True)
                trols.refresh_lists()
                self.write_controlled_one_bit_gate(
                    bits[3], trols, OneBitGates.rot_ax, [sign*theta, 1])

    @staticmethod
    def extended_combi(combi):
        """
        This function takes a combi and it returns its extended combi. The
        extended combi has ints added to fill range(combi[0]+1, combi[1]) (
        and also range(combi[2]+1, combi[3]) for 4 bits). These bits are
        used by the JW dumbbells of the gate which combi represents, and may
        not be used by other non-overlapping gates at the same time.

        Parameters
        ----------
        combi : tuple[int]

        Returns
        -------
        tuple[int]

        """
        le = len(combi)
        assert le >= 2
        lista = list(range(combi[0], combi[1]+1))
        if le == 2:
            pass
        elif le == 3:
            lista += [combi[2]]
        elif le == 4:
            lista += list(range(combi[2], combi[3]+1))
        else:
            assert False
        return tuple(lista)

    def write_diag_1bit_spoke(self, theta_frac=1):
        """
        This function writes a spoke consisting of a sigz rotation at each
        qubit. This spoke is a diagonal matrix.

        Parameters
        ----------
        theta_frac : float
            theta fraction, multiply all thetas by this fraction

        Returns
        -------
        None

        """

        self.write_NOTA("****Begin diag 1 bit spoke", self.data.do_notas)

        self.write_global_phase_fac(self.data.global_theta*theta_frac)

        for k, theta in self.data.diag_bit_to_theta.items():
            self.write_one_bit_gate(
                k, OneBitGates.rot_ax, [theta*theta_frac, 3])

        self.write_NOTA("****End diag 1 bit spoke", self.data.do_notas)

    def write_diag_2bit_spoke(self, bunch, theta_frac=1):
        """
        This function writes a spoke consisting of a bunch of interlaced
        2 qubit gates. This spoke is a diagonal matrix.

        Parameters
        ----------
        bunch : tuple[tuple[int]]
        theta_frac : float
            theta fraction, multiply all thetas by this fraction

        Returns
        -------
        None

        """

        self.write_NOTA("****Begin diag 2bit spoke", self.data.do_notas)

        for combi in bunch:
            trols = Controls.new_knob(self.data.num_orbitals+1, combi[1], True)
            theta = self.data.diag_2bits_to_theta[combi]
            self.write_controlled_one_bit_gate(
                combi[0], trols, OneBitGates.rot_ax, [theta*theta_frac, 3])

        self.write_NOTA("****End diag 2bit spoke", self.data.do_notas)

    def write_nondiag_spoke(self, bunch, theta_frac=1):
        """
        This function writes a spoke (spoke = left endcap, core,
        right endcap) whose core consists of a bunch of interlaced 2,
        3 or 4 qubit gates.

        Parameters
        ----------
        bunch : tuple[tuple[int]]
        theta_frac : float
            theta fraction, multiply all thetas by this fraction

        Returns
        -------
        None

        """

        for combi in bunch:
            le = len(combi)
            if le in [2, 3]:
                jw_range = range(combi[0]+1, combi[1])
                self.write_jw_dumbbell_stair(jw_range, self.data.num_orbitals)
            elif le == 4:
                jw_range = range(combi[0]+1, combi[1])
                self.write_jw_dumbbell_stair(jw_range, self.data.num_orbitals)
                jw_range = range(combi[2]+1, combi[3])
                self.write_jw_dumbbell_stair(jw_range, self.data.num_orbitals)
            else:
                assert False

        for combi in bunch:
            le = len(combi)
            if le in [2, 3]:
                self.write_controlled_sigx_stair(combi[0:1], combi[1])
            elif le == 4:
                self.write_controlled_sigx_stair(combi[0:3], combi[3])
            else:
                assert False

        self.write_NOTA("****Begin spoke core", self.data.do_notas)
        self.write_one_bit_gate(self.data.num_orbitals, OneBitGates.sigx)
        for combi in bunch:
            le = len(combi)
            if le == 2:
                self.write_2bit_ckt(combi,
                    self.data.its_2bits_to_theta[combi]*theta_frac)
            elif le == 3:
                self.write_3bit_ckt(combi,
                    self.data.its_3bits_to_theta[combi]*theta_frac)
            elif le == 4:
                frac_thetas = [xx*theta_frac for
                               xx in self.data.its_4bits_to_3thetas[combi]]
                self.write_4bit_ckt(combi, frac_thetas)
            else:
                assert False
        self.write_NOTA("****End spoke core", self.data.do_notas)

        for combi in bunch:
            le = len(combi)
            if le in [2, 3]:
                self.write_controlled_sigx_stair_herm(combi[0:1], combi[1])
            elif le == 4:
                self.write_controlled_sigx_stair_herm(combi[0:3], combi[3])
            else:
                assert False

        for combi in bunch:
            le = len(combi)
            if le in [2, 3]:
                jw_range = range(combi[0]+1, combi[1])
                self.write_jw_dumbbell_stair_herm(jw_range,
                                                  self.data.num_orbitals)
            elif le == 4:
                jw_range = range(combi[2]+1, combi[3])
                self.write_jw_dumbbell_stair_herm(jw_range,
                                                  self.data.num_orbitals)
                jw_range = range(combi[0]+1, combi[1])
                self.write_jw_dumbbell_stair_herm(jw_range,
                                                  self.data.num_orbitals)
            else:
                assert False

    def write_trotterized_evol_op(self):
        """
        This function writes a Trotter approximation of an evolution
        operator U = exp(-i t H)

        Returns
        -------
        None

        """
        assert self.data.num_orbitals+1 == self.emb.num_bits_bef
        self.loop_count -= 1
        self.write_LOOP(self.loop_count, self.data.num_trot_cycles)
        
        self.write_diag_1bit_spoke()

        ci = CombisInterlacer(n=self.data.num_orbitals,
            combi_list=self.data.diag_2bits_to_theta.keys())
        for bunch in ci.bunch_gen():
            self.write_diag_2bit_spoke(bunch)

        xcombi_to_combi = {}
        xcombi_to_combi.update({self.extended_combi(combi): combi for
                            combi in self.data.its_4bits_to_3thetas})
        xcombi_to_combi.update({self.extended_combi(combi): combi for
                            combi in self.data.its_3bits_to_theta})
        xcombi_to_combi.update({self.extended_combi(combi): combi for
                            combi in self.data.its_2bits_to_theta})

        ci = CombisInterlacer(n=self.data.num_orbitals,
                              combi_list=xcombi_to_combi.keys())
        for xbunch in ci.bunch_gen():
            bunch = tuple([xcombi_to_combi[xcombi] for xcombi in xbunch])
            self.write_nondiag_spoke(bunch)
            
        self.write_NEXT(self.loop_count)

    def write_sec_ord_evol_op(self):
        """
        This function writes a second order Trotter-Suzuki approximation of
        an evolution operator U = exp(-i t H)

        Returns
        -------
        None

        """
        assert self.data.num_orbitals+1 == self.emb.num_bits_bef
        self.loop_count -= 1
        self.write_LOOP(self.loop_count, self.data.num_trot_cycles)

        self.write_diag_1bit_spoke(theta_frac=1/2)

        ci_diag = CombisInterlacer(n=self.data.num_orbitals,
            combi_list=self.data.diag_2bits_to_theta.keys())
        for bunch in ci_diag.bunch_gen():
            self.write_diag_2bit_spoke(bunch, theta_frac=1/2)

        xcombi_to_combi = {}
        xcombi_to_combi.update({self.extended_combi(combi): combi for
                            combi in self.data.its_4bits_to_3thetas})
        xcombi_to_combi.update({self.extended_combi(combi): combi for
                            combi in self.data.its_3bits_to_theta})
        xcombi_to_combi.update({self.extended_combi(combi): combi for
                            combi in self.data.its_2bits_to_theta})

        ci = CombisInterlacer(n=self.data.num_orbitals,
                              combi_list=xcombi_to_combi.keys())
        bunch_ctr = 0
        num_bunches = len(ci.bunch_to_combis)
        for xbunch in ci.bunch_gen():
            bunch_ctr += 1
            if bunch_ctr != num_bunches:
                frac = 1/2
            else:
                frac = 1.0
                self.write_NOTA(
                    "**********Begin central spoke", self.data.do_notas)
            bunch = tuple([xcombi_to_combi[xcombi] for xcombi in xbunch])
            self.write_nondiag_spoke(bunch, theta_frac=frac)

        self.write_NOTA("**********End central spoke", self.data.do_notas)
        bunch_ctr = 0
        for xbunch in ci.reversed_bunch_gen():
            bunch_ctr += 1
            if bunch_ctr != 1:
                bunch = tuple([xcombi_to_combi[xcombi] for xcombi in xbunch])
                self.write_nondiag_spoke(bunch, theta_frac=1/2)

        for bunch in ci_diag.reversed_bunch_gen():
            self.write_diag_2bit_spoke(bunch, theta_frac=1/2)

        self.write_diag_1bit_spoke(theta_frac=1/2)

        self.write_NEXT(self.loop_count)

    def write(self):
        """
        This function writes a Trotter or a second order Trotter-Suzuki
        approximation of an evolution operator U = exp(-i t H)


        Returns
        -------
        None

        """
        if self.data.approx == 1:
            self.write_trotterized_evol_op()
        elif self.data.approx == 2:
            self.write_sec_ord_evol_op()
        else:
            assert False

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


if __name__ == "__main__":
    data = MolEvolOpData(test=True)

    bit_map = list(range(5))
    emb = CktEmbedder(5, 5, bit_map)

    data.approx = 1
    file_prefix = "chem_io_folder/trotter_evol"
    wr = MolEvolOpSEO_writer(
        data,
        do_write=False,
        file_prefix=file_prefix,
        emb=emb)
    wr.write_pow(12)

    data.approx = 2
    file_prefix = "chem_io_folder/trotter_suzuki_evol"
    wr = MolEvolOpSEO_writer(
        data,
        do_write=False,
        file_prefix=file_prefix,
        emb=emb)
    wr.write_pow(12)
