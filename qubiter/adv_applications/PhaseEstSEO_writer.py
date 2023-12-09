from qubiter.SEO_writer import *
from qubiter.FouSEO_writer import *
import copy as cp


class PhaseEstSEO_writer(SEO_writer):
    """
    The quantum phase estimation algorithm or PEA (see

    https://en.wikipedia.org/wiki/Quantum_phase_estimation_algorithm

    for an introduction) was invented by Kitaev in 1995. Since then, it has
    been used for many purposes such as for finding the ground state energy
    of molecules.

    This class is a subclass of SEO_writer. It writes the PEA circuit that
    is given, for example, in the Wikipedia article just mentioned.

    We will divide the qubits used by PEA into 2 disjoint sets and refer to
    those sets as: probe qubits and atom qubits. A matrix U that acts on the
    atom qubits will be referred to as the atom matrix or just the atom.
    First a Hadamard matrix is applied to each probe qubit. Then each probe
    qubit interacts with a different power of the atom matrix U. Finally,
    an inverse quantum Fourier transform is applied to all the probe qubits
    together.

    Even though these probe-atom interactions would not change the state of
    the probe qubits if they were classical particles, the probes do become
    correlated with the atom powers and this shows up after we take the
    inverse Fourier transform.

    Note that this class writes the whole PEA circuit, but it requires as
    input an object of a subclass of the class AtomWriter given below. This
    object will write the atom matrix powers.


    Attributes
    ----------
    atom_wr : AtomWriter
        An object of a subclass of the class AtomWriter given below. This
        object will write the atom matrix powers.
    do_perm : bool
        True if want quantum Fourier Transform circuit to include
        permutation that reverses qbit order
    num_probe_bits : int
        Number of probe qubits.

    """

    def __init__(self, do_write, num_probe_bits, atom_writer,
                 file_prefix, emb, do_perm=True, **kwargs):
        """
        Constructor

        Parameters
        ----------
        do_write : bool
            True if want constructor to write automatically without being
            asked.
        atom_writer : AtomWriter
        num_probe_bits : int
        file_prefix : str
        emb : CktEmbedder
        do_perm : bool
            True if want quantum Fourier Transform circuit to include
            permutation that reverses qbit order

        Returns
        -------

        """
        self.do_perm = do_perm
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        self.num_probe_bits = num_probe_bits

        self.atom_wr = atom_writer
        self.atom_wr.english_out = self.english_out
        self.atom_wr.picture_out = self.picture_out
        self.atom_wr.ZL = self.ZL

        if do_write:
            self.write()

    def write(self):
        """
        Writes the circuit for quantum phase estimation.

        Returns
        -------
        None

        """

        num_qbits = self.emb.num_qbits_bef
        num_atom_bits = num_qbits - self.num_probe_bits
        assert num_atom_bits >= 1, "must have >=1 probe bits"

        # first write the Hadamards
        for k in range(self.num_probe_bits):
            self.write_one_qbit_gate(k, OneQubitGate.had2)

        # next write the probe controlled atoms

        # this pre_emb maps atom -> (atom + probes)
        num_qbits_bef = num_atom_bits
        num_qbits_aft = num_qbits
        bit_map = list(range(self.num_probe_bits, num_qbits))
        pre_emb = CktEmbedder(num_qbits_bef, num_qbits_aft, bit_map)

        for k in range(self.num_probe_bits):
            pre_emb.extra_controls = \
                Controls.new_single_trol(num_qbits, k, True)
            compo_emb = CktEmbedder.composition(self.emb, pre_emb)
            self.atom_wr.emb = compo_emb
            self.atom_wr.write_pow(1 << k)

        # finally write the inverse Fourier transform

        # this pre_emb maps probe bits -> (atom + probes)
        num_qbits_bef = self.num_probe_bits
        num_qbits_aft = num_qbits
        bit_map = list(range(num_qbits_bef))
        pre_emb = CktEmbedder(num_qbits_bef, num_qbits_aft, bit_map)

        compo_emb = CktEmbedder.composition(self.emb, pre_emb)
        fou_writer = FouSEO_writer(
            do_write=False,
            file_prefix='blank',
            emb=compo_emb,
            do_perm=self.do_perm,
            english_out=self.english_out,
            picture_out=self.picture_out,
            ZL=self.ZL)
        fou_writer.write_hermitian()

    def write_hermitian(self):
        """
        Write Hermitian conjugate of circuit written by write().

        Returns
        -------
        None

        """

        num_qbits = self.emb.num_qbits_bef
        num_atom_bits = num_qbits - self.num_probe_bits
        assert num_atom_bits >= 1, "must have >=1 probe bits"

        # first write the Fourier transform

        # this pre_emb maps probe bits -> (atom + probes)
        num_qbits_bef = self.num_probe_bits
        num_qbits_aft = num_qbits
        bit_map = list(range(num_qbits_bef))
        pre_emb = CktEmbedder(num_qbits_bef, num_qbits_aft, bit_map)

        compo_emb = CktEmbedder.composition(self.emb, pre_emb)
        fou_writer = FouSEO_writer(
            do_write=False,
            file_prefix='blank',
            emb=compo_emb,
            do_perm=self.do_perm,
            english_out=self.english_out,
            picture_out=self.picture_out,
            ZL=self.ZL)
        fou_writer.write()

        # next write the probe controlled atoms

        # this pre_emb maps atom -> (atom + probes)
        num_qbits_bef = num_atom_bits
        num_qbits_aft = num_qbits
        bit_map = list(range(self.num_probe_bits, num_qbits))
        pre_emb = CktEmbedder(num_qbits_bef, num_qbits_aft, bit_map)

        for k in reversed(range(self.num_probe_bits)):
            pre_emb.extra_controls = \
                Controls.new_single_trol(num_qbits, k, True)
            compo_emb = CktEmbedder.composition(self.emb, pre_emb)
            self.atom_wr.emb = compo_emb
            self.atom_wr.write_pow_hermitian(1 << k)

        # finally write the Hadamards
        for k in reversed(range(self.num_probe_bits)):
            self.write_one_qbit_gate(k, OneQubitGate.had2)


class AtomWriter(SEO_writer):
    """
    An object of this class or of a subclass thereof is an attribute of
    PhaseEstSEO_writer(). If test=False, this class must be subclassed. If
    test=True, you get an example. In this example, the atom matrix is a
    simple controlled Ry rotation, but more generally, it can be a whole
    circuit.

    Attributes
    ----------

    This class has all the attributes of SEO_writer

    test : bool
        If test=True, the class uses testing parameters. If test=False,
        this becomes an abstract class that must be subclassed.

    """

    def __init__(self, do_write, test=False, file_prefix=None,
                 emb=None,  **kwargs):
        """
        Constructor

        Parameters
        ----------
        do_write : bool
        test :  bool

        file_prefix : str
        emb : CktEmbedder

        Returns
        -------

        """
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        self.test = test
        if do_write:
            self.write()

    def write_pow(self, power):
        """
        Writes circuit for U^power, where U is the atom matrix.

        Parameters
        ----------
        power : int

        Returns
        -------
        None

        """
        if not self.test:
            assert False

        num_qbits = self.emb.num_qbits_bef
        trols = Controls(num_qbits)
        tar_bit_pos = 0
        for k in range(1, num_qbits):
            trols.set_control(k, False)
        trols.refresh_lists()
        self.write_controlled_one_qbit_gate(
            tar_bit_pos,
            trols,
            OneQubitGate.rot_ax,
            [30*power*np.pi/180, 2])

    def write_pow_hermitian(self, power):
        """
        Write Hermitian conjugate of circuit written by write_pow().

        Parameters
        ----------
        power : int

        Returns
        -------
        None

        """
        if not self.test:
            assert False

        num_qbits = self.emb.num_qbits_bef
        trols = Controls(num_qbits)
        tar_bit_pos = 0
        for k in range(1, num_qbits):
            trols.set_control(k, False)
        trols.refresh_lists()
        self.write_controlled_one_qbit_gate(
            tar_bit_pos,
            trols,
            OneQubitGate.rot_ax,
            fun_arg_list=[-30*power*np.pi/180, 2])

    def write(self):
        """
        Same as write_pow(1)

        Returns
        -------
        None

        """
        self.write_pow(1)

    def write_hermitian(self):
        """
        Write Hermitian conjugate of circuit written by write().

        Returns
        -------
        None

        """
        self.write_pow_hermitian(1)


if __name__ == "__main__":
    def main():
        bit_map = list(range(7))
        fin_emb = CktEmbedder(7, 8, bit_map)
        atom_wr = AtomWriter(do_write=False, test=True)
        for ZL in [True, False]:
            wr = PhaseEstSEO_writer(do_write=False,
                                    num_probe_bits=4,
                                    atom_writer=atom_wr,
                                    file_prefix="ph_est_test",
                                    emb=fin_emb,
                                    ZL=ZL)
            wr.write()
            wr.write_NOTA("next write h.c.")
            wr.write_hermitian()
            wr.close_files()
    main()
