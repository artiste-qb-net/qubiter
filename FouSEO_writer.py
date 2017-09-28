from SEO_writer import *


class FouSEO_writer(SEO_writer):
    """
    This class is a subclass of SEO_writer so read that class' docstrings
    first. This class writes English and Picture files for the exact
    compilation, discovered by Coppersmith, of the Discrete Fourier
    Transform. We follow the conventions of quant-ph 0407215,
    "QC Paulinesia" by R.R. Tucci.

    Attributes
    ----------
    do_perm : bool
        True if want circuit to include permutation that reverses qbit order

    """
    def __init__(self, do_write, file_prefix, emb, do_perm=True, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        do_write : bool
            True if want constructor to write automatically without
            being asked.
        emb : CktEmbedder
        do_perm : bool
        True if want circuit to include permutation that reverses qbit order

        kwargs :

        Returns
        -------

        """
        self.do_perm = do_perm
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        if do_write:
            self.write()

    def write(self):
        """
        Writes circuit for quantum Fourier transform.

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef

        # permutation R
        if self.do_perm:
            for r in range(num_bits-1, 0, -1):
                for k in range(r-1, -1, -1):
                    self.write_bit_swap(r, k)

        for k in range(num_bits):
            self.write_one_bit_gate(k, OneBitGates.had2)
            trols = Controls.new_knob(num_bits, k, True)
            for r in range(k+1, num_bits):
                # note r>k
                self.write_controlled_one_bit_gate(
                    r,  # target bit pos
                    trols,
                    OneBitGates.P_1_phase_fac,
                    [np.pi/(1 << (r-k))]
                )

    def write_hermitian(self):
        """
        Write Hermitian conjugate of circuit written by write()

        Returns
        -------
        None

        """
        num_bits = self.emb.num_bits_bef

        for k in range(num_bits-1, -1, -1):
            trols = Controls.new_knob(num_bits, k, True)
            for r in range(num_bits-1, k, -1):
                # note r>k
                self.write_controlled_one_bit_gate(
                    r,  # target bit pos
                    trols,
                    OneBitGates.P_1_phase_fac,
                    [-np.pi/(1 << (r-k))]  # negative of write()
                )
            self.write_one_bit_gate(k, OneBitGates.had2)
        if self.do_perm:
            for r in range(1, num_bits):
                for k in range(r):
                    self.write_bit_swap(r, k)

    @staticmethod
    def fourier_trans_mat(num_rows, herm_conj=False):
        """
        This function returns a numpy array with the quantum Fourier
        transform in it.

        Parameters
        ----------
        num_rows : int
            number of rows
        herm_conj : bool
            If True, changes sign of all phases, which produces the
            hermitian conjugate of the matrix given when this is False

        Returns
        -------
        np.array

        """
        norma = np.sqrt(num_rows)
        sign = 1
        if herm_conj:
            sign = -1
        mat = np.empty((num_rows, num_rows), dtype=complex)
        for r in range(num_rows):
            for s in range(num_rows):
                mat[r, s] = np.exp(1j*sign*2*np.pi*r*s/num_rows)/norma
        return mat

if __name__ == "__main__":
    num_bits_bef = 4
    num_bits_aft = 6
    bit_map = list(range(num_bits_bef))
    emb = CktEmbedder(num_bits_bef, num_bits_aft, bit_map)

    for zf in [True, False]:
        wr = FouSEO_writer(True,
                           'io_folder/fou_test',
                           emb,
                           zero_bit_first=zf)
        wr.write_NOTA('do h.c. next')
        wr.write_hermitian()
        wr.close_files()

    print(FouSEO_writer.fourier_trans_mat(3))

