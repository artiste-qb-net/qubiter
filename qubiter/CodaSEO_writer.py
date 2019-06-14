from qubiter.SEO_writer import *
from shutil import copyfile


class CodaSEO_writer(SEO_writer):
    """
    This class is a child of class `SEO_writer`. The constructor of this
    class accepts as input the paths to a pair of well-formed initial
    English and Picture files. The class makes copies of those 2 initial
    files and opens the copies, called the final English and Picture files,
    in the append mode. Thereafter, the user can write to those final files
    using methods that this class inherits from its parent SEO_writer.
    Cola SEO writer means tail-end SEO writer, which accurately describes
    what this class does.

    Attributes
    ----------
    fin_eng_path : str
        path to final English file that starts life as a verbatim copy of
        the initial English file and is thereafter appended to.

    fin_pic_path : str
        path to final Picture file that starts life as a verbatim copy of
        the initial Picture file and is thereafter appended to.
    """

    def __init__(self, init_file_prefix, fin_file_prefix,
                 num_bits, ZL=True, fin_emb=None):
        """

        Parameters
        ----------
        init_file_prefix : str
        fin_file_prefix : str
        num_bits : int
        ZL : bool
        fin_emb : CktEmbedder
            circuit embedder for the writer of the final English and Picture
            files

        Returns
        -------

        """

        assert init_file_prefix != fin_file_prefix,\
            "intial and final file prefixes can't be the same."

        ZL_str = 'ZL' if ZL else 'ZF'

        def eng_path(file_prefix):
            return file_prefix +\
                   '_' + str(num_bits) + '_eng.txt'

        def pic_path(file_prefix):
            return file_prefix +\
                   '_' + str(num_bits) + '_' + ZL_str + 'pic.txt'

        init_eng_path = eng_path(init_file_prefix)
        fin_eng_path = eng_path(fin_file_prefix)
        self.fin_eng_path = fin_eng_path

        init_pic_path = pic_path(init_file_prefix)
        fin_pic_path = pic_path(fin_file_prefix)
        self.fin_pic_path = fin_pic_path

        # copying files
        try:
            copyfile(init_eng_path, fin_eng_path)
            # print(',,ll', open(fin_eng_path).read())
        except:
            assert False, 'Could not copy file\n' + init_eng_path

        try:
            copyfile(init_pic_path, fin_pic_path)
        except:
            assert False, 'Could not copy file\n' + init_pic_path

        # opening copies of files
        try:
            fin_eng_out = open(fin_eng_path, 'a')
        except:
            assert False, 'Could not open file\n' + fin_eng_path

        try:
            fin_pic_out = open(fin_pic_path, 'a')
        except:
            assert False, 'Could not open file\n' + fin_pic_path

        emb = CktEmbedder(num_bits, num_bits)
        if fin_emb:
            emb = fin_emb

        SEO_writer.__init__(self, fin_file_prefix, emb,
                            ZL=ZL,
                            english_out=fin_eng_out,
                            picture_out=fin_pic_out)

    def write_xy_measurements(self, bit_pos_to_xy_str, write_notas=True):
        """
        This method will append to the final English and Picture files a SEO
        of rotations, one rotation at each key of the input dictionary
        `bit_pos_to_xy_str`. That input dictionary `bit_pos_to_xy_str` maps
        qubit positions (int) to xy str values that are: either 'X' or 'Y'.

        Let sigx, sigy and sigz be the Pauli matrices.

        When key=b and value='X', a rotation exp(i*sigy*pi/4) is applied to
        qubit b. This rotation is chosen because it diagonalizes sigx,
        according to the equation: exp(-i*sigy*pi/4)*sigz*exp(i*sigy*pi/4) =
        sigx.

        When key=b and value='Y', a rotation exp(-i*sigx*pi/4) is applied
        to qubit b. This rotation is chosen because it diagonalizes sigy,
        according to the equation: exp(i*sigx*pi/4)*sigz*exp(-i*sigx*pi/4) =
        sigy.

        These rotations are applied in order to convert an X or Y qubit
        measurement into the standard Z measurements which are along the
        eigenvectors of sigz. There is no need to apply any rotations when Z
        is measured because the operating qubit basis is already the
        eigenvectors of sigz.

        Parameters
        ----------
        bit_pos_to_xy_str : dict[int, str]
        write_notas : bool
            whether to write a NOTA before each rotation explaining it

        Returns
        -------
        None

        """
        for bit_pos, xy_str in bit_pos_to_xy_str.items():
            if xy_str == 'X':
                if write_notas:
                    self.write_NOTA("change |0_X>, |1_X> to |0>, |1>")
                # exp(-i*sigy*pi/4)*sigz*exp(i*sigy*pi/4) = sigx
                self.write_Ry(bit_pos, np.pi/4)
            elif xy_str == 'Y':
                if write_notas:
                    self.write_NOTA("change |0_Y>, |1_Y> to |0>, |1>")
                # exp(i*sigx*pi/4)*sigz*exp(-i*sigx*pi/4) = sigy
                self.write_Rx(bit_pos, -np.pi/4)
            else:
                assert False, "Unsupported qbit measurement. '" + \
                            xy_str + "' Should be either 'X' or 'Y'"

if __name__ == "__main__":
    def main():
        init_file_prefix = 'io_folder/coda_writer_test_init'
        fin_file_prefix = 'io_folder/coda_writer_test_fin'
        num_bits = 4

        wr = CodaSEO_writer(init_file_prefix, fin_file_prefix, num_bits)
        wr.write_NOTA('the coda follows')
        wr.write_H(2)
        wr.write_cnot(0, 1)
        wr.write_xy_measurements({1: 'X', 3: 'Y'})
        wr.close_files()
    main()
