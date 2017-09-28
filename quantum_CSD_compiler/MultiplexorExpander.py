from CktEmbedder import *  # if don't include this, can't find Controls
from EchoingSEO_reader import *
from quantum_CSD_compiler.MultiplexorSEO_writer import *


class MultiplexorExpander(EchoingSEO_reader):
    """
    This class is a child of EchoingSEO_reader. The class reads any
    previously created Qubiter English file and it writes new English &
    Picture files wherein every line of the original English file that
    doesn't start with MP_Y (multiplexor) is echoed faithfully whereas lines
    which do start with MP_Y are expanded via class MultiplexorSEO_writer
    into a sequence in the style specified as input to the class constructor.

    So, to expand English & Picture files that contain DIAG and MP_Y lines
    into just CNOTs and qubit rotations, first use this class and the
    analogous one for DIAG, then use class CktExpander.

    If the input English file has in_file_prefix as file prefix, then the
    output English & Picture files have as file prefix in_file_prefix + '_X1',
    assuming that '_X' + str(k) for some integer k is not already the ending
    of in_file_prefix. If it is, then the ending is changed to '_X' + str(
    k+1).


    Attributes
    ----------


    """

    def __init__(self, file_prefix, num_bits, style, gbit_list=None):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        style : str
        gbit_list : list(int)

        Returns
        -------
        None

        """

        if gbit_list:
            num_gbits = len(gbit_list)
        else:
            num_gbits = 0

        # default embedder and rad_angles
        emb = CktEmbedder(num_bits, num_bits)
        rad_angles = None
        out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
        wr = MultiplexorSEO_writer(out_file_prefix, emb,
            style, rad_angles, num_gbits=num_gbits)

        EchoingSEO_reader.__init__(self, file_prefix, num_bits, wr, style,
                                   gbit_list)

        self.wr.close_files()

    def emb_for_plexor(self, tar_bit_pos, controls):
        """
        This is an internal function used inside the function use_MP_Y().
        The function returns emb, nt, nf to be used as arguments of a
        MultiplexorSEO_writer that will be used to expand the MP_y line
        currently being considered. emb is a circuit embedder, nt is the
        number of T bits and nf is the number of F bits detected in the
        input argument 'controls'.

        Parameters
        ----------
        tar_bit_pos : int
            target bit position of multiplexor currently being considered.
        controls : Controls
            controls of the MP_Y currently being considered.

        Returns
        -------
        CktEmbedder, int, int

        """
        T_bpos = []
        F_bpos = []
        MP_bpos = []
        for bpos, kind in controls.bit_pos_to_kind.items():
            # int is subclass of bool
            # so isinstance(x, int) will be true for x bool too!
            if isinstance(kind, bool):
                if kind:
                    T_bpos.append(bpos)
                else:
                    F_bpos.append(bpos)
            else:
                MP_bpos.append(bpos)
        T_bpos.sort()
        F_bpos.sort()
        MP_bpos.sort()
        if self.gbit_list:
            g_bpos = self.gbit_list.sort()
        else:
            g_bpos = []
        bit_map = T_bpos + F_bpos + MP_bpos + [tar_bit_pos] + g_bpos
        # print("bit_map", bit_map)
        assert len(bit_map) == len(set(bit_map)),\
            "bits used to define multiplexor are not unique"
        assert len(bit_map) <= self.num_bits

        nt = len(T_bpos)
        nf = len(F_bpos)
        emb = CktEmbedder(num_bits, num_bits, bit_map)
        return emb, nt, nf

    def use_MP_Y(self, tar_bit_pos, controls, rad_angles):
        """
        This is an override of a function in the parent class
        EchoingSEO_reader. This is the only use_ function of this class that
        doesn't simply echo its input line. This function does most of its
        work inside the Multiplexor_SEO_writer.write() function that it calls.

        Parameters
        ----------
        tar_bit_pos : int
        controls : Controls
        rad_angles : list[float]

        Returns
        -------
        None

        """

        emb, nt, nf = self.emb_for_plexor(tar_bit_pos, controls)
        self.wr.emb = emb
        self.wr.rad_angles = rad_angles
        self.wr.num_T_trols = nt
        self.wr.num_F_trols = nf
        # style and num_gbits for wr are set by constructor

        self.wr.write()
        # revert to default embedder
        self.wr.emb = CktEmbedder(self.num_bits, self.num_bits)

if __name__ == "__main__":
    num_bits = 6
    file_prefix = "../io_folder/plexor_test_one_line"
    style = 'exact'
    xer = MultiplexorExpander(file_prefix, num_bits, style)
