import copy as cp


class DumbBellNeutralizer:
    """
    We call a dumbbell a controlled sigz, where sigz is the z Pauli matrix.
    sigz(a) = (-1)^{n(a)} so sigz(a)^{n(b)} = (-1)^{n(a)n(b)} where a and b
    are bit labels and n(a) = |1><1| at a = P_1(a). So these dumbbell gates
    are symmetric under exchange of a and b, just like a real dumbbell.
    Since H.sigz.H = sigx, a dumbbell can be easily converted to a CNOT,
    which is just a controlled sigx, CNOT(a, b) =  sigx(a)^{n(b)}. One
    advantage of dumbbells over CNOTs is that dumbbells with different
    endpoints always commute but CNOTs sometimes don't (when the control of
    one CNOT equals the target of the other.).

    In doing phase estimation to find the ground state of molecules,
    one uses Jordan Wigner (JW) tails (strings) to model the fermionic
    degrees of freedom. These JW tails lead to clusters of adjacent JW
    dumbbells.

    The square of a dumbbell equals one. The square of a CNOT does too,
    and Microsoft has a patent (by Wecker) that basically boils down to: We
    claim that we are the first to replace the square of a CNOT by one. I
    would like to forewarn the user of my-chemistry. It's possible that MS
    will claim in the future that this class infringes on their CNOT^2 = 1
    patent.

    This class reads an english file to identify clusters of adjacent
    dumbbells. Then it looks at each cluster and tries to cancel (
    neutralize) pairs of dumbbells with the same endpoints. Then it writes a
    new english file which is identical to the original one, except that all
    neutralized dumbbells are changed from a SIGZ to a NOTA line. The new
    english file has an X added to the name of the original english file,
    so as not to overwrite the original one.

    Attributes
    ----------
    blanked_lines : list[int]
        list of lines that will be blanked out by turning them from SIGZ to
        NOTA
    cluster_to_beg_line : dict[int:int]
        dictionary of cluster number to the number of its beginning line
    cluster_to_end_line : dict[int:int]
        dictionary of cluster number to the number of its ending line
    cur_cluster : int
        current cluster
    cur_line : int
        1-based current line number
    english_in :  _io.TextIOWrapper
        input text file that stores English description of circuit
    file_prefix : str
         beginning of the name of English file being scanned
    inside_cluster : bool
        True if the program is inside a cluster of dumbbells while reading
        input english file the first time
    line_to_dumbbell_min_max_bit : dict[int:list[int]]
        dictionary that maps line number of a dumbbell to a list with two
        entries, the beginning and ending qbit of that dumbbell
    num_bits : int
        number of qbits of whole circuit
    split_line : list[str]
        storage space for a list of strings obtained by splitting a line
    verb : bool
        set this to True if wish verbose output
    zf : bool
        set this to True if input picture file is of the zero-bit-first kind.


    """

    def __init__(self, file_prefix, num_bits,
                 zero_bit_first=False, verbose=False):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
            file must be called file_prefix + '_' + num_bits + "_eng.txt"
        num_bits : int
            total number of qubits of circuit.
        zero_bit_first : bool
        measured_bits : list(int)
            list of bits that have been measured with type 2 measurement and
            haven't been reset to |0> or |1>

        verbose : bool

        Returns
        -------

        """
        self.file_prefix = file_prefix
        self.num_bits = num_bits
        self.zf = zero_bit_first
        self.verb = verbose
        self.english_in = open(
            file_prefix + '_' + str(num_bits) + '_eng.txt', 'rt')

        self.split_line = None

        self.cur_line = 0  # 1 based
        self.cur_cluster = -1
        self.inside_cluster = False
        self.cluster_to_beg_line = {}
        self.cluster_to_end_line = {}

        # self.line_to_offset = {}
        self.line_to_dumbbell_min_max_bit = {}

        self.blanked_lines = []

        while not self.english_in.closed:
            self.scan_next_line()

        self.build_blanked_lines_list()
        self.write_blanked_files()

    def scan_next_line(self):
        """
        Scans one line. Skips over any line that doesn't start with SIGZ

        Parameters
        ----------

        Returns
        -------
        None

        """
        line = self.english_in.readline()

        if not line:
            self.english_in.close()
            return

        self.split_line = line.split()
        line_name = self.split_line[0]
        self.cur_line += 1
        if line_name == 'NOTA':
            pass
        elif line_name == "SIGZ":
            if not self.inside_cluster:
                self.cur_cluster += 1
                self.cluster_to_beg_line[self.cur_cluster] = self.cur_line
            self.inside_cluster = True
            self.cluster_to_end_line[self.cur_cluster] = self.cur_line
            self.scan_SIGZ()
        else:
            if self.inside_cluster:
                self.cluster_to_end_line[self.cur_cluster] = self.cur_line
            self.inside_cluster = False
        if self.verb:
            print('\n', line)
            print('cur_line=', self.cur_line)
            print('cur_cluster=', self.cur_cluster)
            print('inside cluster=', self.inside_cluster)
            print('cluster beg=', self.cluster_to_beg_line)
            print('cluster end=', self.cluster_to_end_line)
            print('line to dumbbell min max',
                  self.line_to_dumbbell_min_max_bit)

    def scan_SIGZ(self):
        """
        Parses line starting with "SIGZ". Stores line number and dumbbell
        endpoints

        Returns
        -------
        None

        """
        # example:
        # SIGZ	AT	1	IF	0T
        tar_bit = int(self.split_line[2])
        trol_bit = int(self.split_line[4][:-1])
        min_bit = min(tar_bit, trol_bit)
        max_bit = max(tar_bit, trol_bit)

        # self.line_to_offset[self.cur_line] = self.english_in.tell()
        self.line_to_dumbbell_min_max_bit[self.cur_line] = [min_bit, max_bit]

    def build_blanked_lines_list(self):
        """
        Builds list of lines that will be blanked, set to NOTA, because
        their dumbbells cancel out.

        Returns
        -------
        None

        """
        if self.verb:
            print('\nbuild blanked line list')
        for clu in self.cluster_to_end_line:
            a = self.cluster_to_beg_line[clu]
            b = self.cluster_to_end_line[clu]
            clu_lines = [k for k in
                         self.line_to_dumbbell_min_max_bit.keys()
                         if a <= k <= b]
            active = cp.copy(clu_lines)
            lea = len(active)
            k1 = 0
            k2 = 1
            while lea >= 2:
                if self.verb:
                    print('k1, k2, active', k1, k2, active)
                line1 = active[k1]
                line2 = active[k2]
                db1 = self.line_to_dumbbell_min_max_bit[line1]
                db2 = self.line_to_dumbbell_min_max_bit[line2]
                if db1 == db2:
                    active.remove(line1)
                    active.remove(line2)
                    lea = len(active)
                    if lea < 2:
                        break
                    else:
                        k1 = 0
                        k2 = 1
                else:
                    if k2 < lea - 1:
                        k2 += 1
                    else:
                        if k1 < lea - 2:
                            k1 += 1
                            k2 = k1 + 1
                        else:
                            break
            if self.verb:
                print('active=', active)
            self.blanked_lines += [k for k in clu_lines if k not in active]
        if self.verb:
            print('blanked lines=', self.blanked_lines)

    def write_blanked_files(self):
        """
        Writes new picture and english files identical to original ones,
        but with lines in blanked lines list set to NOTA. New files have the
        same name as the original ones except an X is added to prevent
        overwriting.

        Returns
        -------
        None

        """

        english_in = open(
            self.file_prefix + '_' + str(self.num_bits) + '_eng.txt', 'rt')
        picture_in = open(
            self.file_prefix + '_' + str(self.num_bits) +
            ('_ZF' if self.zf else '_ZL') + 'pic.txt', 'rt')

        english_out = open(
            self.file_prefix + 'X_' + str(self.num_bits) + '_eng.txt', 'wt')
        picture_out = open(
            self.file_prefix + 'X_' + str(self.num_bits) +
            ('_ZF' if self.zf else '_ZL') + 'pic.txt', 'wt')

        self.cur_line = 0  # 1 based
        while not english_in.closed:
            eng_line = english_in.readline()
            pic_line = picture_in.readline()
            if not eng_line:
                english_in.close()
                picture_in.close()
                english_out.close()
                picture_out.close()
                return

            self.split_line = eng_line.split()
            line_name = self.split_line[0]
            self.cur_line += 1
            if line_name == "SIGZ":
                if self.cur_line in self.blanked_lines:
                    english_out.write('NOTA\tNeutralized dumbbell\n')
                    picture_out.write('NOTA\tNeutralized dumbbell\n')
                else:
                    english_out.write(eng_line)
                    picture_out.write(pic_line)
            else:
                english_out.write(eng_line)
                picture_out.write(pic_line)


if __name__ == "__main__":
    prefix = 'chem_io_folder/db_neutralizer_test'
    neu = DumbBellNeutralizer(prefix, 6, verbose=True)

    prefix = 'chem_io_folder/trotter_evol'
    neu = DumbBellNeutralizer(prefix, 5, verbose=True)
