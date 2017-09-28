

class SEO_pre_reader:
    """
    This class' constructor scans (pre-reads) an English file (a type of txt
    file). It skips all lines except the ones starting with LOOP or NEXT. By
    doing so, it collects information about LOOPs (where they start and end,
    their id number and their number of reps). This class is inherited by
    class SEO_reader which uses its information to handle embedded loops.

    See docstring for class SEO_writer for more info about English files.

    Attributes
    ----------
    english_in : _io.TextIOWrapper
        file object for input text file that stores English description of
        circuit
    file_prefix : str
        beginning of the name of English file being scanned
    loop_queue : list[int]
        a queue of loops labelled by their id number
    loop_to_reps : dict[int, int]
        a dictionary mapping loop number TO total number of repetitions of
        loop
    loop_to_start_line : dict[int, int]
        a dictionary mapping loop number TO loop line + 1
    loop_to_start_offset : dict[int, int]
        a dictionary mapping loop number TO offset of loop's start
    num_bits : int
        number of qubits in whole circuit
    split_line : list[str]
        storage space for a list of strings obtained by splitting a line
    tot_num_lines : int
        number of lines in English file

    """

    def __init__(self, file_prefix, num_bits):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
            file must be called file_prefix + '_' + num_bits + "_eng.txt"
        num_bits : int
            total number of qubits of circuit.

        Returns
        -------

        """
        self.file_prefix = file_prefix
        self.num_bits = num_bits
        self.english_in = open(
            file_prefix + '_' + str(num_bits) + '_eng.txt', 'rt')
        self.split_line = None

        self.tot_num_lines = 0
        self.loop_to_start_offset = {}
        self.loop_to_start_line = {}
        self.loop_to_reps = {}
        self.loop_queue = []

        while not self.english_in.closed:
            self.scan_next_line()

    def scan_next_line(self):
        """
        Scans one line. Skips over any line that doesn't start with LOOP or
        NEXT.

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
        self.tot_num_lines += 1
        if line_name == "LOOP":
            self.scan_LOOP()
        elif line_name == "NEXT":
            self.scan_NEXT()
        else:
            pass

    def scan_LOOP(self):
        """
        Parses and uses line starting with "LOOP".

        Returns
        -------
        None

        """
        # example:
        # LOOP 5 REPS: 2
        loop_num = int(self.split_line[1])
        reps = int(self.split_line[3])
        # print(self.split_line)
        assert loop_num not in self.loop_to_reps.keys(),\
            "this loop number has occurred before"
        self.loop_to_start_offset[loop_num] = self.english_in.tell()
        self.loop_to_start_line[loop_num] = self.tot_num_lines + 1
        self.loop_to_reps[loop_num] = reps
        self.loop_queue += [loop_num]

    def scan_NEXT(self):
        """
        Parses and uses line starting with "NEXT".

        Returns
        -------
        None

        """
        # example:
        # NEXT 5
        loop_num = int(self.split_line[1])
        # print(self.split_line)
        if not self.loop_queue:
            assert False, "unmatched NEXT"
        if loop_num == self.loop_queue[-1]:
            del self.loop_queue[-1]
        else:
            assert False, "improperly nested loops"

if __name__ == "__main__":
    print(5)
