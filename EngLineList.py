from EchoingSEO_reader import *


class EngLineList:
    """
    Eng=English File

    In programs like PyQuil (Rigetti) and Cirq (Google), circuits (aka
    programs) are stored in memory essentially as Python lists of gates.
    Lists of this type can be sliced, combined, etc. They are very
    convenient, more flexible than English files, when one wants to load the
    whole circuit into memory, so as to be able to access any gate of it at
    will. This is convenient, for instance, when doing circuit optimizations
    (i.e., replacing the circuit by an equivalent but hopefully shorter one,
    what IBM qiskit calls "transpiling").

    In this class, we support Qubiter's version of PyQuil's and Cirq's gate
    lists. In Qubiter, we use simply a Python list of the lines, stored as
    strings, of the circuit's English file.

    """
    def __init__(self, line_list, num_bits):
        """
        Costructor

        Returns
        -------

        """
        self.line_list = line_list
        self.num_bits = num_bits

    def print(self):
        """
        prints self.line_list, one item per line

        Returns
        -------
        None

        """
        print('\n'.join(self.line_list))

    @staticmethod
    def eng_file_to_line_list(file_prefix, num_bits):
        """
        This static method reads an English file with file prefix
        `file_prefix` and it returns a list of its line strings.

        Parameters
        ----------
        file_prefix : str
        num_bits : int

        Returns
        -------
        list[str]

        """
        path = file_prefix + '_' + str(num_bits) + '_eng.txt'
        with open(path, 'r') as f:
            line_list = [line.rstrip('\n') for line in f]
        return line_list

    @staticmethod
    def line_list_to_eng_and_pic_files(line_list, file_prefix, num_bits):
        """
        This method does the reverse of eng_file_to_line_list(). It writes
        both an English file and a Picture file with file prefix=file_prefix.

        Parameters
        ----------
        line_list : list[str]
        file_prefix : str
        num_bits : int

        Returns
        -------
        None

        """
        end_str = '_' + str(num_bits) + '_eng.txt'
        with open(file_prefix + end_str, 'w') as f:
            for line in line_list:
                f.write(line + '\n')

        # this writes a Picture file from the English file just created
        EchoingSEO_reader.pic_file_from_eng_file(file_prefix, num_bits)

    def write_eng_and_pic_files(self, file_prefix):
        """
        This method does the same as line_list_to_eng_and_pic_files(),
        except it is called by an object of this class instead of being a
        static method.

        Parameters
        ----------
        file_prefix :

        Returns
        -------
        None

        """
        EngLineList.line_list_to_eng_and_pic_files(
            self.line_list, file_prefix, self.num_bits)


if __name__ == "__main__":
    def main():
        num_bits = 4
        file_prefix = 'io_folder/eng_line_list_test'
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_Rx(2, rads=np.pi)
        wr.write_Rn(3, rads_list=[np.pi/2, np.pi/2, np.pi/2])
        wr.write_cnot(2, 3)
        wr.close_files()

        list1 = EngLineList.eng_file_to_line_list(file_prefix, num_bits)
        list2 = list1 + list1[1:]
        ell = EngLineList(list2, num_bits)
        ell.print()
        ell.write_eng_and_pic_files(file_prefix + '_sum')
    main()
