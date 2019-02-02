from SEO_reader import *
from LoopyPlaceholderManager import *


class LoopFileGenerator(SEO_reader):
    """
    This class is important to you only if you are interested in using loops
    in the English file.

    This class reads an English file (a .txt text file) with nested loops in
    it and it writes a Loop File (a .py python script). Like all
    SEO_reader's, this class has a PlaceholderManager call vars_manager (
    variables manager), but for this class, the vars_manager is of the
    special type LoopyPlaceholderManager.

    Attributes
    ----------
    indentation_li : list[int]
        one item list containing an int that stores the current number
        of indentation spaces in the Loop File being written
    loop_out : _io.TextIOWrapper
        out stream to Loop File being generated.

    """

    def __init__(self, file_prefix, num_bits, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        kwargs : dict

        Returns
        -------

        """
        self.loop_out = open(
            file_prefix + '_' + str(num_bits) + '_loop.py', 'wt')

        # to pass integer by reference, must put inside a list
        self.indentation_li = [0]

        vman = LoopyPlaceholderManager(self.loop_out,
                            self.indentation_li, eval_all_vars=False)

        vman.write_loop_file_prelude()

        # Must set write_log to True or else SEO_reader will abort. However,
        #  we also override do_log in this class to prevent log from being
        # written.
        SEO_reader.__init__(self, file_prefix, num_bits, write_log=True,
                            vars_manager=vman, **kwargs)

        vman.write_loop_file_ending()
        self.loop_out.close()

    def do_log(self):
        """
        We override this to prevent writing of log, because when we do,
        we get a wrong count for the number of elem ops, because this
        SEO_reader is performing a flat reading of the English file (i.e.,
        not jumping back at NEXTs)

        Returns
        -------
        None

        """
        pass

    def use_LOOP(self, loop_num, nreps):
        """
        This method overrides method of same name in parent class SEO_reader.

        Its purpose is to write a "for" line in the Loop File. It also
        increases the indentation in the Loop File by 4 spaces. This gives
        Python style indentation for for-loops in Loop File.

        Parameters
        ----------
        loop_num : int
        nreps : int

        Returns
        -------
        None

        """
        self.loop_out.write(' ' * self.indentation_li[0] +
            'for j' + str(loop_num) + ' in range(' + str(nreps) + '):' +
            '  # line ' + str(self.line_count) + "\n")
        self.indentation_li[0] += 4

    def use_NEXT(self, loop_num):
        """
        This method overrides method of same name in parent class
        SEO_reader. It reduces the indentation in the Loop File by 4 spaces.
        This gives Python style indentation for for-loops in Loop File.

        The use_NEXT in the parent class is responsible for jumping back
        from NEXT line to LOOP line if the total number of reps (nreps)
        hasn't been completed. This override doesn't do that jump. This
        gives a flat reading of the English file, without jumping back at
        NEXTs.

        Parameters
        ----------
        loop_num : int

        Returns
        -------
        None

        """
        self.indentation_li[0] -= 4


if __name__ == "__main__":
    from SEO_writer import *
    from SEO_simulator import *
    from StateVec import *

    def main():
        file_prefix = 'io_folder/loop_gen_test'
        num_bits = 4

        # write the English and Picture files
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_controlled_one_bit_gate(0,
                                         Controls.new_knob(num_bits, 2, False),
                                         OneBitGates.rot_ax,
                                         ['#1', 1])
        wr.write_LOOP(20, nreps=2)
        wr.write_controlled_one_bit_gate(1,
                                         Controls.new_knob(num_bits, 2, False),
                                         OneBitGates.rot_ax,
                                         ['-my_fun1#1#2', 2])
        wr.write_LOOP(10, nreps=4)
        wr.write_controlled_one_bit_gate(2,
                                         Controls.new_knob(num_bits, 3, True),
                                         OneBitGates.rot,
                                         ['-#1*.5', '#2',  '-my_fun3#3'])
        wr.write_NEXT(10)
        wr.write_controlled_one_bit_gate(1,
                                         Controls.new_knob(num_bits, 2, False),
                                         OneBitGates.rot_ax,
                                         ['my_fun1#1#2', 2])
        wr.write_NEXT(20)
        wr.write_controlled_one_bit_gate(0,
                                 Controls.new_knob(num_bits, 2, False),
                                 OneBitGates.rot_ax,
                                 ['#1*.3', 1])
        wr.write_one_bit_gate(1, OneBitGates.rot_ax, ['my_fun#1', 1])
        wr.close_files()

        # write a log
        SEO_reader(file_prefix, num_bits, write_log=True)

        # write a Loop File
        LoopFileGenerator(file_prefix, num_bits)

        # read a Loop xfile and do simulation
        sim = SEO_simulator(file_prefix, num_bits, verbose=False,
                            xfile_num=1)
        print("\n----------------------------------------")
        StateVec.describe_st_vec_dict(sim.cur_st_vec_dict)
    main()


