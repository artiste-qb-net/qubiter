import numpy as np


class PlaceholderManager:
    """
    The word "placeholder" (as used here and in Tensorflow) refers to a
    variable whose evaluation is postponed until a later date. In this
    class, we consider gate angle variables/placeholders.

    An object of this class is owned by SEO_reader and subclasses thereof.
    When a circuit is read, the methods owned by this class are called to
    try to replace the pound (#) prefixed gate variables by values. This
    class owns a dictionary var_num_to_degs that is used to store values of
    gate variables.


    Attributes
    ----------
    no_vars : bool
        will abort if this is True and a variable is detected
    eval_all_vars : bool
        will abort if this is True and a variable can't be evaluated
    num_vars : int
        number of variables
    cur_var_num : int
        latest variable number encountered
    var_num_to_degs : dict[int, float]
        a dict mapping variable numbers to a float for degrees. Used for
        evaluating gate variables in a SEO_reader or its sub classes.

    """

    def __init__(self, no_vars=False, eval_all_vars=True,
                 var_num_to_degs=None):
        """
        Constructor

        Parameters
        ----------
        no_vars : bool
        eval_all_vars : bool
        var_num_to_degs : dict[int, float]

        Returns
        -------


        """
        
        self.no_vars = no_vars
        self.eval_all_vars = eval_all_vars
        self.var_num_to_degs = var_num_to_degs
        self.num_vars = 0
        self.cur_var_num = -1
        
    def degs_str_to_rads(self, degs_str):
        """
        This method takes in a string degs_str which might represent either
        degrees or a variable of the form '#'  + str( int). If degs_str
        starts with the char '#', the method tries to find a value for that
        variable in the self.var_num_to_degs dictionary. If it finds a
        value there, the method returns value*pi/180. If no value can be
        found, the method outputs the input string unchanged.

        Parameters
        ----------
        degs_str : str

        Returns
        -------
        float | str

        """
        if degs_str[0] == '#':
            assert not self.no_vars, 'no circuit variables allowed'
            # remove first char #
            var_num = int(degs_str[1:])
            if var_num > self.cur_var_num:
                self.num_vars += 1
            self.cur_var_num = var_num
            if self.var_num_to_degs:
                key_exists = var_num in self.var_num_to_degs.keys()
                if key_exists:
                    return self.var_num_to_degs[var_num] * np.pi / 180
                else:  # key doesn't exist
                    if self.eval_all_vars:
                        assert False, 'no value for variable #' + str(var_num)
                    else:
                        return degs_str
            else:  # no self.var_num_to_degs dict
                if self.eval_all_vars:
                    assert False, 'no self.var_num_to_degs dict'
                else:
                    return degs_str
        else:
            return float(degs_str)*np.pi/180

    @staticmethod
    def assert_no_vars_in(rads):
        """
        The default value for the flag self.self.eval_all_vars is True. As
        long as this flag is not changed, overrides of SEO_reader:use_
        methods will get rads arguments that contain no pound prefixed
        strings. If, however, this flag is changed to False, then should use
        this test inside the use_ method if that method assumes that all its
        rads arguments are floats, not strings.


        Parameters
        ----------
        rads : float | list[]

        Returns
        -------
        None

        """
        if isinstance(rads, str):
            assert False
        else:
            for x in rads:
                assert not isinstance(x, str)

if __name__ == "__main__":
    from SEO_writer import *
    from SEO_reader import *
    from EchoingSEO_reader import *
    from SEO_simulator import *

    def main():
        # This produces English and Picture files as usual.
        # Note that 3 angles have been entered as None.
        # In the English file, you will see #0, #1, #2 where
        # the numerical values of those angles would have been.
        num_bits = 4
        file_prefix = 'io_folder/gate_vars_test'
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_Rx(2, rads=None)
        wr.write_Rn(3, rads_list=[None, np.pi/2, None])
        wr.write_cnot(2, 3)
        wr.close_files()

        # this produces a log file with the given file prefix,
        # wherein the number of gate variables is reported
        SEO_reader(file_prefix, num_bits, write_log=True)

        # partial substitution, this creates new files
        # with #0=30, #1=60 but #2 still undecided
        vman = PlaceholderManager(eval_all_vars=False,
                                  var_num_to_degs={0: 30, 1: 60})
        wr = SEO_writer(file_prefix + '_eval01', emb)
        EchoingSEO_reader(file_prefix, num_bits, wr,
                          vars_manager=vman)

        # this runs the simulator after substituting
        # #0=30, #1=60, #2=45
        vman = PlaceholderManager(
            var_num_to_degs={0: 30, 1: 60, 2: 45})
        sim = SEO_simulator(file_prefix, num_bits, verbose=True,
                            vars_manager=vman)
        print("\n----------------------------------------")
        StateVec.describe_st_vec_dict(sim.cur_st_vec_dict)

    main()
