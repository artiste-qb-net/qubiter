import numpy as np
from SEO_writer import *


class PlaceholderManager:
    """
    The word "placeholder" (as used here and in Tensorflow) refers to a
    variable whose evaluation is postponed until a later date. In this
    class, we consider gate angle variables/placeholders.

    Variable names must be of form '#3' or '-#3' with 3 replaced by any
    other non-negative int. The 3 is called the variable number.

    An object of this class is owned by SEO_reader and subclasses thereof.
    When a circuit is read, the methods owned by this class are called to
    try to replace the pound (#) prefixed gate variables by values. This
    class owns a dictionary var_num_to_degs that is used to store values of
    gate variables.

    Attributes
    ----------
    eval_all_vars : bool
        will abort if this is True and a variable can't be evaluated
    no_vars : bool
        will abort if this is True and a variable is detected
    var_num_to_degs : dict[int, float]
        a dict mapping variable numbers to a float for degrees. Used for
        evaluating gate variables in a SEO_reader or its child classes.
    var_nums_list : list[int]
        a list of all distinct numbers of the variables encountered

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
        self.var_nums_list = []
        
    def degs_str_to_rads(self, degs_str):
        """
        This method takes in a string degs_str which might represent either
        a str(float) expressing degrees or a legal var name.

        If degs_str is not a legal_var_name, the method assumes (This
        assumption is okay because the file was written by SEO_writer so it
        is always well formed) it's a str( float) and returns float(
        deg_str)*pi/180.

        If degs_str is a legal variable name, the method tries to find a
        value for that variable in the self.var_num_to_degs dictionary. If
        it finds a value there, the method returns value*pi/180. If no value
        can be found, the method outputs the input string unchanged (unless
        eval_all_vars=True in which case the method aborts).

        Parameters
        ----------
        degs_str : str

        Returns
        -------
        float | str

        """
        if SEO_writer.is_legal_var_name(degs_str):
            assert not self.no_vars, 'no circuit variables allowed'
            if degs_str[0] == '#':
                var_num = int(degs_str[1:])
                sign = +1
            else:
                var_num = int(degs_str[2:])
                sign = -1
            if var_num not in self.var_nums_list:
                self.var_nums_list.append(var_num)
            if self.var_num_to_degs:
                key_exists = var_num in self.var_num_to_degs.keys()
                if key_exists:
                    return sign*self.var_num_to_degs[var_num]*np.pi/180
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

if __name__ == "__main__":
    from SEO_writer import *
    from SEO_reader import *
    from EchoingSEO_reader import *
    from SEO_simulator import *

    def main():
        # This produces English and Picture files as usual. Note that 4
        # angles have been entered as '#1', '-#1', '#2', '#3', all legal
        # variable names . In the English file, you will see those legal
        # names where the numerical values of those angles would have been.

        num_bits = 4
        file_prefix = 'io_folder/gate_vars_test'
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_Rx(2, rads=np.pi/7)
        wr.write_Rx(1, rads='#2')
        wr.write_Rn(3, rads_list=['#1', '-#1', '#3'])
        wr.write_cnot(2, 3)
        wr.close_files()

        # this produces a log file with the given file prefix,
        # wherein the a list of distinct gate numbers encountered is reported
        SEO_reader(file_prefix, num_bits, write_log=True)

        # partial substitution, this creates new files
        # with #1=30, #2=60 but #3 still undecided
        vman = PlaceholderManager(eval_all_vars=False,
                                  var_num_to_degs={1: 30, 2: 60})
        wr = SEO_writer(file_prefix + '_eval01', emb)
        EchoingSEO_reader(file_prefix, num_bits, wr,
                          vars_manager=vman)

        # this runs the simulator after substituting
        # #1=30, #2=60, #3=90
        vman = PlaceholderManager(
            var_num_to_degs={1: 30, 2: 60, 3: 90})
        sim = SEO_simulator(file_prefix, num_bits, verbose=True,
                            vars_manager=vman)
        print("\n----------------------------------------")
        StateVec.describe_st_vec_dict(sim.cur_st_vec_dict)

    main()
