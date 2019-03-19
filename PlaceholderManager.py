import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np
else:
    import autograd.numpy as np
# from SEO_writer import *
import utilities_gen as ug
from collections import defaultdict


class PlaceholderManager:
    """
    The word "placeholder" (as used here and in Tensorflow) refers to a
    variable whose evaluation is postponed until a later time. In this
    class, we consider gate angle variables-placeholders.

    Legal variable names are described in the docstring of the function
    is_legal_var_name() in this class.

    An object of this class is owned by SEO_reader and subclasses thereof.
    When a circuit is read, the methods owned by this class are called to
    try to replace the placeholders by values. This class owns a dictionary
    `var_num_to_rads` that is used to store values of all angle variables
    used in all placeholders. The class also owns a dictionary
    `fun_name_to_fun` that is used to store functions used in "functional"
    type placeholders. Together, these 2 dictionaries are used to find a
    float value ("resolve") each placeholder.

    Attributes
    ----------
    all_fun_names : list[str]
        a list of all the distinct function names encountered in circuit
    all_var_nums : list[int]
        a list of all distinct variable numbers encountered in circuit
    eval_all_vars : bool
        will abort if this is True and a variable can't be evaluated
    fun_name_to_fun : dict[str, function]
        dictionary mapping a function name to a function
    fun_name_to_hist : dict[str, list[function|None]]
        dictionary mapping a function name to its history (a list of
        function or None). This is only used when using loops in the English
        file. See LoopFileGenerator and LoopyPlaceholder classes for more
        info about this variable.
    fun_name_to_use_count: dict[str, int]
        dictionary mapping a function name to its current use count (the
        number of times it has been used. Repeat appearances inside a loop
        or nested loops are considered as different uses.)
    no_vars : bool
        will abort if this is True and a variable is detected
    var_num_to_hist : dict[int, list[float|None]]
        dictionary mapping a variable number to its history (a list of float
        or None). This is only used when using loops in the English file. See
        LoopFileGenerator and LoopyPlaceholder classes for more info about
        this variable.
    var_num_to_rads : dict[int, float]
        dictionary mapping a variable number to its value in radians
    var_num_to_use_count : dict[int, int]
        dictionary mapping a variable number to its current use count (the
        number of times it has been used. Repeat appearances inside a loop
        or nested loops are considered as different uses.)

    """

    def __init__(self, no_vars=False, eval_all_vars=True,
                 var_num_to_rads=None, fun_name_to_fun=None):
        """
        Constructor

        Parameters
        ----------
        no_vars : bool
        eval_all_vars : bool
        var_num_to_rads : dict[int, float]
            a dict mapping variable numbers to a float for radians. Used by
            a SEO_reader or its children to evaluate hash prefixed angle
            variables in a placeholder.
        fun_name_to_fun : dict[str, function]
            a dictionary mapping each function name to a function that
            returns a float that stands for radians. Used by a SEO_reader or
            its children to evaluate function names in a functional
            placeholder.

        Returns
        -------


        """
        
        self.no_vars = no_vars
        self.eval_all_vars = eval_all_vars
        self.var_num_to_rads = var_num_to_rads
        self.fun_name_to_fun = fun_name_to_fun

        # will be filled by this class
        self.all_var_nums = []
        self.all_fun_names = []
        self.var_num_to_hist = defaultdict(list)
        self.fun_name_to_hist = defaultdict(list)
        self.var_num_to_use_count = {}
        self.fun_name_to_use_count = {}

    @staticmethod
    def is_legal_var_name(name):
        """
        This method returns True iff name is a legal variable name.

        Legal variable names must be of form `#3` or `-#3` or `#3*.5` or
        `-#3*.5` where 3 can be replaced by any non-negative int, and .5 can
        be replaced by anything that can be an argument of float() without
        throwing an exception. In this example, the 3 that follows the hash
        character is called the variable number

        NEW! (functional placeholder variables)

        Now legal variable names can ALSO be of the form `my_fun#1#2` or
        `-my_fun#1#2`, where

        * the 1 and 2 can be replaced by any non-negative integers and there
        might be any number > 0 of hash variables. Thus, there need not
        always be precisely 2 hash variables as in the example.

        * `my_fun` can be replaced by the name of any function with one or
        more input floats (2 inputs in the example), as long as the first
        character of the function's name is a lower case letter.

        The strings `my_fun#1#2` or `-my_fun#1#2` indicate than one wants to
        use for the angle being replaced, the values of `my_fun(#1, #2)` or
        `-my_fun(#1, #2)`, respectively, where the inputs #1 and #2 are
        floats standing for radians and the output is also a float standing
        for radians.

        Parameters
        ----------
        name : str

        Returns
        -------
        bool

        """
        if not isinstance(name, str) or len(name) < 2:
            return False
        nom = name
        if name[0] == "-":
            nom = name[1:]
        if len(nom) < 2:
            return False

        if nom[0] == "#":
            tokens = nom[1:].split('*')
            if len(tokens) > 2:
                return False
            if not ug.is_non_neg_int(tokens[0]):
                return False
            if len(tokens) == 2:
                try:
                    # print(',,', nom[star_ind+1:])
                    x = float(tokens[1])
                except ValueError:
                    return False

        elif nom[0].islower():
            tokens = nom.split('#')
            if len(tokens) < 2:
                return False
            for t in tokens[1:]:
                if not ug.is_non_neg_int(t):
                    return False
        else:
            return False
        return True

    @staticmethod
    def is_functional_var(var_name):
        """
        Assumes var_name is a legal variable name (doesn't check). Returns
        True iff var_name is a functional variable name

        Parameters
        ----------
        var_name : str

        Returns
        -------
        bool

        """
        assert isinstance(var_name, str) and len(var_name) >= 2
        if var_name[0] == '-':
            if not var_name[1].islower():
                return False
        else:
            if not var_name[0].islower():
                return False
        return True

    @staticmethod
    def get_leg_var_sign(var_name):
        """
        Assumes var_name is a legal variable name (doesn't check). Returns
        sign of variable, either +1 or -1

        Parameters
        ----------
        var_name : str

        Returns
        -------
        int

        """
        return -1 if var_name[0] == '-' else 1

    @staticmethod
    def get_leg_var_var_nums(var_name):
        """
        Assumes var_name is a legal variable name (doesn't check). Returns a
        list of the variable numbers

        Parameters
        ----------
        var_name : str

        Returns
        -------
        list[int]

        """
        is_fun_var = PlaceholderManager.is_functional_var(var_name)

        sign = PlaceholderManager.get_leg_var_sign(var_name)

        # this gives -1 if "*" not found
        star_ind = var_name.find("*")
        if star_ind < 0:
            star_ind = len(var_name)

        # this excludes sign char and scale factor if there is one
        tokens = var_name[(0 if sign == 1 else 1): star_ind].split("#")
        # if var_name starts with # or -#, first token will be
        # empty string so remove it
        if tokens[0] == '':
            tokens = tokens[1:]

        # create list of variable number strings
        # this excludes function name if there is one
        var_num_tokens = tokens
        if is_fun_var:
            var_num_tokens = tokens[1:]
        # print("--", var_num_tokens)
        return [int(t) for t in var_num_tokens]

    @staticmethod
    def get_leg_var_fun_name(var_name):
        """
        Assumes var_name is a legal variable name (doesn't check). Returns
        the function name if there is one or None otherwise.

        Parameters
        ----------
        var_name : str

        Returns
        -------
        str | None

        """
        is_fun_var = PlaceholderManager.is_functional_var(var_name)
        if not is_fun_var:
            return None

        sign = PlaceholderManager.get_leg_var_sign(var_name)

        # this finds index of first "#"
        hash_ind = var_name.find("#")

        return var_name[0 if sign == 1 else 1: hash_ind]

    @staticmethod
    def get_leg_var_scale_fac(var_name):
        """
        Assumes var_name is a legal variable name (doesn't check). If it's a
        functional variable, it returns None. If it's not a functional
        variable, it returns the scale factor, the float after the *,
        if there is a * or a 1 if there is no *

        Parameters
        ----------
        var_name : str

        Returns
        -------
        float | None

        """
        is_fun_var = PlaceholderManager.is_functional_var(var_name)
        if is_fun_var:
            return None

        # this finds index of first "*", returns -1 if there is none
        star_ind = var_name.find("*")

        if star_ind < 0:
            return 1

        return float(var_name[star_ind+1:])

    def degs_str_to_rads(self, degs_str, line_count):
        """
        This method takes in a string degs_str which might represent either
        a str() of a float expressing degrees or a legal var name.

        If degs_str is not a legal_var_name, the method assumes (This
        assumption is okay because the file was written by SEO_writer so it
        is always well formed) it's a str() of a float and returns float(
        deg_str)*pi/180.

        If degs_str is a legal variable name, the method tries to resolve it
        into a rads_value using the 2 dictionaries: self.var_num_to_rads and
        self.fun_name_to_fun. If it finds a rads_value, the method returns
        rads_value. If no rads_value can be found, the method outputs the
        input string unchanged (unless eval_all_vars=True in which case the
        method aborts).

        Parameters
        ----------
        degs_str : str
        line_count : int
            this is the line_count in the English file being read. This is
            not used by this function but is used by overrides of it like
            the one in LoopyPlaceholderManager

        Returns
        -------
        float | str

        """
        if not PlaceholderManager.is_legal_var_name(degs_str):
            try:
                y = float(degs_str)*np.pi/180
            except:
                assert False, 'tried to convert ' + str(degs_str)\
                 + ' to a float'

            return y
        else:  # is a legal variable name
            assert not self.no_vars, 'no circuit variables allowed'

            is_fun_var = PlaceholderManager.is_functional_var(degs_str)
            sign = PlaceholderManager.get_leg_var_sign(degs_str)
            scale_fac = PlaceholderManager.get_leg_var_scale_fac(degs_str)
            fun_name = PlaceholderManager.get_leg_var_fun_name(degs_str)
            token_var_nums = PlaceholderManager.get_leg_var_var_nums(degs_str)

            # store var numbers and update their use_count
            for var_num in token_var_nums:
                if var_num not in self.all_var_nums:
                    self.all_var_nums.append(var_num)
                if var_num not in self.var_num_to_use_count:
                    self.var_num_to_use_count[var_num] = 1
                else:
                    self.var_num_to_use_count[var_num] += 1

            # store function name and update its use_count
            if fun_name:
                if fun_name not in self.all_fun_names:
                    self.all_fun_names.append(fun_name)
                if fun_name not in self.fun_name_to_use_count.keys():
                    self.fun_name_to_use_count[fun_name] = 1
                else:
                    self.fun_name_to_use_count[fun_name] += 1

            if self.var_num_to_rads is None:
                self.var_num_to_rads = {}
            if self.fun_name_to_fun is None:
                self.fun_name_to_fun = {}

            var_num_to_rads = self.var_num_to_rads
            fun_name_to_fun = self.fun_name_to_fun

            # if there's a var_num_to_hist, ignore input value of
            # self.var_num_to_rads. Refill it using the history
            if self.var_num_to_hist:
                for var_num in token_var_nums:
                    var_num_to_rads[var_num] = \
                        self.var_num_to_hist[var_num][
                            self.var_num_to_use_count[var_num] - 1]

            # if there's a fun_name_to_hist, ignore input value of
            # self.fun_name_to_rads. Refill it using the history
            if self.fun_name_to_hist:
                if fun_name:
                    fun_name_to_fun[fun_name] = \
                        self.fun_name_to_hist[fun_name][
                            self.fun_name_to_use_count[fun_name] - 1]

            if var_num_to_rads:
                if not is_fun_var:
                    # if its not a functional placeholder,
                    # there is only one fun var
                    var_num = token_var_nums[0]
                    key_exists = var_num in var_num_to_rads.keys()
                    if key_exists:
                        return sign*var_num_to_rads[var_num]*scale_fac
                    else:  # key doesn't exist
                        if self.eval_all_vars:
                            assert False, "eval_all_vars is True and " +\
                                'var_num_to_rads has no value for ' + \
                                'variable #' + str(var_num)
                        else:
                            return degs_str
                else:  # is_fun_var
                    if fun_name_to_fun:
                        all_var_keys_exist = \
                            all([var_num in var_num_to_rads.keys()
                                          for var_num in token_var_nums])
                        fun_key_exists = \
                            fun_name in fun_name_to_fun.keys()

                        if all_var_keys_exist and fun_key_exists:
                            arg_float_list = [var_num_to_rads[var_num]
                                              for var_num in token_var_nums]
                            fun = fun_name_to_fun[fun_name]
                            return sign*fun(*arg_float_list)

                        else:
                            if self.eval_all_vars:
                                assert False, "eval_all_vars is True and " +\
                                    "can't resolve the function " +\
                                    fun_name + \
                                    " acting on variables " + \
                                    str(token_var_nums)
                            else:
                                return degs_str
                    else:  # no fun_name_to_fun
                        if self.eval_all_vars:
                            assert False, 'eval_all_vars is True and ' \
                                'there is no fun_name_to_fun dictionary' + \
                                " so can't evaluate " + degs_str
                        else:
                            return degs_str

            else:  # no var_num_to_rads
                if self.eval_all_vars:
                    assert False, 'eval_all_vars is True and ' \
                        'there is no var_num_to_rads dictionary' + \
                        " so can't evaluate " + degs_str
                else:
                    return degs_str

    @staticmethod
    def have_resolved_history(hist):
        """
        This function tries to resolve a history. It returns True iff it
        succeeds.

        A history of a variable is a list of all the values it will assume
        all the times it is used in the circuit (each repetition in a loop
        or nested loops is counted as a different use.) For a hash
        placeholder, a history is a list of floats or None's. For a
        functional placeholder, a history is a list of functions or None's.

        Resolving a history successfully means replacing all None's by a
        value that is not a None. This is accomplished by moving in order of
        ascending index on the list and replacing each None by its
        predecessor. The very first use cannot be a None.

        Parameters
        ----------
        hist : list[float | function | None]

        Returns
        -------
        bool

        """
        for time, event in enumerate(hist):
            if event is None:
                if time == 0:
                    return False
                else:
                    hist[time] = hist[time-1]
        return True

    def resolve_all_histories(self):
        """
        This function tries to resolve the histories of all hash and
        functional placeholder variables. It aborts if it fails to resolve
        any of those histories.

        Returns
        -------
        None

        """
        for val_num in self.all_var_nums:
            assert PlaceholderManager.have_resolved_history(
                self.var_num_to_hist[val_num])

        for fun_name in self.all_fun_names:
            assert PlaceholderManager.have_resolved_history(
                self.fun_name_to_hist[fun_name])

if __name__ == "__main__":
    from SEO_writer import *
    from SEO_reader import *
    from EchoingSEO_reader import *
    from SEO_simulator import *

    def main():
        # We begin by writing a simple circuit with 4 qubits. As usual,
        # the following code will write an English and a Picture file in the
        #  io_folder directory. Note that some rotation angles have been
        # entered into the write() Python functions as legal variable names
        # instead of floats. In the English file, you will see those legal
        # names where the numerical values of those angles would have been.

        num_bits = 4
        file_prefix = 'io_folder/placeholder_test'
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(file_prefix, emb)
        wr.write_Rx(2, rads=np.pi/7)
        wr.write_Rx(1, rads='#2*.5')
        wr.write_Rx(1, rads='my_fun1#2')
        wr.write_Rn(3, rads_list=['#1', '-#1*3', '#3'])
        wr.write_Rx(1, rads='-my_fun2#2#1')
        wr.write_cnot(2, 3)
        wr.close_files()

        # Simply by creating an object of the class SEO_reader with the flag
        #  `write_log` set equal to True, you can create a log file which
        # contains
        # (1) a list of distinct variable numbers
        # (2) a list of distinct function names
        # encountered in the English file
        SEO_reader(file_prefix, num_bits, write_log=True)

        def my_fun1(x):
            return x*.5

        def my_fun2(x, y):
            return x + y

        # partial substitution, this creates new files
        # with #1=30, #2=60, 'my_fun1'->my_fun1,
        # but #3  and 'my_fun2' still undecided
        vman = PlaceholderManager(eval_all_vars=False,
                    var_num_to_rads={1: np.pi/6, 2: np.pi/3},
                    fun_name_to_fun={'my_fun1': my_fun1})
        wr = SEO_writer(file_prefix + '_eval01', emb)
        EchoingSEO_reader(file_prefix, num_bits, wr,
                          vars_manager=vman)

        # this runs the simulator after substituting
        # #1=30, #2=60, #3=90, 'my_fun1'->my_fun1, 'my_fun2'->my_fun2
        vman = PlaceholderManager(
            var_num_to_rads={1: np.pi/6, 2: np.pi/3, 3: np.pi/2},
            fun_name_to_fun={'my_fun1': my_fun1, 'my_fun2': my_fun2}
        )
        sim = SEO_simulator(file_prefix, num_bits, verbose=True,
                            vars_manager=vman)
        print("\n----------------------------------------")
        StateVec.describe_st_vec_dict(sim.cur_st_vec_dict)

    main()
