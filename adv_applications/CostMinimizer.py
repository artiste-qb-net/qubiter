class CostMinimizer:
    """
    This is an abstract class embodying the ideal of a class for minimizing
    a real valued cost function `cost_fun` with respect to its N variables
    x_val.

    Attributes
    ----------
    cur_cost : float
        current cost
    cur_targ_cost : float
        current target cost
    cur_x_val : np.ndarray
        current x value
    init_x_val : np.ndarray
        initial x value
    iter_count : int
        iteration count. An iteration is every time the cost function is
        called
    print_hiatus : int
        the current values of x_val and cost are printed when iter_count %
        print_hiatus == 0  iff print_hiatus > 0
    verbose : bool

    """

    def __init__(self, init_x_val, print_hiatus=1, verbose=False):
        """
        Constructor

        Parameters
        ----------
        init_x_val : np.ndarray
        print_hiatus : int
        verbose : bool

        Returns
        -------

        """
        self.init_x_val = init_x_val
        self.print_hiatus = print_hiatus
        self.verbose = verbose
        
        self.cur_x_val = init_x_val
        self.cur_cost = None
        self.cur_targ_cost = None
        self.iter_count = 0

    def broadcast_cost_fun_call(self):
        """
        This method prints current cost and x_val each time cost function is
        called. It also prints targ_cost if there is one.

        Returns
        -------
        None

        """
        if self.print_hiatus < 1:
            return
        x_val_str = ', '.join(['{0:0.6f}'.format(x) for x in self.cur_x_val])
        if self.iter_count % self.print_hiatus == 0:
            s = 'iter=' + str(self.iter_count)
            s += ', cost=' + '{0:0.6f}'.format(self.cur_cost)
            if self.cur_targ_cost is not None:
                s += ', targ_cost=' + '{0:0.6f}'.format(self.cur_targ_cost)
            s += ', x_val=' + x_val_str
            print(s)

    def cost_fun(self, x_val):
        """
        Abstract method. Given x_val, return cost (float).

        Parameters
        ----------
        x_val : np.ndarray

        Returns
        -------
        float

        """
        assert False

    def find_min(self, interface, **kwargs):
        """
        Abstract method. Returns the minimum (float) of the cost function

        Parameters
        ----------
        interface : str
        kwargs: dict

        Returns
        -------

        """
        assert False

if __name__ == "__main__":
    def main():
        print(5)
    main()

