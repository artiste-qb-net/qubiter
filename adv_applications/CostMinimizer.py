class CostMinimizer:
    """
    This is an abstract class embodying the ideal of a class for minimizing
    a real valued cost function `cost_fun` with respect to its N variables
    x_val.

    Attributes
    ----------
    cur_cost : float
        current cost
    cur_x_val : tuple[float]
        current x value
    init_x_val : tuple[float]
        initial x value
    iter_count : int
        iteration count. An iteration is every time the cost function is
        called
    mfun_kwargs : dict
        mimimizer function's key-word arguments
    minimizer_fun : function
        the minimizer function. For example, scipy.optimizer.minimizer
    print_hiatus : int
        the current values of x_val and cost are printed when iter_count %
        print_hiatus == 0  iff print_hiatus > 0
    verbose : bool

    """

    def __init__(self, minimizer_fun, init_x_val, 
                 print_hiatus=0, verbose=False, **mfun_kwargs):
        """
        Constructor

        Parameters
        ----------
        minimizer_fun : function
        init_x_val : tuple[float]
        print_hiatus : int
        verbose : bool
        mfun_kwargs : dict

        Returns
        -------

        """
        self.minimizer_fun = minimizer_fun
        self.init_x_val = init_x_val
        self.print_hiatus = print_hiatus
        if not print_hiatus:
            self.print_hiatus = 0
        self.verbose = verbose
        self.mfun_kwargs = mfun_kwargs
        
        self.cur_x_val = init_x_val
        self.cur_cost = None
        self.iter_count = 0

    def broadcast_cost_fun_call(self):
        """

        Returns
        -------
        None

        """
        if self.print_hiatus < 1:
            return
        if self.iter_count % self.print_hiatus == 0:
            print('iter=', self.iter_count,
                  ', cost=', self.cur_cost,
                  ', x_val=', self.cur_x_val)

    def cost_fun(self, x_val):
        """
        Abstract method. Given x_val, return float for cost.

        Parameters
        ----------
        x_val : tuple[float]

        Returns
        -------
        float

        """
        assert False

    def find_min(self):
        """
        Abstract method. Returns the minimum (float) of the cost function

        Returns
        -------

        """
        assert False

if __name__ == "__main__":
    def main():
        print(5)
    main()

