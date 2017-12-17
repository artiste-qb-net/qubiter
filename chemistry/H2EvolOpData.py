from chemistry.MolEvolOpData import *


class H2EvolOpData(MolEvolOpData):
    """
    This is a subclass of MolEvolOpData with data for H2 molecule.

    Attributes
    ----------

    """

    def __init__(self):
        """
        Constructor

        Returns
        -------

        """
        MolEvolOpData.__init__(self)

    def finish_init(self):
        """
        This function overrides function with same name in parent class.

        Returns
        -------
        None

        """

        self.do_notas = True
        self.approx = 2

        self.num_orbitals = 4
        self.num_trot_cycles = 2
        self.total_time = .5

        # this data needs to be checked. Don't trust it for now
        self.set_hdiag(0, -1.252477)
        self.set_hdiag(1, -1.252477)
        self.set_hdiag(2, -0.475934)
        self.set_hdiag(3, -0.475934)

        self.set_hdiag((0, 1), 2*.674493)
        self.set_hdiag((0, 2), 2*.663472 - 2*.181287)
        self.set_hdiag((0, 3), 2*.663472)
        self.set_hdiag((1, 2), 2*.663472)
        self.set_hdiag((1, 3), 2*.663472 - 2*.181287)
        self.set_hdiag((2, 3), 2*.697397)

        self.set_h((0, 1, 2, 3), 2*.181287, 0)
        self.set_h((0, 1, 2, 3), 0., 2)

if __name__ == "__main__":
    data = H2EvolOpData()
    data.print_params()
