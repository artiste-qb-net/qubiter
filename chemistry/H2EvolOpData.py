from chemistry.MolEvolOpData import *


class H2EvolOpData(MolEvolOpData):
    """
    This is a subclass of MolEvolOpData with data for H2 molecule.

    """

    def __init__(self):
        """
        Constructor

        Returns
        -------

        """
        MolEvolOpData.__init__(self)

    def set_arg_list(self, do_notas=True):
        """
        This function overrides function with same name in parent class.

        Parameters
        ----------
        do_notas : bool

        Returns
        -------
        None

        """

        self.do_notas = do_notas

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
    data.set_arg_list()
    arg_list = data.get_arg_list()
    for k in range(len(arg_list)):
        print("\nk=", k)
        print(arg_list[k])

        # arg_list = [
        #     self.do_notas,
        #     self.num_orbitals,
        #     self.num_trot_cycles,
        #     self.total_time,
        #     self.global_theta,
        #     self.diag_bit_to_theta,
        #     self.diag_2bits_to_theta,
        #     self._2bits_to_theta,
        #     self._3bits_to_theta,
        #     self._4bits_to_3thetas
        # ]
