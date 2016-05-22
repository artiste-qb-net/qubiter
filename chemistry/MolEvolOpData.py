import Utilities as ut


class MolEvolOpData:
    """
    Molecular Evolution Operator Data. This class should be subclassed for
    each particular molecule. For example, H2EvolOpData is a child of this
    class for the H2 molecule. The child class should override the function
    set_arg_list() to enter h_values of the molecule's Hamiltonian H, where

    Hamiltonian H =
    \sum_{k1, k2} h_value(k1, k2)* a^\dag(k1) a(k2)
    +
    \sum_{k1, k2, k3, k4} h_value(k1, k2, k3, k4)*
        a^\dag(k1) a^\dag(k2) a(k3) a(k4)

    where a^\dag(k) and a(k) are the creation and annihilation operators for
    qubit k

    Parameters
    ----------
    global_theta : float
        global phase factor in radians
    diag_bit_to_theta : dict[int: float]
        for diagonal part of evolution operator, dictionary of int to theta
        angle in radians
    diag_2bits_to_theta : dict[list(int): float]
        for diagonal part of evolution operator, dictionary of list of 2
        ints to theta angle in radians
    _2bits_to_theta : dict[list(int): float]
        for non-diagonal part of evolution operator, dictionary of list of 2
        ints to theta angle in radians
    _3bits_to_theta : dict[list(int): float]
        for non-diagonal part of evolution operator, dictionary of list of 3
        ints to theta angle in radians
    _4bits_to_3thetas : dict[list(int): list(float)]
        for non-diagonal part of evolution operator, dictionary of list of 4
        ints to list of 3 theta angles in radians
    num_orbitals : int
        number of qubits = number of orbitals + 1. One ancilla to implement
        JW tails
    num_trot_cycles : int
        number of Trotter cycles
    total_time : float
        total time of evolution, t in exp( - i t H)
    do_notas : bool

    """
    
    def __init__(self):
        """
        Constructor

        Returns
        -------

        """
        self.do_notas = None
        self.num_orbitals = 0
        self.num_trot_cycles = 0
        self.total_time = 0

        # the following parameters are filled by using set_h
        self.global_theta = 0
        self.diag_bit_to_theta = {}
        self.diag_2bits_to_theta = {}
        self._2bits_to_theta = {}
        self._3bits_to_theta = {}
        self._4bits_to_3thetas = {}
        
    def get_arg_list(self, do_notas=True):
        """
        Returns a list of arguments. This list is one of the input
        parameters for the constructor of the class PhaseEstSEO_writer

        Parameters
        ----------
        do_notas : bool

        Returns
        -------
        list[]

        """
        self.do_notas = do_notas
        arg_list = [
            self.do_notas,
            self.num_orbitals,
            self.num_trot_cycles,
            self.total_time,
            self.global_theta,
            self.diag_bit_to_theta,
            self.diag_2bits_to_theta,
            self._2bits_to_theta,
            self._3bits_to_theta,
            self._4bits_to_3thetas
        ]
        return arg_list

    def set_hdiag(self, bits, h_value):
        """
        for diagonal part of Hamiltonian H, use this function inside the
        function set_arg_list() to input a coefficient (h_value) of H.

        diagonal parts of H:
            sum_{k1} h_value(k1)* n(k1)
            - sum_{k1 < k2} h_value(k1, k2)* n(k1)n(k2)
        for qubits k1 and k2, where n(k) = |1><1| at k = P_1(k)

        Parameters
        ----------
        bits : int | list[int]
        h_value : float

        Returns
        -------
        None

        """

        if not bits:
            a = 0
        elif isinstance(bits, int):
            a = 1
        else:
            a = len(bits)

        dt = self.total_time / self.num_trot_cycles
        # this accounts for minus sign in evolution operator exp(-i t H )
        dt = -dt

        if a == 0:
            self.global_theta += h_value

        elif a == 1:
            self.global_theta += h_value * dt / 2

            di = self.diag_bit_to_theta
            inc = -h_value*dt/2
            ut.increment_dict(di, bits, inc)

        elif a == 2:
            assert bits[0] < bits[1]
            self.global_theta += h_value * dt / 4

            di = self.diag_bit_to_theta
            inc = -h_value*dt/4
            ut.increment_dict(di, bits[0], inc)

            di = self.diag_2bits_to_theta
            inc = -h_value*dt/2
            ut.increment_dict(di, bits, inc)

        else:
            assert False

    def set_h(self, bits, h_value, index=None):
        """
        for non-diagonal part of Hamiltonian H, use this function inside the
        function set_arg_list() to input a coefficient (h_value) of H.

        non-diagonal parts of H:
            sum_{k1<k2}  h_value(k1, k2)*a^\dag(k1) a(k2) + h.c.

            sum_{k1<k2, k3 != k1 or k2}
                h_value(k1, k2, k3)* a^\dag(k1) a(k2) n(k3) + h.c.

            index = 0
            sum_{k1<k2<k3<k4}
            h_value(k1, k2, k3, k4)* a^\dag(k1) a^\dag(k2) a(k3) a(k4) + h.c.

            index = 1
            sum_{k1<k2<k3<k4}
            -h_value(k1, k2, k3, k4)* a^\dag(k1) a(k2) a^\dag(k3) a(k4) + h.c.

            index = 2
            sum_{k1<k2<k3<k4}
            h_value(k1, k2, k3, k4)* a^\dag(k1) a(k2) a(k3) a^\dag(k4) + h.c.

        where n(k) = |1><1| at k = P_1(k) and where a^\dag(k) and a(k) are
        the creation and annihilation operators for qubit k


        Parameters
        ----------
        bits : list[int]
        h_value : float
        index : int
            for the case of 4 bit combis, this ranges over 0, 1, 2

        Returns
        -------
        None

        """

        a = len(bits)

        dt = self.total_time / self.num_trot_cycles
        # this accounts for minus sign in evolution operator exp(-i t H )
        dt = -dt

        if a == 2:
            assert bits[0] < bits[1]
            self._2bits_to_theta[bits] = h_value*dt

        elif a == 3:
            assert bits[0] < bits[1]
            self._3bits_to_theta[bits] = h_value*dt

        elif a == 4:
            assert bits[0] < bits[1] < bits[2] < bits[3]
            if bits not in self._4bits_to_3thetas:
                self._4bits_to_3thetas[bits] = [0., 0., 0.]
            self._4bits_to_3thetas[bits][index] = h_value*dt

        else:
            assert False

    def set_arg_list(self, do_notas=True):
        """
        You must override this function in a child class for each molecule.
        See set_test_arg_list() for an example of how to write this function.

        Parameters
        ----------
        do_notas : bool

        Returns
        -------
        None

        """
        assert False

    def set_test_arg_list(self, do_notas=True):
        """
        This function is an example of what set_arg_list() should be like.
        In this function, we use functions set_hdiag( ) and set_h() to enter
        numerical values of coefficients of H

        Hamiltonian H =
        \sum_{k1, k2} h_value(k1, k2)* a^\dag(k1) a(k2)
        +
        \sum_{k1, k2, k3, k4} h_value(k1, k2, k3, k4)*
            a^\dag(k1) a^\dag(k2) a(k3) a(k4)

        where a^\dag(k) and a(k) are the creation and annihilation operators
        for qubit k

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

        self.set_hdiag(0, .5)
        self.set_hdiag(1, .5)
        self.set_hdiag(2, .5)
        self.set_hdiag(3, .5)

        self.set_hdiag((0, 1), .5)
        self.set_hdiag((0, 2), .5)
        self.set_hdiag((0, 3), .5)
        self.set_hdiag((1, 2), .5)
        self.set_hdiag((1, 3), .5)
        self.set_hdiag((2, 3), .5)

        self.set_h((0, 1), .5)
        self.set_h((0, 2), .5)
        self.set_h((0, 3), .5)
        self.set_h((1, 2), .5)
        self.set_h((1, 3), .5)
        self.set_h((2, 3), .5)

        self.set_h((0, 1, 2), .5)
        self.set_h((0, 1, 3), .5)

        self.set_h((0, 1, 2, 3), .5, 0)
        self.set_h((0, 1, 2, 3), .5, 1)
        self.set_h((0, 1, 2, 3), .5, 2)


if __name__ == "__main__":
    data = MolEvolOpData()
    data.set_test_arg_list()
    arg_list = data.get_arg_list()
    for k in range(len(arg_list)):
        print("\nk=", k)
        print(arg_list[k])
