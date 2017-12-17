import Utilities as ut


class MolEvolOpData:
    """
    Molecular Evolution Operator Data. This class should be subclassed for
    each particular molecule. For example, H2EvolOpData is a child of this
    class for the H2 molecule. The child class should override the function
    finish_init() to enter h_values of the molecule's Hamiltonian H, where

    Hamiltonian H =
    \sum_{k1, k2} h_value(k1, k2)* a^\dag(k1) a(k2)
    +
    \sum_{k1, k2, k3, k4} h_value(k1, k2, k3, k4)*
        a^\dag(k1) a^\dag(k2) a(k3) a(k4)

    where a^\dag(k) and a(k) are the creation and annihilation operators for
    qubit k

    Attributes
    ----------
    _2bits_to_theta : dict[list(int): float]
        for non-diagonal part of evolution operator, dictionary of list of 2
        ints to theta angle in radians
    _3bits_to_theta : dict[list(int): float]
        for non-diagonal part of evolution operator, dictionary of list of 3
        ints to theta angle in radians
    _4bits_to_3thetas : dict[list(int): list(float)]
        for non-diagonal part of evolution operator, dictionary of list of 4
        ints to list of 3 theta angles in radians
    approx : int
        1 for Trotter approx,
        2 for second order Trotter-Suzuki approx
    diag_2bits_to_theta : dict[list(int): float]
        for diagonal part of evolution operator, dictionary of list of 2
        ints to theta angle in radians
    diag_bit_to_theta : dict[int: float]
        for diagonal part of evolution operator, dictionary of int to theta
        angle in radians
    do_notas : bool
        setting this to False turns off many NOTAS at once
    global_theta : float
        global phase factor in radians
    num_orbitals : int
        number of qubits = number of orbitals + 1. One ancilla to implement
        JW tails
    num_trot_cycles : int
        number of Trotter cycles
    test : bool
        True if want to load testing data
    total_time : float
        total time of evolution, t in exp( - i t H)

    """
    
    def __init__(self, test=False):
        """
        Constructor

        Returns
        -------

        """
        self.do_notas = False
        self.approx = 2
        self.num_orbitals = 0
        self.num_trot_cycles = 0
        self.total_time = 0
        self.test = test

        # the following parameters are filled by using set_h
        self.global_theta = 0
        self.diag_bit_to_theta = {}
        self.diag_2bits_to_theta = {}
        self._2bits_to_theta = {}
        self._3bits_to_theta = {}
        self._4bits_to_3thetas = {}

        self.finish_init()

    @property
    def its_2bits_to_theta(self):
        return self._2bits_to_theta

    @property
    def its_3bits_to_theta(self):
        return self._3bits_to_theta

    @property
    def its_4bits_to_3thetas(self):
        return self._4bits_to_3thetas

    def set_hdiag(self, bits, h_value):
        """
        for diagonal part of Hamiltonian H, use this function inside the
        function set_arg_list() to input a coefficient (h_value) of H.

        diagonal parts of H:
            sum_{k1} h_value(k1)* n(k1)
            + sum_{k1 < k2} h_value(k1, k2)* n(k1)n(k2)
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
        2-qubits
            sum_{k1<k2}  h_value(k1, k2)*[a^\dag(k1) a(k2) + h.c.]
        3-qubits
        Let k = (k1, k2, k3) next
            sum_{k1<k2} \sum{ k3 != k1 or k2}
            h_value(k)* [a^\dag(k1) a(k2) n(k3) + h.c.]
        4-qubits
            Let k = (k1, k2, k3, k4) next

            index = 0
            sum_{k1<k2<k3<k4}
            h_value(k)*[a^\dag(k1) a^\dag(k2) a(k3) a(k4) + h.c.]

            index = 1
            sum_{k1<k2<k3<k4}
            h_value(k)*[a^\dag(k1) a(k2) a^\dag(k3) a(k4) + h.c.]

            index = 2
            sum_{k1<k2<k3<k4}
            h_value(k)*[a^\dag(k1) a(k2) a(k3) a^\dag(k4) + h.c.]

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

    def finish_init(self):
        """
        By default, self.test is set to False in order to force you to
        override this function in a child class for each molecule. If you
        set self.test to True, you get an example of the proper way of
        overriding this function. A proper overrider function should use
        set_hdiag( ) and set_h() to enter numerical values of coefficients
        of H

        Hamiltonian H =
        \sum_{k1, k2} h_value(k1, k2)* a^\dag(k1) a(k2)
        +
        \sum_{k1, k2, k3, k4} h_value(k1, k2, k3, k4)*
            a^\dag(k1) a^\dag(k2) a(k3) a(k4)

        where a^\dag(k) and a(k) are the creation and annihilation operators
        for qubit k

        Returns
        -------
        None

        """
        if not self.test:
            assert False

        self.do_notas = True
        self.approx = 2

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

    def print_params(self):
        print('do_notas', self.do_notas)
        print('approx=', self.approx)
        print('num_orbitals=', self.num_orbitals)
        print('num_trot_cycles=', self.num_trot_cycles)
        print('total_time=', self.total_time)
        print('test=', self.test)
        print('global_theta=', self.global_theta)
        print('diag_bit_to_theta=', self.diag_bit_to_theta)
        print('diag_2bits_to_theta=', self.diag_2bits_to_theta)
        print('_2bits_to_theta=', self._2bits_to_theta)
        print('_3bits_to_theta=', self._3bits_to_theta)
        print('_4bits_to_3thetas=', self._4bits_to_3thetas)


if __name__ == "__main__":
    data = MolEvolOpData(test=True)
    data.print_params()

