class QbitPlanarLattice:
    """
    This class translates between int and int pair coordinates for qubits on
    a  planar lattice. From a rectangular ascii picture of the planar chip (
    for instance, device_specific.chip_couplings_google._BRISTLECONE_GRID),
    the class finds neighbors of each qubit. Two qubits are neighbors iff
    they are adjacent and have the same column or row. Two qubits are
    neighbors also iff a CNOT is physically allowed with either one of the
    qubits as target and the other as control.

    Attributes
    ----------
    num_bits : int
    qbit_2d_coords : list[tuple[int, int]]
        List of 2d coords for each qubit

    References
    ----------
    1. cirq/google/xmon_device.py
    2. cirq/devices/grid_qubit.py
    3. cirq/google/known_devices.py

    """

    def __init__(self, ascii_pic):
        """
        Constructor

        Parameters
        ----------
        ascii_pic : str
            rectangular ascii picture of planar chip. For instance,
            device_specific.chip_couplings_google._BRISTLECONE_GRID

        Returns
        -------


        """

        self.qbit_2d_coords = []
        lines = ascii_pic.strip().split('\n')
        for row, line in enumerate(lines):
            for col, char in enumerate(line.strip()):
                if char != '-':
                    self.qbit_2d_coords.append((row, col))
        self.num_bits = len(self.qbit_2d_coords)

    def two2one(self, pair):
        """
        Translates from int pair coords to int coord.

        Parameters
        ----------
        pair : tuple[int, int]

        Returns
        -------
        int

        """
        if pair in self.qbit_2d_coords:
            return self.qbit_2d_coords.index(pair)
        else:
            return None

    def one2two(self, index):
        """
        Translates from int coord to int pair coords.

        Parameters
        ----------
        index : int

        Returns
        -------
        tuple[int, int]

        """
        assert 0 <= index < self.num_bits
        return self.qbit_2d_coords[index]

    def is_empty(self, pair):
        """
        True iff pair=(row, col)  position on grid has no qubit.

        Parameters
        ----------
        pair : tuple[int, int]

        Returns
        -------
        bool

        """
        return pair not in self.qbit_2d_coords

    def are_neighbors(self, ind1, ind2):
        """
        Returns True iff qubits with int coords ind1 and ind2 are neighbors.

        Parameters
        ----------
        ind1 : int
        ind2 : int

        Returns
        -------
        bool

        """
        r1, c1 = self.one2two(ind1)
        r2, c2 = self.one2two(ind2)
        return abs(r1 - r2) + abs(c1 - c2) == 1

    def neighbors_of(self, ind):
        """
        Returns list of ints that are int coords of qubits that are
        neighbors of qubit with int coord ind.

        Parameters
        ----------
        ind : int

        Returns
        -------
        list[int]

        """
        r, c = self.one2two(ind)
        nbors = []
        for r1, c1 in [(r, c+1), (r, c-1), (r+1, c), (r-1, c)]:
                ind1 = self.two2one((r1, c1))
                if ind1 is not None:
                    nbors.append(ind1)
        return nbors

    def get_c_to_tars(self):
        """
        Returns a dictionary mapping each qubit j to a list of all
        physically allowed target qubits, assuming j is a control of a CNOT.
        All qubits are specified by their int coords.

        Returns
        -------
        dict[int, list[int]]

        """
        c_to_tars = {}
        for ind in range(self.num_bits):
            c_to_tars[ind] = self.neighbors_of(ind)
        return c_to_tars

if __name__ == "__main__":
    import device_specific.chip_couplings_google as cc
    import pprint as pp

    def main():
        lattice = QbitPlanarLattice(cc.BRISTLECONE_GRID)
        pp.pprint(lattice.qbit_2d_coords)
        print('\nnum_bits=', lattice.num_bits, '\n')
        c_to_tars = lattice.get_c_to_tars()
        pp.pprint(c_to_tars)
    main()







