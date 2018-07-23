import networkx as nx

import device_specific.utilities as dut
from Controls import *
from EchoingSEO_reader import *


class ForbiddenCNotExpander(EchoingSEO_reader):
    """
    Most chips are not fully connected (not all pairs of qubits are
    physically connected). Furthermore, even if two qubits are connected,
    one of them may be disallowed, forbidden, as a target of a CNOT between
    the 2 qubits. This class is designed to circumvent this chip limitation.

    This class is a child of the class EchoingSEO_reader. It is one of
    several expander classes that replace certain single gates by expansions
    (sequences) of other gates.

    The class reads an English file and outputs a new English file and
    corresponding Picture file. The new English file echoes every line of
    the original English file except for those lines which are SIGX with one
    or more controls.

    If this class reads a line which is a SIGX with > 1 controls, it outputs
    an error message. In such a case, you should use the class CGateExpander
    first to expand such gates into single qubit rotations and simple CNOTs
    with a single control.

    If this class reads a line which is a SIGX with a single control,
    it echoes it if such a CNOT is allowed according to the input list
    'c_to_tars'. Otherwise, it replaces the disallowed (a.k.a., forbidden,
    unphysical, between disconnected qubits) CNOT by a sequence of Hadamards
    and allowed, elementary CNOTs.

    Next we explain the expansion used by this class to replace forbidden
    CNOTs.

    Let us denote a CNot with control a and target b by C(a->b)=C(a, b)

    Note that if C(a, b) is forbidden but C(b, a) is allowed, we can express
    the forbidden one in terms of the allowed one and four Hadamard matrices
    using the identity (X is the target SIGX and @ is the True control)

    X---@
    equals
    H   H
    @---X
    H   H

    Note that

    X---+---@
    equals
    X---@   |
    |   X---@
    X---@   |
    |   X---@
    equals
    |   X---@
    X---@   |
    |   X---@
    X---@   |

    One can generalize the previous identity as follows:

    X---+---+---@
    equals
    X---+---@   |
    |   |   X---@
    X---+---@   |
    |   |   X---@
    equals
    |   X---@   |
    X---@   |   |
    |   X---@   |
    X---@   |   |
    |   |   X---@
    X---@   |   |
    |   X---@   |
    X---@   |   |
    |   X---@   |
    |   |   X---@
    equals (cancel two internal CNots)
    |   X---@   |
    X---@   |   |
    |   X---@   |
    |   |   X---@
    |   X---@   |
    X---@   |   |
    |   X---@   |
    |   |   X---@

    One can generalize the previous identity as follows:

    X---+---+---+---@
    equals
    |   |   X---@   |
    |   X---@   |   |
    X---@   |   |   |
    |   X---@   |   |
    |   |   X---@   |
    |   |   |   X---@
    |   |   X---@   |
    |   X---@   |   |
    X---@   |   |   |
    |   X---@   |   |
    |   |   X---@   |
    |   |   |   X---@

    In general, let's define a composite V gate (called V because it looks
    like a V lying on its side) as follows:

    V(0->4) =
    |   |   |   X---@
    |   |   X---@   |
    |   X---@   |   |
    X---@   |   |   |
    |   X---@   |   |
    |   |   X---@   |
    |   |   |   X---@

    Above, 0, 1, 2, 3, 4 can be replaced by any other distinct qubits. Also,
    on can define an analogous V for any number >= 2 of qubits.

    If
    C(0->4)=
    X---+---+---+---@
    then we proved above that

    C(0->4)= V(0->4)V(1->4)

    In fact, we also proved

    C(0->j)= V(0->j)V(1->j) for j = 2, 3, 4, ...

    We like to refer to the last equation as the vv expansion of C(0->j). In
    this class, we expand a forbidden CNot C(trol->targ) using the last
    equation with the qubit positions 0, 1, 2, ..., j in the last equation
    mapped in a 1-1 onto fashion to qubit positions along a path of qubits
    connecting the two qubits trol and targ. The path is found by calling
    the python networkx function that yields the shortest path between two
    nodes of an undirected graph G. We let G be the undirected graph that
    has as edges all pairs of qubits that are coupled according to the input
    `c_to_d`.

    Attributes
    ----------
    c_to_tars : dict[int, list[int]]
        a dictionary mapping j in range(num_bits) to a list, possibly empty,
        of the physically allowed targets of qubit j, when j is the control
        of a CNOT.
    graph : networkx.Graph
        A networkx undirected graph derived from `c_to_tars` by taking all
        items in ForbiddenCNotExpander.get_dir_edges_from_c_to_tars(
        c_to_tars) as edges.

    """

    def __init__(self, file_prefix, num_bits, c_to_tars):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        c_to_tars : dict[int, list[int]]

        Returns
        -------
        None

        """

        self.c_to_tars = c_to_tars

        self.graph = nx.Graph()
        dir_edges = dut.get_dir_edges_from_c_to_tars(c_to_tars)
        self.graph.add_edges_from(dir_edges)
        # print("graph", self.graph.edges())

        out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
        emb = CktEmbedder(num_bits, num_bits)
        wr = SEO_writer(out_file_prefix, emb)

        EchoingSEO_reader.__init__(self, file_prefix, num_bits, wr)

        self.wr.close_files()

    def edge_type(self, x, y):
        """
        Returns 0 if C(x->y) and C(y->x) are both allowed.
        Returns 1 if C(x->y) but not C(y->x) are allowed
        Returns -1 if C(y->x) but not C(x->y) are allowed
        Returns error message if neither is allowed

        Parameters
        ----------
        x : int
        y : int

        Returns
        -------
        int

        """
        pos = y in self.c_to_tars[x]
        neg = x in self.c_to_tars[y]
        if pos and neg:
            ty = 0
        elif pos:
            ty = 1
        elif neg:
            ty = -1
        else:
            assert False
        return ty

    def get_symbolic_vv_expansion(self, trol_pos, targ_pos):
        """
        This function returns a list called `expansion`. The items in
        `expansion` can be either a tuple (int, int) or a tuple (int,
        bool). A tuple (c, t): (int, int) signifies a CNot with control c
        and target t.  A tuple (k, keep): (int, bool) signifies a Hadamard
        matrix at qubit k and whether to keep it or not--not if it cancels
        another Hadamard. This symbolic expansion is supposed to replace a
        CNot trol_pos->targ_pos.

        Parameters
        ----------
        trol_pos : int
        targ_pos : int

        Returns
        -------
        list[tuple(int, int)|tuple(int, bool)]

        """
        path = nx.shortest_path(self.graph, trol_pos, targ_pos)

        expansion = []

        def do_v(start_bit, end_bit):
            for j in range(start_bit, end_bit):
                ty = self.edge_type(path[j], path[j+1])
                # print(path[j], '->', path[j+1], 'ty=', ty)
                if ty in [0, 1]:
                    # CNot(path[j]->path[j+1]) is allowed
                    expansion.append((path[j], path[j+1]))
                else:
                    # CNot(path[j]->path[j+1]) is not allowed so reverse it
                    expansion.extend([
                        (path[j], True),
                        (path[j+1], True),
                        (path[j+1], path[j]),
                        (path[j+1], True),
                        (path[j], True)])
            for j in reversed(range(start_bit, end_bit-1)):
                ty = self.edge_type(path[j], path[j+1])
                # print(path[j], '->', path[j+1], 'ty=', ty)
                if ty in [0, 1]:
                    # CNot(path[j]->path[j+1]) is allowed
                    expansion.append((path[j], path[j+1]))
                else:
                    # CNot(path[j]->path[j+1]) is not allowed so reverse it
                    expansion.extend([
                        (path[j], True),
                        (path[j+1], True),
                        (path[j+1], path[j]),
                        (path[j+1], True),
                        (path[j], True)])

        do_v(1, len(path)-1)  # small v
        do_v(0, len(path)-1)  # large v
        # now cancel out (set keep to false)
        # Hadamard pairs that have "air" in between
        for j, ej in enumerate(expansion):
            if isinstance(ej[1], bool) and ej[1]:
                # print("-----------j,ej", j, ej)
                for k in range(j+1, len(expansion)):
                    ek = expansion[k]
                    # print("k, ek", k, ek)
                    # a bool is an int as far as isinstance() is concerned
                    # but an int is not a bool
                    if isinstance(ek[1], bool):
                        if ek[1] and (ek[0] == ej[0]):
                            # print("b2", k)
                            expansion[j] = (ej[0], False)
                            expansion[k] = (ek[0], False)
                            break
                    else:  # if ek[1] is int
                        if ej[0] in ek:
                            # print("b1", k)
                            break
        return expansion

    def use_SIG(self, axis, tar_bit_pos, controls):
        """
        This function overrides echoing function in parent class
        EchoingSEO_reader.

        Parameters
        ----------
        axis : int
        tar_bit_pos : int
        controls : Controls

        Returns
        -------
        None

        """
        assert len(controls.bit_pos) <= 1, "Found gate with >= 2 controls"
        if len(controls.bit_pos) == 0:
            EchoingSEO_reader.use_SIG(self, axis, tar_bit_pos, controls)
            return
        # number of controls is 1 from here on
        assert axis == 1, "Found sigy, or sigz with one control"
        assert controls.kinds[0], "Found sigx with False control"
        trol_bit_pos = controls.bit_pos[0]
        sym_exp = self.get_symbolic_vv_expansion(trol_bit_pos, tar_bit_pos)
        # print("expansion", trol_bit_pos, tar_bit_pos, sym_exp)
        for x in sym_exp:
            if isinstance(x[1], bool):
                if x[1]:
                    self.wr.write_one_bit_gate(x[0], OneBitGates.had2)
            else:  # x[1] is int
                self.wr.write_cnot(x[0], x[1])

if __name__ == "__main__":
    def main():
        import device_specific.chip_couplings_ibm as ibm
        file_prefix = "../io_folder/forbidden_cnots_ibm"
        print(file_prefix)
        num_bits = 5
        c_to_tars = ibm.ibmqx2_c_to_tars
        ForbiddenCNotExpander(file_prefix, num_bits, c_to_tars)

        file_prefix = "../io_folder/forbidden_cnots1"
        print(file_prefix)
        num_bits = 4
        c_to_tars = {0: [1], 1: [2], 2: [3], 3: []}
        ForbiddenCNotExpander(file_prefix, num_bits, c_to_tars)
    main()



