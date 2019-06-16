import itertools as it

import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso

import qubiter.device_specific.utilities_ds as uds
from qubiter.CktEmbedder import *
from qubiter.EchoingSEO_reader import *


class ChipCouplingsFitter:
    """
    Recall that an undirected graph G is defined as a pair of sets (V, E),
    where V is a set of vertices, and E, the edges, is a set of tuples (c,
    t), where c and t are in V. For undirected graphs, the order of c and t
    is ignored.

    This class reads an English file in order to build a list of the CNots
    used in the file. The CNots in the list are specified as (c, t) tuples.
    Then the class builds an undirected graph, call it GE=(V_GE, E_GE) for
    Graph_English, from that CNots list. Then the class builds another
    graph, call it GP=(V_GP, E_GP) for Graph_Physical, from the c_to_tars
    describing the couplings of a particular chip. Then the class uses the
    python networkx function for finding "graph isomorphisms" between
    undirected graphs. The class attempts to find a map phi() from the
    vertices of GE to the vertices of GP, so that phi(x_GE) is a subset of
    x_GP with x=V, E. If the search for phi() succeeds, then phi induces a
    permutation map, call it bit_map: range(num_bits)->range(num_bits),
    bit_ge->bit_gp, of the qubits of the circuit from which GE was
    assembled. Under this permutation, the CNots of the permuted circuit are
    all allowed** by the chip constraint c_to_tars. In that sense, the map
    bit_map() is a "fit" to the chip couplings.

    **except, some CNots may have to be reversed (control and target
    swapped) using Hadamards. These reversals can be performed with the
    class ForbiddenCNotExpander.

    Attributes
    ----------
    bit_map: list[int]
        bit_map[bit_ge] = bit_gp defines a map from the bits bit_ge
        of the graph GE to the bits bit_gp of the graph GP.

    """
    def __init__(self, file_prefix, num_bits, c_to_tars, verbose=False):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
            file prefix of English file which is to be read to assemble
            a list of CNots used in the file.
        num_bits : int
            Number of qubits used in English file with file prefix
            `file_prefix`. IMP: We assume that c_to_tars refers to a chip with
            num_bits too. Both the English file and the chip must have the
            same number of qubits. This is no loss of generality. As long as
            the English file doesn't mention qubit positions >= num_bits,
            all you have to do to conform is to change the name of the
            English file so that it claims to pertain to num_bits qubits.
        c_to_tars : dict[int, list[int]]
            a dictionary mapping j in range(num_bits) to a list, possibly
            empty, of the physically allowed targets of qubit j, when j is
            the control of a CNOT.

        verbose : bool

        Returns
        -------


        """

        old_cnots = ChipCouplingsFitter.get_cnots_in_file(
            file_prefix, num_bits, verbose)
        self.bit_map = ChipCouplingsFitter.get_bit_map_from_c_to_tars(
            num_bits, old_cnots, c_to_tars, verbose)
        emb = CktEmbedder(num_bits, num_bits, self.bit_map)
        out_file_prefix = SEO_reader.xed_file_prefix(file_prefix)
        wr = SEO_writer(out_file_prefix, emb)
        EchoingSEO_reader(file_prefix, num_bits, wr)

    @staticmethod
    def get_cnots_in_file(file_prefix, num_bits, verbose=False):
        """
        This function reads an English file with file prefix `file_prefix`
        pertaining to a circuit with num_bit many qubits. It returns a tuple
        of the CNots mentioned in the English file. The CNots are specified
        as tuples (c, t), where c is the control bit position and t is the
        target bit position.

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        verbose : bool

        Returns
        -------
        tuple[tuple[int, int]]

        """
        old_cnots = []
        english_in = open(utg.preface(
            file_prefix + '_' + str(num_bits) + '_eng.txt'), 'rt')

        while not english_in.closed:
            line = english_in.readline()
            if not line:
                english_in.close()
                break
            split_line = line.split()
            line_name = split_line[0]
            if line_name == "SIGX":
                # example:
                # SIGX AT 1 IF 3F 2T
                tar_bit_pos = int(split_line[2])
                if len(split_line) > 3:
                    if verbose:
                        print(line[:-1])
                    assert len(split_line) == 5, \
                        "Only SIGX with <= 1 controls are allowed."
                    trol_bit_pos = int(split_line[4][:-1])
                    old_cnots.append((trol_bit_pos, tar_bit_pos))
        if verbose:
            print("old_cnots=", old_cnots)
        return tuple(old_cnots)

    @staticmethod
    def draw_phys_and_eng_graphs(file_prefix, num_bits, c_to_tars):
        """
        Draws the Physical and English undirected graphs. This is useful in
        case you want to try to use human pattern recognition to embed the
        English graph inside the Physical graph.

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        c_to_tars : dict[int, list[int]]

        Returns
        -------
        None

        """

        plt.figure(1)
        GP = nx.Graph()
        dir_edges = uds.get_dir_edges_from_c_to_tars(c_to_tars)
        GP.add_edges_from(dir_edges)

        plt.title('Physical graph')
        nx.draw(GP, with_labels=True, node_color='white')

        plt.figure(2)
        old_cnots = ChipCouplingsFitter.get_cnots_in_file(file_prefix,
                                                          num_bits)
        GE = nx.Graph()
        GE.add_edges_from(old_cnots)
        plt.title('English graph')
        nx.draw(GE, with_labels=True, node_color='white')

        plt.show()

    @staticmethod
    def get_bit_map_from_c_to_tars(
            num_bits, old_cnots, c_to_tars, verbose=False):
        """
        This function has as inputs `old_cnot` describing the CNots in an
        English file for a circuit with num_bits qubits and `c_to_tars`
        describing the couplings of a chip with num_bits qubits. The
        function returns a list `bit_map` of num_bit many ints that
        describes a permutation map mapping range(num_bits)->range(
        num_bits), bit_ge->bit_gp. Under this permutation, the CNots of the
        circuit described by the input English file are mapped into new
        CNots which are all allowed** by the chip couplings constraint
        c_to_tars.

        ** except their targets and controls may have to be reversed.

        Parameters
        ----------
        num_bits : int
        old_cnots : tuple[tuple[int, int]]
        c_to_tars : dict[int, list[int]]
        verbose : bool

        Returns
        -------
        list[int]

        """
        GP = nx.Graph()
        dir_edges = uds.get_dir_edges_from_c_to_tars(c_to_tars)
        GP.add_edges_from(dir_edges)
        if verbose:
            print("GP=", GP.edges())

        GE = nx.Graph()
        GE.add_edges_from(old_cnots)
        if verbose:
            print("GE=", GE.edges())

        ge_num_edges = len(GE.edges())
        g_match = None
        is_iso = False
        for edge_sublist in it.combinations(GP.edges(), ge_num_edges):
            subgraph = nx.Graph()
            subgraph.add_edges_from(edge_sublist)
            if verbose:
                print("subgraph=", subgraph.edges())
            g_match = iso.GraphMatcher(subgraph, GE)
            if g_match.is_isomorphic():
                is_iso = True
                break
        assert is_iso, "Could not find node mapping."
        gp_bit_to_ge_bit = g_match.mapping
        if verbose:
            print('gp_bit_to_ge_bit=', gp_bit_to_ge_bit)
        bit_map = [-1]*num_bits
        gp_bits, ge_bits = zip(*gp_bit_to_ge_bit.items())
        for k, ge_bit in enumerate(ge_bits):
            bit_map[ge_bit] = gp_bits[k]
        undefined_domain = [k for k in range(num_bits) if bit_map[k] == -1]
        complement_of_range = [k for k in range(num_bits) if k not in bit_map]
        for y, x in enumerate(undefined_domain):
            bit_map[x] = complement_of_range[y]
        if verbose:
            print("bit_map=", bit_map)
        return bit_map
            
if __name__ == "__main__":
    def main():
        import qubiter.device_specific.chip_couplings_ibm as ibm
        c_to_tars = ibm.ibmq5YorktownTenerife_c_to_tars
        num_bits = 5
        file_prefix = "couplings_fitter"

        print("control_to_targets=", c_to_tars)

        verbose = True
        fitter = ChipCouplingsFitter(
            file_prefix, num_bits, c_to_tars, verbose=verbose)

        print("Must close or save/close all matplotlib "
              "windows in order to finish execution of script.")
        ChipCouplingsFitter.draw_phys_and_eng_graphs(file_prefix,
                                                     num_bits,
                                                     c_to_tars)
    main()

