from SEO_writer import *
from quantum_CSD_compiler.Node import *
import collections as co


class Tree(SEO_writer):
    """
    This class creates a binary tree of nodes whose cargo is contained in
    the attributes of class Node. This class, being a child of class
    SEO_writer, is also capable of writing English & Picture files. After
    creating a binary tree, it proceeds to use that tree to produce a CS
    decomposition of the unitary matrix init_unitary_mat that is fed into its
    constructor. This CS (cosine-sine) decomp consists of a sequence of
    diagonal unitaries (DIAG lines in English file) and multiplexors (MP_Y
    lines in English file) whose product equals init_unitary_mat.

    If you wish to expand DIAG and MP_Y lines into cnots and single qubit
    rotations, use DiagUnitaryExpander and MultiplexorExpander classes.

    The CS decomposition was a famous decomp of Linear Algebra well before
    quantum computing. It was first applied to quantum computing in the 1999
    paper and accompanying C++ program cited below. Much of the code of the
    original C++ Qubiter has been rewritten in Python for the new pythonic
    Qubiter.

    Let init_unitary_mat be N dimensional, with N = 2^n, where n = number of
    qubits. A general N dimensional unitary matrix has N^2 dofs (real
    degrees of freedom). That's because it has N^2 complex entries, so 2*N^2
    real parameters, but those parameters are subject to N real constraints
    and N(N-1)/2 complex constraints, for a total of N^2 real constraints.
    So 2N^2 real parameters minus N^2 real constraints gives N^2 dofs.

    (a) Each DIAG (MP_Y, resp.) line of the CS decomp of init_unitary_mat
    depends on N (N/2, resp.) angles and there are about N DIAG and N MP_Y
    lines. So the DIAG lines alone have enough dofs, N^2 of them, to cover
    all N^2 dofs of init_unitary_mat. So clearly, there is a lot of
    redundancy in the CS decomp used by Qubiter. But, there is hope: the CS
    decomp is not unique, and it might be possible to choose a CS decomp
    that makes zero many of the angles in the DIAG and MP_Y lines. Some of
    those "compiler optimizations" are considered in references below.

    (b) The CS decomp as used here leads to order N^2 = 2^{2n} cnots and
    qubit rotations so it is impractical for large N. But for small N,
    it can be useful. For large N, it might be possible to discover
    approximations to individual MP_Y and DIAG lines. An approximation of
    this type is considered in MultiplexorExpander.

    Clearly, there is much room for future research to improve (a) and (b).

    References
    ----------
    1. R.R. Tucci, A Rudimentary Quantum Compiler(2cnd Ed.)
    https://arxiv.org/abs/quant-ph/9902062

    2. Qubiter 1.11, a C++ program whose first version was released together
    with Ref.1 above. Qubiter 1.11 is included in the
    quantum_CSD_compiler/LEGACY folder of this newer, pythonic version of Qubiter

    3. R.R. Tucci, Quantum Fast Fourier Transform Viewed as a Special Case
    of Recursive Application of Cosine-Sine Decomposition,
    https://arxiv.org/abs/quant-ph/0411097

    Attributes
    ----------
    global_phase_rads : float
        If arr is the initial unitary matrix fed to the constructor,
        then this equals delta, where arr = exp(i*delta) arr1, where arr1 is
        a special unitary matrix (det(arr1) = 1)
    root_nd : Node
        The root or starting node of the tree. The only node without parents.
        Each node remembers its children, so you only need the root_nd to
        access all other nodes.

    """

    def __init__(self, do_write, file_prefix, emb, init_unitary_mat,
                 verbose=False, **kwargs):
        """
        Constructor

        Parameters
        ----------
        do_write : bool
        file_prefix : str
        emb : CktEmbedder
        init_unitary_mat : np.ndarray
            This is the matrix that is fed to cs_decomp() in root node
            constructor.
        verbose : bool
        kwargs : dict()

        Returns
        -------
        None

        """
        self.verbose = verbose
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        assert UnitaryMat.is_unitary(init_unitary_mat)
        self.global_phase_rads = \
            UnitaryMat.global_phase_rads(init_unitary_mat)
        ph_fac = np.exp(1j*self.global_phase_rads)
        self.root_nd = self.build_tree(init_unitary_mat/ph_fac)
        if do_write:
            self.write()
        
    def build_tree(self, init_unitary_mat):
        """
        This function is called by the constructor to build a tree of
        Node's. It returns the root node of the tree.

        Parameters
        ----------
        init_unitary_mat : np.ndarray

        Returns
        -------
        Node

        """

        nd_ctr = 0

        num_bits = self.emb.num_bits_bef
        num_rows = (1 << num_bits)
        assert init_unitary_mat.shape == (num_rows, num_rows)
        root_nd = Node(nd_ctr, None, None,
                            init_unitary_mat=init_unitary_mat)
        if self.verbose:
            print('building tree------------')
            print(root_nd)
        node_q = co.deque([root_nd])

        # level = level of tree splitting = len(node_q)
        # level = 1 for root node
        # level = num_of_bits+1 for node whose
        # central_mat is list of 1 dim arrays
        level = 1

        while level != 0:
            # since level!=0, cur_nd is not None here
            cur_nd = node_q[0]
            if level == num_bits+1 or cur_nd.is_barren():
                node_q.popleft()
                level -= 1
            else:
                if cur_nd.left_nd is None:
                    nd_ctr += 1
                    next_nd = Node(nd_ctr, cur_nd, 'left')
                    if self.verbose:
                        print(cur_nd, '-left->', next_nd)
                    node_q.appendleft(next_nd)
                    level += 1
                elif cur_nd.right_nd is None:
                    nd_ctr += 1
                    next_nd = Node(nd_ctr, cur_nd, 'right')
                    if self.verbose:
                        print(cur_nd, '-right->', next_nd)
                    node_q.appendleft(next_nd)
                    level += 1
                else:
                    node_q.popleft()
                    level -= 1

        return root_nd

    def write(self):
        """
        This function writes English & Picture files. It visits all the
        Node's of the tree from right to left (this way: <--). It calls
        self.write_node() for each node.

        Returns
        -------
        None

        """

        node_q = co.deque()
        nd = self.root_nd
        if self.verbose:
            print("writing tree------------")
            print(nd)
        while True:
            if nd is not None:
                node_q.appendleft(nd)
                if self.verbose:
                    if nd.right_nd is not None:
                        print(nd, '-right->', nd.right_nd)
                    else:
                        print(nd, '-right->', 'None')
                nd = nd.right_nd
            else:
                # Extract first of the node_q and assign it to nd.
                # Exit while() loop if node_q is empty.
                try:
                    nd = node_q.popleft()
                    self.write_node(nd)
                    if self.verbose:
                        if nd.left_nd is not None:
                            print(nd, '-left->', nd.left_nd)
                        else:
                            print(nd, '-left->', 'None')
                    nd = nd.left_nd
                except:
                    break

    def write_node(self, nd):
        """
        This function is called by self.write() for each node of the tree.
        For a node with level <= num_bits, the function writes an MP_Y line,
        whereas if level = num_bits + 1, it writes a DIAG line.

        Parameters
        ----------
        nd : Node

        Returns
        -------
        None

        """
        if self.verbose:
            self.write_NOTA(str(nd) + "next:")
            print('------start writing ', nd)
        if nd.is_barren():
            self.write_NOTA("barren node")
            return
        num_bits = self.emb.num_bits_bef

        assert 1 <= nd.level <= num_bits+1
        # tar_bit_pos = num_bits - 1 for level=1
        # tar_bit_pos = 0 for level=num_bits
        # tar_bit_pos = -1 for level=num_bits+1
        tar_bit_pos = num_bits - nd.level

        trols = Controls(num_bits)
        if tar_bit_pos >= 0:
            trols.bit_pos_to_kind = {c: c for c in range(tar_bit_pos)}
            for c in range(tar_bit_pos, num_bits-1):
                trols.bit_pos_to_kind[c+1] = c
        else:
            trols.bit_pos_to_kind = {c: c for c in range(num_bits)}
        trols.refresh_lists()

        rad_angles = []
        # central_mats is list of numpy arrays
        for dmat in nd.central_mats:
            rad_angles += list(dmat.flatten())

        # permute arr bit indices
        if 0 <= tar_bit_pos <= num_bits-3:
            # turn rad_angles into equivalent bit indexed tensor
            arr = np.array(rad_angles).reshape([2]*(num_bits-1))
            perm = list(range(tar_bit_pos)) + \
                list(range(tar_bit_pos+1, num_bits-1)) + [tar_bit_pos]
            if self.verbose:
                print("permutation", perm)
            arr.transpose(perm)
            # flatten arr and turn it into a list
            rad_angles = list(arr.flatten())
        # print(rad_angles)

        if self.verbose:
            print("target bit", tar_bit_pos)
            print("controls", trols.bit_pos_to_kind)
            print("rad_angles", rad_angles)
        if tar_bit_pos >= 0:
            self.write_controlled_multiplexor_gate(
                    tar_bit_pos, trols, rad_angles)
        else:
            self.write_controlled_diag_unitary_gate(trols, rad_angles)

if __name__ == "__main__":
    from FouSEO_writer import *
    num_bits = 3
    init_unitary_mat = FouSEO_writer.fourier_trans_mat(1 << num_bits)
    emb = CktEmbedder(num_bits, num_bits)
    file_prefix = '../io_folder/csd_test'
    t = Tree(True, file_prefix, emb, init_unitary_mat, verbose=False)
