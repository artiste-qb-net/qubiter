from quantum_CSD_compiler.UnitaryMat import *
import Utilities as ut


class Node:
    """
    This class carries the "cargo" of each node of a tree. Included in that
    cargo are

    * self's parent and 2 children nodes,
    * its side, level and id,
    * its left, central and right matrix lists.

    This class also performs the very important task of calling within its
    constructor the function UnitaryMat.cs_decomp() which fills the
    node's left, central and right matrix lists.

    Attributes
    ----------
    central_mats : list(np.ndarray)
        Central matrix list returned by call to UnitaryMat.cs_decomp( ). A
        central_mats is a list of dmats. A dmat= D matrix  is numpy array
        containing floats (radian angles).
    left_mats : list(np.ndarray)
        Left matrix list returned by call to UnitaryMat.cs_decomp()
    left_nd : Node
        Node to left of self.
    level : int
        1<= level <= num_bits+1. level = 1 for root node, level =
        num_of_bits+1 for node whose central_mat is list of 1 dim arrays
    nd_id : int
        node id, int assigned by Tree, nd_id=0 for first (root) node created
        by Tree, nd_id=1 for second node created, etc.
    pa_nd : Node
        parent node
    right_mats : list(np.ndarray)
        Right matrix list returned by call to UnitaryMat.cs_decomp()
    right_nd : Node
        Node to right of self.
    side : str
        to which side of its parent does self find itself, either 'right' or
        'left'


    """

    def __init__(self, nd_id,  pa_nd, side, init_unitary_mat=None):
        """
        Constructor

        Parameters
        ----------
        nd_id : int
        pa_nd : Node
        side : str
        init_unitary_mat : np.ndarray
            This is the matrix that is fed to UnitaryMat.cs_decomp() in root
            node constructor. pa_nd and side are ignored if this is not None.
        Returns
        -------
        None

        """

        self.nd_id = nd_id
        self.pa_nd = pa_nd  # pa=parent, nd=node
        self.side = side  # either 'left' or 'right'
        self.level = None
        # "is None" does not work for numpy array,
        if ut.is_arr(init_unitary_mat):
            pa_nd = None
            side = None

        self.left_nd = None
        self.right_nd = None

        # mats = matrices
        # left_mats, central_mats and right_mats are all list(nd.array)
        self.left_mats = None
        self.central_mats = None
        self.right_mats = None

        if ut.is_arr(init_unitary_mat):
            self.level = 1
            [self.left_mats, self.central_mats, self.right_mats] =\
                UnitaryMat.cs_decomp([init_unitary_mat])
            # release memory
            init_unitary_mat = None
        else:
            self.level = pa_nd.level + 1
            in_mats = None
            if side == 'left':
                in_mats = pa_nd.left_mats
                pa_nd.left_nd = self
            elif side == 'right':
                in_mats = pa_nd.right_mats
                pa_nd.right_nd = self
            else:
                assert False
            [self.left_mats, self.central_mats, self.right_mats] =\
                UnitaryMat.cs_decomp(in_mats)
            if ut.is_arr(in_mats):
                # release memory
                in_mats = None

    def is_barren(self):
        """
        Returns True iff node's left, central and right matrix lists are all
        None.

        Returns
        -------
        bool

        """
        return self.left_mats is None and\
            self.central_mats is None and\
            self.right_mats is None

    def make_barren(self):
        """
        Sets node's left, central and right matrix lists to None.

        Returns
        -------
        None

        """
        self.left_mats = None
        self.central_mats = None
        self.right_mats = None

    def __str__(self):
        """
        Gives a readable description of self when self is ordered to print.
        For example, if nd_id = 3, level = 4, node prints as '3(L4)'.

        Returns
        -------
        str

        """
        return str(self.nd_id) + '(L' + str(self.level) + ')'

if __name__ == "__main__":
    print(5)
