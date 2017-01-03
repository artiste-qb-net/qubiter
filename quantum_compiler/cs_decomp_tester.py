import numpy as np
import Utilities as ut
import cuncsd

def is_unitary(arr):
    """
    Returns True iff arr is a numpy array that represents a unitary
    matrix.

    Parameters
    ----------
    arr : np.ndarray

    Returns
    -------
    bool

    """
    num_rows = arr.shape[0]
    assert arr.shape == (num_rows, num_rows)
    err = np.linalg.norm(np.dot(arr, arr.conj().T) -
                         np.identity(num_rows))
    return err < 1e-6


def cs_decomp(unitary_mats):
    """
    This function does a CS (cosine-sine) decomposition (by calling the
    LAPACK function cuncsd.f. The old C++ Qubiter called zggsvd.f
    instead) of each unitary matrix in the list of arrays unitary_mats.
    This function is called by the constructor of the class Node and is
    fundamental for decomposing a unitary matrix into multiplexors and
    diagonal unitaries.

    Parameters
    ----------
    unitary_mats : list(np.ndarray)

    Returns
    -------
    list(np.ndarray), list(np.ndarray), list(np.ndarray)

    """
    block_size = unitary_mats[0].shape[0]
    num_mats = len(unitary_mats)
    for mat in unitary_mats:
        assert mat.shape == (block_size, block_size)

        # In[3]: import cuncsd
        # In[4]: print(cuncsd.cuncsd.__doc__)
        # x11,x12,x21,x22,theta,u1,u2,v1t,v2t,work,rwork,iwork,info =
        # cuncsd(p,x11,x12,x21,x22,[jobu1,jobu2,jobv1t,jobv2t,trans,signs,
        # m,q,ldx11,ldx12,ldx21,ldx22,ldu1,ldu2,ldv1t,ldv2t,lwork,lrwork])
        #
        # Wrapper for ``cuncsd``.
        #
        # Parameters
        # ----------
        # p : input int
        # x11 : input rank-2 array('F') with bounds (p,*)
        # x12 : input rank-2 array('F') with bounds (p,*)
        # x21 : input rank-2 array('F') with bounds (p,*)
        # x22 : input rank-2 array('F') with bounds (p,*)
        #
        # Other Parameters
        # ----------------
        # jobu1 : input string(len=1), optional
        #     Default: 'Y'
        # jobu2 : input string(len=1), optional
        #     Default: 'Y'
        # jobv1t : input string(len=1), optional
        #     Default: 'Y'
        # jobv2t : input string(len=1), optional
        #     Default: 'Y'
        # trans : input string(len=1), optional
        #     Default: 'T'
        # signs : input string(len=1), optional
        #     Default: 'O'
        # m : input int, optional
        #     Default: 2*p
        # q : input int, optional
        #     Default: p
        # ldx11 : input int, optional
        #     Default: p
        # ldx12 : input int, optional
        #     Default: p
        # ldx21 : input int, optional
        #     Default: p
        # ldx22 : input int, optional
        #     Default: p
        # ldu1 : input int, optional
        #     Default: p
        # ldu2 : input int, optional
        #     Default: p
        # ldv1t : input int, optional
        #     Default: p
        # ldv2t : input int, optional
        #     Default: p
        # lwork : input int, optional
        #     Default: -1
        # lrwork : input int, optional
        #     Default: -1
        #
        # Returns
        # -------
        # x11 : rank-2 array('F') with bounds (p,*)
        # x12 : rank-2 array('F') with bounds (p,*)
        # x21 : rank-2 array('F') with bounds (p,*)
        # x22 : rank-2 array('F') with bounds (p,*)
        # theta : rank-1 array('f') with bounds (*)
        # u1 : rank-2 array('F') with bounds (ldu1,*)
        # u2 : rank-2 array('F') with bounds (ldu2,*)
        # v1t : rank-2 array('F') with bounds (ldv1t,*)
        # v2t : rank-2 array('F') with bounds (ldv2t,*)
        # work : rank-1 array('F') with bounds (*)
        # rwork : rank-1 array('f') with bounds (*)
        # iwork : rank-1 array('i') with bounds (*)
        # info : int

    left_mats = []
    central_mats = []
    right_mats = []

    for mat in unitary_mats:
        dim = mat.shape[0]
        assert dim % 2 == 0
        hdim = dim >> 1  # half dimension

        x11 = np.copy(mat[0:hdim, 0:hdim], 'C')
        x12 = np.copy(mat[0:hdim, hdim:dim], 'C')
        x21 = np.copy(mat[hdim:dim, 0:hdim], 'C')
        x22 = np.copy(mat[hdim:dim, hdim:dim], 'C')

        print(x11)
        print(x12)
        print(x21)
        print(x22)

        x11, x12, x21, x22, theta, u1, u2, v1t, v2t,\
            work, rwork, iwork, info =\
                cuncsd.cuncsd(p=hdim, x11=x11, x12=x12, x21=x21, x22=x22 )
        left_mats.append(u1).append(u2)
        central_mats.append(theta)
        right_mats.append(v1t).append(v2t)

    return left_mats, central_mats, right_mats


if __name__ == "__main__":
    from FouSEO_writer import *
    num_bits = 2
    num_rows = 1 << num_bits
    mat = FouSEO_writer.fourier_trans_mat(num_rows)
    print(mat)
    assert is_unitary(mat)
    left_mats, central_mats, right_mats = cs_decomp([mat])
    print(left_mats)
    print(central_mats)
    print(right_mats)