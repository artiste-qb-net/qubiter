import math
import numpy as np
import Utilities as ut
from quantum_CSD_compiler.MultiplexorSEO_writer import *
import cuncsd_sq as csd

from FouSEO_writer import *


class UnitaryMat:
    """
    This class contains only static methods and no constructor. It contains
    some functions associated with unitary matrices.

    Note that any numpy array called arr can be saved and loaded from a text
    file called file_name using

        np.savetxt(file_name, arr)
    and
        arr = np.loadtxt(file_name)
    """

    @staticmethod
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

    @staticmethod
    def global_phase_rads(arr):
        """
        Assuming that arr is a unitary matrix, this function returns delta
        such that arr = exp(i*delta) arr1, where arr1 is a special unitary
        matrix (det(arr1) = 1).

        Parameters
        ----------
        arr : np.ndarray

        Returns
        -------
        float

        """
        det = np.linalg.det(arr)
        c = det.real
        s = det.imag
        num_rows = arr.shape[0]
        return np.arctan2(s, c)/num_rows

    @staticmethod
    def u2_from_params(delta, rot_rads, unit_vec):
        """
        This function returns a numpy array arr = exp(i*delta)exp(
        i*rot_rads*sig_w) such that unit_vec = [wx, wy, wz] is a unit vector
        and sig_w = wx*sigx + wy*sigy + wz*sigz.

        u2_to_params() maps arr -> (delta, rot_rads, unit_vec). This 
        function maps (delta, rot_rads, unit_vec) -> arr so it is the 
        inverse of u2_to_params(). But be careful, because delta, rot_rads 
        are not single valued, they are only valid mod 2*pi. 

        Parameters
        ----------
        delta : float
        rot_rads : float
        unit_vec : list[float]

        Returns
        -------
        np.ndarray

        """
        x_rads = rot_rads*unit_vec[0]
        y_rads = rot_rads*unit_vec[1]
        z_rads = rot_rads*unit_vec[2]
        arr = np.exp(1j*delta)*OneBitGates.rot(x_rads, y_rads, z_rads)
        return arr

    @staticmethod
    def u2_to_params(arr):
        """
        Assuming that arr is a U(2) matrix, this function returns the
        parameters (delta, rot_rads, unit_vec) such that arr = exp(
        i*delta)exp( i*rot_rads*sig_w) such that unit_vec = [wx, wy, wz] is
        a unit vector and sig_w = wx*sigx + wy*sigy + wz*sigz.

        Parameters
        ----------
        arr : np.ndarray

        Returns
        -------
        float, float, list[float]

        """
        delta = UnitaryMat.global_phase_rads(arr)
        arr1 = arr/np.exp(1j*delta)

        cw = arr1[0, 0].real

        wx = arr1[0, 1].imag
        wy = arr1[0, 1].real
        wz = arr1[0, 0].imag

        sw = np.sqrt(wx**2 + wy**2 + wz**2)
        wx /= sw
        wy /= sw
        wz /= sw

        rot_rads = np.arctan2(sw, cw)
        unit_vec = [wx, wy, wz]

        return delta, rot_rads, unit_vec

    @staticmethod
    def u2_zyz_decomp(arr):
        """
        Assuming that arr is a U(2) matrix, this function returns (delta,
        left_rads, center_rads, right_rads) such that

        arr = exp(i*delta)*
        exp(i*left_rads*sigz) exp(i*center_rads*sigy) exp(i*right_rads*sigz)

        If change axes by x->y, y->z, z->x, get xzx decomp with same angles.

        Parameters
        ----------
        arr : np.ndarray

        Returns
        -------
        float, float, float, float

        """
        delta, rot_rads, unit_vec = UnitaryMat.u2_to_params(arr)
        [wx, wy, wz] = unit_vec
        sw = np.sin(rot_rads)
        cw = np.cos(rot_rads)

        theta1 = np.arctan2(wz*sw, cw)
        theta2 = np.arctan2(wx, wy)
        left_rads = (theta1 + theta2)/2
        right_rads = (theta1 - theta2)/2
        center_rads = np.arctan2(sw*np.sqrt(wx**2 + wy**2),
                           np.sqrt(cw**2 + (wz*sw)**2))

        return delta, left_rads, center_rads, right_rads

    @staticmethod
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
        if block_size == 1:
            left_mats = None
            right_mats = None
            vec = np.array([unitary_mats[k][0, 0]
                            for k in range(0, num_mats)])
            vec1 = vec[0]*np.ones((num_mats,))
            if np.linalg.norm(vec-vec1) < 1e-6:
                central_mats = None
            else:
                c_vec = np.real(vec)
                s_vec = np.imag(vec)
                central_mats = np.arctan2(s_vec, c_vec)
        else:
            # Henning Dekant constructed a python wrapper for LAPACK's cuncsd.f
            # via the python application f2py. Thanks Henning!

            # In[2]: import cuncsd
            # In[3]: print(cuncsd.cuncsd.__doc__)
            # x11,x12,x21,x22,theta,u1,u2,v1t,v2t,work,rwork,iwork,info =
            # cuncsd(p,x11,x12,x21,x22,lwork,lrwork,
            # [jobu1,jobu2,jobv1t,jobv2t,trans,signs,m,q,
            # ldx11,ldx12,ldx21,ldx22,ldu1,ldu2,ldv1t,ldv2t,credit])
            #
            # Wrapper for ``cuncsd``.
            #
            # Parameters
            # ----------
            # p : input int
            # x11 : input rank-2 array('F') with bounds (p,p)
            # x12 : input rank-2 array('F') with bounds (p,p)
            # x21 : input rank-2 array('F') with bounds (p,p)
            # x22 : input rank-2 array('F') with bounds (p,p)
            # lwork : input int
            # lrwork : input int
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
            # credit : input int, optional
            #     Default: 0
            #
            # Returns
            # -------
            # x11 : rank-2 array('F') with bounds (p,p)
            # x12 : rank-2 array('F') with bounds (p,p)
            # x21 : rank-2 array('F') with bounds (p,p)
            # x22 : rank-2 array('F') with bounds (p,p)
            # theta : rank-1 array('f') with bounds (p)
            # u1 : rank-2 array('F') with bounds (p,p)
            # u2 : rank-2 array('F') with bounds (p,p)
            # v1t : rank-2 array('F') with bounds (p,p)
            # v2t : rank-2 array('F') with bounds (p,p)
            # work : rank-1 array('F') with bounds (abs(lwork))
            # rwork : rank-1 array('f') with bounds (abs(lrwork))
            # iwork : rank-1 array('i') with bounds (p)
            # info : int

            left_mats = []
            central_mats = []
            right_mats = []

            for mat in unitary_mats:
                dim = mat.shape[0]
                assert dim % 2 == 0
                hdim = dim >> 1  # half dimension

                p = hdim
                # x11 = np.copy(mat[0:hdim, 0:hdim], 'C')
                # x12 = np.copy(mat[0:hdim, hdim:dim], 'C')
                # x21 = np.copy(mat[hdim:dim, 0:hdim], 'C')
                # x22 = np.copy(mat[hdim:dim, hdim:dim], 'C')

                x11 = mat[0:hdim, 0:hdim]
                x12 = mat[0:hdim, hdim:dim]
                x21 = mat[hdim:dim, 0:hdim]
                x22 = mat[hdim:dim, hdim:dim]

                # print('mat\n', mat)
                # print('x11\n', x11)
                # print('x12\n', x12)
                # print('x21\n', x21)
                # print('x22\n', x22)
                x11, x12, x21, x22, theta, u1, u2, v1t, v2t,\
                    work, rwork, iwork, info =\
                        csd.cuncsd(p, x11, x12, x21, x22, lwork=-1, lrwork=-1,
                                      trans='F')
                # print('x11\n', x11)
                # print('x12\n', x12)
                # print('x21\n', x21)
                # print('x22\n', x22)

                lw = math.ceil(work[0].real)
                lrw = math.ceil(rwork[0].real)
                # print("work query:", lw)
                # print("rwork query:", lrw)
                x11, x12, x21, x22, theta, u1, u2, v1t, v2t,\
                    work, rwork, iwork, info =\
                        csd.cuncsd(p, x11, x12, x21, x22, lwork=lw, lrwork=lrw,
                                      trans='F')
                # print('info', info)
                # print('u1 continguous', u1.flags.contiguous)
                # u1 = np.ascontiguousarray(u1)
                # u2 = np.ascontiguousarray(u2)
                # v1t = np.ascontiguousarray(v1t)
                # v2t = np.ascontiguousarray(v2t)
                # print('u1 continguous', u1.flags.contiguous)

                left_mats.append(u1)
                left_mats.append(u2)
                central_mats.append(theta)
                right_mats.append(v1t)
                right_mats.append(v2t)
        return left_mats, central_mats, right_mats


if __name__ == "__main__":
    from FouSEO_writer import *
    unit_vec = np.array([1, 2, 3])
    unit_vec = unit_vec/np.linalg.norm(unit_vec)
    unit_vec = list(unit_vec)
    delta = 3
    rot_rads = .3
    print("delta in=", delta)
    print("rot_rads in=", rot_rads)
    print("unit_vec in=", unit_vec)
    arr_in = UnitaryMat.u2_from_params(delta, rot_rads, unit_vec)
    print("arr_in:\n", arr_in)
    delta1, rot_rads1, unit_vec1 = UnitaryMat.u2_to_params(arr_in)
    print("delta out=", delta1)
    print("rot_rads out=", rot_rads1)
    print("unit_vec out=", unit_vec1)
    arr_out = UnitaryMat.u2_from_params(delta1, rot_rads1, unit_vec1)
    print("arr_out=\n", arr_out)

    delta, left_rads, center_rads, right_rads =\
        UnitaryMat.u2_zyz_decomp(arr_in)
    c = np.cos(center_rads)
    s = np.sin(center_rads)
    a = np.exp(1j*left_rads)
    ah = np.conj(a)
    b = np.exp(1j*right_rads)
    bh = np.conj(b)
    arr = np.empty((2, 2), dtype=complex)
    arr[0, 0] = c*a*b
    arr[0, 1] = s*a*bh
    arr[1, 0] = -s*ah*b
    arr[1, 1] = c*ah*bh
    arr = np.exp(1j*delta)*arr
    print("zyz decomp arr=\n", arr)

    print("\ncs decomp example-------------")
    num_bits = 2
    num_rows = 1 << num_bits
    mat = FouSEO_writer.fourier_trans_mat(num_rows)
    assert UnitaryMat.is_unitary(mat)
    left_mats, central_mats, right_mats = UnitaryMat.cs_decomp([mat])
    # print('left mats\n', left_mats)
    # print('central_mats\n', central_mats)
    # print('right_mats\n', right_mats)
    
    left = np.zeros((num_rows, num_rows),  dtype=complex)
    left[0:num_rows/2, 0:num_rows/2] = left_mats[0]
    left[num_rows/2:num_rows, num_rows/2:num_rows] = left_mats[1]
    # print('left', left)
    
    right = np.zeros((num_rows, num_rows),  dtype=complex)
    right[0:num_rows/2, 0:num_rows/2] = right_mats[0]
    right[num_rows/2:num_rows, num_rows/2:num_rows] = right_mats[1]
    
    center = MultiplexorSEO_writer.mp_mat(central_mats[0])
    # print('center', center)
    mat_prod = np.dot(np.dot(left, center), right)
    # print('mat\n', mat)
    # print('mat_prod\n', mat_prod)
    err = np.linalg.norm(mat - mat_prod)
    print('err=', err)
