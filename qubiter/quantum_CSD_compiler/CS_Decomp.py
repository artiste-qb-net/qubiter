import math
import cuncsd_sq as csd
from qubiter.UnitaryMat import *


class CS_Decomp:

    @staticmethod
    def get_csd(unitary_mats):
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
            vec1 = vec[0] * np.ones((num_mats,))
            if np.linalg.norm(vec - vec1) < 1e-6:
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
                x11, x12, x21, x22, theta, u1, u2, v1t, v2t, \
                work, rwork, iwork, info = \
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
                x11, x12, x21, x22, theta, u1, u2, v1t, v2t, \
                work, rwork, iwork, info = \
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
    from qubiter.FouSEO_writer import *
    from qubiter.quantum_CSD_compiler.MultiplexorSEO_writer import *

    def main():
        print("\ncs decomp example-------------")
        num_bits = 2
        num_rows = 1 << num_bits
        mat = FouSEO_writer.fourier_trans_mat(num_rows)
        assert UnitaryMat.is_unitary(mat)
        left_mats, central_mats, right_mats = CS_Decomp.get_csd([mat])
        # print('left mats\n', left_mats)
        # print('central_mats\n', central_mats)
        # print('right_mats\n', right_mats)

        left = np.zeros((num_rows, num_rows), dtype=complex)
        half_nrows = num_rows // 2
        left[0:half_nrows, 0:half_nrows] = left_mats[0]
        left[half_nrows:num_rows, half_nrows:num_rows] = left_mats[1]
        # print('left', left)

        right = np.zeros((num_rows, num_rows), dtype=complex)
        right[0:half_nrows, 0:half_nrows] = right_mats[0]
        right[half_nrows:num_rows, half_nrows:num_rows] = right_mats[1]

        center = MultiplexorSEO_writer.mp_mat(central_mats[0])
        # print('center', center)
        mat_prod = np.dot(np.dot(left, center), right)
        # print('mat\n', mat)
        # print('mat_prod\n', mat_prod)
        err = np.linalg.norm(mat - mat_prod)
        print('err=', err)

    main()


