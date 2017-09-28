import numpy as np


class OneBitGates:
    """
    This class has no attributes or constructor. It is simply a collection
    of static methods, all of which return a complex 2 by 2 matrix (numpy
    array). In cases where the entries of the matrix are all real,
    an is_quantum bool option is given to choose between a float64 or
    complex128 array.


    Attributes
    ----------

    """

    @staticmethod
    def sigx(is_quantum=True):
        """
        Returns \sigma_x Pauli matrix.

        Parameters
        ----------
        is_quantum : bool

        Returns
        -------
        np.ndarray

        """
        if not is_quantum:
            ty = np.float64
        else:
            ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 1] = 1
        mat[1, 0] = 1
        return mat

    @staticmethod
    def sigy():
        """
        Returns \sigma_y Pauli matrix.

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 1] = -1j
        mat[1, 0] = 1j
        return mat

    @staticmethod
    def sigz(is_quantum=True):
        """
        Returns \sigma_z Pauli matrix.

        Parameters
        ----------
        is_quantum : bool

        Returns
        -------
        np.ndarray

        """
        if not is_quantum:
            ty = np.float64
        else:
            ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 0] = 1
        mat[1, 1] = -1
        return mat

    @staticmethod
    def had2(is_quantum=True):
        """
        Returns 2 dimensional Hadamard matrix (\sigma_x + \sigma_z)/sqrt(2)

        Parameters
        ----------
        is_quantum : bool

        Returns
        -------
        np.ndarray

        """
        if not is_quantum:
            ty = np.float64
        else:
            ty = np.complex128
        x = 1/np.sqrt(2)
        mat = np.zeros([2, 2], dtype=ty) + x
        mat[1, 1] = - mat[1, 1]
        return mat

    @staticmethod
    def P_0(is_quantum=True):
        """
        Returns projection operator P_0 = |0><0| = nbar, where |0> = [1,
        0]^T and |1> = [0, 1]^T, T = transpose

        Parameters
        ----------
        is_quantum : bool

        Returns
        -------
        np.ndarray

        """
        if not is_quantum:
            ty = np.float64
        else:
            ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 0] = 1
        return mat

    @staticmethod
    def P_1(is_quantum=True):
        """
        Returns projection operator P_1 = |1><1| = nbar, where |0> = [1,
        0]^T and |1> = [0, 1]^T, T = transpose

        Parameters
        ----------
        is_quantum : bool

        Returns
        -------
        np.ndarray

        """
        if not is_quantum:
            ty = np.float64
        else:
            ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[1, 1] = 1
        return mat

    @staticmethod
    def P_0_phase_fac(ang_rads):
        """
        Returns

        exp(1j*ang_rads*P_0) = [[x, 0],[0, 1]] with x = exp(1j*ang_rads)

        Parameters
        ----------
        ang_rads : float

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 0] = np.exp(1j*ang_rads)
        mat[1, 1] = 1
        return mat

    @staticmethod
    def P_1_phase_fac(ang_rads):
        """
        Returns

        exp(1j*ang_rads*P_1) = [[1, 0],[0, x]] with x = exp(1j*ang_rads)

        Parameters
        ----------
        ang_rads : float

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[1, 1] = np.exp(1j*ang_rads)
        mat[0, 0] = 1
        return mat

    @staticmethod
    def phase_fac(ang_rads):
        """
        Returns

        exp(1j*ang_rads*I_2) = [[x, 0],[0, x]] with x = exp(1j*ang_rads)


        Parameters
        ----------
        ang_rads : float

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        x = np.exp(1j*ang_rads)
        mat[1, 1] = x
        mat[0, 0] = x
        return mat

    @staticmethod
    def rot(rad_ang_x, rad_ang_y, rad_ang_z):
        """
        Returns

        exp(1j*(rad_ang_x*sig_x + rad_ang_y*sig_y + rad_ang_z*sig_z)

        where rad_ang_x is an angle in radians and sig_x is the x Pauli
        matrix, etc.

        Parameters
        ----------
        rad_ang_x : float
        rad_ang_y : float
        rad_ang_z : float

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        vec = np.array([rad_ang_x, rad_ang_y, rad_ang_z])
        n = np.linalg.norm(vec)  # sqrt(dot(vec, vec.conj))
        if abs(n) < 1e-8:
            mat[0, 0] = 1
            mat[1, 1] = 1
        else:
            nx = rad_ang_x/n
            ny = rad_ang_y/n
            nz = rad_ang_z/n
            c = np.cos(n)
            s = np.sin(n)
            mat[0, 0] = c + 1j*s*nz
            mat[0, 1] = s*ny + 1j*s*nx
            mat[1, 0] = -s*ny + 1j*s*nx
            mat[1, 1] = c - 1j*s*nz
        return mat

    @staticmethod
    def rot_ax(rad_ang, axis):
        """
        Returns

        exp(1j*rad_ang*sig_n)

        where n = x if axis = 1, n = y if axis = 2 and n = z if axis = 3

        Parameters
        ----------
        rad_ang : float
        axis : int

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        c = np.cos(rad_ang)
        s = np.sin(rad_ang)

        if axis == 1:
            mat[0, 0] = c
            mat[0, 1] = 1j*s
            mat[1, 0] = 1j*s
            mat[1, 1] = c
        elif axis == 2:
            mat[0, 0] = c
            mat[0, 1] = s
            mat[1, 0] = -s
            mat[1, 1] = c
        elif axis == 3:
            mat[0, 0] = c + 1j*s
            mat[1, 1] = c - 1j*s
        else:
            assert False, "axis not in [1,2,3]"

        return mat

    @staticmethod
    def mat_S(herm=False):
        """
        Returns

        [[1, 0],[0, x*sign]] where x=exp(j*pi/2)=j

        where sign = 1 if herm=False and sign = -1 if herm=True

        Parameters
        ----------
        herm : bool

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 0] = 1
        if not herm:
            sign = 1
        else:
            sign = -1
        mat[1, 1] = 1j*sign
        return mat

    @staticmethod
    def mat_T(herm=False):
        """
        Returns

        [[1, 0],[0, exp(j*pi/4*sign)]]

        where sign = 1 if herm=False and sign = -1 if herm=True

        Parameters
        ----------
        herm : bool

        Returns
        -------
        np.ndarray

        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        mat[0, 0] = 1
        if not herm:
            sign = 1
        else:
            sign = -1
        mat[1, 1] = np.exp(1j*np.pi/4*sign)
        return mat

    @staticmethod
    def mat_Tdag():
        """
        returns T^\dag

        Returns
        -------
        np.ndarray

        """
        return OneBitGates.mat_T(True)

    @staticmethod
    def mat_Sdag():
        """
        returns S^\dag

        Returns
        -------
        np.ndarray

        """
        return OneBitGates.mat_S(True)

if __name__ == "__main__":
    print('sigx= ', OneBitGates.sigx())
    print('sigy= ', OneBitGates.sigy())
    print('sigz= ', OneBitGates.sigz())

    print('had2= ', OneBitGates.had2())

    print('P_0= ', OneBitGates.P_0())
    print('P_1= ', OneBitGates.P_1())

    print('P_0_phase_fac= ', OneBitGates.P_0_phase_fac(10))
    print('P_1_phase_fac= ', OneBitGates.P_1_phase_fac(10))
    print('phase_fac= ', OneBitGates.phase_fac(10))

    mat = OneBitGates.rot(10, 20, 30)
    print('rot*rot^H= ', np.dot(mat, mat.conj().T))

    mat = OneBitGates.rot_ax(10, 1)
    print('rotx*rotx^H= ', np.dot(mat, mat.conj().T))

    mat = OneBitGates.rot_ax(10, 2)
    print('roty*roty^H= ', np.dot(mat, mat.conj().T))

    mat = OneBitGates.rot_ax(10, 3)
    print('rotz*rotz^H= ', np.dot(mat, mat.conj().T))

    print('mat_S=', OneBitGates.mat_S())
    print('mat_Sdag=', OneBitGates.mat_Sdag())
    print('mat_T=', OneBitGates.mat_T())
    print('mat_Tdag=', OneBitGates.mat_Tdag())