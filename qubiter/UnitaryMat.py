import numpy as np
from qubiter.OneQubitGate import *


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

        params_from_u2() maps arr -> (delta, rot_rads, unit_vec). This
        function maps (delta, rot_rads, unit_vec) -> arr so it is the
        inverse of params_from_u2(). But be careful, because delta, rot_rads
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
        arr = np.exp(1j*delta)*OneQubitGate.rot(x_rads, y_rads, z_rads)
        return arr

    @staticmethod
    def params_from_u2(arr):
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
        delta, rot_rads, unit_vec = UnitaryMat.params_from_u2(arr)
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


if __name__ == "__main__":
    from FouSEO_writer import *

    def main():
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
        delta1, rot_rads1, unit_vec1 = UnitaryMat.params_from_u2(arr_in)
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

    main()

