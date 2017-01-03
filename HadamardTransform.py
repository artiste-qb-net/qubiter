import numpy as np
from BitVector import *


class HadamardTransform:
    """
    This class contains only static methods and no constructor. Its
    functions are related to the Hadamard Transform.

    """

    @staticmethod
    def ht(num_bits, in_arr):
        """
        This function calculates the Hadamard transform of in_arr. Let H be
        the 2 dimensional Hadamard matrix and let ht be the num_bit-fold
        tensor product of H. Then this function returns the matrix product
        ht*in_arr = out_arr. in_arr and out_arr both have the same shape (
        2^num_bits,).

        Parameters
        ----------
        num_bits : int
            The number of bits. The dimension of the Hadamard tranform is
            2^num_bits

        in_arr : np.ndarray
            Input array.

        Returns
        -------
        np.ndarray

        """
        length = (1 << num_bits)
        half_len = (length >> 1)
        assert len(in_arr) == length, \
            "in_arr for Hadamard Transform has wrong length"
        prev_arr = np.zeros(length, dtype=float)
        root2 = np.sqrt(2)
        out_arr = np.copy(in_arr)
        for beta in range(num_bits):
            prev_arr[:] = out_arr[:]
            for k in range(half_len):
                x = prev_arr[2*k]
                y = prev_arr[2*k+1]
                out_arr[k] = (x + y)/root2
                out_arr[half_len+k] = (x - y)/root2
        return out_arr

    @staticmethod
    def hadamard_mat(num_bits, is_quantum=True):
        """
        This function return a numpy array with the num_bits-fold tensor
        product of the 2 dim Hadamard matrix H. If is_quantum=True (False,
        resp.), it returns a complex (real, resp.) array.

        Parameters
        ----------
        num_bits : int
        is_quantum : bool

        Returns
        -------

        """
        num_rows = (1 << num_bits)
        norma = np.sqrt(num_rows)
        if is_quantum:
            ty = complex
        else:
            ty = float
        mat = np.full((num_rows, num_rows),
                      fill_value=1/norma, dtype=ty)
        bvec = BitVector(num_bits, 0)
        for j in range(num_rows):
            for k in range(num_rows):
                bvec.dec_rep = j & k
                if bvec.get_num_T_bits() % 2 == 1:
                    mat[j, k] = - mat[j, k]
        return mat


if __name__ == "__main__":
    num_bits = 3
    length = 1 << num_bits
    in_arr = np.random.rand(length) - 0.5
    out_arr = HadamardTransform.ht(num_bits, in_arr)
    in_arr1 = HadamardTransform.ht(num_bits, out_arr)
    # print("in_arr=", in_arr)
    # print("out_arr=", out_arr)
    err = np.linalg.norm(in_arr-in_arr1)
    print("error=", err)

    print(HadamardTransform.hadamard_mat(2))
