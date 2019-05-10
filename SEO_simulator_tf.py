from SEO_simulator import *
import tensorflow as tf
tf.enable_eager_execution()


class SEO_simulator_tf(SEO_simulator):
    """
    TF= TensorFlow

    This class is a child of SEO_simulator. The class replaces numpy methods
    by TensorFlow (2.0, Eager) methods when those methods perform the state
    vector evolution. This allows those tf methods to be taped, and the tape
    to be run in reverse to calculate gradients.

    At the beginning (in the constructor), the class converts the arrays in
    `self.cur_st_vec_dict` from numpy to tf.

    Certain "interface" methods of the parent class were written assuming
    that the arrays of `self.cur_st_vec_dict` are numpy. Such methods call
    convert_tensors_to_numpy() at their beginning to convert the tf arrays
    to numpy. They also call convert_tensors_to_tf() at their end to
    convert the numpy tensors back to tf, the way they were upon starting
    the method.


    Attributes
    ----------

    """

    def __init__(self, file_prefix, num_bits, init_st_vec=None, **kwargs):
        """
        Constructor

        Parameters
        ----------
        file_prefix : str
        num_bits : int
        init_st_vec : StateVec
        kwargs : dict
            keyword args of SEO_simulator

        Returns
        -------

        """
        SEO_simulator.__init__(self, file_prefix, num_bits,
                               init_st_vec=init_st_vec, **kwargs)
        self.convert_tensors_to_tf()
        SEO_simulator.transpose = tf.transpose
        SEO_simulator.tensordot = tf.tensordot
        self.use_tf = True

    def convert_tensors_to_tf(self):
        """
        This method converts the arrays in `cur_st_vec_dict` from numpy to
        tensorflow.

        Returns
        -------
        None

        """
        for br_key, st_vec in self.cur_st_vec_dict.items():
            if isinstance(st_vec.arr, np.ndarray):
                st_vec.arr = tf.convert_to_tensor(st_vec.arr)

    def convert_tensors_to_numpy(self):
        """
        This method converts the arrays in `cur_st_vec_dict` from tensorflow
        to numpy.

        Returns
        -------
        None

        """
        for br_key, st_vec in self.cur_st_vec_dict.items():
            if isinstance(st_vec.arr, tf.Tensor):
                st_vec.arr = st_vec.arr.numpy()

    def do_tf_ruse(self, br_key, slicex, sub_arr):
        """
        internal function used in evolve_ methods iff tf eager is on. Should
        have same effect as

        self.cur_st_vec_dict[br_key].arr[slicex] = sub_arr

        Parameters
        ----------
        br_key : str
        slicex : tuple
        sub_arr : tf.Tensor

        Returns
        -------
        None

        """
        test = False
        arr = self.cur_st_vec_dict[br_key].arr
        if test:
            arr1 = tf.Tensor(arr)
            # print('arr1 bef', arr1)
            arr1[slicex] = sub_arr
            # print('arr1 aft', arr1)
        on_slicex = tf.fill(arr.shape, False)
        on_slicex[slicex] = True
        bigger_shape = [1]*self.num_bits  # slicex is num_bits long
        k = 0
        # print('wwwww', sub_arr.shape, slicex)
        for bit, kind in enumerate(slicex):
            if kind not in [0, 1]:
                bigger_shape[bit] = sub_arr.shape[k]
                k += 1
        sub_arr = tf.reshape(sub_arr, tuple(bigger_shape))
        arr = \
            arr*tf.cast(tf.math.logical_not(on_slicex), dtype=tf.int32)\
            + sub_arr*tf.cast(on_slicex, dtype=tf.int32)
        self.cur_st_vec_dict[br_key].arr = arr
        if test:
            # print('arr aft', arr)
            print('testing tf simulator')
            assert tf.linalg.norm(arr-arr1) < 1e-6, 'tf sim test failed'


if __name__ == "__main__":
    def main():
        # use test = 0 if want to run all tests at once.
        test = 0
        # test = 3
        if test in [0, 1]:
            # test on circuit for a quantum fourier transform
            # (no loops, no internal measurements)
            sim = SEO_simulator_tf('io_folder/sim_test1', 6, verbose=True)

        if test in [0, 2]:
            # test embedded loops
            sim = SEO_simulator('io_folder/sim_test2', 4, verbose=True)

        if test in [0, 3]:
            # test MEAS branching. Each kind 2 measurement doubles number of
            # branches
            sim = SEO_simulator('io_folder/sim_test3', 4, verbose=True)

    main()


