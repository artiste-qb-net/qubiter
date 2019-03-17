# from pyquil import get_qc
# from pyquil.api._base_connection import ForestConnection
import numpy as np


class Cloud_rigetti:
    """
    This class has no constructor. It consists of static methods that are
    useful when communicating with the Rigetti Cloud.

    """

    # @staticmethod
    # def get_qc(device_name, noisy=False, **kwargs):
    #     """
    #     This method creates ForestConnection object and calls PyQuil get_qc()
    #     method.
    #
    #     Parameters
    #     ----------
    #     device_name : str
    #     noisy : bool
    #     kwargs : dict
    #
    #     Returns
    #     -------
    #     QuantumComputer
    #
    #     """
    #
    #     qvm_url = "http://127.0.0.1:5000"
    #     compiler_url = "http://127.0.0.1:6000"
    #     forest_url = "https://forest-server.qcs.rigetti.com"
    #
    #     con = ForestConnection(
    #         sync_endpoint=qvm_url,
    #         compiler_endpoint=compiler_url,
    #         forest_cloud_endpoint=forest_url)
    #     qc = get_qc(device_name, connection=con,
    #               noisy=noisy, **kwargs)
    #     return qc

    @staticmethod
    def obs_vec_from_bitstrings(bitstrings, num_qbits, bs_is_array):
        """
        This method converts a PyQuil `bitstrings` into a Qubiter
        observation vector which it then returns.

        Qubiter likes to state the results of an experiment repeated
        num_shot times by what it calls an observation vector. An obs vec is
        a 1-dim array, num_shots long, whose entries are integers which are
        the decimal representation of a string, num_qbits long, of zeros and
        ones. This string of zeros and ones gives the state of each qubit in
        the ZL convention

        PyQuil likes to state the results of an experiment repeated num_shot
        times by what it calls bitstrings. A bitstrings is a dict that maps
        qubit number to a 1-dim array, num_shots long, of zeros and ones.
        Here is an example:

        # [1]:
        #
        # bitstrings = qc.run_and_measure(program, trials=10)
        # bitstrings
        #
        # [2]:
        #
        # {0: array([1, 0, 0, 1, 1, 1, 1, 0, 1, 0]),
        #  1: array([1, 0, 0, 1, 1, 1, 1, 0, 1, 0]),
        #  2: array([1, 0, 0, 1, 1, 1, 1, 0, 1, 0])}

        However, qc.run() returns a numpy array of zeros and ones and shape
        (num_shots, num_qbits), formed from the bitstrings dict just
        described. If bs_is_array=False, we assume the input bitstrings is a
        dict, and if True, we assume it is an array.


        Parameters
        ----------
        bitstrings : dict[int, np.ndarray]
        num_qbits : int
        bs_is_array : bool
            stands for: bitstrings is array

        Returns
        -------
        np.ndarray
            shape = (num_shots,)

        """
        if not bs_is_array:
            assert isinstance(bitstrings, dict)
            assert len(bitstrings) == num_qbits,\
                "for num_qbits = " + str(num_qbits) + ' got bitstrings=\n' +\
                str(bitstrings)
            num_shots = bitstrings[0].shape[0]
            bs_array = np.vstack([bitstrings[q] for q in range(num_qbits)])
        else:
            assert isinstance(bitstrings, np.ndarray)
            assert bitstrings.shape[1] == num_qbits
            num_shots = bitstrings.shape[0]
            bs_array = bitstrings.T

        obs_vec = np.zeros((num_shots,), dtype=int)
        for shot in range(num_shots):
            shot_array = bs_array[:, shot]
            s = ''.join([str(shot_array[q])
                         for q in reversed(range(num_qbits))])
            obs_vec[shot] = int(s, 2)
        return obs_vec

if __name__ == "__main__":
    def main1():
        bitstrings = {
            0: np.array([1, 0, 1, 1, 1, 1, 0, 1, 1, 0]),
            1: np.array([1, 1, 0, 0, 1, 0, 1, 0, 1, 0]),
            2: np.array([0, 0, 0, 1, 1, 1, 1, 0, 1, 1])}
        num_qbits = 3
        obs_vec = Cloud_rigetti.obs_vec_from_bitstrings(
            bitstrings, num_qbits, bs_is_array=False)
        print(obs_vec)

    def main2():
        bitstrings = np.vstack([
            [1,1,0], [0,1,0], [1,0,0], [1,0,1], [1,1,1],
            [1,0,1], [0,1,1], [1,0,0], [1,1,1], [0,0,1]
        ])
        num_qbits = 3
        obs_vec = Cloud_rigetti.obs_vec_from_bitstrings(
            bitstrings, num_qbits, bs_is_array=True)
        print(obs_vec)
    main1()
    main2()
