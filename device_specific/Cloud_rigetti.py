# from pyquil import get_qc
# from pyquil.api._base_connection import ForestConnection
import numpy as np


class Cloud_rigetti:
    """
    This class has no constructor. It consists of static methods that are
    useful when communicating with the Rigetti Cloud.

    """

    @staticmethod
    def get_qc(device_name, as_qvm=False, noisy=False, **kwargs):
        """
        This method creates ForestConnection object and calls PyQuil get_qc()
        method.

        Parameters
        ----------
        device_name : str
        as_qvm : bool
        noisy : bool
        kwargs : dict

        Returns
        -------
        QuantumComputer

        """

        qvm_url = "http://127.0.0.1:5000"
        compiler_url = "http://127.0.0.1:6000"
        forest_url = "https://forest-server.qcs.rigetti.com"

        con = ForestConnection(
            sync_endpoint=qvm_url,
            compiler_endpoint=compiler_url,
            forest_cloud_endpoint=forest_url)
        qc = get_qc(device_name, connection=con,
                as_qvm=as_qvm, noisy=noisy, **kwargs)
        return qc

    @staticmethod
    def obs_vec_from_bitstrings(bitstrings, num_qbits):
        """
        This method converts a PyQuil `bitstrings` into a Qubiter
        observation vector which it then returns.

        Qubiter likes to state the results of an experiment repeated
        num_shot times by what I call an observation vector. An obs vec is a
        1-dim array, num_shots long, whose entries are integers which are
        the decimal representation of a string, num_qbits long, of zeros and
        ones. This string of zeros and ones gives the state of each qubit in
        the ZL convention

        The PyQuil authors like to state the results of an experiment
        repeated num_shot times by what they call bitstrings. A bitstrings
        is a dict that maps qubit number to a 1-dim array, num_shots long,
        of zeros and ones. Here is an example:

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

        Parameters
        ----------
        bitstrings : dict[int, np.ndarray]
        num_qbits : int

        Returns
        -------
        np.ndarray
            num_shots = bitstrings[0].shape[0], shape = (num_shots,)

        """
        assert len(bitstrings) == num_qbits
        num_shots = bitstrings[0].shape[0]
        obs_vec = np.zeros((num_shots,), dtype=int)
        bs_array = np.vstack([bitstrings[q] for q in range(num_qbits)])
        for shot in range(num_shots):
            shot_array = bs_array[:, shot]
            s = ''.join([str(shot_array[q])
                         for q in reversed(range(num_qbits))])
            obs_vec[shot] = int(s, 2)
        return obs_vec

if __name__ == "__main__":
    def main():
        bitstrings = {
            0: np.array([1, 0, 1, 1, 1, 1, 0, 1, 1, 0]),
            1: np.array([1, 1, 0, 0, 1, 0, 1, 0, 1, 0]),
            2: np.array([0, 0, 0, 1, 1, 1, 1, 0, 1, 1])}
        num_qbits = 3
        obs_vec = Cloud_rigetti.obs_vec_from_bitstrings(
            bitstrings, num_qbits)
        print(obs_vec)
    main()
