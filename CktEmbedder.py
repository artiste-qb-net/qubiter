

class CktEmbedder:
    """
    This class stores parameters for embedding a quantum circuit within a
    larger circuit (one with more qubits). Note that we will use the words
    bit to mean either cbits or qbits.

    Attributes
    ----------
    num_bits_bef : int
        number of qubits before circuit embedding
    num_bits_aft : int
        number of qubits after circuit embedding
    bit_map : list[int]
        1-1 but not onto map of range(num_bits_bef) into range(num_bits_aft)

    """

    def __init__(self, num_bits_bef, num_bits_aft):
        """
        Constructor

        Parameters
        ----------
        num_bits_bef : int
        num_bits_aft : int

        Returns
        -------

        """
        self.num_bits_bef = num_bits_bef
        self.num_bits_aft = num_bits_aft
        self.bit_map = None

    def get_num_new_bits(self):
        """
        Returns num_bits_aft - num_bits_bef

        Returns
        -------
        int

        """
        return self.num_bits_aft - self.num_bits_bef

    def is_identity(self):
        """
        Returns True (False) if bit_map is None (list)

        Returns
        -------
        bool

        """
        return self.bit_map is None

    def aft(self, bef):
        """
        Returns bit_map[bef] if bit_map is list. If bit_map
        is None, returns input 'bef' untouched.

        Parameters
        ----------
        bef : int

        Returns
        -------
        int

        """
        if self.is_identity():
            return bef
        else:
            return self.bit_map[bef]

    def is_old_bit(self, aft):
        """
        Returns True if aft is in image of bit_map and False otherwise.

        Parameters
        ----------
        aft : int

        Returns
        -------
        bool

        """
        if self.is_identity():
            return True
        else:
            return aft in self.bit_map

    def get_composite_map_of_you_followed_by(self, emb):
        """
        Returns composite map emb(self()).

        Parameters
        ----------
        emb : CktEmbedder

        Returns
        -------
        CktEmbedder

        """
        if self.is_identity():
            return emb
        if emb.is_identity():
            return self
        assert self.num_bits_aft == emb.num_bits_bef,\
            "can't chain embedders"
        new = CktEmbedder(self.num_bits_bef, emb.num_bits_aft)

        new.bit_map = [0]*new.num_bits_bef
        for k in range(new.num_bits_bef):
            new.bit_map[k] = emb.aft(self.aft(k))

        return new

if __name__ == "__main__":
    print(5)
