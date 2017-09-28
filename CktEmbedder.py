from Controls import Controls


class CktEmbedder:
    """
    This class stores parameters for embedding a quantum circuit within a
    larger circuit (one with more qubits). Note that we will use the word
    "bit" to mean either cbits or qbits.

    Attributes
    ----------
    bit_map : list[int]
        1-1 but not onto map of range(num_bits_bef) into range(num_bits_aft)
    extra_controls : Controls
        When embedding controls, these extra ones will be added to the before
        embedding controls
    num_bits_aft : int
        number of qubits after circuit embedding
    num_bits_bef : int
        number of qubits before circuit embedding

    """

    def __init__(self, num_bits_bef, num_bits_aft, bit_map=None):
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
        assert num_bits_bef <= num_bits_aft
        if num_bits_aft > num_bits_bef:
            assert bit_map, "must give a bit_map"
        self.bit_map = bit_map
        self.extra_controls = Controls(num_bits_aft)

    def get_num_new_bits(self):
        """
        Returns num_bits_aft - num_bits_bef

        Returns
        -------
        int

        """
        return self.num_bits_aft - self.num_bits_bef

    def is_identity_map(self):
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
        if self.is_identity_map():
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
        if self.is_identity_map():
            return True
        else:
            return aft in self.bit_map

    @staticmethod
    def composition(emb2, emb1):
        """
        Returns composite map emb2(emb1()). This entails defining a new
        CktEmbedder whose bit_map is the composition of the bit maps of the
        other two. Also, the extra controls of emb2 and emb1 must be merged
        and given to the composition.

        Parameters
        ----------
        emb1 : CktEmbedder
        emb2 : CktEmbedder

        Returns
        -------
        CktEmbedder

        """
        if emb1.is_identity_map():
            return emb2
        if emb2.is_identity_map():
            return emb1
        assert emb1.num_bits_aft == emb2.num_bits_bef,\
            "can't chain embedders"

        bit_map = [emb2.aft(emb1.aft(k)) for
                         k in range(emb1.num_bits_bef)]
        compo_emb = CktEmbedder(
            emb1.num_bits_bef,
            emb2.num_bits_aft,
            bit_map)

        # new_dict is an empty dictionary initially
        new_dict = compo_emb.extra_controls.bit_pos_to_kind
        # add after controls of emb2
        new_dict.update(emb2.extra_controls.bit_pos_to_kind)
        # add after controls of emb1, after mapping them
        if emb1.extra_controls:
            old_dict = emb1.extra_controls.bit_pos_to_kind
            old_dict1 = {emb2.aft(bit): old_dict[bit] for bit in old_dict}
            new_dict.update(old_dict1)
        compo_emb.extra_controls.refresh_lists()

        return compo_emb

if __name__ == "__main__":
    print(5)
