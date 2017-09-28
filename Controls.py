# from CktEmbedder import *

class Controls:
    """
    This class stores a dictionary called self.bit_pos_to_kind containing
    key-value pairs of the form (a control's qubit number: its kind). Kinds
    can be either a bool or a non-negative integer. Kind is True if control
    is P_1 = n = |1><1|. Kind is False if control is P_0 = nbar = |0><0|.
    Kind is a non-negative integer for MP_Y and DIAG controls.

    Attributes
    ----------
    bit_pos : list[int]
        This is the key's half if you unzip bit_pos_to_kind, in decreasing
        order
    bit_pos_to_kind : dict[int, bool|int]
        Dictionary matching control bit position with its kind.
        The domain of the map bit_pos_to_kind is a subset of range(num_bits)
    kinds : list[bool|int]
        this is the value's half if you unzip bit_pos_to_kind.
    num_bits : int
        number of qubits in full quantum circuit

    """
    # combines my JAVA classes Control, TFControls, MultiControls

    def __init__(self, num_bits):
        """
        Constructor

        Parameters
        ----------
        num_bits : int

        Returns
        -------

        """
        self.num_bits = num_bits
        # num_controls = len(bit_pos_to_kind) = len(bit_pos) = len(kinds)
        self.bit_pos_to_kind = {}
        self.bit_pos = []
        self.kinds = []

    @staticmethod
    def copy(old, extra_bits=0):
        """
        Create a copy of self but add extra_bit many uncontrolled bits at the
        end.

        Parameters
        ----------
        old : Controls
        extra_bits : int

        Returns
        -------
        Controls

        """
        new = Controls(old.num_bits + extra_bits)
        new.bit_pos_to_kind = old.bit_pos_to_kind
        new.bit_pos = old.bit_pos
        new.kinds = old.kinds
        return new

    def set_control(self, bit_pos, kind, do_refresh=False):
        """
        Add key-value pair (bit_pos: kind) to self.bit_pos_to_kind dictionary

        Parameters
        ----------
        bit_pos : int
        kind : bool|int
        do_refresh : bool

        Returns
        -------
        None

        """
        assert -1 < bit_pos < self.num_bits, \
            "bit position is out of range"
        self.bit_pos_to_kind[bit_pos] = kind
        if do_refresh:
            self.refresh_lists()

    def refresh_lists(self):
        """
        Replace the 2 lists bit_pos and kinds by unzipping (actually zip()
        does both zip and unzip) the dictionary bit_pos_to_kind.

        Returns
        -------
        None

        """
        # li is list of tuples.
        # sort li items so that in decreasing order of first entry of tuples
        if not self.bit_pos_to_kind:
            self.bit_pos = []
            self.kinds = []
            return
        li = sorted(
            self.bit_pos_to_kind.items(), key=lambda t: t[0], reverse=True)
        self.bit_pos, self.kinds = zip(*li)


    @staticmethod
    def new_knob(num_bits, bit_pos, kind):
        """
        Returns a single control

        Parameters
        ----------
        num_bits : int
        bit_pos : int
        kind : bool|int

        Returns
        -------
        Controls

        """

        new = Controls(num_bits)
        new.set_control(bit_pos, kind)
        new.refresh_lists()
        return new

    def is_control(self, bit_pos):
        """
        True if bit_pos in the keys of bit_pos_to_kind. False otherwise.

        Parameters
        ----------
        bit_pos : int

        Returns
        -------
        bool

        """
        return bit_pos in self.bit_pos

    def get_workspace_bit_pot(self, target_bit_pos):
        """
        Find a bit that is neither the target bit nor one of the controls.
        Make it as big as possible.

        Parameters
        ----------
        target_bit_pos : int

        Returns
        -------
        int

        """
        num_controls = len(self.bit_pos_to_kind)
        # must have room for num_controls + one target + one workspace bit
        assert self.num_bits >= num_controls + 2,\
            "must have more qubits to find a working space bit"

        li = [target_bit_pos] + self.bit_pos
        for k in reversed(range(self.num_bits)):
            if k not in li:
                break
        return k

    def get_num_int_controls(self):
        """
        Find number of controls that are not of boolean kind.

        Returns
        -------
        int

        """
        num = 0
        num_controls = len(self.kinds)
        for k in range(num_controls):
            if not isinstance(self.kinds[k], bool):
                num += 1
        return num

    def new_embedded_self(self, emb):
        """
        In bit_pos_to_kind, replace each key by a new key emb.bit_map[key].
        Also add extra controls carried by emb.

        Parameters
        ----------
        emb : CktEmbedder

        Returns
        -------
        Controls

        """
        if emb.is_identity_map():
            return self
        if not self.bit_pos_to_kind:
            return self

        assert self.num_bits == emb.num_bits_bef,\
            "controls cannot be reset with CktEmbedder"
        new = Controls(emb.num_bits_aft)

        self.refresh_lists()
        bef_num_controls = len(self.bit_pos)
        new.bit_pos_to_kind = {emb.aft(self.bit_pos[c]): self.kinds[c] for
                    c in range(bef_num_controls)}
        # embedding can add new controls
        extra_trols = emb.extra_controls.bit_pos_to_kind
        if extra_trols:
            new.bit_pos_to_kind.update(extra_trols)
        new.refresh_lists()
        return new

if __name__ == "__main__":
    print(5)


