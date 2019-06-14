
class BitVector:
    """
    This class uses an int called dec_rep as a vector of self.len many bits,
    where self.len <= self.max_len. The class wraps some common bitwise
    operations, and some less common ones too (like Gray coding that is
    needed by Qubiter). In some cases, the bitwise manipulation might be
    more succinct than the corresponding function in this wrapper, but the
    wrapper function's name spells out in words what is wanted.

    Attributes
    ----------
    dec_rep : int
        decimal representation, the int whose binary representation carries
        a bit vector of length self.len
    len : int
        the length of the bit vector
    max_len : int
        maximum self.len allowed
    """

    def __init__(self, length, dec_rep):
        """
        Constructor

        Parameters
        ----------
        length : int
        dec_rep : int

        Returns
        -------

        """
        self.len = length
        self.dec_rep = dec_rep
        self.max_len = 16
        assert length <= self.max_len, "bit vector is too long"
        assert length > 0, "bit vector len must be >=1"

    @staticmethod
    def copy(bvec):
        """
        Copy constructor, returns a new BitVector which is a copy of the
        BitVector bvec.

        Parameters
        ----------
        bvec : BitVector

        Returns
        -------
        BitVector

        """
        return BitVector(bvec.len, bvec.dec_rep)

    def bit_is_T(self, bpos):
        """
        Returns True iff bit at position bpos is 1 (True)

        Parameters
        ----------
        bpos : int
            bit position

        Returns
        -------
        bool

        """
        assert bpos < self.len, "bit position is too large"
        mask = (1 << bpos)
        return (self.dec_rep & mask) == mask

    def set_bit_T(self, bpos):
        """
        Sets to 1 (True) the bit of self at position bpos.

        Parameters
        ----------
        bpos : int
            bit position

        Returns
        -------
        None

        """
        assert bpos < self.len, "bit position is too large"
        self.dec_rep |= (1 << bpos)

    def set_bit_F(self, bpos):
        """
        Sets to 0 (False) the bit of self at position bpos.

        Parameters
        ----------
        bpos : int
            bit position

        Returns
        -------
        None

        """
        assert bpos < self.len, "bit position is too large"
        self.dec_rep &= ~(1 << bpos)

    def set_all_bits_T(self):
        """
        Sets to 1 (True) the bits of self at position bpos
        from 0 to len-1 inclusive.

        Returns
        -------
        None

        """
        # example: len = 3, dec_rep becomes 15 = b0111
        self.dec_rep = (1 << self.len + 1) - 1

    def set_all_bits_F(self):
        """
        Sets to 0 (False) the bits of self at positions bpos
        from 0 to len-1 inclusive.

        Returns
        -------
        None

        """
        self.dec_rep = 0

    def get_num_T_bits(self):
        """
        Returns the number of 1 (True) bits at positions bpos
        from 0 to len-1 inclusive.

        Returns
        -------
        int

        """
        count = 0
        for bpos in range(self.len):
            if self.bit_is_T(bpos):
                count += 1
        return count

    def find_T_bit_to_right_of(self, bpos):
        """
        Returns position of 1 (True) bit immediately to the right of
        position bpos. Returns -1 if there is no such bit.

        Parameters
        ----------
        bpos : int
            bit position

        Returns
        -------
        int

        """

        if bpos <= 0:
            return -1

        right_T_bit = bpos
        mask = (1 << right_T_bit)
        found_it = False
        while True:
            right_T_bit -= 1
            mask >>= 1
            found_it = ((self.dec_rep & mask) == mask)
            if right_T_bit == 0 or found_it:
                break
        if found_it:
            return right_T_bit
        else:
            return -1

    def find_T_bit_to_left_of(self, bpos):
        """
        Returns position of 1 (True) bit immediately to the left of position
        bpos. Returns -1 if there is no such bit.

        Parameters
        ----------
        bpos : int
            bit position

        Returns
        -------
        int

        """

        if bpos >= self.len-1:
            return -1

        left_T_bit = bpos
        mask = (1 << left_T_bit)
        found_it = False
        while True:
            left_T_bit += 1
            mask <<= 1
            found_it = ((self.dec_rep & mask) == mask)
            if left_T_bit == self.len-1 or found_it:
                break
        if found_it:
            return left_T_bit
        else:
            return -1

    def find_leftmost_T_bit(self):
        """
        Out of all 1 (True) bits, returns position of the leftmost one of
        those to the the left of position bpos. Returns -1 if there is no
        such bit.

        Returns
        -------
        int

        """
        if self.bit_is_T(self.len-1):
            return self.len-1
        else:
            return self.find_T_bit_to_right_of(self.len - 1)

    def find_rightmost_T_bit(self):
        """
        Out of all 1 (True) bits, returns position of the rightmost one of
        those to the the right of position bpos. Returns -1 if there is no
        such bit.

        Returns
        -------
        int

        """
        if self.bit_is_T(0):
            return 0
        else:
            return self.find_T_bit_to_left_of(0)

    def get_bit_string(self):
        """
        Returns self represented as string of length self.len of ones and
        zeros. If bit_str is the output, [int(x) for x in bit_str] will turn
        result to list of ints.

        Returns
        -------
        str

        """
        bit_str = ''
        for beta in range(self.len-1, -1, -1):
            if self.bit_is_T(beta):
                bit_str += '1'
            else:
                bit_str += '0'
        return bit_str

    def __str__(self):
        """
        Readable representation of self

        Returns
        -------
        str

        """

        return self.get_bit_string() + '=' + str(self.dec_rep)

    @staticmethod
    def new_with_T_on_diff(bvec1, bvec2):
        """
        Given two BitVectors bevc1 and bvec2, this return a BitVector which
        is a bitwise xor (mod 2 sum) of the bits of bvec1 and bvec2.

        Parameters
        ----------
        bvec1 : BitVector
        bvec2 : BitVector

        Returns
        -------
        BitVector

        """
        assert bvec1.len == bvec2.len
        return BitVector(bvec1.len, bvec1.dec_rep ^ bvec2.dec_rep)

    @staticmethod
    def get_lazy_from_normal(bit_len, normal):
        """
        Throughout Qubiter, we will often refer to "Gray Code" as "lazy
        ordering". In lazy ordering with bit_len many bits, one gives a
        sequence of bit vectors of length bit_len, so that two adjacent
        items of the sequence differ by just one bit. For example 000=0,
        100=4, 110=6, 010=2, 011=3, 111=7, 101=5, 001=1. Each element of
        this sequence represented as an int will be called lazy, and each
        int in the sequence 0, 1, 2, 3, 4, 5, 6, 7 will be called normal.
        Normal ordering is usually called dictionary ordering. Normal and
        lazy sequences both start at 0.

        Suppose bit_len = 3. The lazy sequence 000, 100, 110, 010, 011, 111,
        101, 001 is easily obtained from the "standard" lazy sequence 000,
        001, 011, 010, 110, 111, 101, 100 by "reflecting" each sequence
        term. We will use the second sequence because it is more common in
        the literature.

        References
        ----------
        1. Martin Gardener, "Knotted DoughNuts and Other
        Mathematical Entertainments", chapt. 2, "The Binary Gray Code"
        2. "Numerical Recipies in C"
        3. Many books on Discrete Mathematics for CompSci types
        4. On the web, in Eric's Treasure Trove/Math/Gray Codes

        Parameters
        ----------
        bit_len : int
        normal : int
            Function returns the lazy int that corresponds to this normal int.

        Returns
        -------
        int

        """
        lazy_bvec = BitVector(bit_len, normal)
        lazy = lazy_bvec.dec_rep
        if bit_len > 1:
            for m in range(bit_len-2, -1, -1):
                # Look at bpos = m+1, if it's ON, then flip bpos=m.
                #  Remember that ^ is same as a mod2 sum.
                lazy ^= (((normal >> m+1) & 1) << m)
        return lazy

    @staticmethod
    def lazy_advance(old_normal, old_lazy):
        """
        This method takes int "old_lazy" (which corresponds to bit vector
        "old_normal"), and changes it to the next lazy int, "new_lazy" (
        which corresponds to "new_normal").

        example:

        lazy sequence: 000, 001, 011, 010, 110, 111, 101, 100

        old_lazy = 011
        old_normal = 2 = 010

        new_normal = 3 = 011
        mask = (new_normal & ~old_normal) = 011 & 101 = 001
        new_lazy = new_normal ^ mask = 011 ^ 001 = 010


        Parameters
        ----------
        old_normal : int
        old_lazy : int

        Returns
        -------
        int, int

        """
        new_normal = old_normal + 1
        new_lazy = old_lazy ^ (new_normal & ~(new_normal-1))
        return new_normal, new_lazy

if __name__ == "__main__":
    def main():
        for length in [3, 4]:
            print('\nlength=', length)
            print('normal, lazy, lazy in binary:')
            max_normal = (1 << length) - 1
            normal = 0
            lazy = 0
            while normal <= max_normal:
                lazy_bvec = BitVector(length, lazy)
                print(normal, lazy, BitVector.get_bit_string(lazy_bvec))
                normal, lazy = BitVector.lazy_advance(normal, lazy)
    main()

