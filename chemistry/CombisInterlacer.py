import itertools as it
from pprint import pprint


class CombisInterlacer:
    """
    Combis is short for combinations. A combi will be represented by
    a list of distinct ints. The integers in a combi are in range(n).
    Specifying a subn is optional. If subn is specified, where subn < n,
    then all combis contain subn elements. If subn is not specified, combis
    may be of arbitrary length.

    This class stores combis in a dictionary bunch_to_combis such that the
    key is an integer labelling the bunch, and the value is a list of
    combis. The combis in a bunch are interlaced, meaning that they are
    disjoint as sets.

    Attributes
    ----------
    bunch_to_avail : dict[int:bool]
        Bunch number to availability. The availability is an array of n bools
    bunch_to_combis : dict[int:list[list[int]]]
        bunch to combinations
    combi_list : list[list[int]]
        optional input. If specified, the combis in this list are separated
        into bunches, with the combis in each bunch being interlaced.
    n : int
        combis will assume values in range(n)
    num_bunches: int
        number of bunches
    subn : int
        optional input. If specified, it means you want all combis to
        contain subn elements, where subn < n
    verbose : bool
    """

    def __init__(self, n, subn=None, combi_list=None, verbose=False):
        """
        Constructor

        Parameters
        ----------
        n : int
        subn : int
        combi_list : list[list[int]]
        verbose : bool

        Returns
        -------

        """
        self.n = n
        self.subn = subn
        if subn:
            assert subn < n
        self.bunch_to_avail = {0: [True]*n}
        self.bunch_to_combis = {0: []}
        self.num_bunches = 1
        self.combi_list = combi_list
        self.verbose = verbose

        self.generate_bunches()

    def combi_fits_in_this_bunch(self, combi, bun):
        """
        This function tells whether or not combi fits into the bunch
        labelled by the int bun


        Parameters
        ----------
        combi : list[int]
        bun : int

        Returns
        -------
        bool

        """
        it_fits = True
        combi_len = len(combi)
        for k in range(combi_len):
            if not self.bunch_to_avail[bun][combi[k]]:
                it_fits = False
                break
        return it_fits

    def insert_combi(self, combi):
        """
        This function inserts combi into bunches_to_combis. It checks to see
        if combi fits into any of the bunches that already exist. If not,
        it creates a new bunch in bunches_to_combis and places combi inside
        that one.


        Parameters
        ----------
        combi : list[int]

        Returns
        -------
        None

        """
        is_in = False
        if self.verbose:
            print('combi=', combi)
        combi_len = len(combi)
        for bun in range(self.num_bunches):
            if self.verbose:
                print('bun=', bun)
            if self.combi_fits_in_this_bunch(combi, bun):
                for k in range(combi_len):
                    self.bunch_to_avail[bun][combi[k]] = False
                self.bunch_to_combis[bun].append(combi)
                is_in = True
                if self.verbose:
                    pprint(self.bunch_to_avail)
                    pprint(self.bunch_to_combis)
                break
        if not is_in:
            self.num_bunches += 1
            bun = self.num_bunches - 1
            if self.verbose:
                print('bun=', bun)
            self.bunch_to_avail[bun] = [True] * self.n
            for k in range(combi_len):
                self.bunch_to_avail[bun][combi[k]] = False
            self.bunch_to_combis[bun] = [combi]
            if self.verbose:
                pprint(self.bunch_to_avail)
                pprint(self.bunch_to_combis)

    def generate_bunches(self):
        """
        This function generates all bunches in bunch_to_combis. If a
        combi_list is specified, bunches are made out of combis in that
        list. Otherwise, if a combis_list is not specified, the list is
        assumed to be the list of all length subn combinations of n letters.

        Returns
        -------
        None


        """

        if self.combi_list is not None:
            for combi in self.combi_list:
                self.insert_combi(combi)
        elif self.subn is not None:
            for combi in it.combinations(range(self.n), self.subn):
                self.insert_combi(combi)
        else:
            assert False

    def make_bunches_into_tuples(self):
        """
        This function changes bunches from lists into tuples.

        Returns
        -------
        None

        """
        for bun in range(self.num_bunches):
            self.bunch_to_combis[bun] = tuple(self.bunch_to_combis[bun])

    def bunch_gen(self):
        """
        This function is a generator that yields all stored bunches, one at
        a time.

        Returns
        -------
        list[list[int]]

        """

        self.make_bunches_into_tuples()
        x = self.bunch_to_combis
        arr = tuple(x.values())

        num = 0
        while num < len(arr):
            yield arr[num]
            num += 1

    def reversed_bunch_gen(self):
        """
        This function is a generator that yields all stored bunches, one at
        a time. It reverses the order of bunch_gen()

        Returns
        -------
        list[list[int]]

        """

        self.make_bunches_into_tuples()
        x = self.bunch_to_combis
        arr = tuple(x.values())

        num = 0
        tot_num = len(arr)
        while num < tot_num:
            yield arr[tot_num - 1 - num]
            num += 1

if __name__ == "__main__":

    ci = CombisInterlacer(n=12, subn=2)
    x = ci.bunch_to_combis

    pprint(x)
    pprint(list(x.values()))

    ci.make_bunches_into_tuples()

    pprint(tuple(x.values()))

