from pprint import pprint


class InterlacedPairs:

    def __init__(self, nn):
        self.nn = int(nn)
        self.pair_avail_mat = [[False] * nn for k in range(nn)]

        # only pairs with 1 <= pair[0] < pair[1]
        # are available initially
        for row in range(nn):
            for col in range(row+1, nn):
                self.pair_avail_mat[row][col] = True

        # pprint(self.pair_avail_mat)

    def pair_is_available(self, pair):
        return self.pair_avail_mat[pair[0]][pair[1]]

    def pair_len_used_up(self, pair_len):
        all_false = True
        for a in range(self.nn - pair_len):
            pair = [a, a + pair_len]
            if self.pair_is_available(pair):
                all_false = False
                break
        return all_false

    def ints_in_plist(self, plist):
        cum = set()
        # print(plist)
        for pair in plist:
            # print(pair)
            cum |= set(pair)
            # print(cum)
        return cum

    def ints_not_in_plist(self, plist):
        return [k for k in range(self.nn)
                if k not in self.ints_in_plist(plist)]

    def plist_likes_new_pair(self, plist, pair):
        return pair[0] in self.ints_not_in_plist(plist) \
                and pair[1] in self.ints_not_in_plist(plist)

    def use_pair(self, pair):
        self.pair_avail_mat[pair[0]][pair[1]] = False
        return pair

    def get_bunches_for_odd_nn(self, verbose=True):
        rem = self.nn % 2  # rem = remainder
        assert rem == 1
        num_bunches = self.nn - 1 + rem
        num_parallel_pairs = self.nn//2
        bunches = {}
        bunches[0] = [self.use_pair([a, a+1]) for a in range(0, self.nn-1, 2)]
        bunches[num_bunches-1] = [self.use_pair([a+1, a+2])
                                 for a in range(0, self.nn-2, 2)]

        for bun in range(1, num_bunches - 1):
            bunches[bun] = [self.use_pair([0, bun + 1])]

        pprint(bunches)
        for pair_len in range(2, self.nn):
            print('pair_len=', pair_len)
            for bun in reversed(range(1, num_bunches -1)):
                print('bun=', bun)
                for a in reversed(range(self.nn - pair_len)):
                    pair = [a, a + pair_len]
                    cond1 = self.pair_is_available(pair)
                    cond2 = self.plist_likes_new_pair(bunches[bun], pair)
                    print('pair, avail, liked=', pair, cond1, cond2)
                    if cond1 and cond2:
                        bunches[bun].append(self.use_pair(pair))
                        pprint(bunches)
                if self.pair_len_used_up(pair_len):
                    break
        return bunches

    def get_bunches(self):
        rem = self.nn % 2
        num_bunches = self.nn - 1 + rem
        if rem == 1:
            return self.get_bunches_for_odd_nn()
        else:
            odd_case = InterlacedPairs(self.nn - 1)
            bunches = odd_case.get_bunches_for_odd_nn()
            for bun in range(num_bunches):
                hole = odd_case.ints_not_in_plist(bunches[bun])[0]
                bunches[bun].append([hole, self.nn-1])
            return bunches

if __name__ == "__main__":

    x = InterlacedPairs(8).get_bunches()
    pprint(x)
    pprint(list(x.values()))
