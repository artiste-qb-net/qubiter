
def get_dir_edges_from_c_to_tars(c_to_tars):
    """
    returns tuple of all allowed directed edges (c, t) where c control
    and t target.

    Parameters
    ----------
    c_to_tars : dict[int, list[int]]
        a dictionary mapping j in range(num_qbits) to a list, possibly
        empty, of the physically allowed targets of qubit j, when j is
        the control of a CNOT.

    Returns
    -------
    tuple[tuple[int, int]]

    """
    dir_edges = []
    # print(c_to_tars)
    for c, tars in c_to_tars.items():
        # print('****', c, tars)
        for t in tars:
            dir_edges.append((c, t))
    return tuple(dir_edges)
