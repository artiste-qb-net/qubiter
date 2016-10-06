import numpy as np

TOL = 1E-8


def centered_rads(ang_rads):
    """
    takes any real number and returns a number between -pi and pi that is
    equal to the original one mod 2pi

    Parameters
    ----------
    ang_rads : float
        angle in radians

    Returns
    -------
    float

    """
    ang_rads = np.fmod(ang_rads, 2*np.pi)
    if ang_rads > np.pi:
        ang_rads -= 2*np.pi
    return ang_rads


def increment_dict(di, key, inc, initial=0):
    """
    Increments dictionary entry at position 'key' by inc. If item at
    position key does not exist coming in, first creates one with value
    initial, then increments it by inc.

    Parameters
    ----------
    di : dict[]
    key :
    inc :
        increment
    initial :

    Returns
    -------

    """
    if key not in di:
        di[key] = initial
    di[key] += inc


def pp_numpy_arr(arr, omit_zero_amps=False):
    """
    pp=pretty print. Print numpy array as column of (index, array value)
    pairs of the form (i, j, k, ...) arr[i, j, k, ...], so zero bit first.

    Parameters
    ----------
    arr : numpy.array
    omit_zero_amps : bool
        If True, will not list states with zero amplitude

    Returns
    -------
    None

    """
    if not omit_zero_amps:
        for index, x in np.ndenumerate(arr):
            print(index, x)
    else:
        for index, x in np.ndenumerate(arr):
            if np.absolute(x) > 1E-6:
                print(index, x)


