import numpy as np


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


def centered_rads1(ang_rads_list):
    """
    Same as centered_rads() but for list

    Parameters
    ----------
    ang_rads_list : list[float]

    Returns
    -------
    list[float]

    """

    return [centered_rads(ang) for ang in ang_rads_list]


def is_prob_dist(pd):
    """
    Returns True iff pd is a probability distribution

    Parameters
    ----------
    pd : np.ndarray

    Returns
    -------
    bool

    """
    return abs(np.sum(pd)-1) < 1e-5


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


def is_arr(x):
    """
    Returns True iff x is a numpy array.

    Parameters
    ----------
    x :

    Returns
    -------
    bool

    """
    return isinstance(x, np.ndarray)


def is_diag_mat(arr):
    """
    Returns True iff arr is numpy array for diagonal square matrix

    Parameters
    ----------
    arr : np.ndarray

    Returns
    -------
    bool

    """
    assert is_arr(arr)

    num_rows = arr.shape[0]
    assert arr.shape == (num_rows, num_rows)
    # this extracts diagonal v, then
    # creates a diagonal matrix with v as diagonal
    arr1 = np.diag(np.diag(arr))
    return np.linalg.norm(arr - arr1) < 1e-6


def is_const_mat(arr):
    """
    Returns True iff arr is numpy array for constant square matrix.

    Parameters
    ----------
    arr : np.ndarray

    Returns
    -------
    bool

    """

    if not is_diag_mat(arr):
        return False
    num_rows = arr.shape[0]
    arr1 = arr[0, 0]*np.ones((num_rows,))
    arr2 = np.diag(arr)
    return np.linalg.norm(arr1 - arr2) < 1e-6


def log_print(x):
    """
    Prints file name of log_print() call, then file line of log_print()
    call, then x

    Parameters
    ----------
    x : object

    Returns
    -------
    None

    """
    from inspect import getframeinfo, stack
    caller = getframeinfo(stack()[1][0])
    print(caller.filename, "line=", caller.lineno, ":\n", x)

