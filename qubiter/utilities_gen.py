from functools import reduce
import os
import posixpath
import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np
else:
    import autograd.numpy as np

# def centered_rads(ang_rads):
#     """
#     Takes any real number and returns a number between -pi and pi that is
#     equal to the original one mod 2pi.
#
#     Parameters
#     ----------
#     ang_rads : float
#         angle in radians
#
#     Returns
#     -------
#     float
#
#     """
#     ang_rads = np.fmod(ang_rads, 2*np.pi)
#     if ang_rads > np.pi:
#         ang_rads -= 2*np.pi
#     return ang_rads
#
#
# def centered_rads1(ang_rads_list):
#     """
#     Same as centered_rads() but for list
#
#     Parameters
#     ----------
#     ang_rads_list : list[float]
#
#     Returns
#     -------
#     list[float]
#
#     """
#
#     return [centered_rads(ang) for ang in ang_rads_list]


def is_prob_dist(pd):
    """
    Returns True iff pd is a probability distribution.

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
    Increments dictionary entry at position `key` by inc. If item at
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
    Returns True iff arr is numpy array for diagonal square matrix.

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
    call, then x.

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


def all_strings(li):
    """
    Returns True iff all items in list are strings.

    Parameters
    ----------
    li : list

    Returns
    -------
    bool

    """
    return all([isinstance(x, str) for x in li])


def all_floats(li):
    """
    Returns True iff all items in list are floats.

    Parameters
    ----------
    li : list

    Returns
    -------
    bool

    """
    return all([isinstance(x, float) for x in li])


def is_non_neg_int(s):
    """
    Returns True iff string s is a non-negative number.

    Parameters
    ----------
    s : str

    Returns
    -------
    bool

    """
    # s.isdigit() with s equal to
    # '0', '1', '001' gives True,
    # '-1', '1.' gives False
    if not s.isdigit():
        return False
    if s[0] == '0' and len(s) > 1:
        return False
    return True


def scalar_prod(scalars_list):
    """
    This method returns the product of the list of scalars which it has as
    input.

    Parameters
    ----------
    scalars_list : list[int|float|complex] | tuple[int|float|complex]

    Returns
    -------
    complex|float|int

    """
    if len(scalars_list) == 1:
        return scalars_list[0]
    else:
        return reduce(lambda x, y: x*y, scalars_list)


def kron_prod(mat_list):
    """
    This method returns the Kronecker product of the list of matrices which
    it has as input.

    Parameters
    ----------
    mat_list : list[np.ndarray]

    Returns
    -------
    np.ndarray

    """
    num_mats = len(mat_list)
    prod = mat_list[0]
    for k in range(1, num_mats, 1):
        prod = np.kron(prod, mat_list[k])
    return prod


def get_eng_file_rel_path(file_prefix, num_qbits):
    """
    Returns path to English file.

    Returns
    -------
    str

    """

    return file_prefix + '_' + str(num_qbits) +\
           '_eng.txt'


def get_pic_file_rel_path(file_prefix, num_qbits, ZL=True):
    """
    Returns path to Picture file.

    Returns
    -------
    str

    """

    return file_prefix + '_' + str(num_qbits) +\
            ('_ZL' if ZL else '_ZF') + 'pic.txt'


def get_value(kwargs, key_str, default_val=None):
    """
    Returns kwargs[key_str] if there is one. Else it returns default_val if
    there is one. Else aborts.

    Parameters
    ----------
    kwargs : dict[str, float]
    key_str : str
    default_val : float

    Returns
    -------
    float

    """
    if key_str in kwargs:
        return kwargs[key_str]
    elif default_val:
        return default_val
    else:
        assert False, "must pass-in keyword " + key_str +\
            ' in ' + str(kwargs)


def find_path_to_qubiter():
    """
    Returns absolute path to this file.

    Returns
    ------
    str

    """
    from inspect import getsourcefile
    path = getsourcefile(find_path_to_qubiter)
    return posixpath.normpath(path)


def preface(a_str):
    """
    Throughout Qubiter, the term `file_prefix` is used for files. If the
    file_prefix string starts with `_`, then the file is created with
    relative path equal to (hence, it shows up in the current working
    directory)

    `./file_prefix + ending`

    If file_prefix doesn't start with `_`, then the file is created with
    absolute path equal to (hence, it shows up in qubiter's io_folder)

    `absolute_path_to_io_folder/file_prefix + ending`.

    Given a_str, if it doesn't start with `_`, this method returns
    `absolute_path_to_io_folder/a_str`. Otherwise, this method just returns
    a_str.

    Parameters
    ----------
    a_str : str

    Returns
    -------
    str

    """
    if a_str[0] == '_':
        return a_str
    # this is something/qubiter/qubiter/utilities_gen
    path1 = find_path_to_qubiter()
    # this is something/qubiter/qubiter
    path1 = os.path.split(path1)[0]
    # using os.path instead of posixpath led to errors in open(file_path) if
    # file_path had directories with blank spaces (whitespace)
    return posixpath.normpath(path1 + '/io_folder/' + a_str)


if __name__ == "__main__":
    def main():
        print('---------------------')
        print(find_path_to_qubiter())
        print(preface('A/B/test.py'))

    main()
