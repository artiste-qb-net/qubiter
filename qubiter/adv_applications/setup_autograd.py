"""

The purpose of this file is to install autograd and its dependencies and to
provide utility functions that are used when it is used in conjunction with
Qubiter.

When using autograd, one declares np to be the alias to module
numpy.autograd. If another file later declares np to be alias to numpy,
all sorts of error messages start cropping up. What I've done to avoid this
is to change the statements `import numpy as np` in some (not all, just the
ones called while using autograd) files by

import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np
else:
    import autograd.numpy as np

References
----------
1. https://github.com/HIPS/autograd/blob/master/docs/tutorial.md
2. https://github.com/HIPS/autograd/blob/master/docs/updateguide.md

"""
import autograd.numpy as np
from autograd import grad, jacobian
from autograd.extend import primitive, defvjp

import sys
print('np installed?', 'np' in sys.modules)  # False
print('numpy installed?', 'numpy' in sys.modules)  # True
print('autograd.numpy installed?', 'autograd.numpy' in sys.modules)  # True


def sig_all():
    """
    This method returns a numpy array of shape=(2, 3, 3) which contains the
    3 Pauli matrices in it. sigx = sig_all[:, :, 0], sigy = sig_all[:, :,
    1], sigz = sig_all[:, :, 2],

    Returns
    -------
    np.ndarray
        shape = (2, 2, 3)

    """
    sigx = np.array([[0, 1], [1, 0]])
    sigy = np.array([[0, -1j], [1j, 0]])
    sigz = np.array([[1, 0], [0, -1]])
    all_paulis = np.vstack([sigx, sigy, sigz])
    all_paulis = np.reshape(all_paulis, (3, 2, 2)).transpose(1, 2, 0)
    return all_paulis


def u2_alt(*tlist):
    """
    An alternative to OneBitGates.u2(). Both should return identical 2-dim
    matrices for identical arguments.

    Parameters
    ----------
    tlist : list[float]
        tlist = [rads0, rads1, rads2, rads3]

    Returns
    -------
    np.ndarray
        shape = (2, 2)

    """
    assert len(tlist) == 4
    t = np.sqrt(tlist[1]**2 + tlist[2]**2 + tlist[3]**2)
    if abs(t) < 1e-6:
        return np.exp(1j*tlist[0])*np.eye(2)
    tvec = np.array([tlist[1], tlist[2], tlist[3]])/t
    out = np.eye(2)*np.cos(t) + 1j*np.dot(sig_all(), tvec)*np.sin(t)
    return np.exp(1j*tlist[0])*out


def d_u2(dwrt, *tlist):
    """
    tlist is a list of 4 floats, and dwrt (which stands for "derivative with
    respect to") is in range(4). This method returns the analytical (not
    numerical, in terms of closed functions) derivative of u2(*tlist) with
    respect to tlist[dwrt].

    The output of this method has been verified by comparing it to same
    derivatives calculated numerically with autograd.

    Parameters
    ----------
    dwrt : int
    tlist : list[float]

    Returns
    -------
    np.ndarray
        shape = (2, 2)

    """
    assert dwrt in range(4)
    assert len(tlist) == 4
    if dwrt == 0:
        return 1j*u2_alt(*tlist)
    dwrt -= 1
    t = np.sqrt(tlist[1]**2 + tlist[2]**2 + tlist[3]**2)
    if abs(t) < 1e-6:
        # we already know dwrt !=0
        return np.zeros((2, 2), dtype=complex)
    tvec = np.array([tlist[1], tlist[2], tlist[3]])/t
    dotted_vec = tvec*tvec[dwrt]*np.cos(t) +\
                 (np.sin(t)/t)*(-tvec*tvec[dwrt] + np.eye(3)[dwrt, :])
    out = -np.sin(t)*tvec[dwrt]*np.eye(2) +\
          1j*np.dot(sig_all(), dotted_vec)
    return np.exp(1j*tlist[0])*out


def d_auto_u2(dwrt, *tlist):
    """
    Returns the automatic (computed by backprop) derivative of 2-dim matrix
    UnitaryMat.u2_alt. UnitaryMat.u2_alt is an alternative to
    OneBitGates.u2. Both functions return same answer for identical input (
    input is 4 real parameters in tlist).

    Parameters
    ----------
    dwrt : int
        stands for 'derivative with respect to'. int in range(4)
    tlist : list[float]
        len = 4

    Returns
    -------
    np.ndarray
        shape=(2,2)

    """

    def u2r(*tlist1):
        return np.real(u2_alt(*tlist1))

    def u2i(*tlist1):
        return np.imag(u2_alt(*tlist1))
    return jacobian(u2r, dwrt)(*tlist) + 1j*jacobian(u2i, dwrt)(*tlist)


@primitive
def pu2r(*tlist):
    """
    Returns real part of u2, and registers it as being primitive.

    Primitive means that its derivative will be provided in a defvjp (
    def of vector-jacobian-product) so no need for autograd to calculate it
    from the u2 definition.

    Parameters
    ----------
    tlist : list[float]
        len = 4

    Returns
    -------
    np.ndarray
        shape=(2,2)

    """
    return np.real(u2_alt(*tlist))


@primitive
def pu2i(*tlist):
    """
    Returns imaginary part of u2, and registers it as being primitive.

    Primitive means that its derivative will be provided in a defvjp (
    def of vector-jacobian-product) so no need for autograd to calculate it
    from the u2 definition.

    Parameters
    ----------
    tlist : list[float]
        len = 4

    Returns
    -------
    np.ndarray
        shape=(2,2)

    """
    return np.imag(u2_alt(*tlist))


def pu2(*tlist):
    """
    Returns primitive u2 as (primitive real part of u2) + j*(primtive
    imaginary part of u2).

    Parameters
    ----------
    tlist : list[float]
        len = 4

    Returns
    -------
    np.ndarray
        shape=(2,2)

    """
    # print('mmbbvv, pu2', pu2r(*tlist) +1j* pu2r(*tlist))
    return pu2r(*tlist) + 1j*pu2i(*tlist)

defvjp(pu2r,
       # defines vector-jacobian-product of pu2r
       # g.shape == pu2r.shape
       lambda ans, *tlist: lambda g: np.sum(
           g*np.real(d_u2(0, *tlist))),
       lambda ans, *tlist: lambda g: np.sum(
           g*np.real(d_u2(1, *tlist))),
       lambda ans, *tlist: lambda g: np.sum(
           g*np.real(d_u2(2, *tlist))),
       lambda ans, *tlist: lambda g: np.sum(
           g*np.real(d_u2(3, *tlist))),
       argnums=range(4))

defvjp(pu2i,
       # defines vector-jacobian-product of pu2i
       # g.shape == pu2i.shape
       lambda ans, *tlist: lambda g: np.sum(
           g*np.imag(d_u2(0, *tlist))),
       lambda ans, *tlist: lambda g: np.sum(
           g*np.imag(d_u2(1, *tlist))),
       lambda ans, *tlist: lambda g: np.sum(
           g*np.imag(d_u2(2, *tlist))),
       lambda ans, *tlist: lambda g: np.sum(
           g*np.imag(d_u2(3, *tlist))),
       argnums=range(4))


def d_auto_pu2(dwrt, *tlist):
    """
    Returns the automatic derivative of pu2. We have defined things so that
    this derivative is stipulated analytically a priori rather than being
    calculated by autograd from def of u2.

    Parameters
    ----------
    dwrt : int
        stands for 'derivative with respect to'. int in range(4)
    tlist : list[float]
        len = 4

    Returns
    -------
    np.ndarray
        shape=(2,2)

    """
    assert dwrt in range(4)
    return jacobian(pu2r, dwrt)(*tlist) + 1j*jacobian(pu2i, dwrt)(*tlist)

if __name__ == "__main__":
    from qubiter.OneBitGates import *

    def main():
        print("\nu2_alt example-------------")
        ex = np.array([1, 0, 0])
        ey = np.array([0, 1, 0])
        ez = np.array([0, 0, 1])
        all_paulis = sig_all()
        sigx_ = np.dot(all_paulis, ex)
        sigy_ = np.dot(all_paulis, ey)
        sigz_ = np.dot(all_paulis, ez)
        print('sigx_=\n', sigx_)
        print('sigy_=\n', sigy_)
        print('sigz_=\n', sigz_)

        rads_list = [.1, .2, .3, .4]
        err = np.linalg.norm(OneBitGates.u2(*rads_list) -
                       u2_alt(*rads_list))
        print('err=', err)

        tlist = [.3, 1.1, .7, .5]
        for dwrt in range(4):
            print('err=', np.linalg.norm(
                d_auto_u2(dwrt, *tlist) - d_u2(dwrt, *tlist)))
        for dwrt in range(4):
            print('err=', np.linalg.norm(
                d_auto_pu2(dwrt, *tlist) - d_u2(dwrt, *tlist)))
    main()
