"""

The purpose of this file is to install autograd and scripts it depends on
and to provide utility functions that are used in conjunction with it.

When using autograd, one declares np to be the alias to module
numpy.autograd. If another file later declares np to be alias to numpy,
all sorts of error messages start cropping up. What I've done to avoid this
is to change the statements `import numpy as np` in some (not all, just the
ones called while using autograd) files by

import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np

"""
import autograd.numpy as np
from autograd import grad, jacobian
from autograd.extend import primitive, defvjp

from quantum_CSD_compiler.UnitaryMat import UnitaryMat

import sys
print('np', 'np' in sys.modules)  # False
print('numpy', 'numpy' in sys.modules)  # True
print('autograd.numpy', 'autograd.numpy' in sys.modules)  # True


def d_u2(dwrt, *tlist):
    """
    tlist is a list of 4 floats, and dwrt, which stands for "derivative
    with respect to", is in range(4). This method returns the derivative
    of u2(*tlist) with respect to tlist[dwrt].

    The output of this method has been verified by comparing it to same
    derivatives calculated numerically with autograd.

    Parameters
    ----------
    dwrt : int
    tlist : list[float]

    Returns
    -------
    np.ndarray

    """
    assert dwrt in range(4)
    assert len(tlist) == 4
    if dwrt == 0:
        return 1j*UnitaryMat.u2_alt(*tlist)
    dwrt -= 1
    t = np.sqrt(tlist[1]**2 + tlist[2]**2 + tlist[3]**2)
    if abs(t) > 1e-6:
        tvec = np.array([tlist[1], tlist[2], tlist[3]])/t
    else:
        tvec = np.zeros((3,))
    dotted_vec = tvec*tvec[dwrt]*np.cos(t) +\
                 (np.sin(t)/t)*(-tvec*tvec[dwrt] + np.eye(3)[dwrt, :])
    out = -np.sin(t)*tvec[dwrt]*np.eye(2) +\
          1j*np.dot(UnitaryMat.sig_all(), dotted_vec)
    return np.exp(1j*tlist[0])*out


def d_auto_u2(dwrt, *tlist):
    """
    Returns the automatic (numerically computed by backprop) derivative of
    2-dim matrix UnitaryMat.u2_alt (UnitaryMat.u2_alt is an alternative to
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
        return np.real(UnitaryMat.u2_alt(*tlist1))

    def u2i(*tlist1):
        return np.imag(UnitaryMat.u2_alt(*tlist1))
    return jacobian(u2r, dwrt)(*tlist) + 1j*jacobian(u2i, dwrt)(*tlist)


@primitive
def pu2r(*tlist):
    """
    Returns real part of u2, and registers it as being primitive.

    Primitive means that its derivative will be provided in a defvjp (
    def of vector jacobian product) so no need for autograd to calculate it
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
    return np.real(UnitaryMat.u2_alt(*tlist))


@primitive
def pu2i(*tlist):
    """
    Returns imaginary part of u2, and registers it as being primitive.

    Primitive means that its derivative will be provided in a defvjp (
    def of vector jacobian product) so no need for autograd to calculate it
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
    return np.imag(UnitaryMat.u2_alt(*tlist))


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
    return pu2r(*tlist) + 1j*pu2i(*tlist)

defvjp(pu2r,
       # defines vector jacobian product of pu2r
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
       # defines vector jacobian product of pu2i
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
    def main():
        tlist = [.3, 1.1, .7, .5]
        for dwrt in range(4):
            print('err=', np.linalg.norm(
                d_auto_u2(dwrt, *tlist) - d_u2(dwrt, *tlist)))
        for dwrt in range(4):
            print('err=', np.linalg.norm(
                d_auto_pu2(dwrt, *tlist) - d_u2(dwrt, *tlist)))
    main()
