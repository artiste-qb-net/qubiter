import pennylane as qml


def Feynman(rads2, rads1):
    # distinct fun names in functional placeholders=
    # ['my_fun']
    def my_fun(x, y):
        return x + .5*y
    qml.RX(-0.8975979010256552, wires=2)
    qml.RX(rads2*.5*(-2), wires=1)
    def rot(rad_ang_x, rad_ang_y, rad_ang_z):
        """
        Returns
    
        exp(1j*(rad_ang_x*sig_x + rad_ang_y*sig_y + rad_ang_z*sig_z))
    
        where rad_ang_x is an angle in radians and sig_x is the x Pauli
        matrix, etc.
    
        Parameters
        ----------
        rad_ang_x : float
        rad_ang_y : float
        rad_ang_z : float
    
        Returns
        -------
        np.ndarray
    
        """
        ty = np.complex128
        mat = np.zeros([2, 2], dtype=ty)
        vec = np.array([rad_ang_x, rad_ang_y, rad_ang_z])
        n = np.linalg.norm(vec)  # sqrt(dot(vec, vec.conj))
        if abs(n) < 1e-8:
            mat[0, 0] = 1
            mat[1, 1] = 1
        else:
            nx = rad_ang_x/n
            ny = rad_ang_y/n
            nz = rad_ang_z/n
            c = np.cos(n)
            s = np.sin(n)
            mat[0, 0] = c + 1j*s*nz
            mat[0, 1] = s*ny + 1j*s*nx
            mat[1, 0] = -s*ny + 1j*s*nx
            mat[1, 1] = c - 1j*s*nz
        return mat
    qml.QubitUnitary(rot(rads1*(-2), -rads1*3*(-2), rads2*(-2)), wires=3)
    qml.RX(-my_fun(rads2, rads1)*(-2), wires=1)
    qml.CNOT(wires=[2, 3])
    return qml.expval.Hermitian(hamil)
