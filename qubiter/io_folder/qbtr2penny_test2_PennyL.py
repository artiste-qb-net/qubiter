import pennylane as qml


def Feynman(rads2, rads1):
    # distinct fun names in functional placeholders=
    # ['my_fun']
    def my_fun(x, y):
        return x + .5*y
    qml.RX(-0.8975979010256552, wires=2)
    qml.RX(rads2*.5*(-2), wires=1)
    def rot(rad_ang_x, rad_ang_y, rad_ang_z, lib='np'):
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
        lib : str
    
        Returns
        -------
        np.ndarray
    
        """
        if 'autograd.numpy' in sys.modules:
            tlist = [0., rad_ang_x, rad_ang_y, rad_ang_z]
            return pu2(*tlist)
    
        mat_dict = OneBitGates.const_dict(0)
        vec = [rad_ang_x, rad_ang_y, rad_ang_z]
        n = OneBitGates.get_fun(lib, 'sqrt')(vec[0]**2 + vec[1]**2 + vec[2]**2)
        if abs(n) < 1e-8:
            mat_dict['00'] = 1
            mat_dict['11'] = 1
        else:
            nx = rad_ang_x/n
            ny = rad_ang_y/n
            nz = rad_ang_z/n
            c = OneBitGates.get_fun(lib, 'cos')(n)
            s = OneBitGates.get_fun(lib, 'sin')(n)
            mat_dict['00'] = c + 1j*s*nz
            mat_dict['01'] = s*ny + 1j*s*nx
            mat_dict['10'] = -s*ny + 1j*s*nx
            mat_dict['11'] = c - 1j*s*nz
        return OneBitGates.get_mat(lib, mat_dict)
    qml.QubitUnitary(rot(rads1*(-2), -rads1*3*(-2), rads2*(-2)), wires=3)
    qml.RX(-my_fun(rads2, rads1)*(-2), wires=1)
    qml.CNOT(wires=[2, 3])
    return qml.expval.Hermitian(hamil)
