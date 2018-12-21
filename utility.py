import numpy as np
import subprocess
import os.path
from scipy.integrate import trapz, cumtrapz
import pygadgetreader as gad


def create_folder(root, path=None):
    """Creates a folder with the name 'root' either in the current folder if path is None or a specified path.

    Parameters
    ----------
    root : str
        The name of the created folder.
    path : str, optional
        The name of the path of the created folder.
    """
    if path is None:
        if os.path.isdir(root) is False:
            subprocess.call('mkdir ' + root, shell=True)
    else:
        if os.path.isdir(path+root) is False:
            subprocess.call('mkdir ' + path + root, shell=True)


def integrate(x, y, total=False, equal_spacing=False):
    """Finds the integral of y wrt to x using the trapezium rule.

    keyword arguments:
    INPUT:
    x -- axis over the which the integration is performed.
    y -- the function f(x) for which the integration is performed.
    total -- if True will only output the integral over the full range of x.
    equal_spacing -- if true then the spacings between x is assumed to be equal, and the integral is performed as a sum.
                     if false then the spacings can change and is performed using the trapezium rule via a fortran sub-
                     routine.
    INTERNAL:
    _dx -- the difference between each ascending x value of an array. Only valid for equal_spacing=True.
    OUTPUT:
    y_integrand -- integral of the function y.
    """
    if equal_spacing is False:
        if total is True:
            y_integrand = trapz(y, x=x)
        else:
            y_integrand = cumtrapz(y, x=x, initial=0)
    else:
        _dx = x[1] - x[0]
        if total is True:
            y_integrand = trapz(y, dx=_dx)
        else:
            y_integrand = cumtrapz(y, dx=_dx, initial=0)
    return y_integrand


def get_gadget_data(gadget_filename, particles='dm', info='pos'):
    """Outputs the desired gadget data.

    keyword arguments:
    INPUT:
    gadget_filename -- name of the gadget file.
    particles -- particles type: dark matter ('dm'), gas or other types.
    info -- positions or velocities.
    OUTPUT:
    x, y, z -- positions or velocities.
    """
    data = gad.readsnap(gadget_filename, info, particles)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    return x, y, z


def get_fortran_unformatted_data(fortran_filename, data_structure):
    """Outputs data from an unformatted fortran binary data file.

    Parameters
    ----------
    fortran_filename : str
        The filename of the fortran file.
    data_structure : str
        Data structure of the unformatted fortran file.

    Returns
    -------
    data : array_like
        The data from the fortran file.
    """
    f = open(fortran_filename, 'rb')
    dt = np.dtype(data_structure)
    data = np.fromfile(f, dtype=dt, count=-1)
    f.close()
    return data


def check_data_within_box(x, y, z, boxsize):
    check1 = np.where((x < 0.) | (x > boxsize))[0]
    check2 = np.where((y < 0.) | (y > boxsize))[0]
    check3 = np.where((z < 0.) | (z > boxsize))[0]
    while len(check1) != 0 or len(check2) != 0 or len(check3) != 0:
        condition = np.where(x <= 0.)[0]
        x[condition] += boxsize
        condition = np.where(y <= 0.)[0]
        y[condition] += boxsize
        condition = np.where(z <= 0.)[0]
        z[condition] += boxsize
        condition = np.where(x >= boxsize)[0]
        x[condition] -= boxsize
        condition = np.where(y >= boxsize)[0]
        y[condition] -= boxsize
        condition = np.where(z >= boxsize)[0]
        z[condition] -= boxsize
        check1 = np.where((x < 0.) | (x > boxsize))[0]
        check2 = np.where((y < 0.) | (y > boxsize))[0]
        check3 = np.where((z < 0.) | (z > boxsize))[0]
    return x, y, z
