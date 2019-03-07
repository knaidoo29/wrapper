import os.path
import numpy as np
import subprocess
import rockstar_utility
import source_path


class Rockstar:
    """A class for running the Rockstar halo finder.

    Parameters
    ----------
    h0 : float
        Hubble constant / 100.
    omega_m : float
        Matter density.
    omega_l : float
        Lambda density.
    boxsize : float
        Size of simulation box, in Mpc/h.
    nsample : int
        Size of particle grid.
    nmesh : int
        size of the density grid.
    redshift : float
        Redshift of the simulation box.
    path : str
        Folder path of the output and input files.
    root : str
        Particle filename root.
    parts : int
        The particle file parts.
    parallel : int
        Number of cores to use.
    use_gadget : bool, optional
        To use gadget file type or not... currently doesn't work.
    periodic : bool, optional
        Set periodic conditions or not... currently doesn't work.
    min_halo_part : int, optional
        Minimum number of particles for a group considered to be called a halo.
    location : str, optional
        The location of the Rockstar executable, dependent on where this being run.
    """

    def __init__(self):
        self.h0 = None
        self.omega_m = None
        self.omega_l = None
        self.boxsize = None
        self.nsample = None
        self.nmesh = None
        self.redshift = None
        self.path = None
        self.root = None
        self.parts = None
        self.parallel = None
        self.periodic = None
        self.min_halo_part = None

    def setup(self, h0, omega_m, omega_l, boxsize, nsample, nmesh, redshift,
              path, root, parts, parallel, periodic=True, min_halo_part=10, use_gadget=False):
        self.h0 = h0
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.boxsize = boxsize
        self.nsample = nsample
        self.nmesh = nmesh
        self.redshift = redshift
        self.path = path
        self.root = root
        self.parts = parts
        self.parallel = parallel
        self.periodic = periodic
        self.min_halo_part = min_halo_part
        rockstar_utility.create_rockstar_config_file(self.h0, self.omega_m, self.omega_l, self.boxsize,
                                                     self.nsample, self.nmesh, self.redshift, self.path, self.root,
                                                     self.parts, self.parallel, use_gadget=use_gadget, periodic=self.periodic,
                                                     min_halo_part=self.min_halo_part)

    def run(self, location='home'):
        camb_path, mpirun_path, mg_picola_exec, rockstar_exec = source_path.get_src(location)
        #fname_exist = False
        #while fname_exist is False:
        subprocess.call(rockstar_exec + " -c " + self.path+'rockstar_temp.cfg &', shell=True)
        subprocess.call(rockstar_exec + " -c " + self.path+'auto-rockstar.cfg', shell=True)
        #fname = self.path + 'halos_0.0.ascii'
        #fname_exist = os.path.isfile(fname)
