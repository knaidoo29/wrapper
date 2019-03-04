"""mg_picola_wrapper.py contains functions for managing the creation of required files and the subsequent execution of
MG-PICOLA"""

import numpy as np
import subprocess
import picola_utility
import utility as util
import source_path


class MGPicola:

    def __init__(self, n_grid=256, n_mesh=None, box_size=1000., seed=None, processors=2, initial_redshift=10,
                 f_r_0=1e-5, omega_m=0.3175, omega_baryon=0.049, h_0=0.6711, a_s=2.13e-9, n_s=0.9624, sum_m_v=0.,
                 location=None):
        """Defines a number of default set up parameters.

        Parameters
        ----------
        n_grid : int, optional
            The size of the N-body particle grid.
        n_mesh : int, optional
            The size of the fft grid.
        box_size : float, optional
            The size of the N-body box along one axis.
        seed : int, optional
            The seed for the initial conditions, if nothing is given a random number will be generated.
        processors : int, optional
            The number of processes used by the simulations.
        f_r_0 : float, optional
            The value of the f_R(R) at the present time.
        omega_m : float, optional
            The matter energy density.
        omega_baryon : float, optional
            The baryon energy density.
        h_0 : float, optional
            The reduced hubble constant at redshift 0.
        n_s : float, optional
            The primordial spectral index.
        sum_m_v : float, optional
            The sum of neutrino masses.
        location : str, optional
            Default location is splinter.
        """
        self.n_grid = n_grid
        if n_mesh is None:
            self.n_mesh = n_grid
        else:
            self.n_mesh = n_mesh
        self.box_size = box_size
        if seed is None:
            self.seed = np.random.randint(10000, size=1)[0]
        else:
            self.seed = seed
        self.processors = processors
        self.initial_redshift = initial_redshift
        self.f_r_0 = f_r_0
        self.omega_m = omega_m
        self.omega_lambda = 1. - self.omega_m
        self.omega_baryon = omega_baryon
        self.omega_neutrino = round(sum_m_v/(94.1*h_0**2.), 6)
        self.h_0 = h_0
        self.sigma_8 = None
        self.a_s = a_s
        self.n_s = n_s
        self.sum_m_v = sum_m_v
        if sum_m_v == 0.:
            self.omega = omega_m - self.omega_neutrino
            self.omega_cdm = omega_m - omega_baryon - self.omega_neutrino
        else:
            self.omega = omega_m
            self.omega_cdm = omega_m - omega_baryon
        self.path = None
        self.root = None
        if location is None:
            self.location = 'splinter'
        else:
            self.location = location
        camb_path, mpirun_path, mg_picola_exec = source_path.get_src(self.location)
        self.camb_path = camb_path
        self.mpirun_path = mpirun_path
        self.mg_picola_exec = mg_picola_exec

    def set_up(self, path, root, redshift, steps, get_pk=False, get_haloes=False):
        """Sets up the folder structure and required inputs.

        Parameters
        ----------
        path : str
            The path of the folder.
        root : str
            The name of the new folder to be created.
        redshift : list
            The desired output redshifts.
        steps : list
            Steps between each redshift output.
        get_pk : bool, optional
            determines whether to calculate the power spectra.
        get_haloes : bool, optional
            default : False, determines whether a halo catalogue is generated.
        """
        self.path = path
        self.root = root
        self.get_haloes = get_haloes
        util.create_folder(root, path=path + '/')
        util.create_folder('input/', path=path + '/' + root+'/')
        util.create_folder('paramfile/', path=path + '/' + root+'/')
        util.create_folder('output/', path=path + '/' + root+'/')
        picola_utility.create_camb_ini_file(self.path + '/' + self.root + '/input', self.h_0, self.omega_baryon,
                                            self.omega_cdm, self.omega_lambda, self.omega_neutrino, self.a_s, self.n_s)
        picola_utility.mg_picola_camb(self.path + '/' + self.root + '/input', self.camb_path)
        camb_line = 'at z =  0.000 sigma8 (all matter) = '
        with open(path + '/' + root + '/input/temp/camb.log') as search:
            for line in search:
                line = line.rstrip()  # remove '\n' at end of line
                if  camb_line == line[:len(camb_line)]:
                    self.sigma_8 = float(line[len(camb_line):])
        picola_utility.mg_picola_redshifts_files(self.path + '/' + self.root, redshift, steps)
        picola_utility.mg_picola_paramfile(self.path + '/' + self.root, self.f_r_0, self.processors, self.n_grid,
                                           self.n_mesh, self.box_size, self.initial_redshift, self.seed, self.omega,
                                           self.omega_baryon, self.h_0, self.n_s, self.sum_m_v, self.sigma_8,
                                           get_pk=get_pk, get_haloes=self.get_haloes)

    def run(self, parallel_setup='machinefile'):
        paramfile_path = self.path + '/' + self.root + '/paramfile/paramfile.txt'
        if parallel_setup == 'machinefile':
            subprocess.call(self.mpirun_path + ' --machinefile $PBS_NODEFILE -x PATH -x LD_LIBRARY_PATH ' + self.mg_picola_exec + ' ' + paramfile_path, shell=True)
        elif parallel_setup == '-np':
            subprocess.call(self.mpirun_path + ' -np ' + str(self.processors) + ' ' + self.mg_picola_exec + ' ' + paramfile_path, shell=True)
