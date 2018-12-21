import subprocess
import os
import numpy as np
import rustico_utility


class RUSTICO:

    def __init__(self, x, y, z, useprint=True):
        self.x = x
        self.y = y
        self.z = z
        self.path = None
        self.root = None
        self.use_splinter = False
        self.which_triangles = "ALL"
        self.triangle_types = np.array(["ALL", "EQU", "ISO", "SQU"])
        self.useprint = useprint
        self.identifier = 1
        self.ascii_filename = "temp"
        self.boxsize = None
        self.ngrid = None
        self.kf = None
        self.kn = None
        self.dk = None
        self.gridpower = None
        self.particle2grid_types = np.array(["NGC", "CIC", "TSC", "PCS"])
        self.particle2grid = "CIC"
        self.nprocessors = 1

    def setup_path(self, path=None, root=None, processors=None, use_splinter=False, identifier=None):
        if path is None:
            path = os.getcwd() + '/'
        else:
            if path[-1] != '/':
                path += '/'
        if root is None:
            root = "temp"
        if identifier is not None:
            self.identifier = identifier
        self.path = path
        self.root = root
        if processors is not None:
            self.nprocessors = processors
        self.use_splinter = use_splinter
        if self.useprint is True:
            print "\n**** SETTING PATHS ****"
            print "======================="
            print "\nPath = ", self.path
            print "Root = ", self.root
            print "Num Processors = ", self.nprocessors
            print "Splinter = ", self.use_splinter
            print "ASCII File = ", self.ascii_filename + "_" + str(self.identifier)

    def setup_param(self, boxsize, dk, gridpower=6, particle2grid="CIC", triangles="ALL"):
        if self.useprint is True:
            print "\n**** SETTING PARAMS ****"
            print "========================"
            print "\nAre paths defined?"
        if self.path is None or self.root is None:
            if self.useprint is True:
                print "  No"
                print "."*60 + "\n"
            self.setup_path()
            if self.useprint is True:
                print "\n" + "."*60
        else:
            if self.useprint is True:
                print "  Yes"
        if self.useprint is True:
                print "\nWhich Triangles?"
        if triangles != self.which_triangles:
            check = np.where(triangles == self.triangle_types)[0]
            if len(check) != 0:
                self.which_triangles = triangles
            else:
                if self.useprint is True:
                    print "  ERROR: '" + triangles + "' is not supported."
        if self.useprint is True:
            print "  Using '" + self.which_triangles + "' triangles"
        self.gridpower = gridpower
        self.boxsize = boxsize
        self.ngrid = 2.**gridpower
        if self.useprint is True:
            print "\nFFT info"
            print "  Boxsize = " + str(self.boxsize)
            print "  FFT Grid = " + str(self.ngrid)
        # fundamental_frequency
        self.kf = 2 * np.pi / self.boxsize
        # nyquist_frequency
        self.kn = np.pi / (self.boxsize/self.ngrid)
        self.dk = dk
        if self.useprint is True:
            print "\nK Range?"
            print "  Fundamental Frequency, kf = " + str(self.kf)
            print "  Nyquist Frequency, kn = " + str(self.kn)
            print "  dk = " + str(self.dk)
            print "  Num k bins per axis = " + str((self.kn - self.kf)/self.dk)
        if self.useprint is True:
                print "\nWhich particle to density assignment?"
        if particle2grid != self.particle2grid:
            check = np.where(particle2grid == self.particle2grid_types)[0]
            if len(check) != 0:
                self.which_particle2grid = triangles
            else:
                if self.useprint is True:
                    print "  ERROR: '" + particle2grid + "' is not supported."
        if self.useprint is True:
            print "  Using '" + self.particle2grid + "' to assign particles to grid."

    def run(self):
        # Saving temp file with data to be read by RUSTICO
        if self.useprint is True:
            print "\n**** RUNNING ****"
            print "================="
            print "\nSaving temporary ascii file to be read by RUSTICO."
        self.x, self.y, self.z = utility.check_data_within_box(self.x, self.y, self.z, self.boxsize)
        np.savetxt(self.path + self.ascii_filename + '_' + str(self.identifier) + '.txt', zip(self.x, self.y, self.z, np.ones(len(self.x))))
        if self.useprint is True:
            print "  Temporary ascii file = " + self.path + self.ascii_filename + '_' + str(self.identifier) + '.txt'
        # Create parameter file
        rustico_utility.create_rustico_paramfile(self.path, self.root,
                                                 self.path + self.ascii_filename + '_' + str(self.identifier) + '.txt',
                                                 self.boxsize, self.dk, gridpower=self.gridpower, particle2grid=self.particle2grid,
                                                 triangles=self.which_triangles, identity=self.identifier, useprint=self.useprint)
        # splinter or not
        if self.useprint is True:
            print "Finding RUSTICO executable"
        if self.use_splinter is False:
            rustico_executable_path = "/Users/krishna/Programs/Rustico-master"
            rustico_executable = "/Users/krishna/Programs/Rustico-master/rustico.o"
        else:
            rustico_executable_path = "/share/data1/knaidoo/splinter_libraries/build/rustico"
            rustico_executable = "/share/data1/knaidoo/splinter_libraries/build/rustico/rustico.o"
        # run rustico
        if self.useprint is True:
            print "\nSetting up number of processors to use = " + str(self.nprocessors)
        subprocess.call("export OMP_NUM_THREADS=" + str(self.nprocessors), shell=True)
        if self.useprint is True:
            print "\nRunning RUSTICO...\n"
            print "-"*50 + "\n"
        current_path = os.getcwd()
        subprocess.call("cd " + rustico_executable_path, shell=True)
        subprocess.call(rustico_executable + " " + self.path + "rustico_paramfile_" + self.root + "_" + str(self.identifier) + ".txt", shell=True)
        subprocess.call("cd " + current_path, shell=True)
        if self.useprint is True:
            print "\n" + "-" * 22 + " DONE " + "-" * 22 + "\n"
            #print "DONE!"

    def get_bk(self):
        bispectrum_file = self.path + 'Bispectrum_' + self.root + '_' + str(self.identifier) + '.txt'
        data = np.loadtxt(bispectrum_file, unpack=True)
        k1_eff = data[1]
        k2_eff = data[3]
        k3_eff = data[5]
        bk = data[6]
        return k1_eff, k2_eff, k3_eff, bk

    def clean(self, deep_clean=False):
        # delete parameter file and rustico data
        subprocess.call('rm ' + self.path + "rustico_paramfile_" + self.root + "_" + str(self.identifier) + ".txt", shell=True)
        subprocess.call('rm ' + self.path + self.ascii_filename + '_' + str(self.identifier) + '.txt', shell=True)
        if deep_clean is True:
            power_spectrum_file = self.path + 'Power_Spectrum_' + self.root + '_' + str(self.identifier) + '.txt'
            bispectrum_file = self.path + 'Bispectrum_' + self.root + '_' + str(self.identifier) + '.txt'
            temp_grid = os.getcwd() + '/temp' + str(int(2**self.gridpower)) + '.dat'
            subprocess.call('rm ' + power_spectrum_file, shell=True)
            subprocess.call('rm ' + bispectrum_file, shell=True)
            subprocess.call('rm ' + temp_grid, shell=True)
        # reset values
        self.x = None
        self.y = None
        self.z = None
        self.path = None
        self.root = None
        self.use_splinter = False
        self.which_triangles = "ALL"
        self.useprint = True
        self.identifier = 1
        self.ascii_filename = "temp"
        self.boxsize = None
        self.ngrid = None
        self.kf = None
        self.kn = None
        self.dk = None
        self.gridpower = None
        self.particle2grid = "CIC"
        self.nprocessors = 1
