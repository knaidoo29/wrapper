import numpy as np
import subprocess
import rockstar_temp


class Rockstar:

    def __init__(self):
        self.h0 = None
        self.omega_m = None
        self.omega_l = None
        self.boxsize = None
        self.nmesh = None
        self.path = None
        self.root = None
        self.parts = None
        self.parallel = None

    def setup(self, h0, omega_m, omega_l, boxsize, nmesh, path, root, parts, parallel):
        self.h0 = h0
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.boxsize = boxsize
        self.nmesh = nmesh
        self.path = path
        self.root = root
        self.parts = parts
        self.parallel = parallel
        rockstar_temp.create_rockstar_config_file(self.h0, self.omega_m, self.omega_l,
                                                  self.boxsize, self.nmesh, self.path,
                                                  self.root, self.parts, self.parallel)

    def run(self, location='macbook'):
        if location == 'macbook':
            rockstar_exec = './Users/krishna/Programs/rockstar/rockstar'
        subprocess.call(rockstar_exec + " -c " + self.path+'/rockstar_temp.cfg', shell=True)
