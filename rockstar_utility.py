import numpy as np
import subprocess

def create_rockstar_config_file(h0, omega_m, omega_l, boxsize, nmesh, path, root, parts, parallel,
                                periodic=True):
    subprocess.call('touch '+path+'/rockstar_temp.cfg', shell=True)
    rockstar_temp = open(path+'/rockstar_temp.cfg', 'w')
    rockstar_temp.write("FILE_FORMAT = \"GADGET2\" \n")
    rockstar_temp.write("PARTICLE_MASS = 0\n")
    rockstar_temp.write("SCALE_NOW = 1\n")
    rockstar_temp.write("h0 = " + str(h0) + "\n")
    rockstar_temp.write("Ol = " + str(omega_l) + "\n")
    rockstar_temp.write("Om = " + str(omega_m) + "\n")
    rockstar_temp.write("GADGET_LENGTH_CONVERSION = 1\n")
    rockstar_temp.write("GADGET_MASS_CONVERSION = 1e+10\n")
    if periodic is False:
        rockstar_temp.write("PERIODIC=0")
    else:
        pass
    rockstar_temp.write("FORCE_RES = " + str(boxsize/(2.*nmesh)) + "\n")
    rockstar_temp.write("PARALLEL_IO=1\n")
    rockstar_temp.write("INBASE=\"" + path + "\"\n")
    rockstar_temp.write("FILENAME=\"" + root + ".<block>\"\n")
    rockstar_temp.write("NUM_SNAPS=1\n")
    rockstar_temp.write("NUM_BLOCKS=" + str(parts) + "\n")
    rockstar_temp.write("NUM_WRITERS=" + str(parallel) + "\n")
    rockstar_temp.close()
