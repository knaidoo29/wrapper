import numpy as np
import subprocess


def get_gadget_particles(root, parts):
    """Opens an ascii gadget particle file.

    Parameters
    ----------
    root : str
        Particle filename root.
    parts : int
        The particle file parts.

    Returns
    -------
    x, y, z : array_like
        3D coordinate positions.
    vx, vy, vz : array_like
        The euclidean components of the velocity of each particles.
    """
    for i in range(0, parts):
        filename = root + '.' + str(i)
        _data = np.loadtxt(filename, unpack=True, skiprows=1)
        if i == 0:
            data = _data
        else:
            data = np.concatenate([data.T, _data.T]).T
    x, y, z, vx, vy, vz = data[0,:], data[1,:], data[2,:], data[3,:], data[4,:], data[5,:]
    return x, y, z, vx, vy, vz


def get_particle_mass(omega_m, boxsize, nsample):
    """Calculates the particle mass in a simulation.

    Parameters
    ----------
    omega_m : float
        Matter density.
    boxsize : float
        Size of simulation box, in Mpc/h.
    nsample : int
        Size of particle grid.

    Returns
    -------
    part_mass : float
        Dark matter particle mass.
    """
    G = 6.67408e-11
    Mpc = 3.085678e22
    H0 = 100.*1000./Mpc
    rho_m0 = ((3.*H0**2.)/(8.*np.pi*G))*omega_m
    M_total = rho_m0*(boxsize*Mpc)**3.
    part_mass = M_total/(nsample**3.)
    part_mass /= 1.989e30
    return part_mass


def create_rockstar_config_file(h0, omega_m, omega_l, boxsize, nsample, nmesh, redshift,
                                path, root, parts, parallel, use_gadget=True, periodic=True,
                                min_halo_part=10):
    """Create the rockstar configuration file.

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
    """
    subprocess.call('touch '+path+'rockstar_temp.cfg', shell=True)
    rockstar_temp = open(path+'rockstar_temp.cfg', 'w')
    if use_gadget is True:
        rockstar_temp.write("FILE_FORMAT = \"GADGET2\" \n")
    else:
        rockstar_temp.write("FILE_FORMAT = \"ASCII\" \n")
    rockstar_temp.write("PARTICLE_MASS = "+str(get_particle_mass(omega_m, boxsize, nsample))+"\n")
    rockstar_temp.write("SCALE_NOW = " + str(1./(1.+redshift)) + "\n")
    rockstar_temp.write("h0 = " + str(h0) + "\n")
    rockstar_temp.write("Ol = " + str(omega_l) + "\n")
    rockstar_temp.write("Om = " + str(omega_m) + "\n")
    if periodic is False:
        rockstar_temp.write("PERIODIC=0\n")
    else:
        pass
    if use_gadget is True:
        rockstar_temp.write("GADGET_LENGTH_CONVERSION = 1\n")
        rockstar_temp.write("GADGET_MASS_CONVERSION = 1e+10\n")
    rockstar_temp.write("BOX_SIZE = " + str(boxsize) + "\n")
    rockstar_temp.write("FORCE_RES = " + str(boxsize/(2.*nmesh)) + "\n")
    rockstar_temp.write("TOTAL_PARTICLES = " + str(nsample**3) + "\n")
    rockstar_temp.write("IGNORE_PARTICLE_IDS = 1\n")
    rockstar_temp.write("PARALLEL_IO=1\n")
    if use_gadget is False:
        rockstar_temp.write("INBASE=\"" + path + "\"\n")
        rockstar_temp.write("FILENAME=\"" + root + "\"\n")
    else:
        rockstar_temp.write("INBASE=\"" + path + "\"\n")
        rockstar_temp.write("FILENAME=\"" + root + ".<block>\"\n")
        rockstar_temp.write("NUM_SNAPS=1\n")
        rockstar_temp.write("NUM_BLOCKS=" + str(parts) + "\n")
    rockstar_temp.write("OUTBASE=\"" + path + "\"\n")
    rockstar_temp.write("MIN_HALO_OUTPUT_SIZE = "+str(min_halo_part)+"\n")
    rockstar_temp.write("NUM_WRITERS=" + str(parallel) + "\n")
    rockstar_temp.write("FORK_READERS_FROM_WRITERS = 1\n")
    rockstar_temp.write("FORK_PROCESSORS_PER_MACHINE =" + str(parallel) + "\n")
    rockstar_temp.close()
