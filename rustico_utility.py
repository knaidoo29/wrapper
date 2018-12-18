import subprocess
import os
import numpy as np

def create_rustico_paramfile(path, root, ascii_file, boxsize, dk, output_path=None,
                             gridpower=6, particle2grid="CIC", triangles="ALL",
                             identity=1, useprint=True):
    if output_path is None:
        output_path = path
    if useprint is True:
        print "\n------------------------------------"
        print "** Writing RUSTICO Parameter File **"
        print "\n  Path = ", path
        print "  Paramfile = ", "rustico_paramfile_" + root + "_" + str(identity) + ".txt"
        print "  RUSTICO outputs:"
        print "  Path = ", output_path
        print "  Root = ", root + "_" + str(identity)
    subprocess.call('touch ' + path + "rustico_paramfile_" + root + "_" + str(identity) + ".txt", shell=True)
    rustico_paramfile = open(path + "rustico_paramfile_" + root + "_" + str(identity) + ".txt", 'w')
    rustico_paramfile.write("#Main parameters\n")
    rustico_paramfile.write("#Type of Box (periodic/cutsky): periodic\n")
    rustico_paramfile.write("#Type of file (ascii/gadget): ascii\n")
    rustico_paramfile.write("#Number of gadget files: 1\n")
    rustico_paramfile.write("#RSD distorsion on gadget periodic box (yes/no): no\n")
    rustico_paramfile.write("#Size of the Box (double/double): 0.0 +" + str(boxsize) + "\n")
    rustico_paramfile.write("#Type of Computation (DS/FFT): FFT\n")
    rustico_paramfile.write("#Binning for the Power Spectrum (linear/log10): linear\n")
    rustico_paramfile.write("#Size of the bin for the power spectrum (double): 0.01\n")
    rustico_paramfile.write("#k-range for computation (double/double): 0 1.\n")
    rustico_paramfile.write("\n")
    rustico_paramfile.write("#Bispectrum parameters\n")
    rustico_paramfile.write("#Do Bispectrum (yes/no): yes\n")
    rustico_paramfile.write("#Do Multigrid (yes/no): no\n")
    rustico_paramfile.write("#Triangle Shapes (ALL/EQU/ISO/SQU): " + triangles + "\n")
    rustico_paramfile.write("#Size of the bin for the bispectrum (double): " + str(dk) + "\n")
    rustico_paramfile.write("#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): FFT\n")
    rustico_paramfile.write("#Write triangles in each bin(yes/no): no\n")
    rustico_paramfile.write("#Path for triangles in each bin: ./power_spectra/triangles\n")
    rustico_paramfile.write("\n")
    rustico_paramfile.write("#Read inout parameters\n")
    rustico_paramfile.write("#Path of data: " + ascii_file + "\n")
    rustico_paramfile.write("#Path of randoms: nothing\n")
    if output_path[-1] == '/':
        output_path = output_path[:-1]
    rustico_paramfile.write("#Path of output: " + output_path + "\n")
    rustico_paramfile.write("#Identifier of output: " + root + "_" + str(identity) + "\n")
    rustico_paramfile.write("#Write header: yes\n")
    rustico_paramfile.write("\n")
    rustico_paramfile.write("#FFT parameters\n")
    rustico_paramfile.write("#Number of Grid Cells power (int): " + str(gridpower) + "\n")
    rustico_paramfile.write("#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): " + str(particle2grid) + "\n")
    rustico_paramfile.write("#Type of Yamamoto (GridCenter/GridAverage): GridCenter\n")
    rustico_paramfile.write("#Number of interlacing steps (int): 2\n")
    rustico_paramfile.write("#Do Grid Correction? (yes/no): yes\n")
    rustico_paramfile.write("\n")
    rustico_paramfile.write("#Cutsky parameters\n")
    rustico_paramfile.write("#Redshift Range (double/double): 0.8 2.2\n")
    rustico_paramfile.write("#Omega matter value (double): 0.31\n")
    rustico_paramfile.write("#Area effective value in deg^2 (double): 1001.2525\n")
    rustico_paramfile.write("#Hexadecapole as (L4/L2L2): L2L2\n")
    rustico_paramfile.write("#Compute Normalization as (area/density): density\n")
    rustico_paramfile.write("#Compute Normalization using (randoms/data): data\n")
    rustico_paramfile.write("#Compute Shot noise as (double): 1.0\n")
    rustico_paramfile.close()
    if useprint is True:
        print '\n  Parameter file written!'
        print "------------------------------------\n"
