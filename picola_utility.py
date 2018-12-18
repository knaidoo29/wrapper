"""camb_wrapper.py contains functions for constructing camb outputs."""

import numpy as np
import subprocess
import glob
import utility as util


def create_camb_ini_file(path, h_0, omega_baryon, omega_cdm, omega_lambda, omega_neutrino, a_s, n_s):
    """Creates a CAMB paramater file.

    Parameters
    ----------
    path : str
        Path of camb parameter file.
    h_0 : float
        Reduced hubble parameter.
    omega_baryon : float
        Baryon energy density.
    omega_cdm : float
        Cold dark matter energy desity.
    omega_lambda : float
        Dark energy energy density.
    omega_neutrino : float
        Neutrino energy density.
    a_s : float
        Amplitude of scalar fluctuations.
    n_s : float
        Primordial scalar index.

    Returns
    -------
    Writes a CAMB ini file.
    """
    subprocess.call('touch '+path+'/camb_temp.ini', shell=True)
    camb_ini = open(path+'/camb_temp.ini', 'w')
    camb_ini.write('get_scalar_cls     = T\n')
    camb_ini.write('get_vector_cls     = F\n')
    camb_ini.write('get_tensor_cls     = F\n')
    camb_ini.write('get_transfer       = T\n')
    camb_ini.write('do_lensing         = T\n')
    camb_ini.write('do_nonlinear       = 0\n')
    camb_ini.write('l_max_scalar       = 2200\n')
    camb_ini.write('l_max_tensor       = 1500\n')
    camb_ini.write('k_eta_max_tensor   = 3000\n')
    camb_ini.write('use_physical         = F\n')
    camb_ini.write('hubble               = ' + str(round(100.*h_0, 2)) + '\n')
    camb_ini.write('w                    = -1\n')
    camb_ini.write('omega_baryon         = ' + str(round(omega_baryon, 3)) + '\n')
    camb_ini.write('omega_cdm            = ' + str(round(omega_cdm, 4)) + '\n')
    camb_ini.write('omega_lambda         = ' + str(round(omega_lambda, 4)) + '\n')
    camb_ini.write('omega_neutrino       = ' + str(round(omega_neutrino, 6))+ '\n')
    camb_ini.write('massive_nu_approx    = 1\n')
    camb_ini.write('massless_neutrinos   = 2.046\n')
    camb_ini.write('massive_neutrinos    = 1\n')
    camb_ini.write('nu_mass_eigenstates  = 1\n')
    camb_ini.write('share_delta_neff     = T\n')
    camb_ini.write('nu_mass_fractions    = 1\n')
    camb_ini.write('nu_mass_degeneracies = \n')
    camb_ini.write('transfer_power_var   = 8\n')
    camb_ini.write('transfer_high_precision     = T\n')
    camb_ini.write('transfer_kmax               = 100\n')
    camb_ini.write('transfer_k_per_logint       = 5\n')
    camb_ini.write('transfer_interp_matterpower = T\n')
    camb_ini.write('initial_power_num         = 1\n')
    camb_ini.write('pivot_scalar              = 0.05\n')
    camb_ini.write('pivot_tensor              = 0.05\n')
    camb_ini.write('scalar_amp(1)             = ' + str(round(a_s*10.**9., 4)) + 'e-9\n')
    camb_ini.write('scalar_spectral_index(1)  = ' + str(round(n_s, 4)) + '\n')
    camb_ini.write('scalar_nrun(1)            = 0\n')
    camb_ini.write('scalar_nrunrun(1)         = 0\n')
    camb_ini.write('tensor_spectral_index(1)  = 0\n')
    camb_ini.write('tensor_nrun(1)            = 0\n')
    camb_ini.write('tensor_parameterization   = 1\n')
    camb_ini.write('initial_ratio(1)          = 1\n')
    camb_ini.write('cs2_lam         = 1\n')
    camb_ini.write('temp_cmb        = 2.7255\n')
    camb_ini.write('helium_fraction = 0.24\n')
    camb_ini.write('reionization             = T\n')
    camb_ini.write('re_use_optical_depth     = T\n')
    camb_ini.write('re_optical_depth         = 0.09\n')
    camb_ini.write('re_redshift              = 11\n')
    camb_ini.write('re_delta_redshift        = 1.5\n')
    camb_ini.write('re_ionization_frac       = -1\n')
    camb_ini.write('re_helium_redshift       = 3.5\n')
    camb_ini.write('re_helium_delta_redshift = 0.5\n')
    camb_ini.write('RECFAST_fudge            = 1.14\n')
    camb_ini.write('RECFAST_fudge_He         = 0.86\n')
    camb_ini.write('RECFAST_Heswitch         = 6\n')
    camb_ini.write('RECFAST_Hswitch          = T\n')
    camb_ini.write('initial_condition = 1\n')
    camb_ini.write('initial_vector    = -1 0 0 0 0\n')
    camb_ini.write('vector_mode       = 0\n')
    camb_ini.write('COBE_normalize    = F\n')
    camb_ini.write('CMB_outputscale   = 7.42835025e12\n')
    camb_ini.write('scalar_output_file         = scalCls.dat\n')
    camb_ini.write('vector_output_file         = vecCls.dat\n')
    camb_ini.write('tensor_output_file         = tensCls.dat\n')
    camb_ini.write('total_output_file          = totCls.dat\n')
    camb_ini.write('lensed_output_file         = lensedCls.dat\n')
    camb_ini.write('lensed_total_output_file   = lensedtotCls.dat\n')
    camb_ini.write('lens_potential_output_file = lenspotentialCls.dat\n')
    camb_ini.write('FITS_filename              = scalCls.fits\n')
    camb_ini.write('do_lensing_bispectrum         = F\n')
    camb_ini.write('do_primordial_bispectrum      = F\n')
    camb_ini.write('bispectrum_nfields            = 1\n')
    camb_ini.write('bispectrum_slice_base_L       = 0\n')
    camb_ini.write('bispectrum_ndelta             = 3\n')
    camb_ini.write('bispectrum_delta(1)           = 0\n')
    camb_ini.write('bispectrum_delta(2)           = 2\n')
    camb_ini.write('bispectrum_delta(3)           = 4\n')
    camb_ini.write('bispectrum_do_fisher          = F\n')
    camb_ini.write('bispectrum_fisher_noise       = 0\n')
    camb_ini.write('bispectrum_fisher_noise_pol   = 0\n')
    camb_ini.write('bispectrum_fisher_fwhm_arcmin = 7\n')
    camb_ini.write('bispectrum_full_output_file   =\n')
    camb_ini.write('bispectrum_full_output_sparse = F\n')
    camb_ini.write('bispectrum_export_alpha_beta  = F\n')
    camb_ini.write('feedback_level      = 1\n')
    camb_ini.write('output_file_headers = T\n')
    camb_ini.write('derived_parameters  = T\n')
    camb_ini.write('lensing_method      = 1\n')
    camb_ini.write('accurate_BB         = F\n')
    camb_ini.write('accurate_polarization   = T\n')
    camb_ini.write('accurate_reionization   = T\n')
    camb_ini.write('do_tensor_neutrinos     = T\n')
    camb_ini.write('do_late_rad_truncation  = T\n')
    camb_ini.write('halofit_version         =\n')
    camb_ini.write('number_of_threads       = 0\n')
    camb_ini.write('high_accuracy_default   = T\n')
    camb_ini.write('accuracy_boost          = 1\n')
    camb_ini.write('l_accuracy_boost        = 1\n')
    camb_ini.write('l_sample_boost          = 1\n')
    camb_ini.close()
    print 'CAMB parameter file written.'


def mg_picola_camb(path, camb_path):
    """Creates an MG-PICOLA CAMB generation file.

    Parameters
    ----------
    path : str
        Path of camb parameter file.
    camb_path : str
        Location of camb executable.

    Returns
    -------
    Writes a generating script to get CAMB outputs.
    """
    subprocess.call('touch '+path+'/gen_camb.sh', shell=True)
    gen_camb_file = open(path+'/gen_camb.sh', 'w')
    gen_camb_file.write('#!/bin/bash\n')
    gen_camb_file.write('suffix="temp"\n')
    gen_camb_file.write('output_folder="' + path + '/${suffix}"\n')
    gen_camb_file.write('camb_template="' + camb_path + 'HighLExtrapTemplate_lenspotentialCls.dat"\n')
    gen_camb_file.write('camb_executable="' + camb_path + 'camb"\n')
    gen_camb_file.write('camb_parameter_file="' + path + '/camb_${suffix}.ini"\n')
    gen_camb_file.write('camb_output_root="${suffix}"\n')
    gen_camb_file.write('camb_input_filename="camb_${camb_output_root}.ini"\n')
    gen_camb_file.write('picola_output_filename="picola_transfer_info_${suffix}.txt"\n')
    gen_camb_file.write('number_redshifts="50"\n')
    gen_camb_file.write('z_min="0.0"\n')
    gen_camb_file.write('z_max="50.0"\n')
    gen_camb_file.write('if [[ ! -e $output_folder ]]; then\n')
    gen_camb_file.write('  mkdir $output_folder\n')
    gen_camb_file.write('fi\n')
    gen_camb_file.write('cd $output_folder\n')
    gen_camb_file.write('if [[ -e $camb_executable ]]; then\n')
    gen_camb_file.write('  cp $camb_executable camb\n')
    gen_camb_file.write('else\n')
    gen_camb_file.write('  echo "Error: cannot find the camb executable [$camb_template]"\n')
    gen_camb_file.write('  exit\n')
    gen_camb_file.write('fi\n')
    gen_camb_file.write('if [[ -e $camb_template ]]; then\n')
    gen_camb_file.write('  cp $camb_template HighLExtrapTemplate_lenspotentialCls.dat\n')
    gen_camb_file.write('else\n')
    gen_camb_file.write('  echo "Error: cannot find the camb highL template [$camb_template]"\n')
    gen_camb_file.write('  exit\n')
    gen_camb_file.write('fi\n')
    gen_camb_file.write('if [[ -e $camb_parameter_file ]]; then\n')
    gen_camb_file.write('  camb_param="$(cat $camb_parameter_file)"\n')
    gen_camb_file.write('else\n')
    gen_camb_file.write('  echo "Error: cannot find the camb parameter file [$camb_parameter_file]"\n')
    gen_camb_file.write('  exit\n')
    gen_camb_file.write('fi\n')
    gen_camb_file.write('tmp_file="output_root              = $camb_output_root\\ntransfer_num_redshifts   = $number_redshifts"\n')
    gen_camb_file.write('picola_file="$(printf "$output_folder $number_redshifts")"\n')
    gen_camb_file.write('for i in $(seq 1 $number_redshifts)\n')
    gen_camb_file.write('do\n')
    gen_camb_file.write("""  cur_z=$(echo $z_min $z_max $i $number_redshifts | awk '{printf "%8.3f\\n", exp( log($1+1) + (log($2+1) - log($1+1)) * ($4 - $3) / ($4-1.0) ) - 1.0 }' )\n""")
    gen_camb_file.write("""  cur_z=$(echo $cur_z | awk '{$1=$1;print}')\n""")
    gen_camb_file.write('  tmp_file=$(printf "$tmp_file\ntransfer_redshift($i)    = $cur_z")\n')
    gen_camb_file.write('  tmp_file=$(printf "$tmp_file\ntransfer_filename($i)    = transfer_z${cur_z}.dat")\n')
    gen_camb_file.write('  tmp_file=$(printf "$tmp_file\ntransfer_matterpower($i) = matterpower_z${cur_z}.dat")\n')
    gen_camb_file.write("""  cur_z=$(echo $z_min $z_max $i $number_redshifts | awk '{printf "%8.3f\\n", exp( log($1+1) + (log($2+1) - log($1+1)) * ($3-1.0) / ($4-1.0) ) - 1.0 }' )\n""")
    gen_camb_file.write("""  cur_z=$(echo $cur_z | awk '{$1=$1;print}')\n""")
    gen_camb_file.write('  picola_file=$(printf "$picola_file\n${camb_output_root}_transfer_z${cur_z}.dat  $cur_z")\n')
    gen_camb_file.write('done\n')
    gen_camb_file.write('printf "$picola_file" > $picola_output_filename\n')
    gen_camb_file.write('printf "$camb_param\n $tmp_file" > $camb_input_filename\n')
    gen_camb_file.write('./camb $camb_input_filename > camb.log\n')
    gen_camb_file.write('rm ' + path + '/temp/camb\n')
    gen_camb_file.write('rm ' + path + '/temp/HighLExtrapTemplate_lenspotentialCls.dat\n')
    gen_camb_file.close()
    subprocess.call('bash '+ path + '/gen_camb.sh', shell=True)


def mg_picola_redshifts_files(path, redshift, steps):
    """Constructs the redshift input file.

    Parameters
    ----------
    path : str
        Folder path.
    redshift : list
        Output redshift list.
    steps : list
        Steps needed for each step.

    Returns
    -------
    Returns a redshift output file for MG-PICOLA.
    """
    subprocess.call('touch ' + path + '/input/output_redshifts.dat', shell=True)
    redshift_file = open(path + '/input/output_redshifts.dat', 'w')
    for i in range(0, len(redshift)):
        redshift_file.write(str(round(redshift[i], 2)) + ', ' + str(steps[i]) + '\n')
    redshift_file.close()


def mg_picola_paramfile(path, f_r_0, processors, n_grid, n_mesh, box_size, initial_redshift, seed, omega, omega_baryon, h_0,
                        n_s, sum_m_v, sigma_8, get_pk=False, get_haloes=False):
    """Constructs the MG-PICOLA parameter file.

    Parameters
    ----------
    path : str
        Folder path.
    f_r_0 : float
        The F_r constant at z=0.
    processors : int
        Number of processors for the function to be run in.
    n_grid : int
        The grid size of the simulation.
    n_mesh : int
        The fft grid size.
    box_size : float
        The size of the box of the simulations.
    initial_redshift : float
        The initial redshift of the simulation.
    seed : int
        The seed for creating the initial conditions of the simulations.
    omega : float
        The energy density of dark matter and baryons.
    omega_baryon : float
        The energy density of baryons.
    h_0 : float
        The value of hubble constant at redshift z=0.
    n_s : float
        The spectral index value.
    sum_m_v : float
        The sum of neutrino masses.
    sigma_8 : float
        sigma_8.
    get_pk : bool, optional
        Ouputs the power spectra of the data.
    get_haloes : bool, optional
        default : False, determines whether a halo catalogue is generated.

    Returns
    -------
    Outputs the parameter file for MG-PICOLA.
    """
    subprocess.call('touch '+path+'/paramfile/paramfile.txt', shell=True)
    paramfile = open(path+'/paramfile/paramfile.txt', 'w')
    paramfile.write('% ====================================== %\n')
    paramfile.write('% = Generated MG-PICOLA parameter file = %\n')
    paramfile.write('% ====================================== %\n')
    paramfile.write('\n')
    paramfile.write('% Modified gravity parameters\n')
    paramfile.write('modified_gravity_active  0\n')
    paramfile.write('fofr0                    ' + str(round(f_r_0*10.**5., 2)) + 'e-5\n')
    paramfile.write('nfofr                    1.0\n')
    paramfile.write('include_screening        1\n')
    paramfile.write('use_lcdm_growth_factors  0\n')
    paramfile.write('input_pofk_is_for_lcdm   1\n')
    paramfile.write('input_sigma8_is_for_lcdm 1\n')
    paramfile.write('inverted_initial_condition 0\n')
    paramfile.write('amplitude_fixed_initial_condiiton 0\n')
    paramfile.write('\n')
    paramfile.write('% Simulation outputs\n')
    paramfile.write('OutputDir                   ' + path + '/output\n')
    paramfile.write('FileBase                    particles\n')
    paramfile.write('OutputRedshiftFile          ' + path + '/input/output_redshifts.dat\n')
    paramfile.write('NumFilesWrittenInParallel   1\n')
    paramfile.write('\n')
    paramfile.write('% Simulation specifications\n')
    paramfile.write('UseCOLA          1\n')
    paramfile.write('Buffer           1.5\n')
    paramfile.write('Nmesh            ' + str(n_mesh) + '\n')
    paramfile.write('Nsample          ' + str(n_grid) + '\n')
    paramfile.write('Box              ' + str(box_size) + '\n')
    paramfile.write('Init_Redshift    ' + str(initial_redshift) + '\n')
    paramfile.write('Seed             ' + str(seed) + '\n')
    paramfile.write('SphereMode       0\n')
    paramfile.write('\n')
    paramfile.write('% Initial conditions\n')
    paramfile.write('WhichSpectrum    1\n')
    paramfile.write('WhichTransfer    0\n')
    paramfile.write('FileWithInputSpectrum -\n')
    paramfile.write('FileWithInputTransfer -\n')
    paramfile.write('\n')
    paramfile.write('% Parameters for massive neutrinos\n')
    paramfile.write('nu_FilenameTransferInfofile  ' + path + '/input/temp/picola_transfer_info_temp.txt\n')
    paramfile.write('nu_include_massive_neutrinos  1\n')
    paramfile.write('nu_SumMassNuEV                ' + str(round(sum_m_v, 4)) + '\n')
    paramfile.write('\n')
    paramfile.write('% Cosmological parameters\n')
    paramfile.write('Omega            ' + str(round(omega, 4)) + '\n')
    paramfile.write('OmegaBaryon      ' + str(round(omega_baryon, 3)) + '\n')
    paramfile.write('HubbleParam      ' + str(round(h_0, 4)) + '\n')
    paramfile.write('Sigma8           ' + str(round(sigma_8, 4)) + '\n')
    paramfile.write('PrimordialIndex  ' + str(round(n_s, 4)) + '\n')
    paramfile.write('\n')
    paramfile.write('% Units\n')
    paramfile.write('UnitLength_in_cm                3.085678e24\n')
    paramfile.write('UnitMass_in_g                   1.989e43\n')
    paramfile.write('UnitVelocity_in_cm_per_s        1e5\n')
    paramfile.write('InputSpectrum_UnitLength_in_cm  3.085678e24\n')
    paramfile.write('\n')
    paramfile.write('% Power-spectrum calculation\n')
    if get_pk is True:
        paramfile.write('pofk_compute_every_step   1\n')
    else:
        paramfile.write('pofk_compute_every_step   0\n')
    paramfile.write('pofk_compute_rsd_pofk     2\n')
    paramfile.write('pofk_nbins                0\n')
    paramfile.write('pofk_bintype              1\n')
    paramfile.write('pofk_subtract_shotnoise   1\n')
    paramfile.write('pofk_kmin                 0\n')
    paramfile.write('pofk_kmax                 0\n')
    paramfile.write('\n')
    paramfile.write('% Halofinding\n')
    if get_haloes is True:
        paramfile.write('mm_run_matchmaker 1\n')
    else:
        paramfile.write('mm_run_matchmaker 0\n')
    paramfile.write('mm_output_format  0\n')
    paramfile.write('mm_min_npart_halo 10\n')
    paramfile.write('mm_linking_length 0.2\n')
    paramfile.write('mm_dx_extra_mpc   10.0\n')
    paramfile.write('mm_output_pernode 1\n')
    paramfile.close()


def calculate_sigma_8(k, pk):
    """Calulates sigma_8 from an input power spectrum.

    Parameters
    ----------
    k, pk : array_like
        The linear power spectrum and corresponding k value.

    Returns
    -------
    sigma_8 : float
        The calculated value of sigma_8.
    """
    w = 3.*(np.sin(k*8.)-(k*8.*np.cos(k*8.)))/(k*8.)**3.
    integrand = k*k*pk*w*w
    sigma_8_2 = (1./(2.*np.pi**2.))*util.integrate(k, integrand, total=True)
    sigma_8 = sigma_8_2**0.5
    return sigma_8


def get_particle_filenames(path, root, redshift):
    """Gets the particle filenames.

    Parameters
    ----------
    path : str
        Path of simulation.
    root : str
        The root of the simulation run.
    redshift : array_like
        The redshift of the snapshots.

    Returns
    -------
    filename_index : int
        Integer index of the particle outputs matched to the redshift.
    particle_filename : str
        The particle output filenames.
    """
    # gets the particle file names
    path2output = path + '/' + root + '/output/'
    particle_file_full_names = (glob.glob(path2output + 'particles_z*.0'))
    particle_filename = []
    for i in range(0, len(particle_file_full_names)):
        particle_filename.append(particle_file_full_names[i][len(path2output):])
    # converts file names
    redshift_file2float = []
    for i in range(0, len(particle_filename)):
        particle_file = particle_filename[i][11:-2]
        redshift_file = float(particle_file[:particle_file.find('p')] + '.' + particle_file[particle_file.find('p')+1:])
        redshift_file2float.append(redshift_file)
    redshift_file2float = np.array(redshift_file2float)
    redshift_map2file = np.copy(redshift_file2float)
    for i in range(0, len(redshift)):
        condition = np.where(redshift_file2float == redshift[i])[0]
        if redshift[i] == 0. and len(condition) == 0:
            condition = np.where(redshift == -1.999)[0]
        elif len(condition) == 0:
            condition = np.argmin(abs(redshift_file2float - redshift[i]))
        else:
            pass
        redshift_map2file[i] = condition
    filename_index = np.ndarray.flatten(np.array(redshift_map2file)).astype('int')
    return filename_index, particle_filename
