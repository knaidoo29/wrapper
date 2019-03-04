
def get_src(location):
    if location == 'splinter':
        camb_path = '/share/data1/knaidoo/splinter_libraries/src/CAMB/'
        mpirun_path = '/opt/openmpi/bin/mpirun'
        mg_picola_exec = '/share/data1/knaidoo/splinter_libraries/src/MG-PICOLA/MG_PICOLA_FOFRNU'
    elif location == 'home':
        camb_path = '/Users/krishna/Programs/CAMB-0.1.7/'
        mpirun_path = '/Users/krishna/Programs/build/openmpi/bin/mpirun'
        mg_picola_exec = '/Users/krishna/Programs/MG-PICOLA-PUBLIC-master/MG_PICOLA_HALOES_FOFRNU'
    else:
        print('Source location is undefined.')
    return camb_path, mpirun_path, mg_picola_exec
