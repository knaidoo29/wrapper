import numpy as np
import wrapper as wrap

path = '/share/data1/knaidoo/splinter_libraries/pylib/wrapper/test'
root = 'splinter_test'
redshift = [0., 0.5, 1., 2.]
steps = [5, 5, 5, 5]

mgp = wrap.MGPicola(n_grid=64, processors=16)
mgp.set_up(path, root, redshift, steps, get_pk=True, get_haloes=True)
mgp.run()#parallel_setup='-np')
