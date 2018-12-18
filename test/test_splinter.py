import numpy as np
import wrapper as wrap

path = '/share/data1/knaidoo/splinter_libraries/pylib/wrapper/test'
root = 'splinter_test'
redshift = [0.5]
steps = [20]

mgp = wrap.MGPicola(n_grid=256, n_mesh=3*256)
mgp.set_up(path, root, redshift, steps, get_pk=True, get_haloes=True)
mgp.run()
