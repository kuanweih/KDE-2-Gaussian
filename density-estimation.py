import numpy as np


# files names
COORDFILE = 'stars-coord-fornax.npy'    # input stars coords
INFOFILE = 'stars-coord-attr.npy'    # input center info
SIGNI_FILE = 'significance'    # output significance file
MESHFILE = 'meshgrids'    # output mesh grids
PARAMETER_TXT = 'param-den.txt'    # parma file txt


# NUM_GRID, SIGMA1, SIGMA2 from parameter file
with open(PARAMETER_TXT, 'r') as param_file:
    exec(param_file.readline())
    exec(param_file.readline())
    exec(param_file.readline())


def gaussian(x, y, s):
    g = np.exp(- 0.5 * (x**2 + y**2) / s**2) / (2. * np.pi * s**2)
    return g


def star_density(x, y, star_x, star_y, s):
    od = np.sum(np.array([gaussian(x - star_x[i], y - star_y[i], s)
                          for i in range(len(star_x))]), axis=0)
    return od


def significance(x, y, s1, s2, star_x, star_y):
    ss = (star_density(x, y, star_x, star_y, s1) -
          star_density(x, y, star_x, star_y, s2))
    ss /= np.sqrt(star_density(x, y, star_x, star_y, s2))
    ss *= np.sqrt(4. * np.pi) * s1
    return ss


# load ra and dec
coords = np.load(COORDFILE)
infos = np.load(INFOFILE)

ra_center = infos[0]
dec_center = infos[1]
width_mesh = infos[2]


# create mesh
x = np.linspace(ra_center - 0.5 * width_mesh,
                ra_center + 0.5 * width_mesh, num=NUM_GRID, endpoint=True)
y = np.linspace(dec_center - 0.5 * width_mesh,
                dec_center + 0.5 * width_mesh, num=NUM_GRID, endpoint=True)
xx, yy = np.meshgrid(x, y, sparse=True)  # TODO: what does sparse mean?


ss = significance(xx, yy, SIGMA1, SIGMA2, coords[0], coords[1])
# print(ss)
np.save(SIGNI_FILE, ss)
np.save(MESHFILE, np.array([x, y]))


# end of code
