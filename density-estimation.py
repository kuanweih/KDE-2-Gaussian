import numpy as np


# files names
COORDFILE = 'stars-coord.npy'    # input stars coords
INFOFILE = 'stars-coord-attr.npy'    # input center info
SIGNI_FILE = 'significance'    # output significance file
MESHFILE = 'meshgrids'    # output mesh grids
PARAMETER_TXT = 'param-den.txt'    # parma file txt

KERNEL_BG = 'poisson'


# NUM_GRID, SIGMA1, SIGMA2 from parameter file
# TODO update to the one to read all lines in param file as in get....py
with open(PARAMETER_TXT, 'r') as param_file:
    exec(param_file.readline())
    exec(param_file.readline())
    exec(param_file.readline())


def gaussian(x, y, s):
    """ 2d Gaussian with a width of s """
    g = np.exp(- 0.5 * (x**2 + y**2) / s**2) / (2. * np.pi * s**2)
    return g


def od_gaussian(x, y, star_x, star_y, s):
    """ overdensity on a mesh with a 2d Gaussian filter
        x, y: mesh arrays, star_x, star_y: position of stars """
    od = np.sum(np.array([gaussian(x - star_x[i], y - star_y[i], s)
                          for i in range(len(star_x))]), axis=0)
    return od


def distance2(x_arr, y_arr, x_cen, y_cen):
    """ 2d distance """
    d2 = (x_arr - x_cen)**2 + (y_arr - y_cen)**2
    return d2


# TODO Need to figure out what scales to feed in...
def od_poisson(x, y, star_x, star_y, s1, s2):
    """ overdensity on a mesh with a Poisson CDF
        x, y: mesh arrays, star_x, star_y: position of stars
        s1, s2: inner and outer scales """
    from scipy.stats import poisson
    n_inner = np.sum(np.array([(distance2(x, y, star_x[i], star_y[i]) < s1**2)
                               for i in range(len(star_x))]), axis=0)

    """ v1: number between s1**2 and s2**2
        a lot of zeros... """
    #
    # n_outer = np.sum(np.array([(s1**2 < distance2(x, y, star_x[i], star_y[i])) *
    #                            (distance2(x, y, star_x[i], star_y[i]) < s2**2)
    #                            for i in range(len(star_x))]), axis=0)

    """ v2: number between (s2**2-s1**2) and s2**2
        significance too large? """
    n_outer = np.sum(np.array([((s2**2 - s1**2) < distance2(x, y, star_x[i], star_y[i])) *
                               (distance2(x, y, star_x[i], star_y[i]) < s2**2)
                               for i in range(len(star_x))]), axis=0)
    od = poisson.cdf(n_inner, n_outer)
    return od


def significance(x, y, s1, s2, star_x, star_y, kernel_bg='gaussian'):
    """ get significance on a mesh with a Poisson CDF
        x, y: mesh arrays, star_x, star_y: position of stars
        s1, s2: target and background scales
        kernel_bg: background kernel: 'gaussian' or 'poisson' """
    if kernel_bg == 'gaussian':
        ss = (od_gaussian(x, y, star_x, star_y, s1) -
              od_gaussian(x, y, star_x, star_y, s2))
        ss /= np.sqrt(od_gaussian(x, y, star_x, star_y, s2))
    elif kernel_bg == 'poisson':
        ss = (od_gaussian(x, y, star_x, star_y, s1) -
              od_poisson(x, y, star_x, star_y, s1, s2))
        ss /= np.sqrt(od_poisson(x, y, star_x, star_y, s1, s2))

        # TODO test what output is from func: od_poisson
        # ss = od_poisson(x, y, star_x, star_y, s1, s2)

    else:
        print('wrong kernel :(')
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


# ss = significance(xx, yy, SIGMA1, SIGMA2, coords[0], coords[1])
ss = significance(xx, yy, SIGMA1, SIGMA2,
                  coords[0], coords[1], kernel_bg=KERNEL_BG)
print(ss)
np.save(SIGNI_FILE, ss)
np.save(MESHFILE, np.array([x, y]))


# end of code
