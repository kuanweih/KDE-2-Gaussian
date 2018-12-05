import numpy as np
from param_den import *


# def gaussian(x, y, s):
#     """
#     2d Gaussian with a width of s
#     """
#     g = np.exp(- 0.5 * (x**2 + y**2) / s**2) / (2. * np.pi * s**2)
#     return g


# def od_gaussian(x, y, star_x, star_y, s):
#     """
#     overdensity on a mesh with a 2d Gaussian filter
#     x, y: mesh arrays. star_x, star_y: position of stars
#     """
#     od = np.sum(np.array([gaussian(x - star_x[i], y - star_y[i], s)
#                           for i in range(len(star_x))]), axis=0)
#     return od


def distance2(x_arr, y_arr, x_cen, y_cen):
    """
    2d distance square
    """
    d2 = (x_arr - x_cen)**2 + (y_arr - y_cen)**2
    return d2


def sig_poisson(x, y, s1, s2, star_x, star_y, r12):
    """
    get z-score as significance using inverse survival function of Poisson.
    x, y: mesh arrays. star_x, star_y: position of stars
    s1, s2: inner and outer scales
    r12: ratio of area between target and background
    """
    r = np.sqrt(s2**2 + r12 * s1**2)    # outer radius
    n_inner = np.sum(np.array([(distance2(x, y, star_x[i], star_y[i]) < s1**2)
                               for i in range(len(star_x))]), axis=0)
    n_outer = np.sum(np.array([(s2**2 < distance2(x, y, star_x[i], star_y[i])) *
                               (distance2(x, y, star_x[i], star_y[i]) < r**2)
                               for i in range(len(star_x))]), axis=0)
    lambda_poisson = n_outer / r12    # estimated background count
    sig = (n_inner - lambda_poisson) / np.sqrt(lambda_poisson)    # z score
    return sig


def sig_2_gaussian(x, y, s1, s2, star_x, star_y):
    """
    get significance using 2 Gaussian kernels.
    x, y: mesh arrays
    star_x, star_y: position of stars
    s1, s2: target and background scales
    """
    from scipy.ndimage import gaussian_filter
    hist2d, x, y = np.histogram2d(star_x, star_y, bins=(x, y))
    od_1 = gaussian_filter(hist2d, sigma=s1)
    od_2 = gaussian_filter(hist2d, sigma=s2)
    # od_1 = od_gaussian(x, y, star_x, star_y, s1)
    # od_2 = od_gaussian(x, y, star_x, star_y, s2)
    sig = (od_1 - od_2) / np.sqrt(od_2 / (4. * np.pi * s1**2))
    return sig


def get_grid_coord(center, width_mesh):
    """
    get grid coordinates according to the center position and width of the mesh
    """
    coord = np.linspace(center - 0.5 * width_mesh,
                        center + 0.5 * width_mesh, num=NUM_GRID, endpoint=True)
    return coord


# def create_mesh(ra_center, dec_center, width_mesh):
#     """
#     create meshgrid according to grid coordinates by np.meshgrid
#     """
#     x = get_grid_coord(ra_center, width_mesh)
#     y = get_grid_coord(dec_center, width_mesh)
#     return np.meshgrid(x, y, sparse=True)


def main():
    """
    calculate kernel density estimation
    """
    # files names
    COORDFILE = 'stars-coord.npy'    # input stars coords
    INFOFILE = 'stars-coord-attr.npy'    # input center info
    SIGNI_FILE = 'significance'    # output significance file
    MESHFILE = 'meshgrids'    # output mesh grids

    # load ra and dec
    coords = np.load(COORDFILE)
    infos = np.load(INFOFILE)

    ra_center = infos[0]
    dec_center = infos[1]
    width_mesh = infos[2]  # This is actually the radius when querying

    # create mesh
    # xx, yy = create_mesh(ra_center, dec_center, width_mesh) #TODO del when done
    xx = get_grid_coord(ra_center, width_mesh)
    yy = get_grid_coord(dec_center, width_mesh)
    print('There %d grids on a side.' % NUM_GRID)
    print('Dectection scale is %0.4f degree' % SIGMA1)
    print('Background scale is %0.4f degree' % SIGMA2)

    # get significance
    if KERNEL_BG == 'gaussian':
        print('We are using 2-Gaussian kernels to estimate the density.')
        sig = sig_2_gaussian(xx, yy, SIGMA1, SIGMA2, coords[0], coords[1])
    elif KERNEL_BG == 'poisson':
        print('We are using Poisson statistics to estimate the density.')
        print('Background area = %0.1f detection area.' % RATIO_AREA_TG_BG)
        meshgrid = np.meshgrid(xx, yy, sparse=True)  # TODO check after done
        # xx, yy = create_mesh(ra_center, dec_center, width_mesh)
        sig = sig_poisson(meshgrid[0], meshgrid[1], SIGMA1, SIGMA2,
                          coords[0], coords[1], RATIO_AREA_TG_BG)
    else:
        print('wrong kernel :(')

    np.save(SIGNI_FILE, sig)
    np.save(MESHFILE, np.array([xx, yy]))

    print('Yeah! Done with density estimation! :)')


if __name__ == '__main__':
    main()
