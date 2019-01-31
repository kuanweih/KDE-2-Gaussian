import numpy as np
from get_ra_dec import *
from param import *


if __name__ == '__main__':
    """
    part 1: query data from database
    """
    print('Begin part 1: querying data \n')
    # files names
    FILENAME = 'queried-data'    # output data file
    INFOFILE = 'queried-data-info'    # info of data file
    SIGNI_FILE = 'significance'    # output significance file
    MESHFILE = 'meshgrids'    # output mesh grids

    # query position data from DATABASE
    catalog_str = """
                  ra, dec, parallax, pmra, pmdec,
                  phot_g_mean_mag, astrometric_excess_noise
                  """
    datas = sql_get(catalog_str)
    mask = [True] * len(datas[0])
    print('We are querying data from {}'.format(DATABASE))
    print('Centered at (%0.3f, %0.3f) within %0.2f degree' % (RA, DEC, RADIUS))
    print('with {} < phot_g_mean_mag < {}\n'.format(G_MAG_MIN, G_MAG_MAX))

    # apply selections
    if ASTRO_EX_NOISE_G_MAG_CUT:
        mask = astro_ex_noise_gmag_cut(datas, mask)
        print('    Applying astro_ex_noise and g_mag cut.')
    if REMOVE_PM_NAN:
        mask = remove_pm_nan(datas, mask)
        print('    Removing nan pm.')
    if PM_CUT:
        mask = pm_cut(datas, mask)
        print('    Selecting pmra and pmdec within {} std.'.format(PM_CUT_STD))
    if REMOVE_PARALLAX_NAN:
        mask = remove_parallax_nan(datas, mask)
        print('    Removing nan parallax.')
    if PARALLAX_CUT:
        mask = parallax_cut(datas, mask)
        print('    Selecting parallax within {} std.'.format(PARALLAX_CUT_STD))

    datas = [data[mask] for data in datas]

    # output
    np.save(FILENAME, np.array(datas))
    np.save(INFOFILE, np.array([RA, DEC, RADIUS, len(datas[0])]))

    print('\nYeah! Done with querying data from {}!\n'.format(DATABASE))


    """
    part 2: calculate kernel density estimation
    """
    print('Begin part 2: KDE \n')
    ra, dec = datas[0], datas[1]

    # create mesh
    width_mesh = 2. * RADIUS
    num_grid = round(width_mesh / PIXEL_SIZE)
    x_mesh = get_grid_coord(RA, width_mesh)
    y_mesh = get_grid_coord(DEC, width_mesh)
    print('Create a mesh with width = %0.2f deg' % width_mesh)
    print('There %d grids on a side.' % num_grid)
    print('Pixel size of each grid = %0.5f deg' % PIXEL_SIZE)
    print('Dectection scale is %0.4f degree' % SIGMA1)
    print('Background scale is %0.4f degree' % SIGMA2)
    print(' ')

    # get significance
    if KERNEL_BG == 'gaussian':
        print('We are using 2-Gaussian kernels to estimate the density...')
        sig = sig_2_gaussian(x_mesh, y_mesh, SIGMA1, SIGMA2, ra, dec)
    elif KERNEL_BG == 'poisson':
        print('We are using Poisson statistics to estimate the density...')
        print('Background area = %0.1f detection area.' % DR_FROM_S2)
        meshgrid = np.meshgrid(x_mesh, y_mesh, sparse=True)
        sig = sig_poisson(meshgrid[0], meshgrid[1],
                          SIGMA1, SIGMA2, ra, dec, DR_FROM_S2)
    else:
        print('wrong kernel :(')

    np.save(SIGNI_FILE, sig)
    np.save(MESHFILE, np.array([x_mesh, y_mesh]))

    print('Yeah! Done with density estimation! :)')
