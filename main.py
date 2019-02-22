import numpy as np
from classMWSatellite import *
from classKDE_MWSatellite import *
from param import *
from kw_wsdb import *


def create_dir(dir_name):
    """
    create directory with a name 'dir_name'
    """
    import os
    import errno
    if not os.path.exists(os.path.dirname(dir_name)):
        try:
            os.makedirs(os.path.dirname(dir_name))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


if __name__ == '__main__':

    # TODO finish the output part for all the details
    dir_name = "{}_G{}_{}".format(NAME, G_MAG_MIN, G_MAG_MAX)
    create_dir(dir_name)

    f= open("stdout.txt".format(dir_name),"w+")


    """ file name: FornaxG17_21/gaussian/w0.25/p0.001/s0.004s0.05s1.00sth3 """

    # Satellite = MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR)
    Satellite = KDE_MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR,
                                PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    f.write(Satellite)



    #TODO database account and password here



    Satellite.sql_get(HOST, USER, PASSWORD)
    Satellite.mask_cut("phot_g_mean_mag", G_MAG_MIN, G_MAG_MAX)
    Satellite.mask_g_mag_astro_noise_cut()
    Satellite.compound_significance()








    # # files names
    # FILENAME = 'queried-data'    # output data file
    # INFOFILE = 'queried-data-info'    # info of data file
    # SIGNI_FILE = 'significance'    # output significance file
    # MESHFILE = 'meshgrids'    # output mesh grids
    # # output
    # np.save(FILENAME, np.array(datas))
    # np.save(INFOFILE, np.array([RA, DEC, RADIUS, len(datas[0])]))
    #
    # print('\nYeah! Done with querying data from {}!\n'.format(DATABASE))

    #
    # # get significance
    # if KERNEL_BG == 'gaussian':
    #     print('We are using 2-Gaussian kernels to estimate the density...')
    #     s1_grid = SIGMA1 / PIXEL_SIZE
    #     s2_grid = SIGMA2 / PIXEL_SIZE
    #     sig = sig_2_gaussian(x_mesh, y_mesh, s1_grid, s2_grid, ra, dec)
    # elif KERNEL_BG == 'poisson':
    #     print('We are using Poisson statistics to estimate the density...')
    #     print('Background area = %0.1f detection area.' % DR_FROM_S2)
    #     meshgrid = np.meshgrid(x_mesh, y_mesh, sparse=True)
    #     sig = sig_poisson(meshgrid[0], meshgrid[1],
    #                       SIGMA1, SIGMA2, ra, dec, DR_FROM_S2)
    # else:
    #     print('wrong kernel :(')
    #
    # np.save(SIGNI_FILE, sig)
    # np.save(MESHFILE, np.array([x_mesh, y_mesh]))

    print('Yeah! Done with density estimation! :)')
