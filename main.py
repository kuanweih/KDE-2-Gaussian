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
    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


if __name__ == '__main__':
    """ create result directory """
    dir_name = "results_{}".format(KERNEL_BG)
    create_dir(dir_name)
    dir_name = "{}/{}/G{}_{}/w{}_p{}".format(dir_name, NAME,
                                             G_MAG_MIN, G_MAG_MAX,
                                             WIDTH, PIXEL_SIZE)
    dir_name = "{}/s{}s{}s{}sth{}".format(dir_name,
                                          SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    create_dir(dir_name)

    """ open text file for dumping imformation """
    f= open("{}/stdout.txt".format(dir_name),"w+")

    Satellite = KDE_MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR,
                                PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    f.write(Satellite.__str__())


    #TODO database account and password here



    f.write("\n\nUsing sqlutilpy.get() to query data...\n")
    Satellite.sql_get(HOST, USER, PASSWORD)
    n_source = len(Satellite.datas[Satellite.catalog_list[0]])
    f.write("{} sources are queried \n\n".format(n_source))

    f.write("--> Cut: {} < {} < {}\n".format(G_MAG_MIN, "phot_g_mean_mag", G_MAG_MAX))
    Satellite.mask_cut("phot_g_mean_mag", G_MAG_MIN, G_MAG_MAX)
    n_source = len(Satellite.datas[Satellite.catalog_list[0]])
    f.write("--> {} sources left \n\n".format(n_source))

    f.write("--> Cut: astrometric_excess_noise and phot_g_mean_mag\n")
    Satellite.mask_g_mag_astro_noise_cut()
    n_source = len(Satellite.datas[Satellite.catalog_list[0]])
    f.write("--> {} sources left \n\n".format(n_source))


    Satellite.compound_significance()



    f.close()







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
