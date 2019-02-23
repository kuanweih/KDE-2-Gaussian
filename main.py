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
    """
    create result directory
    """
    dir_name = "results-{}".format(KERNEL_BG)
    dir_name = "{}/{}/G{}-{}".format(dir_name, NAME, G_MAG_MIN, G_MAG_MAX)
    dir_name = "{}/w{}-lp{}".format(dir_name, WIDTH, PIXEL_SIZE)
    dir_name = "{}/s{}s{}s{}sth{}".format(dir_name, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    create_dir(dir_name)

    """
    open text file for dumping imformation
    """
    f= open("{}/stdout.txt".format(dir_name),"w+")



    Satellite = KDE_MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR,
                                PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    f.write(Satellite.__str__())

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


    np.save("{}/{}".format(dir_name, FILE_STAR), Satellite.datas)
    np.save("{}/{}".format(dir_name, FILE_SIG), Satellite.sig_gaussian)
    np.save("{}/{}".format(dir_name, FILE_MESH),
            np.array([Satellite.x_mesh, Satellite.y_mesh]))


    f.close()




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


    print('Yeah! Done with density estimation! :)')
