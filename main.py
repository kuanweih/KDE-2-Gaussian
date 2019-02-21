import numpy as np
from classMWSatellite import *
from classKDE_MWSatellite import *
from param import *
from kw_wsdb import *


if __name__ == '__main__':
    print('Begin part 1: querying data \n')
    # Satellite = MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR)
    Satellite = KDE_MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR,
                                PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    print(Satellite)

    #TODO database account and password here




    Satellite.sql_get(HOST, USER, PASSWORD)
    print(Satellite.datas["phot_g_mean_mag"].shape)
    Satellite.mask_cut("phot_g_mean_mag", G_MAG_MIN, G_MAG_MAX)
    print(Satellite.datas["phot_g_mean_mag"].shape)
    Satellite.mask_g_mag_astro_noise_cut()
    print(Satellite.datas["phot_g_mean_mag"].shape)
    # print(Satellite.overdensity(SIGMA2))
    print(Satellite.compound_significance().max())
    print(Satellite.compound_significance())


    # create mesh
    # x_mesh = get_grid_coord(RA, WIDTH, NUM_GRID)
    # y_mesh = get_grid_coord(DEC, WIDTH, NUM_GRID)
    # print('Create a mesh with width = %0.2f deg' % WIDTH)
    # print('There %d grids on a side.' % NUM_GRID)
    # print('Pixel size of each grid = %0.5f deg' % PIXEL_SIZE)
    # print('Dectection scale is %0.4f degree' % SIGMA1)
    # print('Background scale is %0.4f degree' % SIGMA2)
    # print(' ')




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
    #
    # """
    # part 2: calculate kernel density estimation
    # """
    # print('Begin part 2: KDE \n')
    # ra, dec = datas[0], datas[1]
    #
    # # create mesh
    # WIDTH = RADIUS
    # NUM_GRID = round(WIDTH / PIXEL_SIZE)
    # x_mesh = get_grid_coord(RA, WIDTH, NUM_GRID)
    # y_mesh = get_grid_coord(DEC, WIDTH, NUM_GRID)
    # print('Create a mesh with width = %0.2f deg' % WIDTH)
    # print('There %d grids on a side.' % NUM_GRID)
    # print('Pixel size of each grid = %0.5f deg' % PIXEL_SIZE)
    # print('Dectection scale is %0.4f degree' % SIGMA1)
    # print('Background scale is %0.4f degree' % SIGMA2)
    # print(' ')
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
