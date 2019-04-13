import numpy as np
import time
from classMWSatellite import *
from classKDE_MWSatellite import *
from param import *
from kw_wsdb import *
from plotting import *


def create_dir(dir_name: str):
    """ create directory with a name 'dir_name' """
    import os
    import errno
    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def get_dir_name() -> str:
    """ get the name of results directory """
    dir_name = "results/{}-G{}-{}".format(NAME, G_MAG_MIN, G_MAG_MAX)
    dir_name = "{}-w{}-lp{}".format(dir_name, WIDTH, PIXEL_SIZE)
    dir_name = "{}-gc{}s{}s{}s{}sth{}".format(dir_name, GC_SIZE, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    return  dir_name


def n_source(Satellite: KDE_MWSatellite) -> int:
    return  len(Satellite.datas[Satellite.catalog_list[0]])



if __name__ == '__main__':
    # create result directory
    dir_name = get_dir_name()
    create_dir(dir_name)

    f= open("{}/stdout.txt".format(dir_name), "w+")    # dumping log information

    # create a KDE_MWSatellite object
    Satellite = KDE_MWSatellite(NAME, RA, DEC, WIDTH, DATABASE, CATALOG_STR,
                                PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3, SIGMA_TH, FACTOR_FROM_SIGMA2)
    f.write(Satellite.__str__())

    # query data
    Satellite.sql_get(HOST, USER, PASSWORD)
    f.write("\n\nUsing sqlutilpy.get() to query data...\n")
    f.write("{} sources are queried \n\n".format(n_source(Satellite)))

    # G band cut
    Satellite.mask_cut("phot_g_mean_mag", G_MAG_MIN, G_MAG_MAX)
    f.write("--> Cut: {} < {} < {}\n".format(G_MAG_MIN, "phot_g_mean_mag", G_MAG_MAX))
    f.write("--> {} sources left \n\n".format(n_source(Satellite)))

    # astrometric_excess_noise cut
    Satellite.mask_g_mag_astro_noise_cut()
    f.write("--> Cut: astrometric_excess_noise and phot_g_mean_mag\n")
    f.write("--> {} sources left \n\n".format(n_source(Satellite)))

    # get significance of gaussian
    t0 = time.time()
    Satellite.compound_significance()
    f.write("took %s seconds to calculate significance of Gaussian\n\n" % (time.time() - t0))

    # get significance of gaussian
    t0 = time.time()
    Satellite.get_sig_poisson()
    f.write("took %s seconds to calculate significance of Poisson\n\n" % (time.time() - t0))

    # append 'significance' and 'is_inside' to datas
    Satellite.append_sig_to_data()

    # save queried data, significance, mesh coordinates
    np.save("{}/{}".format(dir_name, FILE_STAR), Satellite.datas)
    np.save("{}/{}".format(dir_name, FILE_SIG_GAUSSIAN), Satellite.sig_gaussian)
    np.save("{}/{}".format(dir_name, FILE_SIG_POISSON), Satellite.sig_poisson)
    np.save("{}/{}".format(dir_name, FILE_MESH), np.array([Satellite.x_mesh, Satellite.y_mesh]))
    f.write("saved output npy files\n\n")

    # pm selection
    Satellite.get_pm_mean_std_inside()
    f.write("Starting source selection based on proper motion\n\n")


    if IS_PM_ERROR_CUT:
        np.seterr(divide='ignore', invalid='ignore')

        Satellite.mask_pm_error_cut(N_ERRORBAR)
        f.write("--> Cut: pm_mean within pm +- {} * pm_error \n".format(N_ERRORBAR))
        f.write("--> {} sources left \n\n".format(n_source(Satellite)))

        f.write("calculating significance pm_mean within pm +- pm_error\n")

        # get significance of gaussian
        t0 = time.time()
        Satellite.compound_significance()
        f.write("took %s seconds to calculate significance of Gaussian\n\n" % (time.time() - t0))

        # get significance of gaussian
        t0 = time.time()
        Satellite.get_sig_poisson()
        f.write("took %s seconds to calculate significance of Poisson\n\n" % (time.time() - t0))

        np.save("{}/{}-pm_error{}".format(dir_name, FILE_STAR, N_ERRORBAR), Satellite.datas)
        np.save("{}/{}-pm_error{}".format(dir_name, FILE_SIG_GAUSSIAN, N_ERRORBAR), Satellite.sig_gaussian)
        np.save("{}/{}-pm_error{}".format(dir_name, FILE_SIG_POISSON, N_ERRORBAR), Satellite.sig_poisson)
        f.write("saved output npy files\n\n")


    # visualize searching results
    plot_dir = "plots"
    fig_name = dir_name.replace("results/", "")
    create_dir(plot_dir)
    visualize_4_panel(dir_name, "{}/{}".format(plot_dir, fig_name), N_ERRORBAR, "gaussian")
    visualize_4_panel(dir_name, "{}/{}".format(plot_dir, fig_name), N_ERRORBAR, "poisson")


    f.write("we are finished :) \n\n")
    f.close()
