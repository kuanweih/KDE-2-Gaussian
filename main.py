import numpy as np
import time
from src.classMWSatellite import *
from src.classKDE_MWSatellite import *
from param.param import *
from param.kw_wsdb import *
from src.plotting import *
from src.peaks import *
from src.tools import *



def get_dir_name() -> str:
    """ Get the name of results directory """
    dir_name = "results/{}-G{}-{}".format(NAME, G_MAG_MIN, G_MAG_MAX)
    dir_name = "{}-w{}-lp{}".format(dir_name, WIDTH, PIXEL_SIZE)
    dir_name = "{}-gc{}s{}s{}s{}sth{}".format(dir_name, GC_SIZE,
                                              SIGMA1, SIGMA2, SIGMA3, SIGMA_TH)
    return  dir_name


def n_source(Satellite: KDE_MWSatellite) -> int:
    """ Calculter number of stars in queried dictionary """
    return  len(Satellite.datas[Satellite.catalog_list[0]])



if __name__ == '__main__':
    # create result directory
    dir_name = get_dir_name()
    create_dir(dir_name)

    f= open("{}/stdout.txt".format(dir_name), "w+")    # dumping log information

    # create a KDE_MWSatellite object
    Satellite = KDE_MWSatellite(NAME, RA, DEC, DISTANCE, WIDTH, DATABASE,
                                CATALOG_STR, PIXEL_SIZE, SIGMA1, SIGMA2,
                                SIGMA3, SIGMA_TH, R_HALFLIGHT)
    f.write(Satellite.__str__())

    # query data
    Satellite.sql_get(HOST, USER, PASSWORD)
    f.write("\n\nUsing sqlutilpy.get() to query data...\n")
    f.write("{} sources are queried \n\n".format(n_source(Satellite)))

    # G band cut
    Satellite.mask_cut("phot_g_mean_mag", G_MAG_MIN, G_MAG_MAX)
    f.write("--> Cut: {} < {} < {}\n".format(G_MAG_MIN,
                                             "phot_g_mean_mag", G_MAG_MAX))
    f.write("--> {} sources left \n\n".format(n_source(Satellite)))

    # astrometric_excess_noise cut
    Satellite.mask_g_mag_astro_noise_cut()
    f.write("--> Cut: astrometric_excess_noise and phot_g_mean_mag\n")
    f.write("--> {} sources left \n\n".format(n_source(Satellite)))

    # get significance of gaussian
    t0 = time.time()
    Satellite.compound_sig_gaussian()
    f.write("took %s sec to calculate Gaussian sig\n\n" % (time.time() - t0))

    # get significance of gaussian
    t0 = time.time()
    Satellite.compound_sig_poisson()
    f.write("took %s sec to calculate Poisson sig\n\n" % (time.time() - t0))

    # append 'significance' and 'is_inside' to datas
    Satellite.append_sig_to_data()

    # save queried data, significance, mesh coordinates
    np.save("{}/{}".format(dir_name, FILE_STAR), Satellite.datas)
    np.save("{}/{}".format(dir_name, FILE_SIG_GAUSSIAN), Satellite.sig_gaussian)
    np.save("{}/{}".format(dir_name, FILE_SIG_POISSON), Satellite.sig_poisson)
    np.save("{}/{}".format(dir_name, FILE_MESH), np.array([Satellite.x_mesh,
                                                           Satellite.y_mesh]))
    f.write("saved output npy files\n\n")

    # pm selection
    Satellite.get_pm_mean_std_inside()
    f.write("Starting source selection based on proper motion\n\n")


    if IS_PM_ERROR_CUT:
        np.seterr(divide='ignore', invalid='ignore')

        Satellite.mask_pm_error_cut(N_ERRORBAR)
        f.write("--> Cut: pm_mean within pm +- {} * pm_error\n".format(
                N_ERRORBAR))
        f.write("--> {} sources left \n\n".format(n_source(Satellite)))

        f.write("calculating significance pm_mean within pm +- pm_error\n")

        # get significance of gaussian
        t0 = time.time()
        Satellite.compound_sig_gaussian()
        f.write("took %s sec to calculate Gaussian sig\n\n" %(time.time() - t0))

        # get significance of gaussian
        t0 = time.time()
        Satellite.compound_sig_poisson()
        f.write("took %s sec to calculate Poisson sig\n\n" % (time.time() - t0))

        np.save("{}/{}-pm_error{}".format(dir_name, FILE_STAR,
                                          N_ERRORBAR), Satellite.datas)
        np.save("{}/{}-pm_error{}".format(dir_name, FILE_SIG_GAUSSIAN,
                                          N_ERRORBAR), Satellite.sig_gaussian)
        np.save("{}/{}-pm_error{}".format(dir_name, FILE_SIG_POISSON,
                                          N_ERRORBAR), Satellite.sig_poisson)
        f.write("saved output npy files\n\n")


    # visualize searching results
    fig_name = dir_name.replace("results/", "")
    plot_dir = "plots"
    create_dir("plots")

    visual_dir = "{}/visual".format(plot_dir)
    create_dir(visual_dir)
    # visualize_4_panel(dir_name, "{}/{}".format(visual_dir, fig_name),
    #                   N_ERRORBAR, "gaussian")
    visualize_4_panel(dir_name, "{}/{}".format(visual_dir, fig_name),
                      N_ERRORBAR, "poisson")

    hist_dir = "{}/hist".format(plot_dir)
    create_dir(hist_dir)
    # hist_2_panel(dir_name, "{}/{}".format(hist_dir, fig_name),
    #              N_ERRORBAR, "gaussian")
    hist_2_panel(dir_name, "{}/{}".format(hist_dir, fig_name),
                 N_ERRORBAR, "poisson")

    peaks_dir = "peaks"
    create_dir(peaks_dir)

    star_dir = "{}/stars".format(peaks_dir)
    create_dir(star_dir)
    # summarize_peaks_star_csv(dir_name, "{}/{}".format(star_dir, fig_name),
    #                          N_ERRORBAR, "gaussian")
    summarize_peaks_star_csv(dir_name, "{}/{}".format(star_dir, fig_name),
                             N_ERRORBAR, "poisson")

    pixel_dir = "{}/pixels".format(peaks_dir)
    create_dir(pixel_dir)
    # summarize_peaks_pixel_csv(dir_name, "{}/{}".format(pixel_dir, fig_name),
    #                           N_ERRORBAR, "gaussian", RA, DEC, WIDTH)
    summarize_peaks_pixel_csv(dir_name, "{}/{}".format(pixel_dir, fig_name),
                              N_ERRORBAR, "poisson", RA, DEC, WIDTH)

    f.write("we are finished :) \n\n")
    f.close()
