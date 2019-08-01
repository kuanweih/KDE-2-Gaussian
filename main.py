import numpy as np

from src.classPatchMWSatellite import *
from src.classKDE_MWSatellite import *
from param.param import *
from param.kw_wsdb import *
from src.plotting import *
from src.peaks import *
from src.tools import *

np.seterr(divide='ignore', invalid='ignore')



def print_sep_line():
    print('------------------------------------------------------------ \n')


def get_dir_name() -> str:
    """ Get the name of results directory """
    if DATABASE == 'gaia_dr2.gaia_source':
        dir_name = "results/gaia_{}-G{}-{}".format(NAME, G_MAG_MIN, G_MAG_MAX)
    else:
        dir_name = "results/{}_{}".format(DATABASE_SHORT, NAME)
    dir_name = "{}-w{}-lp{}".format(dir_name, WIDTH, PIXEL_SIZE)
    dir_name = "{}-gc{}s{}s{}s{}".format(dir_name, GC_SIZE, SIGMA1, SIGMA2, SIGMA3)
    return  dir_name


def gaia_patch_gmag_cut_astro_noise_cut(patch: PatchMWSatellite):
    """ Query data and then apply G band cut and astro noise cut.
    This function just calls a few methods from PatchMWSatellite
    to modify Patch object.
    """
    patch.mask_cut("phot_g_mean_mag", G_MAG_MIN, G_MAG_MAX)    # G band cut
    patch.mask_g_mag_astro_noise_cut()    # astrometric_excess_noise cut



if __name__ == '__main__':
    # create result directory
    dir_name = get_dir_name()
    create_dir(dir_name)

    if DATABASE == 'gaia_dr2.gaia_source':
        print('Creating a Dwarf object within rh = %0.4f deg for proper motion selection: \n'
              % R_HALFLIGHT)    # hard coding of the size for pm selection.
        Dwarf = PatchMWSatellite(NAME, RA_DWARF, DEC_DWARF,
                                 DISTANCE, R_HALFLIGHT, DATABASE, CATALOG_STR)
        print(Dwarf.__str__())

        Dwarf.sql_get(HOST, USER, PASSWORD)    # query data
        gaia_patch_gmag_cut_astro_noise_cut(Dwarf)

        pmra_mean = np.nanmean(Dwarf.datas['pmra'])
        pmra_std = np.nanstd(Dwarf.datas['pmra'])
        pmdec_mean = np.nanmean(Dwarf.datas['pmdec'])
        pmdec_std = np.nanstd(Dwarf.datas['pmdec'])
        n_source_rh = Dwarf.n_source()    # number of stars are withing rh

        print('The proper motions of the dwarf are:')
        print('    pmra: mean = %0.4f, std = %0.4f' %(pmra_mean, pmra_std))
        print('    pmdec: mean = %0.4f, std = %0.4f \n' %(pmdec_mean, pmdec_std))

        del Dwarf    # free the memory though it might not be necessary
        print('Removed the Dwarf object since it is no longer needed. \n')


        print_sep_line()


    print('Creating a Patch object for main KDE calcuation: \n')
    Patch = PatchMWSatellite(NAME, RA, DEC, DISTANCE, WIDTH, DATABASE, CATALOG_STR)
    print(Patch.__str__())

    Patch.sql_get(HOST, USER, PASSWORD)    # query data
    if DATABASE == 'gaia_dr2.gaia_source':
        gaia_patch_gmag_cut_astro_noise_cut(Patch)

    Patch.append_is_inside(RA_DWARF, DEC_DWARF, R_HALFLIGHT)    # TODO add a factor here


    print_sep_line()


    print('Creating a KDEPatch object and start the KDE calcuation: \n')
    KDEPatch = KDE_MWSatellite(RA, DEC, WIDTH, PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3)
    print(KDEPatch.__str__())

    KDEPatch.np_hist2d(Patch.datas['ra'], Patch.datas['dec'])
    KDEPatch.is_inside_2d(RA_DWARF, DEC_DWARF, R_HALFLIGHT)    # TODO add a factor here

    KDEPatch.compound_sig_gaussian()
    KDEPatch.compound_sig_poisson(R_HALFLIGHT)

    Patch.append_sig_to_data(KDEPatch.x_mesh, KDEPatch.y_mesh,
                             KDEPatch.sig_gaussian, KDEPatch.sig_poisson)

    print('Saving datas, sigs, meshgrids ...')
    np.save('{}/{}'.format(dir_name, FILE_STAR), Patch.datas)
    np.save("{}/{}".format(dir_name, FILE_SIG_GAUSSIAN), KDEPatch.sig_gaussian)
    np.save("{}/{}".format(dir_name, FILE_SIG_POISSON), KDEPatch.sig_poisson)
    np.save("{}/{}".format(dir_name, FILE_MESH), np.array([KDEPatch.x_mesh,
                                                           KDEPatch.y_mesh]))
    print('Done =) \n')


    print_sep_line()


    if IS_PM_ERROR_CUT:
        print('Selecting proper motion within %d sigma: \n' % N_ERRORBAR)
        pmramax = pmra_mean + N_ERRORBAR * pmra_std
        pmramin = pmra_mean - N_ERRORBAR * pmra_std
        pmdecmax = pmdec_mean + N_ERRORBAR * pmdec_std
        pmdecmin = pmdec_mean - N_ERRORBAR * pmdec_std

        # proper motion selection
        Patch.mask_cut("pmra", pmramin, pmramax)
        Patch.mask_cut("pmdec", pmdecmin, pmdecmax)

        print('Creating a KDEPatch object and start the KDE calcuation: \n')
        KDEPatch = KDE_MWSatellite(RA, DEC, WIDTH, PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3)
        print(KDEPatch.__str__())

        KDEPatch.np_hist2d(Patch.datas['ra'], Patch.datas['dec'])
        KDEPatch.is_inside_2d(RA_DWARF, DEC_DWARF, R_HALFLIGHT)    # TODO add a factor here

        KDEPatch.compound_sig_gaussian()
        KDEPatch.compound_sig_poisson(R_HALFLIGHT)

        Patch.append_sig_to_data(KDEPatch.x_mesh, KDEPatch.y_mesh,
                                 KDEPatch.sig_gaussian, KDEPatch.sig_poisson)

        print('Saving datas, sigs, meshgrids ...')
        np.save('{}/{}-pm{}std'.format(dir_name, FILE_STAR, N_ERRORBAR), Patch.datas)
        np.save("{}/{}-pm{}std".format(dir_name,
                FILE_SIG_GAUSSIAN, N_ERRORBAR), KDEPatch.sig_gaussian)
        np.save("{}/{}-pm{}std".format(dir_name,
                FILE_SIG_POISSON, N_ERRORBAR), KDEPatch.sig_poisson)
        np.save("{}/{}-pm{}std".format(dir_name, FILE_MESH,
                N_ERRORBAR), np.array([KDEPatch.x_mesh, KDEPatch.y_mesh]))
        print('Done =) \n')

    print('Finished KDE calculation. \n')


    print_sep_line()


    print('Generating plots and tables ... \n')


    # TODO OVERLAPPING FILE NAME TO BE FIXED!!!


    # visualize searching results
    fig_name = dir_name.replace("results/", "")
    plot_dir = "plots"
    create_dir("plots")

    visual_dir = "{}/visual".format(plot_dir)
    create_dir(visual_dir)
    # visualize_4_panel(dir_name, "{}/{}".format(visual_dir, fig_name),
    #                   N_ERRORBAR, "gaussian")
    if DATABASE == 'gaia_dr2.gaia_source':
        visualize_4_panel(dir_name, "{}/{}".format(visual_dir, fig_name), N_ERRORBAR, "poisson")
    else:
        visualize_2_panel(dir_name, "{}/{}".format(visual_dir, fig_name), "poisson")


    hist_dir = "{}/hist".format(plot_dir)
    create_dir(hist_dir)
    # hist_2_panel(dir_name, "{}/{}".format(hist_dir, fig_name),
    #              N_ERRORBAR, "gaussian")
    if DATABASE == 'gaia_dr2.gaia_source':
        hist_4_panel(dir_name, "{}/{}".format(hist_dir, fig_name), N_ERRORBAR, "poisson")
    else:
        hist_2_panel(dir_name, "{}/{}".format(hist_dir, fig_name), "poisson")

    peaks_dir = "peaks"
    create_dir(peaks_dir)

    star_dir = "{}/stars".format(peaks_dir)
    create_dir(star_dir)
    # summarize_peaks_star_csv(dir_name, "{}/{}".format(star_dir, fig_name),
    #                          N_ERRORBAR, "gaussian")
    if DATABASE == 'gaia_dr2.gaia_source':
        summarize_peaks_star_csv_gaia(dir_name,
            "{}/{}".format(star_dir, fig_name), N_ERRORBAR, "poisson")
    else:
        summarize_peaks_star_csv(dir_name,
            "{}/{}".format(star_dir, fig_name), "poisson")

    pixel_dir = "{}/pixels".format(peaks_dir)
    create_dir(pixel_dir)
    # summarize_peaks_pixel_csv(dir_name, "{}/{}".format(pixel_dir, fig_name),
    #                           N_ERRORBAR, "gaussian", RA, DEC, WIDTH)
    if DATABASE == 'gaia_dr2.gaia_source':
        summarize_peaks_pixel_csv_gaia(dir_name, "{}/{}".format(pixel_dir, fig_name),
                              N_ERRORBAR, "poisson", RA, DEC, WIDTH)
    else:
        summarize_peaks_pixel_csv(dir_name, "{}/{}".format(pixel_dir, fig_name),
                              "poisson", RA, DEC, WIDTH)

    print("Done. \n")
    print("We are finished :) \n")
