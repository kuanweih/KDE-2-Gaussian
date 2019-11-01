import numpy as np

from typing import Tuple
from wsdb import HOST, USER, PASSWORD

from src.param import *
from src.tools import create_dir, print_sep_line
from src.classKDE_MWSatellite import KDE_MWSatellite
from src.classPatchMWSatellite import PatchMWSatellite
from src.plotting import visualize_2_panel, hist_2_panel
from src.peaks import summarize_peaks_star_csv, summarize_peaks_pixel_csv


np.seterr(divide='ignore', invalid='ignore')



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


def execute_kde_routine(patch: PatchMWSatellite, kdepatch: KDE_MWSatellite):
    kdepatch.np_hist2d(patch.datas['ra'], patch.datas['dec'])
    kdepatch.add_masks_on_pixels(RA_DWARF, DEC_DWARF, R_HALFLIGHT)
    kdepatch.compound_sig_gaussian()
    kdepatch.compound_sig_poisson()
    patch.append_sig_to_data(kdepatch.x_mesh, kdepatch.y_mesh,
                             kdepatch.sig_gaussian, kdepatch.sig_poisson)



if __name__ == '__main__':
    # create result directory
    dir_name = get_dir_name()
    create_dir(dir_name)


    print('Creating a Patch object for main KDE calcuation: \n')
    Patch = PatchMWSatellite(NAME, RA, DEC, DISTANCE, WIDTH, DATABASE, CATALOG_STR)
    print(Patch.__str__())

    Patch.sql_get(HOST, USER, PASSWORD)    # query data

    # cuts based on surveys
    if DATABASE == 'gaia_dr2.gaia_source':
        gaia_patch_gmag_cut_astro_noise_cut(Patch)
    elif DATABASE == 'panstarrs_dr1.stackobjectthin':
        Patch.mask_panstarrs_stargalaxy_sep()

    Patch.append_is_inside(RA_DWARF, DEC_DWARF, R_HALFLIGHT)    # TODO add a factor here

    print_sep_line()

    if IS_PM_ERROR_CUT:
        Patch.mask_pm_error(PMRA_DWARF, PMDEC_DWARF, N_ERRORBAR)

    print('Creating a KDEPatch object and start the KDE calcuation: \n')
    KDEPatch = KDE_MWSatellite(RA, DEC, WIDTH, PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3, R_HALFLIGHT)
    print(KDEPatch.__str__())

    execute_kde_routine(Patch, KDEPatch)

    print('Saving datas, sigs, meshgrids ...')
    np.save('{}/{}'.format(dir_name, FILE_STAR), Patch.datas)
    np.save("{}/{}".format(dir_name, FILE_SIG_GAUSSIAN), KDEPatch.sig_gaussian)
    np.save("{}/{}".format(dir_name, FILE_SIG_POISSON), KDEPatch.sig_poisson)
    _meshs = np.array([KDEPatch.x_mesh, KDEPatch.y_mesh])
    np.save("{}/{}".format(dir_name, FILE_MESH), _meshs)
    print('Done =) \n')
    print('Finished KDE calculation. \n')
    print_sep_line()
    print('Generating plots and tables ... \n')

    # visualize searching results
    fig_name = dir_name.replace("results/", "")

    plot_dir = "plots"
    peaks_dir = "peaks"
    create_dir(plot_dir)
    create_dir(peaks_dir)

    visual_dir = "{}/visual".format(plot_dir)
    hist_dir = "{}/hist".format(plot_dir)
    create_dir(visual_dir)
    create_dir(hist_dir)

    star_dir = "{}/stars".format(peaks_dir)
    pixel_dir = "{}/pixels".format(peaks_dir)
    create_dir(star_dir)
    create_dir(pixel_dir)

    _kernels = ['poisson']    # ['gaussian', 'poisson']

    for k in _kernels:
        visualize_2_panel(dir_name, "{}/{}".format(visual_dir, fig_name), k)
        hist_2_panel(dir_name, "{}/{}".format(hist_dir, fig_name), k)

    _name_star = "{}/{}".format(star_dir, fig_name)
    _name_pixel = "{}/{}".format(pixel_dir, fig_name)
    summarize_peaks_star_csv(dir_name, _name_star, 'poisson')
    summarize_peaks_pixel_csv(dir_name, _name_pixel, 'poisson', RA, DEC, WIDTH)

    print("Done. \n")
    print("We are finished :) \n")
