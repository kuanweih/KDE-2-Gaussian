import numpy as np

from src.classPatchMWSatellite import *
from src.classKDE_MWSatellite import *
from param.param import *
from param.wsdb import *
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


def calc_pm_dwarf_gaia(dwarf: PatchMWSatellite) -> Tuple[float]:
    """ Calculate the proper motions of a dwarf.

    : dwarf : PatchMWSatellite of the dwarf with a size of rh
    : return : mean and std of pmra and pmdec
    """
    dwarf.sql_get(HOST, USER, PASSWORD)    # query data
    gaia_patch_gmag_cut_astro_noise_cut(dwarf)

    pmra_mean = np.nanmean(dwarf.datas['pmra'])
    pmra_std = np.nanstd(dwarf.datas['pmra'])
    pmdec_mean = np.nanmean(dwarf.datas['pmdec'])
    pmdec_std = np.nanstd(dwarf.datas['pmdec'])

    print('The proper motions of the dwarf are:')
    print('    pmra: mean = %0.4f, std = %0.4f' %(pmra_mean, pmra_std))
    print('    pmdec: mean = %0.4f, std = %0.4f \n' %(pmdec_mean, pmdec_std))

    return  pmra_mean, pmra_std, pmdec_mean, pmdec_std


def execute_kde_routine(patch: PatchMWSatellite, kdepatch: KDE_MWSatellite):
    kdepatch.np_hist2d(patch.datas['ra'], patch.datas['dec'])
    kdepatch.is_inside_2d(RA_DWARF, DEC_DWARF, R_HALFLIGHT)    # TODO add a factor here
    kdepatch.compound_sig_gaussian()
    kdepatch.compound_sig_poisson(R_HALFLIGHT)
    patch.append_sig_to_data(kdepatch.x_mesh, kdepatch.y_mesh,
                             kdepatch.sig_gaussian, kdepatch.sig_poisson)


def apply_pm_cut_gaia(patch: PatchMWSatellite, n_error: float):
    """ Create a patch object for the dwarf so as to get pm limits.
    Then use the limits to cut the data.

    : patch : patch object
    : n_error : N_ERRORBAR
    """
    _str1 = 'Creating a Dwarf object within'
    _str2 = 'for proper motion selection: \n'
    print('{} rh = %0.4f deg {}'.format(_str1, _str2) % R_HALFLIGHT)
    # hard coding of the size for pm selection.
    Dwarf = PatchMWSatellite(NAME, RA_DWARF, DEC_DWARF,
                             DISTANCE, R_HALFLIGHT, DATABASE, CATALOG_STR)
    print(Dwarf.__str__())
    pmra_mean, pmra_std, pmdec_mean, pmdec_std = calc_pm_dwarf_gaia(Dwarf)

    del Dwarf    # free the memory though it might not be necessary
    print('Removed the Dwarf object since it is no longer needed. \n')
    print_sep_line()
    print('Selecting proper motion within %d sigma: \n' % n_error)

    pmramax = pmra_mean + n_error * pmra_std
    pmramin = pmra_mean - n_error * pmra_std
    pmdecmax = pmdec_mean + n_error * pmdec_std
    pmdecmin = pmdec_mean - n_error * pmdec_std

    patch.mask_cut("pmra", pmramin, pmramax)
    patch.mask_cut("pmdec", pmdecmin, pmdecmax)
    print_sep_line()



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
        apply_pm_cut_gaia(Patch, N_ERRORBAR)

    print('Creating a KDEPatch object and start the KDE calcuation: \n')
    KDEPatch = KDE_MWSatellite(RA, DEC, WIDTH, PIXEL_SIZE, SIGMA1, SIGMA2, SIGMA3)
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
    create_dir("plots")
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
