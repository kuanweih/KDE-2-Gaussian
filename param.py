""" Parameter file for KDE detector """


""" default (manual) target parameters """
NAME = 'Fornax'    # name of the dwarf
RA = 39.997      # ra of target (in deg)
DEC = -34.551    # dec of target (in deg)
WIDTH = 0.25     # map width when querying data (in deg)

PIXEL_SIZE = 0.001    # 1d pixel size in deg
SIGMA1 = 0.004    # searching scale in deg
SIGMA2 = 0.02    # background scale (smaller) in deg
SIGMA3 = 1.00    # background scale (larger) in deg

GC_SIZE = 10    # size of target globular clusters (pc)
SIGMA_TH = 1    # sigma threshold to define inside or outside


""" data base and catalog """
DATABASE = 'gaia_dr2.gaia_source'
CATALOG_STR = """
              ra, dec, parallax, pmra, pmdec, pmra_error, pmdec_error,
              phot_g_mean_mag, astrometric_excess_noise
              """


""" g-band cut """
G_MAG_MIN = 17
G_MAG_MAX = 22    # fainter cut at G=22 for Gaia DR2


""" pm cut based on std of the dwarf """
IS_PM_CUT_STD = False
if IS_PM_CUT_STD:
    PM_IN_STD = [3, 2, 1]    # must be in decresing order


""" pm cut based on pm_error """
IS_PM_ERROR_CUT = True
if IS_PM_ERROR_CUT:
    N_ERRORBAR = 2


""" Gaussian or Poisson """
KERNEL_BG = 'gaussian'    # background distribution: default 'gaussian'
if KERNEL_BG == 'poisson':
    DR_FROM_S2 = 5.    # delta distance outside from sigma2 in degree


""" output file name """
FILE_STAR = 'queried-data'    # output data file
FILE_SIG = 'significance'    # output significance file
FILE_MESH = 'meshgrids'    # output mesh grids


""" parse arguments from McConnachie list """
IS_FROM_McConnachie = False
if IS_FROM_McConnachie:
    import argparse
    import numpy as np

    parser = argparse.ArgumentParser(description='Set parameters for a specific dwarf')
    parser.add_argument('--name_dwarf', type=str, help='A dwarf name from McConnachie list')
    parser.add_argument('--gc_size_pc', type=int, help='Size of globular clusters: 10~100 pc')
    parser.add_argument('--scale_sigma2', type=float, nargs='?', const=1, default=1.,
                        help='sigma2 = scale_sigma2 * sigma2')
    args = parser.parse_args()

    NAME = args.name_dwarf    # name of the dwarf

    path_dwarfs = "dwarfs-McConnachie/dwarfs-McConnachie.npy"
    dwarfs_dict = np.load(path_dwarfs).item()

    mask = dwarfs_dict["GalaxyName"] == NAME

    for key, val in dwarfs_dict.items():
        dwarfs_dict[key] = val[mask]

    keys_need = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

    if dwarfs_dict["GalaxyName"][0] != NAME:
        print("Cannot find %s in GalaxyName" %NAME) # TODO RaiseError?

    RA = dwarfs_dict["RA_deg"][0]
    DEC = dwarfs_dict["Dec_deg"][0]
    WIDTH = 8. * dwarfs_dict["rh(arcmins)"][0] / 60. # TODO: round? or including e?
    WIDTH = float("{0:.4f}".format(WIDTH))

    SIGMA1 = GC_SIZE / dwarfs_dict["Distance_pc"][0] * 180. / np.pi
    SIGMA1 = float("{0:.4f}".format(SIGMA1))

    SIGMA2 = float("{0:.4f}".format(0.1 * dwarfs_dict["rh(arcmins)"][0] / 60.))
    SIGMA2 *= args.scale_sigma2

    SIGMA3 = 0.5 * WIDTH
    PIXEL_SIZE = 0.25 * SIGMA1

    GC_SIZE = args.gc_size_pc







#
#
# """ astro_ex_noise cut in eq. 1 of Koposov et al 2017: MNRAS 470"""
# ASTRO_EX_NOISE_G_MAG_CUT = True    # True:on, Flase:off
# if ASTRO_EX_NOISE_G_MAG_CUT:
#     ASTRO_EX_NOISE_G_MAG_MIN = None
#     ASTRO_EX_NOISE_G_MAG_MAX = 0
#
#
# """ remove nan proper motion stars """
# REMOVE_PM_NAN = False    # True:on, Flase:off
#
#
# """ proper motion cut for pmra and pmdec """
# PM_CUT = False    # True:on, Flase:off
# if PM_CUT:
#     PM_CUT_STD = 3.
#
#
# """ remove nan parallax stars """
# REMOVE_PARALLAX_NAN = False    # True:on, Flase:off
#
#
# """ parallax cut """
# PARALLAX_CUT = False    # True:on, Flase:off
# if PARALLAX_CUT:
#     PARALLAX_CUT_STD = 3.
#
#

#
#
