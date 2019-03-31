


NAME = 'Fornax'
RA = 39.997      # ra of target (in deg)
DEC = -34.551    # dec of target (in deg)
WIDTH = 0.25     # map width when querying data (in deg)

PIXEL_SIZE = 0.001    # 1d pixel size in deg
SIGMA1 = 0.004    # searching scale in deg
SIGMA2 = 0.02    # background scale (smaller) in deg
SIGMA3 = 1.00    # background scale (larger) in deg

SIGMA_TH = 1    # sigma threshold to define inside or outside


""" data base and catalog """
DATABASE = 'gaia_dr2.gaia_source'
CATALOG_STR = """
              ra, dec, parallax, pmra, pmdec, pmra_error, pmdec_error,
              phot_g_mean_mag, astrometric_excess_noise
              """


""" g-band cut """
G_MAG_MIN = 17
G_MAG_MAX = 22


""" pm cut based on std of the dwarf """
IS_PM_CUT_STD = False
if IS_PM_CUT_STD:
    PM_IN_STD = [3, 2, 1]    # must be in decresing order


""" pm cut based on pm_error """
IS_PM_ERROR_CUT = False


""" Gaussian or Poisson """
KERNEL_BG = 'gaussian'    # background distribution: default 'gaussian'
if KERNEL_BG == 'poisson':
    DR_FROM_S2 = 5.    # delta distance outside from sigma2 in degree


""" output file name """
FILE_STAR = 'queried-data'    # output data file
FILE_SIG = 'significance'    # output significance file
FILE_MESH = 'meshgrids'    # output mesh grids



""" test how to read the dwarf list """
import numpy as np

path_dwarfs = "dwarfs-McConnachie/dwarfs-McConnachie.npy"
dwarfs_dict = np.load(path_dwarfs).item()

mask = dwarfs_dict["GalaxyName"] == NAME

for key, val in dwarfs_dict.items():
    dwarfs_dict[key] = val[mask]

keys_need = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

for key in keys_need:
    print(dwarfs_dict[key])






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
