NAME = 'Fornax'
RA = 39.997      # ra of target (in deg)
DEC = -34.551    # dec of target (in deg)
WIDTH = 0.25     # map width when querying data (in deg)
DATABASE = 'gaia_dr2.gaia_source'
CATALOG_STR = """
              ra, dec, parallax, pmra, pmdec,
              phot_g_mean_mag, astrometric_excess_noise
              """

""" g-band cut"""
G_MAG_MIN = 17
G_MAG_MAX = 21

PIXEL_SIZE = 0.001    # 1d pixel size in deg
SIGMA1 = 0.004    # searching scale in deg
SIGMA2 = 0.05    # background scale (smaller) in deg
SIGMA3 = 1.00    # background scale (larger) in deg
SIGMA_TH = 3    # sigma threshold to define inside or outside



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
# KERNEL_BG = 'poisson'    # background distribution: default 'gaussian'
# if KERNEL_BG == 'poisson':
#     DR_FROM_S2 = 5.    # delta distance outside from sigma2 in degree
