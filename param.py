RA = 39.997      # ra of target (in deg)
DEC = -34.551    # dec of target (in deg)
RADIUS = 0.25                  # radius when querying data (in deg)
DATABASE = 'gaia_dr2.gaia_source'

""" g-band cut"""
G_MAG_MIN = 17
G_MAG_MAX = 21


""" astro_ex_noise cut in eq. 1 of Koposov et al 2017: MNRAS 470"""
ASTRO_EX_NOISE_G_MAG_CUT = True    # True:on, Flase:off
if ASTRO_EX_NOISE_G_MAG_CUT:
    ASTRO_EX_NOISE_G_MAG_MIN = None
    ASTRO_EX_NOISE_G_MAG_MAX = 0


""" remove nan proper motion stars """
REMOVE_PM_NAN = False    # True:on, Flase:off


""" proper motion cut for pmra and pmdec """
PM_CUT = False    # True:on, Flase:off
if PM_CUT:
    PM_CUT_STD = 3.


""" remove nan parallax stars """
REMOVE_PARALLAX_NAN = False    # True:on, Flase:off


""" parallax cut """
PARALLAX_CUT = False    # True:on, Flase:off
if PARALLAX_CUT:
    PARALLAX_CUT_STD = 3.


PIXEL_SIZE = 0.002    # pixel size 1d in deg
SIGMA1 = 0.01    # searching scale (smaller) in degree
SIGMA2 = 0.05    # background scale (larger) in degree


KERNEL_BG = 'poisson'    # background distribution: default 'gaussian'
if KERNEL_BG == 'poisson':
    DR_FROM_S2 = 5.    # delta distance outside from sigma2 in degree
