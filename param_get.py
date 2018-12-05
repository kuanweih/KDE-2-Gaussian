RA = 39.99708333333333      # ra of target (in deg)
DEC = -34.55083333333334    # dec of target (in deg)
RADIUS = 0.25                  # radius when querying data (in deg)
DATABASE = 'gaia_dr2.gaia_source'


""" remove nan proper motion stars """
REMOVE_PM_NAN = False    # True:on, Flase:off


""" proper motion cut for pmra and pmdec """
PM_CUT = False    # True:on, Flase:off
if PM_CUT:
    PM_CUT_STD = 1.


""" g-band cut"""
G_MAG_CUT = True    # True:on, Flase:off
if G_MAG_CUT:
    G_MAG_MIN = 17
    G_MAG_MAX = 21


""" astro_ex_noise cut in eq. 1 of Koposov et al 2017: MNRAS 470"""
ASTRO_EX_NOISE_G_MAG_CUT = True    # True:on, Flase:off
if ASTRO_EX_NOISE_G_MAG_CUT:
    ASTRO_EX_NOISE_G_MAG_MIN = None
    ASTRO_EX_NOISE_G_MAG_MAX = 0
