RA = 39.99708333333333      # ra of target (in deg)
DEC = -34.55083333333334    # dec of target (in deg)
RADIUS = 3                  # radius when querying data (in deg)
DATABASE = 'gaia_dr2.gaia_source'
CATALOGS = 'ra, dec, phot_g_mean_mag, astrometric_excess_noise'

""" g-band cut"""
G_MAG_CUT = 0    # 1:on 0:off
G_MAG_MIN = 17
G_MAG_MAX = 21

""" astro_ex_noise cut in eq. 1 of Koposov et al 2017: MNRAS 470"""
ASTRO_EX_NOISE_G_MAG_CUT = 0    # 1:on 0:off
ASTRO_EX_NOISE_G_MAG_MIN = None
ASTRO_EX_NOISE_G_MAG_MAX = 0
