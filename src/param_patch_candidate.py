""" parameters about patches """
# preprocess.py
PATCH_DIST = 0.5
N_PATCH_MAX = 7

WIDTH = 1.5    # width of a patch in deg = PATCH_DIST + 2 * sigma3



""" parameters about candidates """
# hips_image.py
NSTAR_MIN = 5    # threshold of min stars in inner aperture
WIDTH_FAC = 10    # width of image = width_fac * sigma1
valid_width = PATCH_DIST    # deal with the boundary


# summary.py
s_above = 5    # significance threshold
res_image = 1000    # image resolution of hips

hips_surveys = ['CDS/P/DES-DR1/g',
                'CDS/P/DECaLS/DR5/color',
                'CDS/P/DSS2/color',
                'CDS/P/SDSS9/color',
                'CDS/P/PanSTARRS/DR1/color-z-zg-g',
                'CDS/P/2MASS/color',
                ]
