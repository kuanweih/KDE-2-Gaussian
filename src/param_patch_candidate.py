""" parameters about patches """
# preprocess.py
PATCH_DIST = 0.9
N_PATCH_MAX = 4

WIDTH = 1    # width of a patch in deg



""" parameters about candidates """
# hips_image.py
NSTAR_MIN = 10    # plotting image if the candidate contains more than 10 stars
WIDTH_FAC = 10    # width of image = width_fac * sigma1
valid_width = 2. * PATCH_DIST - WIDTH    # deal with the boundary

# candidates.py
OUTPUT_PATH = 'plots/candidates'


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
