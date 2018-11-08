NUM_GRID = 100    # number of meshgrids
SIGMA1 = 0.01    # searching scale (smaller) in degree
SIGMA2 = 0.05    # background scale (larger) in degree


KERNEL_BG = 'poisson'    # background distribution: default 'gaussian'
if KERNEL_BG == 'poisson':
    RATIO_AREA_TG_BG = 5.    # ratio of area between target and background
