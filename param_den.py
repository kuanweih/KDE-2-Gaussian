NUM_GRID = 200    # number of meshgrids
SIGMA1 = 0.01    # searching scale (smaller) in degree
SIGMA2 = 0.05    # background scale (larger) in degree


KERNEL_BG = 'poisson'    # background distribution: default 'gaussian'
if KERNEL_BG == 'poisson':
    DR_FROM_S2 = 5.    # delta distance outside from sigma2 in degree
