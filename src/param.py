""" Parameter file for KDE detector """


""" which dwarf list to be used, if none, then use the default manual one """
IS_DWARF_LIST = False    # use joint list
IS_DWARF_SPLIT_LIST = True    # use joint-split list


""" data base and catalog """
DATABASE = 'gaia_dr2.gaia_source'
# DATABASE = 'panstarrs_dr1.stackobjectthin'

if DATABASE == 'gaia_dr2.gaia_source':
    DATABASE_SHORT = 'gaia'
    CATALOG_STR = """ ra, dec, pmra, pmdec, bp_rp,
                      phot_g_mean_mag, astrometric_excess_noise """
elif DATABASE == 'panstarrs_dr1.stackobjectthin':
    DATABASE_SHORT = 'panstarrs'
    CATALOG_STR = 'ra, dec, rpsfmag, rkronmag'
else:
    raise ValueError('Wrong database name: %s is not valid.' % DATABASE)


""" default (manual) target parameters """
NAME = 'Fornax'    # name of the dwarf
RA = 39.997      # ra of the center of a patch (in deg)
DEC = -34.449    # dec of the center of a patch (in deg)
WIDTH = 0.25     # map width when querying data (in deg)
DISTANCE = 140000.    # distance in pc

RA_DWARF = RA      # ra of the dwarf (in deg)
DEC_DWARF = DEC    # dec of the dwarf (in deg)

PIXEL_SIZE = 0.001    # 1d pixel size in deg
SIGMA1 = 0.004    # searching scale in deg
SIGMA2 = 0.02    # background scale (smaller) in deg
SIGMA3 = 1.00    # background scale (larger) in deg

GC_SIZE = 10    # size of target globular clusters (pc)
R_HALFLIGHT = 0.28    # half light radius in deg


""" gaia cuts """
if DATABASE == 'gaia_dr2.gaia_source':
    # g-band cut
    G_MAG_MIN = 17
    G_MAG_MAX = 22    # fainter cut at G=22 for Gaia DR2

    # pm cut based on pm_error
    IS_PM_ERROR_CUT = False
    if IS_PM_ERROR_CUT:
        N_ERRORBAR = 5


""" panstarrs cuts """
if DATABASE == 'panstarrs_dr1.stackobjectthin':
    # pm cut based on pm_error
    IS_PM_ERROR_CUT = False


""" output file name """
FILE_STAR = 'queried-data'    # output data file
FILE_SIG_GAUSSIAN = 'sig_gaussian'    # output significance file
FILE_SIG_POISSON = 'sig_poisson'    # output significance file
FILE_MESH = 'meshgrids'    # output mesh grids


""" parse arguments from the joint-split dwarf list or the joint dwarf list """
if IS_DWARF_SPLIT_LIST or IS_DWARF_LIST:
    import argparse
    import numpy as np

    parser = argparse.ArgumentParser(description='Set parameters for a specific dwarf')
    parser.add_argument('--name_dwarf', type=str, help='A dwarf name from McConnachie list')
    parser.add_argument('--gc_size_pc', type=int, help='Size of globular clusters: e.g. 1~10 pc')
    parser.add_argument('--scale_sigma2', type=float, nargs='?', const=1,
                        default=1., help='sigma2 = scale_sigma2 * sigma2')
    args = parser.parse_args()

    if IS_DWARF_SPLIT_LIST:
        path_dwarfs = "dwarfs/dwarfs-joint-split.npy"
    elif IS_DWARF_LIST:
        path_dwarfs = "dwarfs/dwarfs-joint.npy"
    else:
        raise ValueError('Wrong list boolean')

    dwarfs_dict = np.load(path_dwarfs).item()

    NAME = args.name_dwarf    # name of the dwarf
    mask = dwarfs_dict["GalaxyName"] == NAME

    for key, val in dwarfs_dict.items():
        dwarfs_dict[key] = val[mask]

    if dwarfs_dict["GalaxyName"][0] != NAME:
        raise ValueError("Cannot find %s in GalaxyName" %NAME)

    RA = dwarfs_dict["RA_deg"][0]    # RA of the center of a patch
    DEC = dwarfs_dict["Dec_deg"][0]    # Dec of the center of a patch
    RA_DWARF = dwarfs_dict["RA_dwarf_deg"][0]      # ra of the dwarf (in deg)
    DEC_DWARF = dwarfs_dict["Dec_dwarf_deg"][0]    # dec of the dwarf (in deg)
    R_HALFLIGHT = dwarfs_dict["rh(arcmins)"][0] / 60.    # rh in deg

    if IS_DWARF_SPLIT_LIST:
        from src.param_patch_candidate import WIDTH
    elif IS_DWARF_LIST:
        WIDTH = float('%0.4f' %(8. * R_HALFLIGHT))
    else:
        raise ValueError('Wrong list boolean')

    GC_SIZE = args.gc_size_pc
    DISTANCE = float('%0.4f' %(dwarfs_dict["Distance_pc"][0]))

    SIGMA1 = GC_SIZE / DISTANCE * 180. / np.pi
    SIGMA1 = float('%0.4f' %(SIGMA1))

    SIGMA2 = float('%0.4f' %(0.2 * R_HALFLIGHT * GC_SIZE / 10.))    # factor 0.2 is made up
    SIGMA2 *= args.scale_sigma2

    SIGMA3 = 0.5    # always use 0.5 deg as outer kernel
    PIXEL_SIZE = 0.25 * SIGMA1


if __name__ == '__main__':
    if WIDTH > 5:
        print(NAME)
