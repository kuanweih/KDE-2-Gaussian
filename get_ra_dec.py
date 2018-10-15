import numpy as np
import sqlutilpy
from param_get import *
from kw_wsdb import *


def con_arr_astro_ex_noise_g_mag(astro_ex_noise, g_mag):
    """
    simplify eq 1 in Koposov et al 2017 (MNRAS 470) into exponential form
    """
    return astro_ex_noise - 10.**(0.15 * (g_mag - 15.) + 0.25)


# TODO to be replaced by returning mask_cut
def cut_datas(datas, con_arr, min_val, max_val):
    """
    apply condition to the queried data
    """
    if min_val != None:
        if max_val != None:
            mask = (min_val < con_arr) & (con_arr < max_val)
        else:
            mask = (min_val < con_arr)
    elif max_val != None:
        mask = (con_arr < max_val)
    else:
        print('Oops, no input minimum or maximum value here.')
    return [data[mask] for data in datas]


def mask_cut(datas, con_arr, min_val, max_val):
    """
    apply condition to the queried data
    """
    if min_val != None:
        if max_val != None:
            mask = (min_val < con_arr) & (con_arr < max_val)
        else:
            mask = (min_val < con_arr)
    elif max_val != None:
        mask = (con_arr < max_val)
    else:
        print('Oops, no input minimum or maximum value here.')
    return mask


def sql_get(catalog):
    """
    query 'catalog' from database using sqlutilpy.get()
    """
    query_str = 'select {} from {} where q3c_radial_query(ra, dec, {}, {}, {})'.format(
                catalog, DATABASE, RA, DEC, RADIUS)
    return sqlutilpy.get(query_str, host=HOST, user=USER, password=PASSWORD)


def main():
    # files names
    FILENAME = 'stars-coord'    # output file name
    INFOFILE = 'stars-coord-attr'    # info of the map

    # query data from DATABASE
    ra, dec = sql_get('ra, dec')
    datas = np.array([ra, dec])
    mask = [True] * len(ra)

    if G_MAG_CUT:
        g_mag = sql_get('phot_g_mean_mag')
        datas = np.concatenate((datas, [g_mag]), axis=0)
        mask = mask & mask_cut(datas, g_mag, G_MAG_MIN, G_MAG_MAX)

    if ASTRO_EX_NOISE_G_MAG_CUT:
        am_ex_noise = sql_get('astrometric_excess_noise')
        am_ex_no_g_mag = con_arr_astro_ex_noise_g_mag(am_ex_noise, g_mag)
        datas = np.concatenate((datas, [am_ex_noise]), axis=0)
        datas = mask & mask_cut(datas, am_ex_no_g_mag,
                                ASTRO_EX_NOISE_G_MAG_MIN,
                                ASTRO_EX_NOISE_G_MAG_MAX)

    datas = [data[mask] for data in datas]

    # output
    np.save(FILENAME, np.array([ra, dec]))
    np.save(INFOFILE, np.array([RA, DEC, RADIUS, len(ra)]))

    print('Yeah! Done with getting ra and dec!')


if __name__ == '__main__':
    main()
