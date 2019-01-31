import numpy as np
import sqlutilpy
from param import *
from kw_wsdb import *


def con_arr_astro_ex_noise_g_mag(astro_ex_noise, g_mag):
    """
    simplify eq 1 in Koposov et al 2017 (MNRAS 470) into exponential form
    """
    # return astro_ex_noise - 10.**(0.15 * (g_mag - 15.) + 0.25)
    return astro_ex_noise - np.exp(1.5 + 0.3 * (g_mag - 18.))


def mask_cut(con_arr, min_val, max_val):
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
    query_str = """
                select {} from {}
                where q3c_radial_query(ra, dec, {}, {}, {})
                      and {} < phot_g_mean_mag and phot_g_mean_mag < {}
                """.format(catalog, DATABASE,
                           RA, DEC, RADIUS,
                           G_MAG_MIN, G_MAG_MAX)
    return sqlutilpy.get(query_str, host=HOST, user=USER, password=PASSWORD)


def astro_ex_noise_gmag_cut(datas, mask):
    g_mag, am_noise = datas[5], datas[6]
    am_ex_no_g_mag = con_arr_astro_ex_noise_g_mag(am_noise, g_mag)
    return mask & mask_cut(am_ex_no_g_mag,
                           ASTRO_EX_NOISE_G_MAG_MIN,
                           ASTRO_EX_NOISE_G_MAG_MAX)


def remove_pm_nan(datas, mask):
    pmra, pmdec = datas[3], datas[4]
    return mask & (~np.isnan(pmra)) & (~np.isnan(pmdec))


def pm_cut(datas, mask):
    pmra, pmdec = datas[3], datas[4]

    pmra_mean = np.mean(pmra[~np.isnan(pmra)])
    pmra_std = np.std(pmra[~np.isnan(pmra)])
    pmdec_mean = np.mean(pmdec[~np.isnan(pmdec)])
    pmdec_std = np.std(pmdec[~np.isnan(pmdec)])

    pmra_min = pmra_mean - PM_CUT_STD * pmra_std
    pmra_max = pmra_mean + PM_CUT_STD * pmra_std
    pmdec_min = pmdec_mean - PM_CUT_STD * pmdec_std
    pmdec_max = pmdec_mean + PM_CUT_STD * pmdec_std

    return mask & mask_cut(pmra, pmra_min,
                           pmra_max) & mask_cut(pmdec, pmdec_min, pmdec_max)


def remove_parallax_nan(datas, mask):
    parallax = datas[2]
    return mask & (~np.isnan(parallax))


def parallax_cut(datas, mask):
    plx = datas[2]    # parallax

    plx_mean = np.mean(plx[~np.isnan(plx)])
    plx_std = np.std(plx[~np.isnan(plx)])

    plx_min = plx_mean - PM_CUT_STD * plx_std
    plx_max = plx_mean + PM_CUT_STD * plx_std

    return mask & mask_cut(plx, plx_min, plx_max)
