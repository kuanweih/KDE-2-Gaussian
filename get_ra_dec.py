import numpy as np
import sqlutilpy
from param import *
from kw_wsdb import *


class MWSatellite(object):
    def __init__(self, name_sat, ra_sat, dec_sat, width,
                 database, catalog_str):
        """
        Milky Way (MW) Satellite object:
        name_sat: name of the satellite, e.g. Fornax
        ra_sat: ra of the satellite in deg
        dec_sat: dec of the satellite in deg
        width: width of the square area when querying data in deg
        database: database to be queried
        catalog_str: a string of catalogs for querying
        """
        self.name_sat = name_sat
        self.ra_sat = ra_sat
        self.dec_sat = dec_sat
        self.width = width
        self.database = database
        self.catalog_str = catalog_str
        self.catalog_list = self.catalog_str.replace("\n", "").replace(" ", "").split(",")
        self.datas = {}

    def __str__(self):
        str1 = "This is a MW Satellite object:\n"
        str2 = "    name = {}\n".format(self.name_sat)
        str3 = "    ra = {}\n    dec = {}\n".format(self.ra_sat, self.dec_sat)
        str4 = "    map width = {} deg\n".format(self.width)
        str5 = "    database = {}\n".format(self.database)
        str = "{}{}{}{}{}".format(str1, str2, str3, str4, str5)
        return str

    def sql_get(self, host, user, password):
        """
        query 'catalog_str' from 'database' using sqlutilpy.get()
        """
        ra_min = self.ra_sat - 0.5 * self.width
        ra_max = self.ra_sat + 0.5 * self.width
        dec_min = self.dec_sat - 0.5 * self.width
        dec_max = self.dec_sat + 0.5 * self.width

        query_str = """
                    select {} from {}
                    where {} < ra and ra < {} and {} < dec and dec < {}
                    """.format(self.catalog_str, self.database,
                               ra_min, ra_max, dec_min, dec_max)

        """ use sqlutilpy.get() to query data """
        datas = sqlutilpy.get(query_str,
                              host=host, user=user, password=password)

        """ create 'datas' dic to store queried data """
        for i, catalog in enumerate(self.catalog_list):
            self.datas[catalog] = datas[i]





""" test code """
catalog_str = """
              ra, dec, parallax, pmra, pmdec,
              phot_g_mean_mag, astrometric_excess_noise
              """
#TODO: CHANGE RADIUS INTO WIDTH
Fornax = MWSatellite('Fornax', RA, DEC, RADIUS, DATABASE, catalog_str)
print(Fornax)
print(Fornax.datas["astrometric_excess_noise"].shape)




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

    plx_min = plx_mean - PARALLAX_CUT_STD * plx_std
    plx_max = plx_mean + PARALLAX_CUT_STD * plx_std

    return mask & mask_cut(plx, plx_min, plx_max)
