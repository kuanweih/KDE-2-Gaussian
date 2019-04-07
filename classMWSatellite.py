import numpy as np
import sqlutilpy


class MWSatellite(object):
    def __init__(self, name_sat: str, ra_sat: float, dec_sat: float,
                 width: float, database: str, catalog_str: str):
        """ Milky Way (MW) Satellite object:
        name_sat: name of the satellite, e.g. Fornax
        ra_sat: RA of the satellite in deg
        dec_sat: Dec of the satellite in deg
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
        str1 = "This is a MWSatellite object:\n"
        str2 = "    name = {}\n".format(self.name_sat)
        str3 = "    ra = {}\n    dec = {}\n".format(self.ra_sat, self.dec_sat)
        str4 = "    map width = {} deg\n".format(self.width)
        str5 = "    database = {}\n".format(self.database)
        str = "{}{}{}{}{}".format(str1, str2, str3, str4, str5)
        return str

    def sql_get(self, host: str, user: str, password: str):
        """ query 'catalog_str' from 'database' using sqlutilpy.get() """
        ra_min = self.ra_sat - 0.5 * self.width
        ra_max = self.ra_sat + 0.5 * self.width
        dec_min = self.dec_sat - 0.5 * self.width
        dec_max = self.dec_sat + 0.5 * self.width

        query_str = """
                    select {} from {}
                    where {} < ra and ra < {} and {} < dec and dec < {}
                    """.format(self.catalog_str, self.database,
                               ra_min, ra_max, dec_min, dec_max)

        datas = sqlutilpy.get(query_str,
                              host=host, user=user, password=password)

        # update 'datas' dic to store queried data
        for i, catalog in enumerate(self.catalog_list):
            self.datas[catalog] = datas[i]

    def cut_datas(self, mask):
        """ cut datas based on the mask
        mask: numpy array
        """
        for key, column in self.datas.items():
            self.datas[key] = column[mask]

    def mask_cut(self, catalog: str, min_val: float, max_val: float):
        """ cut the data with a min and a max value """
        maskleft = min_val < self.datas[catalog]
        maskright = self.datas[catalog] < max_val
        mask = maskleft & maskright
        self.cut_datas(mask)

    def mask_g_mag_astro_noise_cut(self):
        """ hard code the astrometric_excess_noise and phot_g_mean_mag cut """
        noise = self.datas["astrometric_excess_noise"]
        g_mag = self.datas["phot_g_mean_mag"]
        maskleft = (g_mag <= 18.) & (noise < np.exp(1.5))
        maskright = (18. < g_mag) & (noise < np.exp(1.5 + 0.3 * (g_mag - 18.)))
        mask = maskleft | maskright
        self.cut_datas(mask)

    def mask_pm_error_cut(self, n_err: float):
        """ hard code the pm_error cut on pmra and pmdec
        n_err: within n_err of error bars for pm_error
        """
        pmra_mean = self.pm_inside["pmra_mean"]
        pmdec_mean = self.pm_inside["pmdec_mean"]

        n = self.datas["is_inside"].sum()
        pmra_err = self.pm_inside["pmra_std"] / n
        pmdec_err = self.pm_inside["pmdec_std"] / n

        maskleft = self.datas["pmra"] - n_err * self.datas["pmra_error"] < pmra_mean + pmra_err
        maskright = pmra_mean - pmra_err < self.datas["pmra"] + n_err * self.datas["pmra_error"]
        mask = maskleft & maskright

        maskleft = self.datas["pmdec"] - n_err * self.datas["pmdec_error"] < pmdec_mean + pmdec_err
        maskright = pmdec_mean - pmdec_err < self.datas["pmdec"] + n_err * self.datas["pmdec_error"]
        mask = maskleft & maskright & mask

        self.cut_datas(mask)
