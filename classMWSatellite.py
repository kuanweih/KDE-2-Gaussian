import numpy as np
import sqlutilpy


class MWSatellite(object):
    def __init__(self, name_sat, ra_sat, dec_sat, width, database, catalog_str):
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
        str1 = "This is a MWSatellite object:\n"
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

        datas = sqlutilpy.get(query_str,
                              host=host, user=user, password=password)

        """ update 'datas' dic to store queried data """
        for i, catalog in enumerate(self.catalog_list):
            self.datas[catalog] = datas[i]

    def cut_datas(self, mask):
        """
        cut datas
        """
        for key, column in self.datas.items():
            self.datas[key] = column[mask]

    def mask_cut(self, catalog, min_val, max_val):
        """
        cut the data with a min and a max value
        """
        maskleft = min_val < self.datas[catalog]
        maskright = self.datas[catalog] < max_val
        mask = maskleft & maskright
        self.cut_datas(mask)

    def mask_g_mag_astro_noise_cut(self):
        """
        hard code the astrometric_excess_noise and phot_g_mean_mag cut
        """
        noise = self.datas["astrometric_excess_noise"]
        g_mag = self.datas["phot_g_mean_mag"]
        maskleft = (g_mag <= 18.) & (noise < np.exp(1.5))
        maskright = (18. < g_mag) & (noise < np.exp(1.5 + 0.3 * (g_mag - 18.)))
        mask = maskleft | maskright
        self.cut_datas(mask)

    def mask_pm_error_cut(self):
        """
        hard code the pm_error cut on pmra and pmdec
        """
        pmra_mean = self.pm_inside["pmra_mean"]
        pmdec_mean = self.pm_inside["pmdec_mean"]

        #  "pmdec_mean":pmdec_mean, "pmdec_std":pmdec_std}
        # self.datas["is_inside"]

        maskleft = self.datas["pmra"] - self.datas["pmra_error"] < pmra_mean
        maskright = pmra_mean < self.datas["pmra"] + self.datas["pmra_error"]
        mask = maskleft & maskright

        maskleft = self.datas["pmdec"] - self.datas["pmdec_error"] < pmdec_mean
        maskright = pmdec_mean < self.datas["pmdec"] + self.datas["pmdec_error"]
        mask = maskleft & maskright & mask

        self.cut_datas(mask)
