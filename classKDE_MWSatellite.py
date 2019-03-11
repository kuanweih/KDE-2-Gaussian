import numpy as np
from param import *
from classMWSatellite import *
from scipy.ndimage import gaussian_filter


class KDE_MWSatellite(MWSatellite):
    def __init__(self, name_sat, ra_sat, dec_sat, width, database, catalog_str,
                 pixel_size, sigma1, sigma2, sigma3, sigma_th):
        """
        Kernel Density Estimation on a MWSatellite object:
        pixel_size: size of pixel in deg
        sigma1: target kernel size in deg: GCs
        sigma2: smaller background kernel size in deg: inside the satellite
        sigma3: larger background kernel size in deg: outside the satellite
        """
        MWSatellite.__init__(self, name_sat, ra_sat, dec_sat, width,
                             database, catalog_str)
        self.pixel_size = pixel_size
        self.num_grid = round(self.width / self.pixel_size)
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3
        self.sigma_th = sigma_th

        self.x_mesh = self.grid_coord(self.ra_sat)
        self.y_mesh = self.grid_coord(self.dec_sat)

    def __str__(self):
        str1 = "This is a KDE_MWSatellite object inherited from\n----"
        str2 = "    pixel size = {}\n".format(self.pixel_size)
        str3 = "    number of grids = {}\n".format(self.num_grid)
        str4 = "    sigma1 = {} deg\n".format(self.sigma1)
        str5 = "    sigma2 inside the satellite = {} deg\n".format(self.sigma2)
        str6 = "    sigma2 outside the satellite = {} deg\n".format(self.sigma3)
        str = "{}\n{}----\n{}{}{}{}{}".format(str1, MWSatellite.__str__(self),
                                              str2, str3, str4, str5, str6)
        return str

    def grid_coord(self, center):
        """
        get grid coordinates according to the center position and the
        width of the mesh
        """
        return np.linspace(center - 0.5 * self.width,
                           center + 0.5 * self.width,
                           num=self.num_grid, endpoint=True)

    def overdensity(self, sigma):
        """
        convolved overdensity map with Gaussian kernel size sigma
        """
        hist2d, x, y = np.histogram2d(self.datas["dec"],
                                      self.datas["ra"],
                                      bins=(self.y_mesh, self.x_mesh))
        s_grid = sigma / self.pixel_size
        od = gaussian_filter(hist2d, s_grid, mode='constant')
        mask2d = np.ones(hist2d.shape)
        od /= gaussian_filter(mask2d, s_grid, mode='constant')
        return od

    def significance(self, sigma1, sigma2):
        """
        get significance map using 2 kernels:
        sigma1: inner kernel
        sigma2: outer kernel
        """
        od_1 = self.overdensity(sigma1)
        od_2 = self.overdensity(sigma2)
        s1 = sigma1 / self.pixel_size
        sig = (od_1 - od_2) / np.sqrt(od_2 / (4. * np.pi * s1**2))
        return sig

    def compound_significance(self):
        """
        significance s12 inside (s23>sigma_th) and s13 outside (s23<sigma_th)
        """
        s12 = self.significance(self.sigma1, self.sigma2)
        s13 = self.significance(self.sigma1, self.sigma3)
        s23 = self.significance(self.sigma2, self.sigma3)
        mask_in = s23 > self.sigma_th    # mask for inside
        sig = s12 * mask_in + s13 * (~mask_in)
        self.is_inside = mask_in
        self.sig_gaussian = sig

    def append_sig_to_data(self):
        """
        append significance of each star to the datas
        """
        n_source = len(self.datas["ra"])

        pixel_size_x = np.max(np.diff(self.x_mesh))
        pixel_size_y = np.max(np.diff(self.y_mesh))

        id_xs = (self.datas["ra"] - self.x_mesh[0]) / pixel_size_x
        id_ys = (self.datas["dec"] - self.y_mesh[0]) / pixel_size_y

        sig_stars = []
        is_insides = []

        for i in range(n_source):
            id_x = int(id_xs[i])
            id_y = int(id_ys[i])
            sig_stars.append(self.sig_gaussian[id_y][id_x])
            is_insides.append(self.is_inside[id_y][id_x])

        self.datas["significance"] = np.array(sig_stars)
        self.datas["is_inside"] = np.array(is_insides)













#TODO: add poisson statistics method into the class

# def distance2(x_arr, y_arr, x_cen, y_cen):
#     """
#     2d distance square
#     """
#     d2 = (x_arr - x_cen)**2 + (y_arr - y_cen)**2
#     return d2
#
#
# def sig_poisson(x, y, s1, s2, star_x, star_y, dr_s2):
#     """
#     get z-score as significance using inverse survival function of Poisson.
#     x, y: mesh arrays
#     star_x, star_y: position of stars
#     s1, s2: inner and outer scales
#     dr_s2: r_out = s2 + dr_s2
#     """
#     r = s2 + dr_s2    # outer radius
#     n_inner = np.sum(np.array([(distance2(x, y, star_x[i], star_y[i]) < s1**2)
#                                for i in range(len(star_x))]), axis=0)
#     n_outer = np.sum(np.array([(s2**2 < distance2(x, y, star_x[i], star_y[i])) *
#                                (distance2(x, y, star_x[i], star_y[i]) < r**2)
#                                for i in range(len(star_x))]), axis=0)
#     r12 = s1**2 / (r**2 - s2**2)    # area ratio = inner / outer
#     lambda_poisson = n_outer * r12    # estimated background count
#     sig = (n_inner - lambda_poisson) / np.sqrt(lambda_poisson)    # z score
#     return sig
#
