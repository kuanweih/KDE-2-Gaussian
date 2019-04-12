import numpy as np
from param import *
from classMWSatellite import *
from scipy.ndimage import gaussian_filter, convolve
from scipy.special import erfinv
from scipy.stats import poisson
from typing import Tuple


class KDE_MWSatellite(MWSatellite):
    def __init__(self, name_sat: str, ra_sat: float, dec_sat: float, width: float, database: str,
                 catalog_str: str, pixel_size: float, sigma1: float, sigma2: float, sigma3: float, sigma_th: int):
        """ Kernel Density Estimation on a MWSatellite object:
        pixel_size: size of pixel in deg
        sigma1: target kernel size in deg: GCs
        sigma2: smaller background kernel size in deg: inside the satellite
        sigma3: larger background kernel size in deg: outside the satellite
        """
        MWSatellite.__init__(self, name_sat, ra_sat, dec_sat, width, database, catalog_str)
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
        str = "{}\n{}----\n{}{}{}{}{}".format(str1, MWSatellite.__str__(self), str2, str3, str4, str5, str6)
        return  str

    def grid_coord(self, center: float) -> np.ndarray:
        """ get grid coordinates according to the center position and the width of the mesh """
        return  np.linspace(center - 0.5 * self.width, center + 0.5 * self.width, num=self.num_grid, endpoint=True)

    def np_hist2d(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ get histogram 2d for the star distribution on the mesh
        return: hist2d, x_coor, y_coor
        """
        return  np.histogram2d(self.datas["dec"], self.datas["ra"], bins=(self.y_mesh, self.x_mesh))

    def overdensity(self, sigma: float) -> np.ndarray:
        """ convolved overdensity map with Gaussian kernel size sigma """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size
        od = gaussian_filter(hist2d, s_grid, mode='constant')
        mask2d = np.ones(hist2d.shape)
        od /= gaussian_filter(mask2d, s_grid, mode='constant')
        return  od

    def significance(self, sigma1: float, sigma2: float) -> np.ndarray:
        """ get significance map using 2 kernels:
        sigma1: inner kernel
        sigma2: outer kernel
        """
        od_1 = self.overdensity(sigma1)
        od_2 = self.overdensity(sigma2)
        s1 = sigma1 / self.pixel_size
        sig = (od_1 - od_2) / np.sqrt(od_2 / (4. * np.pi * s1**2))
        return  sig

    def compound_significance(self):
        """ significance s12 inside (s23>sigma_th) and s13 outside (s23<sigma_th) """
        s12 = self.significance(self.sigma1, self.sigma2)
        s13 = self.significance(self.sigma1, self.sigma3)
        s23 = self.significance(self.sigma2, self.sigma3)
        mask_in = s23 > self.sigma_th    # mask for inside
        sig = s12 * mask_in + s13 * (~mask_in)
        self.is_inside = mask_in
        self.sig_gaussian = sig

    def append_sig_to_data(self):
        """ append significance of each star to the datas """
        n_source = len(self.datas["ra"])

        pixel_size_x = np.max(np.diff(self.x_mesh))
        pixel_size_y = np.max(np.diff(self.y_mesh))

        id_xs = (self.datas["ra"] - self.x_mesh[0]) / pixel_size_x
        id_ys = (self.datas["dec"] - self.y_mesh[0]) / pixel_size_y

        sig_stars = []
        is_insides = []

        for i in range(n_source):
            id_x, id_y = int(id_xs[i]), int(id_ys[i])
            sig_stars.append(self.sig_gaussian[id_y][id_x])
            is_insides.append(self.is_inside[id_y][id_x])

        self.datas["significance"] = np.array(sig_stars)
        self.datas["is_inside"] = np.array(is_insides)

    def get_pm_mean_std_inside(self):
        """ calculate the mean and std for the stars in the targeted dwarf """
        pmra = self.datas["pmra"]
        pmdec = self.datas["pmdec"]
        is_inside = self.datas["is_inside"]

        pmra = pmra[is_inside & ~np.isnan(pmra)]
        pmdec = pmdec[is_inside & ~np.isnan(pmdec)]

        pmra_mean = np.mean(pmra)
        pmdec_mean = np.mean(pmdec)
        pmra_std = np.std(pmra)
        pmdec_std = np.std(pmdec)

        self.pm_inside = {"pmra_mean":pmra_mean, "pmra_std":pmra_std,
                          "pmdec_mean":pmdec_mean, "pmdec_std":pmdec_std}

    def z_score_poisson(self, lambda_poisson: np.ndarray, x: np.ndarray) -> np.ndarray:
        """ z = sqrt(2) * erfinv(1 - 2 * sf(x, lambda)), where sf (survival function) = 1 - CDF
        lambda_poisson: lambda of poisson, here means the background expected number count
        x: number count of observed events, here means number of stars in the inner aperture
        """
        return  np.sqrt(2.) * erfinv(1. - 2. * poisson.sf(x, lambda_poisson))

    # def distance2(self, x_mesh: np.ndarray, y_mesh: np.ndarray, x_stars: np.ndarray, y_stars: np.ndarray) -> np.ndarray:
    #     """ 2d distance square between pixels and stars """
    #     dx = x_mesh[:, None] - x_stars
    #     dy = y_mesh[:, None] - y_stars
    #     return  dx**2 + dy**2

    def circular_kernel(self, radius: int) -> np.ndarray:
        """ a normalized circular kernel with radius according to sigma """
        kernel = np.zeros((2 * radius + 1, 2 * radius + 1))
        y,x = np.ogrid[- radius:radius + 1, - radius:radius + 1]
        mask = x**2 + y**2 <= radius**2
        kernel[mask] = 1
        return  kernel / np.sum(kernel)

    def poisson_inner_number_count(self, sigma: float) -> np.ndarray:
        """ calculate inner number count of stars within sigma """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size
        kernel = self.circular_kernel(round(s_grid))
        conv = convolve(hist2d, kernel, mode='constant')
        mask2d = np.ones(hist2d.shape)
        conv /= convolve(mask2d, kernel, mode='constant')
        return  conv


    def sig_poisson(self, sigma1: float, sigma2: float, factor_sigma2: float):
        n_inner = self.poisson_inner_number_count(self, sigma1)
        print(n_inner)




        # x, y = np.meshgrid(self.x_mesh, self.y_mesh)
        # print(x)
        # print(y)
        # dx = x[:, :, None] - self.datas["ra"]
        # print(dx.shape)
        # print(dx)
        # hist2d, x, y = np.histogram2d(self.datas["dec"], self.datas["ra"], bins=(self.y_mesh, self.x_mesh))


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
    #
    #     r12 = s1**2 / (r**2 - s2**2)    # area ratio = inner / outer
    #     lambda_poisson = n_outer * r12    # estimated background count
    #     sig = (n_inner - lambda_poisson) / np.sqrt(lambda_poisson)    # z score
    #     return sig







# >>> x1 = np.array([1, 2, 3, 4, 5])
# >>> x2 = np.array([1.3, 3.5, 2.8])
# >>> dx = x1[:,None] - x2
# >>> dx
# array([[-0.3, -2.5, -1.8],
#        [ 0.7, -1.5, -0.8],
#        [ 1.7, -0.5,  0.2],
#        [ 2.7,  0.5,  1.2],
#        [ 3.7,  1.5,  2.2]])
# >>> dx<1
# array([[ True,  True,  True],
#        [ True,  True,  True],
#        [False,  True,  True],
#        [False,  True, False],
#        [False, False, False]])
# >>> np.sum(dx<1, axis=1)
# array([3, 3, 2, 1, 0])





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
