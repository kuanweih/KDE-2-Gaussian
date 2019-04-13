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

    def circular_kernel(self, radius: int) -> np.ndarray:
        """ calculate the circular kernel with radius according to sigma.
            this kernel is not normalized because we are using it to sum the number count """
        kernel = np.zeros((2 * radius + 1, 2 * radius + 1))
        y,x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
        mask = x**2 + y**2 <= radius**2
        kernel[mask] = 1
        return  kernel

    def poisson_inner_number_count(self, sigma: float) -> Tuple[np.ndarray, float]:
        """ calculate inner number count of stars within sigma """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size
        kernel = self.circular_kernel(round(s_grid))
        norm_kernel = np.sum(kernel)
        conv = convolve(hist2d, kernel / norm_kernel, mode='constant')
        mask2d = np.ones(hist2d.shape)
        conv /= convolve(mask2d, kernel / norm_kernel, mode='constant')
        return  conv * norm_kernel, norm_kernel

    def poisson_outer_expected_background(self, sigma: float, factor_sigma: float) -> Tuple[np.ndarray, float]:
        """ calculate expected outer number count of stars for background estimation,
            which will be using to calculate 'lambda' for poisson """
        hist2d, _, _ = self.np_hist2d()
        # note that in and out here mean sigma2 and outer radius.
        s_grid_in = round(sigma / self.pixel_size)
        s_grid_out = round(factor_sigma * sigma / self.pixel_size)
        kernel_out = self.circular_kernel(s_grid_out)
        ds_pad = s_grid_out - s_grid_in
        kernel_in_pad = np.pad(self.circular_kernel(s_grid_in), ds_pad, 'constant', constant_values=0)

        kernel = kernel_out - kernel_in_pad
        norm_kernel = np.sum(kernel)

        conv = convolve(hist2d, kernel / norm_kernel, mode='constant')
        mask2d = np.ones(hist2d.shape)
        conv /= convolve(mask2d, kernel / norm_kernel, mode='constant')
        return  conv * norm_kernel, norm_kernel

    def sig_poisson(self, sigma1: float, sigma2: float, factor_sigma2: float):
        n_inner, area_inner = self.poisson_inner_number_count(sigma1)
        n_outer, area_outer = self.poisson_outer_expected_background(sigma2, factor_sigma2)
        ratio = area_inner / area_outer    # area ratio = inner / outer
        lambda_poisson = n_outer * ratio    # estimated background count in inner area
        sig = self.z_score_poisson(lambda_poisson, n_inner)
        self.sig_poisson = sig
