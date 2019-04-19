import numpy as np
from param import *
from classMWSatellite import *
from scipy.special import erfcinv
from scipy.stats import poisson
from typing import Tuple
from scipy.signal import fftconvolve, gaussian



class KDE_MWSatellite(MWSatellite):
    def __init__(self, name_sat: str, ra_sat: float, dec_sat: float, width: float, database: str,
                 catalog_str: str, pixel_size: float, sigma1: float, sigma2: float, sigma3: float,
                 sigma_th: int, factor_from_sigma2: float):
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
        self.factor_from_sigma2 = factor_from_sigma2

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

    def fftconvolve_boundary_adjust(self, hist2d: np.ndarray, kernel: np.ndarray) -> np.ndarray:
        """ use scipy signal fftconvolve to calculate the convolved map. edge effect is also taken care of.
            using fftconvolve will yeild some negative elements while they are actuaclly 0 when using
            convolve. therefore, conv.clip(min=0) is applied for the result.

        Parameters
        ----------
        hist2d : 2d histogram of the source distribution
        kernel : kernel matrix

        Returns
        -------
        np.ndarray : convolved map with non negative values
        """
        conv = fftconvolve(hist2d, kernel, mode='same')
        mask2d = np.ones(hist2d.shape)
        conv /= fftconvolve(mask2d, kernel, mode='same')
        return  conv.clip(min=0)

    def overdensity(self, sigma: float) -> np.ndarray:
        """ convolved overdensity map with Gaussian kernel size sigma """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size

        truncate = 5
        kernel_grid = int(truncate * s_grid)
        kernel = np.outer(gaussian(kernel_grid, s_grid), gaussian(kernel_grid, s_grid))
        kernel /= 2. * np.pi * s_grid**2
        return  self.fftconvolve_boundary_adjust(hist2d, kernel)

    def get_sig_gaussian(self, od_1: np.ndarray, od_2: np.ndarray, sigma1: float, sigma2: float) -> np.ndarray:
        """ get significance map using 2-gaussian kernel

        Parameters
        ----------
        od_1 : overdensity with sigma1
        od_2 : overdensity with sigma2
        sigma1 : inner kernel
        sigma2 : outer kernel

        Returns
        -------
        np.ndarray : significance from 2-gaussian kernel density estimation
        """
        s1 = sigma1 / self.pixel_size
        sig = (od_1 - od_2) * np.sqrt(4. * np.pi * s1**2)
        sig = np.divide(sig, np.sqrt(od_2), out=np.zeros_like(sig), where=od_2!=0)   # force 0 / 0 = 0
        return  sig

    def compound_significance(self):
        """ significance s12 inside (s23>sigma_th) and s13 outside (s23<sigma_th) """
        od_1 = self.overdensity(self.sigma1)
        od_2 = self.overdensity(self.sigma2)
        od_3 = self.overdensity(self.sigma3)

        s12 = self.get_sig_gaussian(od_1, od_2, self.sigma1, self.sigma2)
        s13 = self.get_sig_gaussian(od_1, od_3, self.sigma1, self.sigma3)
        s23 = self.get_sig_gaussian(od_2, od_3, self.sigma2, self.sigma3)
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

        is_insides = []
        sig_gaussian_stars = []
        sig_poisson_stars = []

        for i in range(n_source):
            id_x, id_y = int(id_xs[i]), int(id_ys[i])
            is_insides.append(self.is_inside[id_y][id_x])
            sig_gaussian_stars.append(self.sig_gaussian[id_y][id_x])
            sig_poisson_stars.append(self.sig_poisson[id_y][id_x])

        self.datas["is_inside"] = np.array(is_insides)
        self.datas["sig_gaussian"] = np.array(sig_gaussian_stars)
        self.datas["sig_poisson"] = np.array(sig_poisson_stars)

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
            calculate the z-score of the tail probability of poisson via N(0, 1)

        Parameters
        ----------
        lambda_poisson : lambda of poisson, here means the background expected number count
        x : number count of observed events, here means number of stars in the inner aperture

        Returns
        -------
        np.ndarray : z-score of poisson map
        """
        return  np.sqrt(2.) * erfcinv(2. * poisson.sf(x, lambda_poisson))

    def circular_kernel(self, radius: int) -> np.ndarray:
        """ calculate the circular kernel with radius according to sigma.
            this kernel is not normalized because we are using it to sum the number count """
        kernel = np.zeros((2 * radius + 1, 2 * radius + 1))
        y,x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
        mask = x**2 + y**2 <= radius**2
        kernel[mask] = 1
        return  kernel

    def poisson_inner_number_count(self, sigma: float) -> Tuple[np.ndarray, float]:
        """ calculate inner number count of stars within sigma,
            still using convolve so that it is cylindrical symmetic for targeted objects """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size
        kernel = self.circular_kernel(round(s_grid))
        norm_kernel = np.sum(kernel)
        conv = self.fftconvolve_boundary_adjust(hist2d, kernel / norm_kernel)
        return  conv * norm_kernel, norm_kernel

    def poisson_outer_expected_background(self, sigma: float, factor_sigma: float) -> Tuple[np.ndarray, float]:
        """ [cylindrical but comsuming large amount of memory]
            calculate expected outer number count of stars for background estimation,
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
        conv = self.fftconvolve_boundary_adjust(hist2d, kernel / norm_kernel)
        return  conv * norm_kernel, norm_kernel

    def get_sig_poisson(self):
        """ calculate significance on each pixel based on poisson statistics """
        n_inner, area_inner = self.poisson_inner_number_count(self.sigma1)
        n_outer, area_outer = self.poisson_outer_expected_background(self.sigma2, self.factor_from_sigma2)
        ratio = area_inner / area_outer    # area ratio = inner / outer
        lambda_poisson = n_outer * ratio    # estimated background count in inner area
        sig = self.z_score_poisson(lambda_poisson, n_inner)
        self.sig_poisson = sig















        #
