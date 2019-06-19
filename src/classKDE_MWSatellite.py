import numpy as np
from param.param import *
from src.classMWSatellite import *
from scipy.special import erfcinv
from scipy.stats import poisson
from typing import Tuple
from scipy.signal import fftconvolve, gaussian



class KDE_MWSatellite(MWSatellite):
    def __init__(self, name_sat: str, ra_sat: float, dec_sat: float,
                 dist: float, width: float, database: str, catalog_str: str,
                 pixel_size: float, sigma1: float, sigma2: float,
                 sigma3: float, sigma_th: int, rh: float):
        """ Kernel Density Estimation on a MWSatellite object

        : pixel_size : size of pixel in deg
        : sigma1 : target kernel size in deg: GCs
        : sigma2 : smaller background kernel size in deg: inside the satellite
        : sigma3 : larger background kernel size in deg: outside the satellite
        : sigma_th : sigma threshold to define inside or outside
        : rh : half-light radius of the satellite in deg
        """
        MWSatellite.__init__(self, name_sat, ra_sat, dec_sat,
                             dist, width, database, catalog_str)

        self.pixel_size = pixel_size
        self.num_grid = round(self.width / self.pixel_size)
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3
        self.sigma_th = sigma_th
        self.rh = rh

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
        return  str

    def grid_coord(self, center: float) -> np.ndarray:
        """ Get grid coordinates according to the center position and the width
        of the mesh """
        return  np.linspace(center - 0.5 * self.width,
                            center + 0.5 * self.width,
                            num=self.num_grid, endpoint=True)

    def np_hist2d(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Get histogram 2d for the star distribution on the mesh

        : return : hist2d, x_coor, y_coor
        """
        return  np.histogram2d(self.datas["dec"], self.datas["ra"],
                               bins=(self.y_mesh, self.x_mesh))

    def fftconvolve_boundary_adjust(self, hist2d: np.ndarray,
                                    kernel: np.ndarray) -> np.ndarray:
        """ Use scipy signal fftconvolve to calculate the convolved map.
        Edge effect is also taken care of. Using fftconvolve will yeild
        some negative elements while they are actuaclly 0 when using
        convolve. Also, there will be some positive noises which need to
        be taken care of. Therefore, conv[conv < 1e-15] = 0 is applied
        to get rid of the false divisions.

        : hist2d : 2d histogram of the source distribution
        : kernel : kernel matrix

        : return : convolved map with non negative values
        """
        conv = fftconvolve(hist2d, kernel, mode='same')
        mask2d = np.ones(hist2d.shape)
        conv /= fftconvolve(mask2d, kernel, mode='same')

        conv[conv < 1e-15] = 0.    # rounding the noise < 1e-15
        return  conv

    def overdensity(self, sigma: float) -> np.ndarray:
        """ Convolved overdensity map with Gaussian kernel size sigma """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size

        truncate = 5
        kernel_grid = int(truncate * s_grid)
        kernel = np.outer(gaussian(kernel_grid, s_grid),
                          gaussian(kernel_grid, s_grid))
        kernel /= 2. * np.pi * s_grid**2
        return  self.fftconvolve_boundary_adjust(hist2d, kernel)

    def get_sig_gaussian(self, od_1: np.ndarray, od_2: np.ndarray,
                         sigma1: float, sigma2: float) -> np.ndarray:
        """ Get significance map using 2-gaussian kernel

        : od_1 : overdensity with sigma1
        : od_2 : overdensity with sigma2
        : sigma1 : inner kernel
        : sigma2 : outer kernel

        : return : significance from 2-gaussian kernel density estimation
        """
        s1 = sigma1 / self.pixel_size
        sig = (od_1 - od_2) * np.sqrt(4. * np.pi * s1**2)
        sig = np.divide(sig, np.sqrt(od_2),
                        out=np.zeros_like(sig), where=od_2!=0)  # force 0/0 = 0
        return  sig

    def compound_sig_gaussian(self):
        """ Compound the Gaussian significance map: s12 inside (s23 > sigma_th)
        and s13 outside (s23 < sigma_th) """
        od_1 = self.overdensity(self.sigma1)
        od_2 = self.overdensity(self.sigma2)
        od_3 = self.overdensity(self.sigma3)

        s12 = self.get_sig_gaussian(od_1, od_2, self.sigma1, self.sigma2)
        s13 = self.get_sig_gaussian(od_1, od_3, self.sigma1, self.sigma3)
        s23 = self.get_sig_gaussian(od_2, od_3, self.sigma2, self.sigma3)

        self.is_inside = s23 > self.sigma_th    # mask for inside
        self.sig_gaussian = s12 * self.is_inside + s13 * (~self.is_inside)

    # TODO separate this method into 3 parts so that I can free
    # sig_gaussian and delete it for memory concern
    def append_sig_to_data(self):
        """ Append significance of each star to the datas """
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
        """ Calculate the mean and std for the stars in the targeted dwarf """
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

    def mask_pm(self, pm: np.ndarray, pm_error: np.ndarray,
                pm_mean: float, pm_err: float, n_err: float) -> np.ndarray:
        """ Calculate the mask based on the pm selection.

        : pm : pm array: pmra or pmdec
        : pm_error : pm_error array: pmra or pmdec
        : pm_mean : mean pm inside the dwarf
        : pm_err : pm error on the mean: pm_std / number of stars in the dwarf
        : n_err : within n_err of error bars for pm_error

        : return : pm mask
        """
        maskleft = pm - n_err * pm_error < pm_mean + pm_err
        maskright = pm_mean - pm_err < pm + n_err * pm_error
        return  maskleft & maskright

    def mask_pm_error_cut(self, n_err: float):
        """ Hard code the pm_error cut on pmra and pmdec.

        : n_err : within n_err of error bars for pm_error
        """
        pmra_mean = self.pm_inside["pmra_mean"]
        pmdec_mean = self.pm_inside["pmdec_mean"]

        n = self.datas["is_inside"].sum()
        pmra_err = self.pm_inside["pmra_std"] / n
        pmdec_err = self.pm_inside["pmdec_std"] / n

        mk_pmra = self.mask_pm(self.datas["pmra"], self.datas["pmra_error"],
                               pmra_mean, pmra_err, n_err)
        mk_pmdec = self.mask_pm(self.datas["pmdec"], self.datas["pmdec_error"],
                                pmdec_mean, pmdec_err, n_err)
        self.cut_datas(mk_pmra & mk_pmdec)

    def z_score_poisson(self, lamb: np.ndarray, x: np.ndarray) -> np.ndarray:
        """ Calculate the z-score of the tail probability of poisson via N(0, 1)
        according to z = sqrt(2) * erfinv(1 - 2 * sf(x, lambda)), where sf
        (survival function) = 1 - CDF.

        : lamb : expected background number count from outer aperture (lambda)
        : x : number count of observed stars in the inner aperture

        : return : z-score of poisson map
        """
        return  np.sqrt(2.) * erfcinv(2. * poisson.sf(x, lamb))

    def circular_kernel(self, radius: int) -> np.ndarray:
        """ Calculate the circular kernel with radius according to sigma. This
        kernel is not normalized because we are using it to sum over the number
        count """
        kernel = np.zeros((2 * radius + 1, 2 * radius + 1))
        y,x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
        mask = x**2 + y**2 <= radius**2
        kernel[mask] = 1
        return  kernel

    def poisson_inner_number_count(self,
                                   sigma: float) -> Tuple[np.ndarray, float]:
        """ Calculate inner number count of stars within the area of radius of
            sigma, using convolution.

        : return : convolved map (number count of the inner aperture)
        : return : number of pixels of the inner aperture
        """
        hist2d, _, _ = self.np_hist2d()
        s_grid = sigma / self.pixel_size
        kernel = self.circular_kernel(round(s_grid))
        norm_kernel = np.sum(kernel)
        conv = self.fftconvolve_boundary_adjust(hist2d, kernel / norm_kernel)
        return  conv * norm_kernel, norm_kernel


    def poisson_outer_expected_background(self, sigma_in: float,
            sigma_out: float) -> Tuple[np.ndarray, float]:
        """ Calculate expected backgound number count of stars based on the
        outer aperture, which will be using to calculate 'lambda' for poisson.

        : sigma_in : inner radius of the outer aperture
        : sigma_out : outer radius of the outer aperture

        : return : convolved map (backgound estimation)
        : return : number of pixels of the outer aperture
        """
        hist2d, _, _ = self.np_hist2d()
        s_grid_in = round(sigma_in / self.pixel_size)
        s_grid_out = round(sigma_out / self.pixel_size)
        kernel_out = self.circular_kernel(s_grid_out)
        ds_pad = s_grid_out - s_grid_in
        kernel_in_pad = np.pad(self.circular_kernel(s_grid_in),
                               ds_pad, 'constant', constant_values=0)

        kernel = kernel_out - kernel_in_pad
        norm_kernel = np.sum(kernel)
        conv = self.fftconvolve_boundary_adjust(hist2d, kernel / norm_kernel)
        return  conv * norm_kernel, norm_kernel


    def get_lambda_poisson(self, n_o: np.ndarray, area_o: np.ndarray,
                           area_i: np.ndarray) -> np.ndarray:
        """ Calculate lambda as the estimated background number count.

        : n_o : number of sources from outer aperture
        : area_o : area of outer aperture
        : area_i : area of inner aperture

        : return : lambda (estimated background count)
        """
        lambda_poisson = n_o * area_i / area_o
        return  lambda_poisson

    def compound_sig_poisson(self):
        """ Compound the Poisson significance map: s12 inside (s23 > sigma_th)
        and s13 outside (s23 < sigma_th) """
        # factors using for outer aperture
        f_in2out = 2.    # r_i = f_in2out * s1
        f_out2out = 0.68    # r_o = f_out2out * rh
        rh_th = 0.1    # threshold of minimum half-light radius in deg

        # inner aperture
        n_inner, area_inner = self.poisson_inner_number_count(self.sigma1)

        # outer aperture outside of the dwarf
        r_i = f_in2out * self.sigma1
        r_o = self.sigma3
        n_outer, area_outer = self.poisson_outer_expected_background(r_i, r_o)
        lambda_out = self.get_lambda_poisson(n_outer, area_outer, area_inner)

        # outer aperture inside the dwarf
        if self.rh < rh_th:
            lambda_in = lambda_out
        else:
            r_o = f_out2out * float(self.rh)
            n_outer, area_outer = self.poisson_outer_expected_background(r_i, r_o)
            lambda_in = self.get_lambda_poisson(n_outer, area_outer, area_inner)

        s12 = self.z_score_poisson(lambda_in, n_inner)
        s13 = self.z_score_poisson(lambda_out, n_inner)

        self.sig_poisson = s12 * self.is_inside + s13 * (~self.is_inside)
