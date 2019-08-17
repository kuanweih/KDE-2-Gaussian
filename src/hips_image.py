import matplotlib
matplotlib.use('Agg')

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from typing import List
from src.tools import dist2
from astropy.coordinates import SkyCoord
from hips import WCSGeometry, make_sky_image
from src.param_patch_candidate import NSTAR_MIN, WIDTH_FAC



def multiprocessing_plot_hips_sky_image(
        name_df: List, label_df: List, hips_df: List, ra_df: List, dec_df: List,
        sigma1_df: List, sig_p_df: List, path: str, res: int, id_: int):
    """ Re-arange the order of arguments such that the only iterable 'id_'
    is at the last one. With this arangement, we can use multiprocessing
    under the help of partial from the functools.

    : name_df : list of system names
    : label_df : list of labels
    : hips_df : list of hips_surveyss
    : ra_df : list of average ra of pixels
    : dec_df : list of average dec of pixels
    : sigma1_df : list of sigma1 of images
    : sig_p_df : list of sig_poisson of pixels
    : path : output dir path
    : res : resolution of the image (number of pixels for the image)
    : id_ : iterable of all the target clusters
    """
    plot_hips_sky_image(ra_df[id_], dec_df[id_], sig_p_df[id_], sigma1_df[id_],
                        hips_df, path, name_df[id_], label_df[id_], res)


def plot_hips_sky_image(ra: float, dec: float, sig_p: float, sigma1: float,
        hips_surveys: List, outpath: str, name: str, label: int, res: int):
    """ Plot sky image using hips

    : ra : ra of the pixel
    : dec : dec of the pixel
    : sig_p : sig_poisson of the pixel
    : sigma1 : sigma1 of the map
    : hips_surveys : list of surveys
    : outpath : output dir path
    : name : name of the system (dwarf and more info)
    : label : label of a cluster pixels
    : res : resolution of the image (number of pixels for the image)
    """
    width = WIDTH_FAC * sigma1

    # Compute the sky image
    geometry = WCSGeometry.create(
        skydir=SkyCoord(ra, dec, unit='deg', frame='icrs'), width=res,
        height=res, fov="%f deg" %width, coordsys='icrs', projection='AIT')

    name_split = name.split('-')
    short_name = name_split[0]
    data = np.load('results/{}'.format(name.replace(
        '-poisson' or '-gaussian', '/queried-data.npy'))).item()

    ra_data, dec_data = data['ra'], data['dec']
    pmra_data, pmdec_data = data['pmra'], data['pmdec']
    bp_rp_data, g_mag_data = data['bp_rp'], data['phot_g_mean_mag']

    data = None    # free memory

    ra_min = ra - 0.5 * width
    ra_max = ra + 0.5 * width
    dec_min = dec - 0.5 * width
    dec_max = dec + 0.5 * width

    mask1 = (ra_min < ra_data) & (ra_data < ra_max)
    mask2 = (dec_min < dec_data) & (dec_data < dec_max)
    mask = mask1 & mask2

    ra_data, dec_data = ra_data[mask], dec_data[mask]
    pmra_data, pmdec_data = pmra_data[mask], pmdec_data[mask]
    bp_rp_data, g_mag_data = bp_rp_data[mask], g_mag_data[mask]

    is_in = dist2(ra_data, dec_data, ra, dec) <= sigma1 ** 2
    n_star_in = len(ra_data[is_in])

    if n_star_in < NSTAR_MIN:
        _s1 = 'skipping image for %s ' % short_name
        _s2 = 'because there are only %d stars in the kernel' % n_star_in
        print(st1 + _st2)
        return    # skip plotting image with fewer than NSTAR_MIN stars

    print('plotting image for %s' %short_name)
    sns.set(style="white", color_codes=True, font_scale=1)

    if 'gaia' in short_name:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        name_hue = ['inside' if i else 'outside' for i in is_in]
        ms = 20

        circle = plt.Circle((ra, dec), sigma1,
                            color='k', fill=False, lw=2, alpha=0.6)
        axes[0, 0].add_artist(circle)

        sns.scatterplot(x=ra_data, y=dec_data, s=ms, hue=name_hue, ax=axes[0, 0])
        sns.scatterplot(x=bp_rp_data, y=g_mag_data, s=ms, hue=name_hue, ax=axes[0, 1])
        sns.scatterplot(x=pmra_data, y=pmdec_data, s=ms, hue=name_hue, ax=axes[0, 2])

        sns.scatterplot(x=ra_data[is_in], y=dec_data[is_in], s=ms, ax=axes[0, 0])
        sns.scatterplot(x=bp_rp_data[is_in], y=g_mag_data[is_in], s=ms, ax=axes[0, 1])
        sns.scatterplot(x=pmra_data[is_in], y=pmdec_data[is_in], s=ms,  ax=axes[0, 2])


        axes[0, 0].set_title('stellar distribution')
        axes[0, 0].set_xlabel('ra (deg)')
        axes[0, 0].set_ylabel('dec (deg)')
        axes[0, 0].set_xlim([ra_min, ra_max])
        axes[0, 0].set_ylim([dec_min, dec_max])
        axes[0, 0].invert_xaxis()

        axes[0, 1].set_title('color mag diagram (CMD)')
        axes[0, 1].invert_yaxis()
        axes[0, 1].set_xlabel('Bp-Rp (mag)')
        axes[0, 1].set_ylabel('G (mag)')

        axes[0, 2].set_title('proper motions')
        axes[0, 2].set_xlabel('pmra (mas/yr)')
        axes[0, 2].set_ylabel('pmdec (mas/yr)')

        for u in range(3):
            for v in range(2):
                _asp = np.diff(axes[0, u].get_xlim())[0]
                _asp /= np.diff(axes[0, u].get_ylim())[0]
                axes[1, u].set_aspect(np.abs(_asp))

        cnt = 0    # counter for how many images have been plotted
        for hips_survey in hips_surveys:
            try:
                result = make_sky_image(
                            geometry=geometry, hips_survey=hips_survey,
                            tile_format='jpg', progress_bar=False)
                axes[1, cnt].imshow(result.image, origin='lower')
                axes[1, cnt].set_title(hips_survey)
                cnt += 1
            except:
                pass

            if cnt == 3:
                break

        for u in range(3):
            axes[1, u].tick_params(axis='both', which='both',
                                   labelleft=False, labelbottom=False)
    else:
        fig, axes = plt.subplots(1, 4, figsize=(15, 5))

        # plot sources
        axes[0].set_title(short_name)
        sns.scatterplot(ra_data, dec_data, ax=axes[0])
        axes[0].set_xlim([ra_min, ra_max])
        axes[0].set_ylim([dec_min, dec_max])
        _asp = np.diff(axes[0].get_xlim())[0] / np.diff(axes[0].get_ylim())[0]
        axes[0].set_aspect(_asp)
        axes[0].set_xlim(axes[0].set_xlim([ra_min, ra_max])[::-1])    # flipping

        # plot circle of the size of sigma1
        circle = plt.Circle((ra, dec), sigma1, color='orange', fill=False, lw=2)
        axes[0].add_artist(circle)

        cnt = 0    # counter for how many images have been plotted
        for hips_survey in hips_surveys:
            try:
                result = make_sky_image(
                            geometry=geometry, hips_survey=hips_survey,
                            tile_format='jpg', progress_bar=False)
                cnt += 1
                axes[cnt].imshow(result.image, origin='lower')
                axes[cnt].set_title(hips_survey)
            except:
                pass

            if cnt == 3:
                break

        for i in range(4):
            axes[i].tick_params(axis='both', which='both',
                                labelleft=False, labelbottom=False)


    _st1 = '{}-label{}-sigp%0.2f-'.format(short_name, label) % sig_p
    _st2 = '-width{}'.format(width)
    fig.suptitle('{}(%0.4f,%0.4f){}'.format(_st1, _st2) %(ra, dec))

    _str = '{}/{}-target{}-'.format(outpath, name, label)
    _str = '{}ra%0.4f-dec%0.4f-sigp%0.2f.jpg'.format(_str) % (ra, dec, sig_p)
    plt.savefig(_str, bbox_inches='tight', dpi=300)
