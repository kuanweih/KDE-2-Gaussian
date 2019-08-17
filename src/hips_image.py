import matplotlib
matplotlib.use('Agg')

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from typing import List
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

    data = None    # free memory

    ra_min = ra - 0.5 * width
    ra_max = ra + 0.5 * width
    dec_min = dec - 0.5 * width
    dec_max = dec + 0.5 * width

    mask1 = (ra_min < ra_data) & (ra_data < ra_max)
    mask2 = (dec_min < dec_data) & (dec_data < dec_max)
    mask = mask1 & mask2

    ra_data = ra_data[mask]
    dec_data = dec_data[mask]

    if len(ra_data) < NSTAR_MIN:
        print('skipping image for %s' %short_name)
        return    # skip plotting image with fewer than NSTAR_MIN stars

    print('plotting image for %s' %short_name)
    sns.set(style="white", color_codes=True, font_scale=1)
    fig, axes = plt.subplots(1, 4, figsize=(15, 5))

    _st1 = '{}-label{}-sigp%0.2f-'.format(short_name, label) % sig_p
    _st2 = '-width{}'.format(width)
    fig.suptitle('{}(%0.4f,%0.4f){}'.format(_st1, _st2) %(ra, dec), y=0.9)

    # plot GAIA sources
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
            result = make_sky_image(geometry=geometry, hips_survey=hips_survey,
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

    # axes = plt.subplot(projection=geometry.wcs)

    _str = '{}/{}-target{}-'.format(outpath, name, label)
    _str = '{}ra%0.4f-dec%0.4f.jpg'.format(_str) % (ra, dec)
    plt.savefig(_str, bbox_inches='tight', dpi=300)
