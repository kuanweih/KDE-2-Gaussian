import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from hips import WCSGeometry, make_sky_image
from typing import List


def multiprocessing_plot_hips_sky_image(name_df: List, label_df: List,
                                        hips_df: List, ra_df: List,
                                        dec_df: List, width_df: List,
                                        path: str, res: int, id_: int):
    """ Re-arange the order of arguments such that the only iterable 'id_'
    is at the last one. With this arangement, we can use multiprocessing
    under the help of partial from the functools.

    : name_df : list of system names
    : label_df : list of labels
    : hips_df : list of hips_surveys
    : ra_df : list of average ra of pixels
    : dec_df : list of average dec of pixels
    : width_df : list of width of images
    : path : output dir path
    : res : resolution of the image (number of pixels for the image)
    : id_ : iterable of all the target clusters
    """
    plot_hips_sky_image(ra_df[id_], dec_df[id_], width_df[id_],
                        hips_df[id_], path, name_df[id_], label_df[id_], res)



def plot_hips_sky_image(ra: float, dec: float, width: float, hips_survey: str,
                        outpath: str, name: str, label: int, res: int):
    """ Plot sky image using hips

    : ra : ra of the pixel
    : dec : dec of the pixel
    : width : width of the map
    : hips_survey : which survey for the image
    : outpath : output dir path
    : name : name of the system (dwarf and more info)
    : label : label of a cluster pixels
    : res : resolution of the image (number of pixels for the image)
    """
    # Compute the sky image
    geometry = WCSGeometry.create(skydir=SkyCoord(ra, dec, unit='deg',
                                                  frame='galactic'),
                                  width=res, height=res, fov="%f deg" %width,
                                  coordsys='galactic', projection='AIT')

    result = make_sky_image(geometry=geometry,
                            hips_survey=hips_survey, tile_format='jpg')

    # Draw the sky image
    ax = plt.subplot(projection=geometry.wcs)
    ax.imshow(result.image, origin='lower')
    plt.savefig("{}/{}-target{}-{}.jpg".format(outpath, name, label,
                                               hips_survey.replace('/', '-')),
                bbox_inches='tight', dpi=300)
