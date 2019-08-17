import matplotlib
matplotlib.use('Agg')

import glob
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from wsdb import HOST, USER, PASSWORD
from src.tools import create_dir, print_sep_line
from src.classPatchMWSatellite import PatchMWSatellite
from src.param_patch_candidate import OUTPUT_PATH



def get_datas_gaia_patch(name: str) -> PatchMWSatellite:
    """ Create a PatchMWSatellite object based on a target of candidates.
    TODO: not implemented pm cuts yet. Perhaps change to use results files
    instead of querying new data.

    : name : file name of images of candidates
    : returns : PatchMWSatellite
    """

    database = 'gaia_dr2.gaia_source'
    cats = """ ra, dec, pmra, pmdec, bp_rp,
               phot_g_mean_mag, astrometric_excess_noise """

    if 'gaia' not in name:
        raise ValueError('Not a gaia name. Not implemented yet.')

    dwarf_name = name.split('gc')[0].split('-')[0]
    name = name.split('gc')[1].split('target')

    sigma1 = name[0].replace('-poisson-', '').replace('-gaussian-', '')
    gc_size = float(sigma1.split('s')[0])
    sigma1 = float(sigma1.split('s')[1])

    coor = name[1].split('-dec')
    dec = float(coor[1])
    ra = float(coor[0].split('-ra')[1])
    target = int(coor[0].split('-ra')[0])

    name = dwarf_name + '-target%d' % target
    dist = gc_size / sigma1 * 180. / np.pi
    width = 10. * sigma1    # TODO the factor is made up

    print('Creating a Patch object for candidate %s: \n' % name)
    Patch = PatchMWSatellite(name, ra, dec, dist, width, database, cats)
    print(Patch.__str__())
    Patch.sql_get(HOST, USER, PASSWORD)
    Patch.mask_cut("phot_g_mean_mag", 17, 22)    # TODO hard code G band cut
    Patch.mask_g_mag_astro_noise_cut()
    Patch.append_is_inside(ra, dec, sigma1)

    return  Patch


def plot_gaia_stars_cmd_pm(patch: PatchMWSatellite):
    sns.set(style="white", color_codes=True, font_scale=1)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.3)

    _title = patch.name_sat + '-dist%dkpc' % int(patch.dist / 1e3)
    _title = _title + '-width%0.4fdeg' % patch.width
    _title = _title + '-ra%0.4f' % patch.ra_sat + '-dec%0.4f' % patch.dec_sat
    fig.suptitle(_title)

    bp_rp = patch.datas['bp_rp']
    is_inside = patch.datas['is_inside']
    g_mag = patch.datas['phot_g_mean_mag']
    ra, dec = patch.datas['ra'], patch.datas['dec']
    pmra, pmdec = patch.datas['pmra'], patch.datas['pmdec']

    name_hue = ['inside' if i else 'outside' for i in is_inside]

    sns.scatterplot(x=ra, y=dec, hue=name_hue, ax=axes[0])
    sns.scatterplot(x=bp_rp, y=g_mag, hue=name_hue, ax=axes[1])
    sns.scatterplot(x=pmra, y=pmdec, hue=name_hue, ax=axes[2])

    sns.scatterplot(x=ra[is_inside], y=dec[is_inside], ax=axes[0])
    sns.scatterplot(x=bp_rp[is_inside], y=g_mag[is_inside], ax=axes[1])
    sns.scatterplot(x=pmra[is_inside], y=pmdec[is_inside],  ax=axes[2])

    fig.autofmt_xdate()
    axes[0].invert_xaxis()
    axes[1].invert_yaxis()
    axes[0].set_title('stellar distribution')
    axes[1].set_title('color mag diagram (CMD)')
    axes[2].set_title('proper motions')

    axes[0].set_xlabel('ra (deg)')
    axes[0].set_ylabel('dec (deg)')
    axes[1].set_xlabel('Bp-Rp (mag)')
    axes[1].set_ylabel('G (mag)')
    axes[2].set_xlabel('pmra (mas/yr)')
    axes[2].set_ylabel('pmdec (mas/yr)')

    _filename = '{}/{}.png'.format(OUTPUT_PATH, _title)
    plt.savefig(_filename, bbox_inches='tight', dpi=300)



if __name__ == '__main__':

    paths = glob.glob('images/*.jpg')
    print('We are now going to plot details of %d candidates \n' % len(paths))

    names = [path.replace('images/', '') for path in paths]
    names = [name.replace('.jpg', '') for name in names]

    print('Creating %s' % OUTPUT_PATH)
    create_dir(OUTPUT_PATH)

    for name in names:
        print_sep_line()
        patch = get_datas_gaia_patch(name)
        plot_gaia_stars_cmd_pm(patch)

    print('We are finished :)')
