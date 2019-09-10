import numpy as np
import pandas as pd

from typing import Dict
from src.tools import get_dic_list_npy, dist2
from src.param_patch_candidate import PATCH_DIST, N_PATCH_MAX



def calc_r_200(m_200: float = 1e9) -> float:
    """ Calculate R_200 based on M_200

    : m_200 : mass within R_200 in unit of Msun
    : returns : R_200 in unit of pc
    """
    delta_c = 200.
    rho_c = 140. * 1e-9    # Msun / pc^3
    return  np.power(3. * m_200 / 4. / np.pi / delta_c / rho_c, 1. / 3.)


def calc_width(r_200: float, dist: float, factor: int = 1) -> float:
    """ Calculate the width of the map based on distance and R_200

    : r_200 : R_200 in unit of pc
    : dist : distance of the dwaf from us in unit of pc
    : factor : decide how many times the width as to R_200
    : returns : width of the map in unit of degree
    """
    return  factor * r_200 / dist * 360. / 2. / np.pi


def plot_rh_sigma(dict_joint: Dict):
    """ Plotting half-light radius and different sigmas

    : dict_joint : dictionary of the dwarf list
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    rh_deg = dict_joint["rh(arcmins)"] / 60.
    s1 = 10. / dict_joint["Distance_pc"] * 180. / np.pi
    r2 = 2. * s1
    r3out = 0.5    # sigma 3: 0.5 deg
    r3in = 0.7 * rh_deg    # need to check the final value of 0.7

    masksort = np.argsort(dict_joint["GalaxyName"])

    sns.set(style="whitegrid", font_scale=1)
    f, ax = plt.subplots(figsize=(8, 12))
    sns.set_color_codes("muted")

    sns.barplot(x=rh_deg[masksort],
                y=np.sort(dict_joint["GalaxyName"]),
                label="rh", color="b", alpha=0.5)
    sns.barplot(x=r3in[masksort],
                y=np.sort(dict_joint["GalaxyName"]),
                label="r3in", color="r", alpha=0.5)
    sns.barplot(x=r2[masksort],
                y=np.sort(dict_joint["GalaxyName"]),
                label="r2", color="g", alpha=0.5)
    sns.barplot(x=10. * s1[masksort],
                y=np.sort(dict_joint["GalaxyName"]),
                label="10 s1", color="orange", alpha=0.5)

    ax.legend(ncol=2, loc="lower right", frameon=True)
    ax.set(xlabel="angular size [deg]")
    sns.despine(left=True, bottom=True)

    plt.savefig("dwarfs/rh_sigma.png", bbox_inches='tight', dpi=200)


def savetxt_gaia_pm_nstar(dict_joint):

    from wsdb import HOST, USER, PASSWORD
    from src.classPatchMWSatellite import PatchMWSatellite

    pmra, pmdec, nrh, nbg = [], [], [], []
    for i, name in enumerate(dict_joint['GalaxyName']):
        ra, dec = dict_joint['RA_deg'][i], dict_joint['Dec_deg'][i]
        dist = dict_joint['Distance_pc'][i]
        rh = dict_joint['rh(arcmins)'][i] / 60.
        width = 4. * rh
        database = 'gaia_dr2.gaia_source'
        cat_str = """ ra, dec, pmra, pmdec,
                      phot_g_mean_mag, astrometric_excess_noise """

        Patch = PatchMWSatellite(name, ra, dec, dist, width, database, cat_str)
        Patch.sql_get(HOST, USER, PASSWORD)
        Patch.mask_cut("phot_g_mean_mag", 17, 22)
        Patch.mask_g_mag_astro_noise_cut()

        mask2rh = dist2(
            Patch.datas['ra'], Patch.datas['dec'], ra, dec) <= (2.*rh)**2
        Patch.cut_datas(mask2rh)
        _n2rh = len(Patch.datas['pmra'])

        maskrh = dist2(Patch.datas['ra'], Patch.datas['dec'], ra, dec) <= rh**2
        Patch.cut_datas(maskrh)
        _nrh = len(Patch.datas['pmra'])
        _pmra = np.nanmean(Patch.datas['pmra'])
        _pmdec = np.nanmean(Patch.datas['pmdec'])

        nrh.append(_nrh)
        pmra.append(_pmra)
        pmdec.append(_pmdec)

        _den = (_n2rh - _nrh) / np.pi / 3. / rh**2
        nbg.append(_den)    # Nbg / deg^2

    dict_joint['Nstar_rh'] = np.array(nrh)
    dict_joint['pmra_rh'] = np.array(pmra)
    dict_joint['pmdec_rh'] = np.array(pmdec)
    dict_joint['Nbg_per_deg2'] = np.array(nbg)
    dict_joint['rh_deg'] = dict_joint['rh(arcmins)'] / 60.
    dict_joint['Nstar_per_rhdeg3'] = dict_joint['Nstar_rh'] / dict_joint['rh_deg']**3

    df = pd.DataFrame(dict_joint)
    df.to_csv('dwarfs/dwarfs_detail.csv')



if __name__ == '__main__':
    """ Concatenate the McConnachie list and my list into a joint one """
    path_mcconnachie = "dwarfs/McConnachie/dwarfs-McConnachie.npy"
    path_more = "dwarfs/more/dwarfs-more.npy"

    quantitys = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

    dict1 = get_dic_list_npy(path_mcconnachie, quantitys)
    dict2 = get_dic_list_npy(path_more, quantitys)

    dict_joint = {q: np.concatenate((dict1[q], dict2[q])) for q in quantitys}

    np.save("dwarfs/dwarfs-joint", dict_joint)
    np.savetxt("dwarfs/dwarfs-names.txt", np.sort(dict_joint["GalaxyName"]), fmt="%s")


    plot_rh_sigma(dict_joint)

    savetxt_gaia_pm_nstar(dict_joint)

    """ Expend the joint list by splitting the original map """
    r_200 = calc_r_200()    # R_200 in pc

    name_split = []
    ra_split, dec_split = [], []
    dist_split, rh_split = [], []
    ra_dwarf, dec_dwarf = [], []

    n_dwarf = len(dict_joint['GalaxyName'])
    for k in range(n_dwarf):
        name = dict_joint['GalaxyName'][k]
        ra = dict_joint['RA_deg'][k]
        dec = dict_joint['Dec_deg'][k]
        dist = dict_joint['Distance_pc'][k]
        rh = dict_joint['rh(arcmins)'][k]

        width = calc_width(r_200, dist)
        nmax_patch = np.ceil(0.5 * (width / PATCH_DIST - 1.))
        nmax_patch = int(min(nmax_patch, N_PATCH_MAX))

        id_ = 0
        for i in range(-nmax_patch, nmax_patch + 1):
            for j in range(-nmax_patch, nmax_patch + 1):
                if (i**2 + j**2) <= nmax_patch**2:
                    xi = ra + i * PATCH_DIST
                    yj = dec + j * PATCH_DIST
                    id_ += 1
                    name_split.append('{}=%d'.format(name) %id_)
                    ra_split.append(xi)
                    dec_split.append(yj)
                    dist_split.append(dist)
                    rh_split.append(rh)
                    ra_dwarf.append(ra)
                    dec_dwarf.append(dec)

    list_split = [name_split, ra_split, dec_split, dist_split, rh_split,
                  ra_dwarf, dec_dwarf]

    quantitys.append("RA_dwarf_deg")
    quantitys.append("Dec_dwarf_deg")

    dict_joint = {q: np.array(list_split[i]) for i, q in enumerate(quantitys)}

    sorted_name = np.sort(dict_joint["GalaxyName"])

    np.save("dwarfs/dwarfs-joint-split", dict_joint)
    np.savetxt("dwarfs/dwarfs-names-split.txt", sorted_name, fmt="%s")
