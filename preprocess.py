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



def expand_joint_dict_split(df: pd.DataFrame, is_pm: bool = False):
    """ Expand joint list to a splitted one.

    : df : DataFrame of csv file for the dwarfs
    : is_pm : specify if it is the original csv or the one with proper motions
    """
    quantitys = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

    dict = {}
    for column in df:
        dict[column] = df[column].values

    r_200 = calc_r_200()    # R_200 in pc

    name_split = []
    ra_split, dec_split = [], []
    dist_split, rh_split = [], []
    ra_dwarf, dec_dwarf = [], []

    if is_pm:
        pmra_dwarf, pmdec_dwarf = [], []

    n_dwarf = len(dict['GalaxyName'])
    for k in range(n_dwarf):
        name = dict['GalaxyName'][k]
        ra = dict['RA_deg'][k]
        dec = dict['Dec_deg'][k]
        dist = dict['Distance_pc'][k]
        rh = dict['rh(arcmins)'][k]
        if is_pm:
            pmra = dict['pmra'][k]
            pmdec = dict['pmdec'][k]

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
                    if is_pm:
                        pmra_dwarf.append(pmra)
                        pmdec_dwarf.append(pmdec)

    list_split = [name_split, ra_split, dec_split, dist_split, rh_split,
                  ra_dwarf, dec_dwarf]
    if is_pm:
        list_split.append(pmra_dwarf)
        list_split.append(pmdec_dwarf)


    quantitys.append("RA_dwarf_deg")
    quantitys.append("Dec_dwarf_deg")
    if is_pm:
        quantitys.append('pmra_dwarf')
        quantitys.append('pmdec_dwarf')

    dict = {q: np.array(list_split[i]) for i, q in enumerate(quantitys)}

    sorted_name = np.sort(dict["GalaxyName"])

    if is_pm:
        np.save("dwarfs/dwarfs-joint-split-pm", dict)
        np.savetxt("dwarfs/dwarfs-names-split-pm.txt", sorted_name, fmt="%s")
    else:
        np.save("dwarfs/dwarfs-joint-split", dict)
        np.savetxt("dwarfs/dwarfs-names-split.txt", sorted_name, fmt="%s")



def expand_joint_dict(df: pd.DataFrame, is_pm: bool = False):
    """ Expand joint list to a splitted one.

    : df : DataFrame of csv file for the dwarfs
    : is_pm : specify if it is the original csv or the one with proper motions
    """
    df['RA_dwarf_deg'], df['Dec_dwarf_deg'] = df['RA_deg'], df['Dec_deg']

    quantitys = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]
    quantitys.append("RA_dwarf_deg")
    quantitys.append("Dec_dwarf_deg")
    if is_pm:
        quantitys.append('pmra')
        quantitys.append('pmdec')

    dict = {q: df[q].values for q in quantitys}
    if is_pm:
        dict['pmra_dwarf'] = dict.pop('pmra')
        dict['pmdec_dwarf'] = dict.pop('pmdec')

    sorted_name = np.sort(dict["GalaxyName"])

    if is_pm:
        np.save("dwarfs/dwarfs-joint-pm", dict)
        np.savetxt("dwarfs/dwarfs-names-pm.txt", sorted_name, fmt="%s")
    else:
        np.save("dwarfs/dwarfs-joint", dict)
        np.savetxt("dwarfs/dwarfs-names.txt", sorted_name, fmt="%s")



if __name__ == '__main__':

    df_ori = pd.read_csv('dwarfs/ori-dwarfs.csv').drop(['Unnamed: 0'], axis=1)
    df_ori_pm = pd.read_csv('dwarfs/ori-dwarfs-pms.csv')

    expand_joint_dict(df_ori)
    expand_joint_dict(df_ori_pm, is_pm=True)

    expand_joint_dict_split(df_ori)
    expand_joint_dict_split(df_ori_pm, is_pm=True)
