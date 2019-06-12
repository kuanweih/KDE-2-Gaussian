import numpy as np

PATCH_DIST = 0.9
N_PATCH_MAX = 4


def get_dic_list(path: str, quantitys: str):
    """ Select dict of a dwaft list based on keys as quantitys

    : path : path of the dict npy file
    : quantitys : target keys
    : return : needed dict
    """
    dwarfs_dict = np.load(path).item()
    dict_need = {q: dwarfs_dict[q] for q in quantitys}
    return  dict_need


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



if __name__ == '__main__':
    """ Concatenate the McConnachie list and my list into a joint one """
    path_mcconnachie = "McConnachie/dwarfs-McConnachie.npy"
    path_more = "more/dwarfs-more.npy"

    quantitys = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

    dict1 = get_dic_list(path_mcconnachie, quantitys)
    dict2 = get_dic_list(path_more, quantitys)

    dict_joint = {q: np.concatenate((dict1[q], dict2[q])) for q in quantitys}

    np.save("dwarfs-joint", dict_joint)
    np.savetxt("dwarfs-names.txt", np.sort(dict_joint["GalaxyName"]), fmt="%s")


    """ Expend the joint list by splitting the original map """
    r_200 = calc_r_200()    # R_200 in pc

    name_split = []
    ra_split, dec_split = [], []
    dist_split, rh_split = [], []

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

    list_split = [name_split, ra_split, dec_split, dist_split, rh_split]
    dict_joint = {q: np.array(list_split[i]) for i, q in enumerate(quantitys)}

    np.save("dwarfs-joint-split", dict_joint)
    np.savetxt("dwarfs-names-split.txt", np.sort(dict_joint["GalaxyName"]), fmt="%s")
