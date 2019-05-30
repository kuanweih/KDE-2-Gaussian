import numpy as np



def get_dic_list(path: str, quantitys: str):
    """ Select dict of a dwaft list based on keys as quantitys

    : path : path of the dict npy file
    : quantitys : target keys

    : return : needed dict
    """
    dwarfs_dict = np.load(path).item()
    dict_need = {q: dwarfs_dict[q] for q in quantitys}
    return  dict_need



if __name__ == '__main__':

    path_mcconnachie = "McConnachie/dwarfs-McConnachie.npy"
    path_more = "more/dwarfs-more.npy"

    quantitys = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

    dict1 = get_dic_list(path_mcconnachie, quantitys)
    dict2 = get_dic_list(path_more, quantitys)

    dict_joint = {q: np.concatenate((dict1[q], dict2[q])) for q in quantitys}

    np.save("dwarfs-joint", dict_joint)
    np.savetxt("dwarfs-names.txt", np.sort(dict_joint["GalaxyName"]), fmt="%s")
