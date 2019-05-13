import numpy as np
from typing import List


def get_galaxys_with_quantitys(quantitys: List[str]) -> List[str]:
    """ get galaxys line by line with quantitys """
    is_print = False
    galaxys = []
    f = open('McConnachie.txt', 'r')
    for line in f:
        if is_print:
            list_line = []
            for i in range(len(quantitys)):
                quantity = line[idstarts[i]:idstarts[i] + idoffsets[i]]
                quantity = quantity.replace("  ", " ").replace("  ", " ")
                list_line.append(quantity)
            galaxys.append(list_line)
        if "0123456789" in line:
            is_print = True
    return  galaxys


def split_into_3_float(input_arr) -> List[float]:
    """ split input numpy array into a list with 3 float-elements

    : input_arr : numpy array
    """
    output = []
    for ii in input_arr:
        ii = ii.split(" ")
        is_remove = True
        while is_remove:
            if "" in ii:
                ii.remove("")
            else:
                is_remove = False
        output.append([float(i) for i in ii])
    return  output


def get_ra_deg(hour: float, minute: float, second: float) -> float:
    """ convert RA from (hr, min, sec) -> degree """
    return  (hour + minute / 60. + second / 3600.) * 15.


def get_dec_deg(degree: float, minute: float, second: float) -> float:
    """ convert Dec from (deg, arcmin, arcsec) -> degree """
    deg_sign = 1. if (degree + np.absolute(degree)) > 0 else - 1.
    return  deg_sign * (np.absolute(degree) + minute / 60. + second / 3600.)


def dist_modulus_to_dist(dm: float, dm_errr: float, dm_errl: float) -> float:
    """ convert distance modulus -> distance (pc)
        currently not taking care of error bars though
        
    : dm : distance modulus
    : dm_errr : dm error +
    : dm_errl : dm error -
    """
    return  10.**(dm / 5. + 1.)



if __name__ == '__main__':
    # get name_cols from txt file
    f = open('McConnachie.txt', 'r')
    for line in f:
        if "GalaxyName" in line:
            name_cols = line
    f.close()

    # targeted quantitys
    quantitys = ["GalaxyName         ",
                 "RA         ",
                 "Dec       ",
                 "(m-M)o          ",
                 "e=1-b/a        ",
                 "rh(arcmins)        "]

    idstarts = [name_cols.find(q) for q in quantitys]
    idoffsets = [len(q) for q in quantitys]

    galaxys = get_galaxys_with_quantitys(quantitys)

    keys = [q.replace(" ", "") for q in quantitys]

    dwarfs = {}

    for i, q in enumerate(keys):
        dwarfs[q] = np.array(galaxys)[:, i]

    for key, value in dwarfs.items():
        if key in keys[1:]:
            dwarfs[key] = np.array(split_into_3_float(value))

    dwarfs["GalaxyName"] = np.array(list(map(lambda x: x.replace(" ", ""), dwarfs["GalaxyName"])))
    dwarfs["RA_deg"] = np.array(list(map(lambda x: get_ra_deg(x[0], x[1], x[2]), dwarfs["RA"])))
    dwarfs["Dec_deg"] = np.array(list(map(lambda x: get_dec_deg(x[0], x[1], x[2]), dwarfs["Dec"])))
    dwarfs["Distance_pc"] = np.array(list(map(lambda x: dist_modulus_to_dist(x[0], x[1], x[2]), dwarfs["(m-M)o"])))
    dwarfs["rh(arcmins)"] = dwarfs["rh(arcmins)"][:, 0]

    mask = dwarfs["Distance_pc"] < 600000    # skip dwarfs too far
    mask = mask & (dwarfs["rh(arcmins)"] < 90)    # skip dwarfs with rh = 99.99

    np.save("dwarfs-McConnachie", dwarfs)
    np.savetxt("dwarfs-names.txt", dwarfs["GalaxyName"][mask], fmt="%s")
    np.savetxt("skip_dwarfs.txt", dwarfs["GalaxyName"][~mask], fmt="%s")
