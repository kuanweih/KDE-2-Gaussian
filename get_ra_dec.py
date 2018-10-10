import numpy as np
import sqlutilpy


def con_arr_astro_ex_noise_g_mag(astro_ex_noise, g_mag):
    """
    simplify eq 1 in Koposov et al 2017 (MNRAS 470) into exponential form
    """
    return astro_ex_noise - 10.**(0.15 * (g_mag - 15.) + 0.25)


def cut_datas(datas, con_arr, min_val, max_val):
    """
    apply condition to the queried data
    """
    if min_val != None:
        if max_val != None:
            mask = (min_val < con_arr) & (con_arr < max_val)
        else:
            mask = (min_val < con_arr)
    elif max_val != None:
        mask = (con_arr < max_val)
    else:
        print('Oops, no input minimum or maximum value here.')
    return [data[mask] for data in datas]


def main():
    # files names
    PERSONAL_FILE = 'kw-wsdb.txt'    # text file of input personal info
    # PARAMETER_TXT = 'param-get.txt'    # text file for parameters of database
    PARAMETER_FILE = 'param-get.py'    # text file for parameters of database
    FILENAME = 'stars-coord'    # output file name
    INFOFILE = 'stars-coord-attr'    # info of the map

    # HOST, USER, PASSWORD from personal txt file
    with open(PERSONAL_FILE, 'r') as info_file:
        for line in info_file:
            exec(line)

    # RA, DEC, RADIUS, DATABASE, CATALOGS etc... from parameter file
    with open(PARAMETER_FILE, 'r') as param_file:
        for line in param_file:
            exec(line)

    # query data from DATABASE
    query_str = 'select {} from {} where q3c_radial_query(ra, dec, {}, {}, {})'.format(
                CATALOGS, DATABASE, RA, DEC, RADIUS)
    datas = sqlutilpy.get(query_str, host=HOST, user=USER, password=PASSWORD)

    if G_MAG_CUT:
        datas = cut_datas(datas, datas[2], G_MAG_MIN, G_MAG_MAX)

    if ASTRO_EX_NOISE_G_MAG_CUT:
        datas = cut_datas(datas, con_arr_astro_ex_noise_g_mag(datas[3], datas[2]),
                          ASTRO_EX_NOISE_G_MAG_MIN, ASTRO_EX_NOISE_G_MAG_MAX)

    # TODO change index to correspond with the CATALOGS string
    ra = datas[0]
    dec = datas[1]

    # output
    np.save(FILENAME, np.array([ra, dec]))
    np.save(INFOFILE, np.array([RA, DEC, RADIUS, len(ra)]))

    print('Yeah! Done!')


if __name__ == '__main__':
    main()
