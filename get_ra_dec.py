import numpy as np
import sqlutilpy


# files names
PERSONAL_TXT = 'kw-wsdb.txt'    # text file of input personal info
PARAMETER_TXT = 'param-get.txt'    # text file for parameters of database
FILENAME = 'stars-coord'    # output file name
INFOFILE = 'stars-coord-attr'    # info of the map


# TODO: to be moved to paramfile
G_MAG_CUT = 0  # condition for g band mag cut
G_MAG_MIN = 17
G_MAG_MAX = 21
G_MAG_ASTRO_EX_NOISE_CUT = 1    # condition for g_mag dependent cut on noise


def cut_datas(datas, con_arr, min_val, max_val):
    """ apply condition to the queried data"""
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


# HOST, USER, PASSWORD from personal txt file
with open(PERSONAL_TXT, 'r') as info_file:
    for line in info_file:
        exec(line)


# RA, DEC, RADIUS, DATABASE, CATALOGS from parameter file
with open(PARAMETER_TXT, 'r') as param_file:
    for line in param_file:
        exec(line)


# query data from DATABASE
query_str = 'select {} from {} where q3c_radial_query(ra, dec, {}, {}, {})'.format(
            CATALOGS, DATABASE, RA, DEC, RADIUS)
datas = sqlutilpy.get(query_str, host=HOST, user=USER, password=PASSWORD)


if G_MAG_CUT == 1:
    datas = cut_datas(datas, datas[2], G_MAG_MIN, G_MAG_MAX)


# TODO change index to correspond with the CATALOGS string
ra = datas[0]
dec = datas[1]

# output
np.save(FILENAME, np.array([ra, dec]))
np.save(INFOFILE, np.array([RA, DEC, RADIUS, len(ra)]))
