import numpy as np
import sqlutilpy


# files names
PERSONAL_TXT = 'kw-wsdb.txt'    # text file of input personal info
PARAMETER_TXT = 'param-get.txt'    # text file for parameters of database
FILENAME = 'stars-coord-fornax'    # output file name
INFOFILE = 'stars-coord-attr'    # info of the map


# HOST, USER, PASSWORD from personal txt file
with open(PERSONAL_TXT, 'r') as info_file:
    for line in info_file:
        exec(line.readline())


# RA, DEC, RADIUS, DATABASE, CATALOGS from parameter file
with open(PARAMETER_TXT, 'r') as param_file:
    exec(param_file.readline())
    exec(param_file.readline())
    exec(param_file.readline())
    exec(param_file.readline())
    exec(param_file.readline())


# query data from DATABASE
query_str = 'select {} from {} where q3c_radial_query(ra, dec, {}, {}, {})'.format(
            CATALOGS, DATABASE, RA, DEC, RADIUS)
query_data = sqlutilpy.get(query_str, host=HOST, user=USER, password=PASSWORD)

ra = query_data[0]
dec = query_data[1]

# output
np.save(FILENAME, np.array([ra, dec]))
np.save(INFOFILE, np.array([RA, DEC, RADIUS, len(ra)]))
