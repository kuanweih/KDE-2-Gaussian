import numpy as np



def split_into_3_float(input_arr):
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
    return output


def get_ra_deg(hour, minute, second):
    return (hour + minute / 60. + second / 3600.) * 15.


def get_dec_deg(degree, minute, second):
    deg_mag = np.absolute(degree)
    deg_sign = degree / deg_mag
    return  deg_sign * (deg_mag + minute / 60. + second / 3600.)




f = open('McConnachie.txt', 'r')
for line in f:
    if "GalaxyName" in line:
        name_cols = line
f.close()


quantitys = ["GalaxyName         ",
             "RA         ",
             "Dec       ",
             "(m-M)o          ",
             "e=1-b/a        ",
             "rh(arcmins)        "]

idstarts = [name_cols.find(q) for q in quantitys]
idoffsets = [len(q) for q in quantitys]


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


keys = [q.replace(" ", "") for q in quantitys]

dwarfs = {}

for i, q in enumerate(keys):
    dwarfs[q] = np.array(galaxys)[:, i]

for key, value in dwarfs.items():
    if key in keys[1:]:
        dwarfs[key] = split_into_3_float(value)

dwarfs["RA_deg"] = np.array(list(map(lambda x: get_ra_deg(x[0], x[1], x[2]),
                                     dwarfs["RA"])))
dwarfs["Dec_deg"] = np.array(list(map(lambda x: get_dec_deg(x[0], x[1], x[2]),
                                      dwarfs["Dec"])))


np.save("dwarfs-McConnachie", dwarfs)

# name_cols = name_cols.replace("   ", " ").replace("  ", " ")
# name_cols = name_cols.replace("  ", " ").replace("  ", " ").replace("\n", "")
# name_cols = name_cols.split(" ")
#
# print(len(name_cols))
# print(name_cols)

# print(galaxys)
