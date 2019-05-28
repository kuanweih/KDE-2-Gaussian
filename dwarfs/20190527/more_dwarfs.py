import numpy as np



if __name__ == '__main__':

    dwarfs = []

    quantitys = ["GalaxyName", "RA_deg", "Dec_deg", "Distance_pc", "rh(arcmins)"]

    # https://arxiv.org/pdf/1601.07178.pdf
    dwarfs.append(['Crater2', 177.31, -18.413, 117500., 31.2])

    # https://iopscience.iop.org/article/10.1088/0004-637X/805/2/130/pdf
    dwarfs.append(['Reticulum2', 53.9256, -54.0492, 30000., 3.64])
    dwarfs.append(['Eridanus2', 56.0878, -43.5332, 380000., 1.53])
    dwarfs.append(['Horologium1', 43.8820, -54.1188, 79000., 1.31])
    dwarfs.append(['Pictoris1', 70.9475, -50.2830, 114000., 0.88])
    dwarfs.append(['Phoenix2',354.9975, -54.4060, 83000., 1.09])
    dwarfs.append(['Indus1', 317.2044, -51.1656, 100000., 1.26])
    dwarfs.append(['Grus1', 344.1765, -50.1633, 120000., 1.77])
    dwarfs.append(['Eridanus3', 35.6897, -52.2837, 87000., 0.54])
    dwarfs.append(['Tucana2', 342.9796, -58.5689, 57000., 9.83])

    # https://arxiv.org/pdf/1811.04082.pdf
    dwarfs.append(['Antlia2', 143.8868, -36.7673, 129400., 75.6])

    # https://arxiv.org/pdf/1801.07279.pdf
    dwarfs.append(['Carina2', 114.1066, -57.9991, 36200., 8.69])
    dwarfs.append(['Carina3', 114.6298, -57.8997, 27800., 3.75])

    # http://cdsads.u-strasbg.fr/abs/2018PASJ...70S..18H
    dwarfs.append(['Carina2', 31.331, -4.270, 251000., 1.23])

    # https://arxiv.org/pdf/1605.05338.pdf
    dwarfs.append(['Aquarius2', 338.4813, -9.3274, 107900., 5.1])

    # https://iopscience.iop.org/article/10.3847/0004-637X/832/1/21/pdf
    dwarfs.append(['Virgo1', 180.04, -0.68, 87000., 1.5])

    dwarfs = np.array(dwarfs)

    dwarfs_dic = {}
    for i, q in enumerate(quantitys):
        dwarfs_dic[q] = dwarfs[:, i] if i==0 else dwarfs[:, i].astype(float)

    np.save("dwarfs-more", dwarfs_dic)
    np.savetxt("dwarfs-names.txt", dwarfs_dic["GalaxyName"], fmt="%s")











#
