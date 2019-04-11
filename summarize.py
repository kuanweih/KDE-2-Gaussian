
import numpy as np
import glob


if __name__ == '__main__':
    paths = glob.glob('results-gaussian/*')

    # sigs = [np.load("{}/significance.npy".format(path)) for path in paths]
    # sig_pms = [np.load("{}/significance-pm_error5.npy".format(path)) for path in paths]

    # names = [path.split("-")[1].replace("gaussian/", "") for path in paths]
    # print(names)

    for path in paths:
        s_above = 5
        name = path.split("-")[1].replace("gaussian/", "")
        sig = np.load("{}/significance.npy".format(path))
        sig_pm = np.load("{}/significance-pm_error5.npy".format(path))

        sig_peaks = np.sum(sig > s_above)
        sig_pm_peaks = np.sum(sig > s_above)

        print(name, sig_peaks, sig_pm_peaks)
