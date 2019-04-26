
import numpy as np
import glob
import pandas as pd


def create_dir(dir_name: str):
    """ create directory with a name 'dir_name' """
    import os
    import errno
    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

# MODE = "create_csv"
MODE = "read_csv"


if __name__ == '__main__':

    if MODE == "create_csv":
        paths = glob.glob('results/*')

        names, gc_sizes, n_stars, n_star_pms = [], [], [], []
        sig_g_peaks, sig_g_pm_peaks, sig_p_peaks, sig_p_pm_peaks = [], [], [], []

        for path in paths:
            s_above = 5
            name = path.split("-")[0].replace("results/", "")
            gc_size = int(path.split("-")[5].split("s")[0].replace("gc", ""))

            n_star = len(np.load("{}/queried-data.npy".format(path)).item()['ra'])
            n_star_pm = len(np.load("{}/queried-data-pm_error5.npy".format(path)).item()['ra'])


            sig_g = np.load("{}/sig_gaussian.npy".format(path))
            sig_g_pm = np.load("{}/sig_gaussian-pm_error5.npy".format(path))

            sig_p = np.load("{}/sig_poisson.npy".format(path))
            sig_p_pm = np.load("{}/sig_poisson-pm_error5.npy".format(path))

            sig_g_peak = np.sum(sig_g > s_above)
            sig_g_pm_peak = np.sum(sig_g_pm > s_above)

            sig_p_peak = np.sum(sig_p > s_above)
            sig_p_pm_peak = np.sum(sig_p_pm > s_above)

            names.append(name)
            gc_sizes.append(gc_size)
            n_stars.append(n_star)
            n_star_pms.append(n_star_pm)
            sig_g_peaks.append(sig_g_peak)
            sig_g_pm_peaks.append(sig_g_pm_peak)
            sig_p_peaks.append(sig_p_peak)
            sig_p_pm_peaks.append(sig_p_pm_peak)

        searching_table = {}
        searching_table["name"] = names
        searching_table["gc_size"] = gc_sizes
        searching_table["n_star"] = n_stars
        searching_table["n_star_pm"] = n_star_pms
        searching_table["sig_g_peak"] = sig_g_peaks
        searching_table["sig_g_pm_peak"] = sig_g_pm_peaks
        searching_table["sig_p_peak"] = sig_p_peaks
        searching_table["sig_p_pm_peak"] = sig_p_pm_peaks

        df = pd.DataFrame(data=searching_table)
        df = df[["name", "gc_size", "n_star", "n_star_pm", "sig_g_peak", "sig_g_pm_peak", "sig_p_peak", "sig_p_pm_peak"]]

        df.to_csv("summary.csv", index=False)



    if MODE == "read_csv":
        # create dir for searching results
        path = "tables"
        create_dir(path)

        df = pd.read_csv("summary.csv")

        pixel_th = 12
        n_star_th = 300
        gc_sizes = [6, 7, 8, 9, 10]

        target_dwarfs = []

        for gc_size in gc_sizes:
            mask = df["sig_p_peak"] > pixel_th
            mask = (df["sig_p_pm_peak"] > pixel_th) | mask

            mask = (df["n_star"] > n_star_th) & mask
            mask = (df["n_star_pm"] > n_star_th) & mask

            mask = (df["gc_size"] == gc_size) & mask

            df[mask].to_csv("{}/gc{}.csv".format(path, gc_size), index=False)

            for name in df[mask]["name"]:
                if name not in target_dwarfs:
                    target_dwarfs.append(name)

        for target_dwarf in target_dwarfs:
            print(target_dwarf)

        for target_dwarf in target_dwarfs:
            mask = df["sig_p_peak"] > pixel_th
            mask = (df["sig_p_pm_peak"] > pixel_th) | mask

            mask = (df["n_star"] > n_star_th) & mask
            mask = (df["n_star_pm"] > n_star_th) & mask

            mask = (df["name"] == target_dwarf) & mask

            df[mask].to_csv("{}/{}.csv".format(path, target_dwarf), index=False)



        # print(df[mask])















    #
