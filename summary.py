import numpy as np
import glob
import pandas as pd
import multiprocessing

from functools import partial
from src.tools import create_dir, df_concat
from src.hips_image import multiprocessing_plot_hips_sky_image


""" parameters """
s_above = 5    # significance threshold

hips_surveys = ['CDS/P/DES-DR1/g',
                'CDS/P/DECaLS/DR5/color',
                'CDS/P/DSS2/color',
                'CDS/P/SDSS9/color',
                'CDS/P/PanSTARRS/DR1/color-z-zg-g',
                'CDS/P/2MASS/color',
                ]

res_image = 1000    # image resolution of hips
width_fac = 10.    # width of image = width_fac * sigma1

num_workers = multiprocessing.cpu_count()
# num_workers = int(0.5 * num_workers)



if __name__ == '__main__':
    output_file = 'summary'
    create_dir(output_file)

    # print('Concatenating all peaks/stars/csv files...\n')
    # paths = glob.glob('peaks/stars/*')
    # df_con = df_concat(paths)
    # df_con.to_csv("{}/all_stars_peaks.csv".format(output_file), index=False)
    # print('Done :) \n')

    print('Concatenating all peaks/pixels/csv files\n')
    paths = glob.glob('peaks/pixels/*')
    df_con = df_concat(paths)
    df_con.to_csv("{}/all_pixels_peaks.csv".format(output_file), index=False)
    print('Done :) \n')


    print('Generating summary.csv: \n')
    paths = glob.glob('results/*')

    names, gc_sizes, n_stars = [], [], []
    sig_g_peaks, sig_p_peaks = [], []

    for path in paths:
        name = path.split("-")[0].replace("results/", "")
        if 'gaia' in name:
            gc_size = int(path.split("-")[5].split("s")[0].replace("gc", ""))
        else:
            gc_size = int(path.split("-")[3].split("s")[0].replace("gc", ""))

        try:
            n_star = len(np.load("{}/queried-data.npy".format(path)).item()['ra'])
            sig_g = np.load("{}/sig_gaussian.npy".format(path))
            sig_p = np.load("{}/sig_poisson.npy".format(path))
        except:
            print('    skipping Error csv: %s' %name)
            continue

        names.append(name)
        gc_sizes.append(gc_size)
        n_stars.append(n_star)
        sig_g_peaks.append(np.sum(sig_g > s_above))
        sig_p_peaks.append(np.sum(sig_p > s_above))

    searching_table = {}
    searching_table["name"] = names
    searching_table["gc_size"] = gc_sizes
    searching_table["n_star"] = n_stars
    searching_table["sig_g_peak"] = sig_g_peaks
    searching_table["sig_p_peak"] = sig_p_peaks

    df = pd.DataFrame(data=searching_table)
    df = df[["name", "gc_size", "n_star", "sig_g_peak", "sig_p_peak"]]
    df.to_csv("{}/summary.csv".format(output_file), index=False)
    print('Done :) \n')


    print('Creating images via hips: \n')
    path = "images"
    create_dir(path)

    name_df, label_df, ra_df, dec_df, width_df = [], [], [], [], []
    df = pd.read_csv("{}/all_pixels_peaks.csv".format(output_file))

    name_set = set(df['name'])
    for name in name_set:
        if "gaussian" in name:
            continue    # don't plot images of gaussian kernel
        df_name = df.loc[df['name'] == name]

        label_set = set(df_name['label'])
        for label in label_set:
            df_name_label = df_name.loc[df_name['label'] == label]
            ra = np.mean(df_name_label['ra'])
            dec = np.mean(df_name_label['dec'])
            name_list = name.split("-")
            if 'gaia' in name:
                width = width_fac * float(name_list[5].split('s')[1])
            else:
                width = width_fac * float(name_list[3].split('s')[1])

            name_df.append(name)
            label_df.append(label)
            ra_df.append(ra)
            dec_df.append(dec)
            width_df.append(width)

    df = None    # free memory

    num_target = len(name_df)
    iterable = np.arange(num_target)
    print('There are %d candidates\n' %num_target)

    pool = multiprocessing.Pool(num_workers)
    func = partial(multiprocessing_plot_hips_sky_image, name_df, label_df,
                   hips_surveys, ra_df, dec_df, width_df, path, res_image)
    pool.map(func, iterable)
    pool.close()
    pool.join()

    print("\nWe are finished :) \n")
