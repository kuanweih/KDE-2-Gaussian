import numpy as np
import pandas as pd

from scipy.ndimage import label as snlabel



def summarize_peaks_star_csv(path: str, outfile: str, n_error: float,
                             kernel: str, s_above=5):
    """ Plotting star distribution (left panels) and density maps (right
    panels). Top row: all stars. Bottom row: pm selection.

    : path : path of all results files
    : outfile : result file of the dwarf
    : n_error : N_ERRORBAR for pm selection
    : kernel : 'gaussian' or 'poisson'
    : s_above : significance threshold, default value = 5
    """
    datas = [np.load('{}/queried-data.npy'.format(path)).item(),
             np.load('{}/queried-data-pm_error{}.npy'.format(path,
                                                             n_error)).item()]

    ras = [data["ra"] for data in datas]
    decs = [data["dec"] for data in datas]
    sigs = [data["sig_{}".format(kernel)] for data in datas]

    # masks for stars with sig > s_above
    masks = [sig > s_above for sig in sigs]

    for v in range(2):
        assert (v < 2), ("wrong index in the for loop for all stars and pm selection. ")

        ra_peaks = ras[v][masks[v]]
        dec_peaks = decs[v][masks[v]]
        sig_peaks = sigs[v][masks[v]]

        n_star = len(sig_peaks)

        peaks_table = {}

        if v==0:    # all stars
            name = "{}-{}".format(outfile, kernel)
        else:    # pm selection
            name = "{}-{}-pm".format(outfile, kernel)

        peaks_table["name"] = np.array([name.replace("peaks/stars/", "")] * n_star)

        peaks_table["ra"] = ra_peaks
        peaks_table["dec"] = dec_peaks
        peaks_table["sig"] = sig_peaks

        df = pd.DataFrame(data=peaks_table)
        df = df[["name", "sig", "ra", "dec"]]

        df.to_csv("{}.csv".format(name), index=False)


def summarize_peaks_pixel_csv(path: str, outfile: str, n_error: float,
                              kernel: str, sat_ra: float, sat_dec: float,
                              width: float, s_above=5):
    """ Plotting star distribution (left panels) and density maps (right
    panels). Top row: all stars. Bottom row: pm selection.

    : path : path of all results files
    : outfile : result file of the dwarf
    : n_error : N_ERRORBAR for pm selection
    : kernel : 'gaussian' or 'poisson'
    : sat_ra : ra of the dwarf
    : sat_dec : dec of the dwarf
    : width : width of the map
    : s_above : significance threshold, default value = 5
    """
    x, y = np.load('{}/meshgrids.npy'.format(path))
    meshgrids = np.meshgrid(x, y)

    sigs = [np.load('{}/sig_{}.npy'.format(path, kernel)),
            np.load('{}/sig_{}-pm_error{}.npy'.format(path, kernel, n_error))]

    masks = [sig > s_above for sig in sigs]

    for v in range(2):
        assert (v < 2), ("wrong index in the for loop for all stars and pm selection. ")

        if v==0:    # all stars
            name = "{}-{}".format(outfile, kernel)
        else:    # pm selection
            name = "{}-{}-pm".format(outfile, kernel)

        labeled_array, num_features = snlabel(masks[v])

        label_s, x_peaks, y_peaks = [], [], []

        for nf in range(1, num_features + 1):    # ignore pixels with 0
            # mask_id = labeled_array == nf
            peak_xid, peak_yid = np.where(labeled_array == nf)
            label_ = [nf] * len(peak_xid)
            label_s.append(label_)
            x_peaks.append(x[peak_xid])
            y_peaks.append(y[peak_yid])

        print(np.array(label_s).flatten())

        pixel_peaks_table = {}

        pixel_peaks_table["name"] = np.array([name.replace("peaks/pixels/", "")] * len(label_s))

        pixel_peaks_table["label"] = label_s
        pixel_peaks_table["ra"] = x_peaks
        pixel_peaks_table["dec"] = y_peaks

        df = pd.DataFrame(data=pixel_peaks_table)
        df = df[["name", "label", "ra", "dec"]]

        df.to_csv("{}.csv".format(name), index=False)




























    #
