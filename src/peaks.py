import numpy as np
import pandas as pd

from scipy.ndimage import label as snlabel



def summarize_peaks_star_csv(path: str, outfile: str, kernel: str, s_above=5):
    """ Plotting star distribution (left panels) and density maps (right
    panels). (Others)

    : path : path of all results files
    : outfile : result file of the dwarf
    : kernel : 'gaussian' or 'poisson'
    : s_above : significance threshold, default value = 5
    """
    data = np.load('{}/queried-data.npy'.format(path)).item()

    ra, dec = data["ra"], data["dec"]
    sig = data["sig_{}".format(kernel)]

    # masks for stars with sig > s_above
    mask = sig > s_above

    ra_peaks = ra[mask]
    dec_peaks = dec[mask]
    sig_peaks = sig[mask]
    n_star = len(sig_peaks)

    peaks_table = {}

    name = "{}-{}".format(outfile, kernel)
    peaks_table["name"] = np.array([name.replace("peaks/stars/", "")] * n_star)

    peaks_table["ra"] = ra_peaks
    peaks_table["dec"] = dec_peaks
    peaks_table["sig"] = sig_peaks

    df = pd.DataFrame(data=peaks_table)
    df = df[["name", "sig", "ra", "dec"]]

    df.to_csv("{}.csv".format(name), index=False)



def summarize_peaks_pixel_csv(path: str, outfile: str,
                              kernel: str, sat_ra: float, sat_dec: float,
                              width: float, s_above=5):
    """ Plotting star distribution (left panels) and density maps (right
    panels). (Others)

    : path : path of all results files
    : outfile : result file of the dwarf
    : kernel : 'gaussian' or 'poisson'
    : sat_ra : ra of the dwarf
    : sat_dec : dec of the dwarf
    : width : width of the map
    : s_above : significance threshold, default value = 5
    """
    x, y = np.load('{}/meshgrids.npy'.format(path))
    sig = np.load('{}/sig_{}.npy'.format(path, kernel))

    mask = sig > s_above

    name = "{}-{}".format(outfile, kernel)

    labeled_array, num_features = snlabel(mask, structure=np.ones((3, 3)))

    if num_features == 1:    # skip the case without any peaks
        return  None

    label_s, x_peaks, y_peaks = [], [], []
    for nf in range(1, num_features + 1):    # ignore pixels with 0
        peak_yid, peak_xid = np.where(labeled_array == nf)    # x <-> y
        label_ = [nf] * len(peak_xid)
        label_s.append(label_)
        x_peaks.append(x[peak_xid])
        y_peaks.append(y[peak_yid])

    if num_features > 1:
        label_s = np.concatenate(label_s)
        x_peaks = np.concatenate(x_peaks)
        y_peaks = np.concatenate(y_peaks)
    else:
        label_s = np.array(label_s)
        x_peaks = np.array(x_peaks)
        y_peaks = np.array(y_peaks)

    pixel_peaks_table = {}

    pixel_peaks_table["name"] = np.array([name.replace("peaks/pixels/", "")] * len(label_s))
    pixel_peaks_table["label"] = label_s
    pixel_peaks_table["ra"] = x_peaks
    pixel_peaks_table["dec"] = y_peaks

    df = pd.DataFrame(data=pixel_peaks_table)
    df = df[["name", "label", "ra", "dec"]]

    df.to_csv("{}.csv".format(name), index=False)


























    #
