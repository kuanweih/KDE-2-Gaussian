import numpy as np
import pandas as pd



def summarize_peaks_csv(path: str, outfile: str, n_error: float, kernel: str,
                        s_above=5):
    """ Plotting star distribution (left panels) and density maps (right
    panels). Top row: all stars. Bottom row: pm selection.

    : path : path of all results files
    : outfile : result file of the dwarf
    : n_error : N_ERRORBAR for pm selection
    : kernel : 'gaussian' or 'poisson'
    : s_above : significance threshold, default value = 5
    """
    # _name = path.replace("results/", "")

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

        peaks_table["name"] = np.array([name.replace("peaks/", "")] * n_star)

        peaks_table["ra"] = ra_peaks
        peaks_table["dec"] = dec_peaks
        peaks_table["sig"] = sig_peaks

        df = pd.DataFrame(data=peaks_table)
        df = df[["name", "sig", "ra", "dec"]]

        df.to_csv("{}.csv".format(name), index=False)






















    #
