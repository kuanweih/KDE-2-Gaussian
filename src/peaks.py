import numpy as np



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
    datas = [np.load('{}/queried-data.npy'.format(path)).item(),
             np.load('{}/queried-data-pm_error{}.npy'.format(path,
                                                             n_error)).item()]

    ras = [data["ra"] for data in datas]
    decs = [data["dec"] for data in datas]
    sigs = [data["sig_{}".format(kernel)] for data in datas]

    # masks for stars with sig > s_above
    masks = [sigs > s_above for data in datas]

    for v in range(2):
        ra_peaks = ras[v][masks[v]]
        dec_peaks = decs[v][masks[v]]
        sig_peaks = sigs[v][masks[v]]

        peaks_table = {}
        peaks_table["ra"] = ra_peaks
        peaks_table["dec"] = dec_peaks
        peaks_table["ra"] = sig_peaks

        df = pd.DataFrame(data=searching_table)
        df = df[["name", "gc_size", "n_star", "n_star_pm", "sig_g_peak", "sig_g_pm_peak", "sig_p_peak", "sig_p_pm_peak"]]

        df.to_csv("summary.csv", index=False)




















    #
