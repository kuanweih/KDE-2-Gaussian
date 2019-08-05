import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as stats

from param.param import *


_st1 = '{}  GC={}pc'.format(NAME, GC_SIZE)
_st2 = 'd={}kpc  w={}$^\circ$'.format(round(DISTANCE / 1e3), WIDTH)
_st3 = 's1={}$^\circ$  s2={}$^\circ$'.format(SIGMA1, SIGMA2)
SUB_TITLE = "{}  {}  {}".format(_st1, _st2, _st3)


def visualize_2_panel(path: str, outfile: str, kernel: str, s_above=5):
    """ Plotting star distribution (left panels) and density maps (right
    panels). (Others)

    : path : path of the result file
    : outfile : where to output the plot
    : kernel : 'gaussian' or 'poisson'
    : s_above : significance threshold, default value = 5
    """
    sns.set(style="white", color_codes=True, font_scale=1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle(SUB_TITLE, y=0.93)

    x, y = np.load('{}/meshgrids.npy'.format(path))    # coordinates
    sig = np.load('{}/sig_{}.npy'.format(path, kernel))
    data = np.load('{}/queried-data.npy'.format(path)).item()

    ra, dec, n_star = data["ra"], data["dec"], len(data["ra"])
    extent = [x.min(), x.max(), y.min(), y.max()]    # arg extent for imshow

    is_peak = sig > s_above
    mask = data["sig_{}".format(kernel)] > s_above

    axes[0].plot(ra, dec, '.', c='deepskyblue', ms=0.5, alpha=0.5)
    axes[0].plot(ra[mask], dec[mask], '.', c='orange', ms=1)
    axes[0].set_title('%d stars' % n_star)
    axes[1].set_title('%s:  sig > %0.1f= %d pixels' % (kernel, s_above, np.sum(is_peak)))

    for u in range(2):
        axes[u].imshow(is_peak, cmap='copper', vmin=-0.01, vmax=1.01,
                       extent=extent, origin='lower')
        axes[u].tick_params(axis='both', which='both',
                            labelleft=False, labelbottom=False)
        axes[u].set_xlim(axes[u].set_xlim()[::-1])    # flipping

    _filename = "{}-{}.png".format(outfile, kernel)
    plt.savefig(_filename, bbox_inches='tight', dpi=300)



def hist_2_panel(path: str, outfile: str, kernel: str, s_above=5):
    """ Plotting histograms (left panels) and normalized histograms (right
    panels). (Others)

    : path : path of the result file
    : outfile : where to output the plot
    : kernel : 'gaussian' or 'poisson'
    : s_above : significance threshold, default value = 5
    """
    sns.set(style="white", color_codes=True, font_scale=1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle(SUB_TITLE, y=0.97)
    plt.subplots_adjust(wspace=0.3)

    sig = np.load('{}/sig_{}.npy'.format(path, kernel))
    is_peak = sig > s_above
    sig_finite_flat = sig[np.isfinite(sig)].flatten()

    bins = 20

    axes[0].hist(sig_finite_flat, bins=bins)
    axes[0].set_title('%e pixels' % len(sig_finite_flat))
    axes[0].set_ylabel('number of pixels')

    axes[1].hist(sig_finite_flat, bins=bins, density=True)
    axes[1].set_title('sig > %0.1f= %d pixels' % (s_above, np.sum(is_peak)))
    axes[1].set_ylabel('normalized density of pixels')

    # norm distribution
    mu, variance = 0., 1.
    sigma = np.sqrt(variance)
    xmin, xmax = sig_finite_flat.min(), sig_finite_flat.max()
    x = np.linspace(mu + xmin * sigma, mu + xmax * sigma, 100)
    axes[1].plot(x, stats.norm.pdf(x, mu, sigma), lw=3)
    axes[1].set_xlim([xmin, 10])
    axes[1].set_ylim([1e-10, 1])

    for u in range(2):
        axes[u].set_xlabel('significance in %s' % kernel)
        axes[u].set_yscale('log')

    _filename = "{}-{}.png".format(outfile, kernel)
    plt.savefig(_filename, bbox_inches='tight', dpi=100)


if __name__ == '__main__':
    """ test plotting """
    from main import get_dir_name
    path_dir = get_dir_name()
    visualize_4_panel(path_dir, "test_g.png", N_ERRORBAR, "gaussian")
    visualize_4_panel(path_dir, "test_p.png", N_ERRORBAR, "poisson")
    # hist_2_panel(path_dir, "test_p_hist.png", N_ERRORBAR, "poisson")
