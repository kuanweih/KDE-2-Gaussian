import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from param import *


def visualize_4_panel(path: str, outfile: str, n_error: float, kernel: str, s_above=5):
    """
    kernel: 'gaussian' or 'poisson'
    """
    sns.set(style="white", color_codes=True, font_scale=2)
    fig, axes = plt.subplots(2, 2, figsize=(18, 18))
    fig.suptitle("{}    GC={}pc    {}\nwidth={}deg    s1={}deg    s2={}deg".format(
                 NAME, GC_SIZE, kernel, WIDTH, SIGMA1, SIGMA2), y=0.93)
    plt.subplots_adjust(wspace=0, hspace=0.1)

    x, y = np.load('{}/meshgrids.npy'.format(path))    # coordinates

    sigs = [np.load('{}/sig_{}.npy'.format(path, kernel)),
            np.load('{}/sig_{}-pm_error{}.npy'.format(path, kernel, n_error))]

    datas = [np.load('{}/queried-data.npy'.format(path)).item(),
             np.load('{}/queried-data-pm_error{}.npy'.format(path, n_error)).item()]

    ras = [data["ra"] for data in datas]
    decs = [data["dec"] for data in datas]

    masks = [data["sig_{}".format(kernel)] > s_above for data in datas]    # masks for stars with sig > s_above

    n_stars = [len(data["ra"]) for data in datas]

    extent = [x.min(), x.max(), y.min(), y.max()]    # arg extent for imshow

    for v in range(2):
        axes[v, 0].imshow(sigs[v] > s_above, cmap='copper', vmin=-0.01, vmax=1.01, extent=extent, origin='lower')
        axes[v, 0].plot(ras[v], decs[v], '.', c='deepskyblue', markersize=0.5, alpha=0.5)
        axes[v, 0].plot(ras[v][masks[v]], decs[v][masks[v]], '.', c='orange', markersize=0.5)
        axes[v, 0].set_title('all: {} stars'.format(n_stars[v]) if v==0 else 'pm: {} stars'.format(n_stars[v]))

        axes[v, 1].imshow(sigs[v], cmap='RdBu_r', vmin=0, vmax=8, extent=extent, origin='lower')
        axes[v, 1].set_title('sig > {}: {} pixels'.format(s_above, np.sum(sigs[v] > s_above)))

        for u in range(2):
            axes[v, u].tick_params(axis='both', which='both', labelleft=False, labelbottom=False)

    plt.savefig("{}-{}.png".format(outfile, kernel), bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    """ test plotting """
    from  main  import  get_dir_name
    visualize_4_panel(get_dir_name(), "test.png", N_ERRORBAR)
