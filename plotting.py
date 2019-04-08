import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from matplotlib import ticker
import numpy as np
import seaborn as sns
from param import *



def visualize_4_panel(path: str, outfile: str, n_error: float, s_above=5):
    """
    path:
    outfile:
    """
    sns.set(style="white", color_codes=True, font_scale=2)
    fig, axes = plt.subplots(2, 2, figsize=(18, 18))
    fig.suptitle("{}    GC={}pc    width={}deg    s1={}    s2={}".format(
                 NAME, GC_SIZE, WIDTH, SIGMA1, SIGMA2), y=0.92)
    plt.subplots_adjust(wspace=0, hspace=0.05)

    x, y = np.load('{}/meshgrids.npy'.format(path))    # coordinates

    sigs = [np.load('{}/significance.npy'.format(path)),
            np.load('{}/significance-pm_error{}.npy'.format(path, n_error))]

    datas = [np.load('{}/queried-data.npy'.format(path)).item(),
             np.load('{}/queried-data-pm_error{}.npy'.format(path, n_error)).item()]

    ras = [data["ra"] for data in datas]
    decs = [data["dec"] for data in datas]

    extent = [x.min(), x.max(), y.min(), y.max()]    # arg extent for imshow

    for v in range(2):
        axes[v, 0].imshow(sigs[v] > s_above, cmap='copper', vmin=-0.01, vmax=1.01, extent=extent, origin='lower')
        axes[v, 0].plot(ras[v], decs[v], '.', c='deepskyblue', markersize=0.5, alpha=0.5)

        axes[v, 1].imshow(sigs[v], cmap='RdBu_r', vmin=-5, vmax=5, extent=extent, origin='lower')

        for u in range(2):
            axes[v, u].tick_params(axis='both', which='both', labelleft=False, labelbottom=False)

    plt.savefig(outfile, bbox_inches='tight', dpi=300)



        # axes[i, j].imshow(sig>above, cmap='RdBu_r', vmin=-0.2, vmax=1.1,
        #                   extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
        #
        # axes[i, j].set_ylabel('{} deg'.format(width))
        # axes[i, j].tick_params(axis='both', which='both', labelleft=False, labelbottom=False)
        # axes[i, j].set_xlabel('sig > {}    s1 = {}    s2 = {}'.format(above, s1, s2))
        # if pm==0:
        #     axes[i, j].set_title('{}    {}'.format(dwarf_name, g_band))
        # elif pm==3:
        #     axes[i, j].set_title('{}    {}    pm_error'.format(dwarf_name, g_band))
        # else:
        #     axes[i, j].set_title('{}    {}    pm<{}'.format(dwarf_name, g_band, pm))




    # # Set up axes
    # ax.set_xticklabels([''] + source_sentence_str, rotation=90)
    # ax.set_yticklabels([''] + target_sentence_str)
    #
    # # Show label at every tick
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    #
    # plt.savefig(outfile)
    #
    # plt.close()


# path = "../KDE-Detector/results-gaussian/Fornax/G17-22/w2-lp0.001/s0.004s0.02s1.0sth1/"
#
# print(path.split("/"))
# print(path.split("/")[6].split("s")[1])
#
#     def plot_sig4(path, above=5):
#     dwarf_name = path.split("/")[3]
#     g_band = path.split("/")[4]
#     width = path.split("/")[5].split("-")[0][1:]
#     s1 = path.split("/")[6].split("s")[1]
#     s2 = path.split("/")[6].split("s")[2]





    # savefig('{}-{}-w{}-s{}-s.png'.format(dwarf_name, g_band, width, s1, s2),
    #         bbox_inches='tight', dpi=300)
