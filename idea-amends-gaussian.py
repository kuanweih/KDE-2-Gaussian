import numpy as np
from scipy.ndimage import gaussian_filter


xedges = np.linspace(0, 10, num=5)
yedges = np.linspace(0, 10, num=5)

x = np.random.normal(3, 1, 100)  # stars' ra
y = np.random.normal(7, 1, 100)  # stars' dec
H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
S1 = gaussian_filter(H, sigma=0.2)
S2 = gaussian_filter(H, sigma=1.2)

print(S1 - S2)
