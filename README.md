# KDE-Detector
Kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K). There are two choices of statistics:
- 2 Gaussian kernel convolution
- Poisson distribution (currently unavailable)


# How to use it:
1. get access to the database and enter information in `kw_wsdb.py`
2. setup paramters in param.py
3. ```python  main.py``` to execute the code

P.S. If running on clusters, a pbs script is provided also.
