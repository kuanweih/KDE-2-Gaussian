# KDE-2-Gaussian
Two Gaussians kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K) and plus another estimation for background overdensity with Poisson distribution. 


# How to use it: 

1. setup paramters in param_get.py and param_den.py
2. run get_ra_dec.py to get numpy array by querying data from database
3. run density_estimation.py to calculate overdensity and significance 

P.S. If running on clusters, a pbs script is provided also. 
