# KDE-Detector
Kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K). There are two choices of statistics: 
- 2 Gaussian kernel convolution 
- Poisson distribution 


# How to use it: 
1. setup paramters in param_get.py and param_den.py
2. ```python get_ra_dec.py``` to get numpy array by querying data from database
3. ```python density_estimation.py``` to calculate overdensity and significance 

P.S. If running on clusters, a pbs script is provided also. 
