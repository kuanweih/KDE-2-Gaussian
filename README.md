# KDE-Detector
Kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K). There are two choices of statistics:
- 2 Gaussian kernel convolution
- Poisson distribution (currently unavailable)


# How to use it:
1. Get access to the database and enter information in `kw_wsdb.py`
2. Setup paramters in param.py
3. Execute the code:
    - `python  main.py  --name_dwarf  "Fornax"`
    - `python  main.py  --name_dwarf  "Fornax"  --scale_sigma2  2.0`
4. If running on a cluster, a slurm job script is provided:
    - `sbatch  slurm-coma.sh`
