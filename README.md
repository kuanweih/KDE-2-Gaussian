# KDE-Detector
Kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K). There are two kinds of statistics:
- 2 Gaussian kernel convolution
- Poisson distribution


# Parameters and Database
Before starting density estimation, one should set up `param/param.py` and `param/wsdb.py` first. `param/param.py` contains all the parameters for density estimation and `param/wsdb.py` contains necessary information for database.


# How to use it:
1. Clean the work directory: `bash bashtools/clean.sh`
2. Get access to the database and enter information in `param/wsdb.py`
3. Set up parameters in `param/param.py`
4. Calculate density estimation:
    - `python  main.py`<br>
    (if using manual mode)
    - `python  main.py  --dwarf_list  "McConnachie"  --name_dwarf  "Fornax"  --gc_size_pc  10`<br>
    (if using a dwarf list: "McConnachie" or "more_dwarfs")
5. If running on a cluster, a slurm job script is provided:
    - `sbatch  bashtools/slurm-slurm.sh`
6. Summarize searching result with `python summary.py`
