# KDE-Detector
Kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K). There are two kinds of statistics:
- 2 Gaussian kernel convolution
- Poisson distribution


# Parameters and Database
Before starting density estimation, one should set up `param/param.py` and `param/wsdb.py` first. `param/param.py` contains all the parameters for density estimation and `param/wsdb.py` contains necessary information for database.


# How to use it:
0. Preprocess a dwarf list (optional): <br>
   `bash  bashtools/preprocessing_dwarflist.sh`
1. Clean the work directory: `bash bashtools/clean.sh`
2. Get access to the database and enter information in `param/wsdb.py`
3. Set up parameters in `param/param.py`, especially the following:
    - `IS_DWARF_LIST = False`    # use joint list
    - `IS_DWARF_SPLIT_LIST = True`    # use joint-split list
4. Calculate density estimation:
    - `python  -W  ignore  main.py`<br>
    (if using manual mode)
    - `python  -W  ignore  main.py  --name_dwarf  "Fornax"  --gc_size_pc  10`<br>
    (if using the joint or joint-split dwarf list: names can be find in the txt files in `dwarfs/`)
5. If running on a cluster, a slurm job script is provided:
    - `bash  bashtools/slurm-slurm.sh`
6. Summarize searching result with `python  -W  ignore  summary.py`
