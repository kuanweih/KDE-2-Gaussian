# KDE-Detector
Kernel density estimation for 2D observational data, based on [Koposov et al 2008](http://cdsads.u-strasbg.fr/abs/2008ApJ...686..279K). There are two kinds of statistics:
- 2 Gaussian kernel convolution
- Poisson distribution


# Parameters and Database
Before starting density estimation, one shall first set up:
- `src/param.py`: parameters for density estimation
- `src/param_patch_candidate.py`: parameters for preprocessing and summary
- `wsdb.py`: permission for wsdb.


# How to use it:
0. Clean the work directory: `bash bashtools/clean.sh`
1. Preprocess a dwarf list (optional): <br>
   `python  preprocess.py`
2. Get access to the database and enter information in `wsdb.py`
3. Set up parameters in `src/param.py`, especially the following:
    - `IS_DWARF_LIST = False`    # use joint list
    - `IS_DWARF_SPLIT_LIST = True`    # use joint-split list
4. Calculate density estimation:
    - `python  -W  ignore  main.py`<br>
    (if using manual mode)
    - `python  -W  ignore  main.py  --name_dwarf  "Fornax"  --gc_size_pc  10`<br>
    (if using the joint or joint-split dwarf list: names can be find in the txt files in `dwarfs/`)
5. Summarize searching result with `python  -W  ignore  summary.py`
6. If running step 4 and 5 on a cluster, slurm job scripts are provided:
    - `bash  bashtools/slurm-slurm.sh`    # make sure the right input txt
    - `sbatch  bashtools/slurm-summary.sh`    # make sure all the KDE searches are done and then run this command
