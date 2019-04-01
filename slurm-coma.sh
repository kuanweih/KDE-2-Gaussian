#!/bin/bash
#SBATCH --export=NONE
#SBATCH --partition=long
#SBATCH --ntasks=20
#SBATCH --time=72:00:00

source  activate  mypython3    # use py3 on coma (don't need this if already using py3)

python  main.py  --name_dwarf  "Fornax"
python  main.py  --name_dwarf  "Fornax"  --gc_size_pc  20  --scale_sigma2  2.0
