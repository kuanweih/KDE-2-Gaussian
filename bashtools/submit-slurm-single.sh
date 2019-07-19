#!/bin/bash
#SBATCH  --export=NONE
#SBATCH  --partition=long
#SBATCH  --constraint=intel_e5_v4
#SBATCH  --ntasks-per-node=16
#SBATCH  --time=24:00:00
#SBATCH  --job-name=dwarf


# use py3 on coma (don't need this if already using py3)
source  activate  mypython3

name_dwarf="$1"
gc_size="$2"

python  -W  ignore  main.py  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size
