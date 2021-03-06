#!/bin/bash


# use py3 on coma (don't need this if already using py3)
source  activate  mypython3

python  preprocess.py

# main search
gc_size_pcs="10"
# gc_size_pcs="5  10"

# input="dwarfs/dwarfs-names-split.txt"
input="dwarfs/dwarfs-names-split-pm.txt"


while IFS= read -r name_dwarf;  do
  for gc_size in $gc_size_pcs;  do
    sbatch  bashtools/submit-slurm-single.sh  $name_dwarf  $gc_size
  done  # for
done < "$input"  # while
