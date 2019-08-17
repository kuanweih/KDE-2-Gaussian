#!/bin/bash


# use py3 on coma (don't need this if already using py3)
source  activate  mypython3


bash  bashtools/preprocessing_dwarflist.sh


# main search
gc_size_pcs="10"
# gc_size_pcs="5  10"

# input="dwarfs/dwarfs-names.txt"
input="dwarfs/dwarfs-names-split.txt"


while IFS= read -r name_dwarf;  do
  for gc_size in $gc_size_pcs;  do
    sbatch  bashtools/submit-slurm-single.sh  $name_dwarf  $gc_size
  done  # for
done < "$input"  # while
