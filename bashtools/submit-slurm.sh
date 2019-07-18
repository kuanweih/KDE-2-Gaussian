#!/bin/bash


# use py3 on coma (don't need this if already using py3)
source  activate  mypython3


# pre-processing dwarf list
cd  ../dwarfs/more
python  more_dwarfs.py
cd  ../McConnachie
python  convert_McConnachie.py
cd  ..
python  joint_list.py
cd  ../bashtools


# main search
# gc_size_pcs="10"
gc_size_pcs="5  10"

# input="dwarfs/dwarfs-names.txt"
input="../dwarfs/dwarfs-names-split.txt"


while IFS= read -r name_dwarf;  do
  for gc_size in $gc_size_pcs;  do
    sbatch  submit-slurm-single.sh  $name_dwarf  $gc_size
  done  # for
done < "$input"  # while



# summary and images
# python  -W  ignore  summary.py
