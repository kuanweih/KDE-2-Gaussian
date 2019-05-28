#!/bin/bash
#SBATCH --export=NONE
#SBATCH --partition=long
#SBATCH --ntasks=20
#SBATCH --time=72:00:00


# use py3 on coma (don't need this if already using py3)
source  activate  mypython3


gc_size_pcs="10"
# gc_size_pcs="5  10"


input="dwarfs/McConnachie/dwarfs-names.txt"
while IFS= read -r name_dwarf;  do
  for gc_size in $gc_size_pcs;  do
    python  main.py  --dwarf_list  "McConnachie"  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size
  done  # for
done < "$input"  # while


input="dwarfs/20190527/dwarfs-names.txt"
while IFS= read -r name_dwarf;  do
  for gc_size in $gc_size_pcs;  do
    python  main.py  --dwarf_list  "more_dwarfs"  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size
  done  # for
done < "$input"  # while


python  summary.py
