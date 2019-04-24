#!/bin/bash
#SBATCH --export=NONE
#SBATCH --partition=long
#SBATCH --ntasks=20
#SBATCH --time=72:00:00



# name_dwarfs="Fornax  UrsaMajorII  UrsaMinor"
# name_dwarfs="*Eridanus3  Eridanus2  Hercules"


gc_size_pcs="10"
# gc_size_pcs="5  10"

input="dwarfs-McConnachie/test-dwarfs-names.txt"

while IFS= read -r name_dwarf;  do
  for gc_size in $gc_size_pcs;  do
    python  main.py  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size
  done  # for
done < "$input"  # while
