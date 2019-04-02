#!/bin/bash
#SBATCH --export=NONE
#SBATCH --partition=long
#SBATCH --ntasks=20
#SBATCH --time=72:00:00

source  activate  mypython3    # use py3 on coma (don't need this if already using py3)


name_dwarfs="Fornax  UrsaMajorII  UrsaMinor"
# name_dwarfs="LeoIV  LeoV  LeoII"
# name_dwarfs="LeoI  Eridanus2  Hercules"
# name_dwarfs="Carina  Sculptor  *Eridanus3"


gc_size_pcs='1  5  10'

# python  main.py  --name_dwarf  $name_dwarf

for name_dwarf in $name_dwarfs;  do
for gc_size in $gc_size_pcs;  do

python  main.py  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size
python  main.py  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size  --scale_sigma2  2.0

done
done
