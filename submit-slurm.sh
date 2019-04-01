#!/bin/bash
#SBATCH --export=NONE
#SBATCH --partition=long
#SBATCH --ntasks=20
#SBATCH --time=72:00:00

source  activate  mypython3    # use py3 on coma (don't need this if already using py3)

name_dwarf="Fornax"
# name_dwarf="UrsaMajorII"
# name_dwarf="UrsaMinor"
# name_dwarf="LeoIV"
# name_dwarf="LeoV"
# name_dwarf="LeoII"
# name_dwarf="LeoI"
# name_dwarf="*Eridanus3"


gc_size_pcs='10 20 30 40 50'

# python  main.py  --name_dwarf  $name_dwarf

for gc_size in $gc_size_pcs;  do

python  main.py  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size
python  main.py  --name_dwarf  $name_dwarf  --gc_size_pc  $gc_size  --scale_sigma2  2.0

done
