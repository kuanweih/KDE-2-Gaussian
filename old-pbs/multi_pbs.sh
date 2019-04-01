#!/bin/bash


job_path="/Users/kwhuang/bash-multi-pbs/Fornax"
source_path="/Users/kwhuang/bash-multi-pbs/KDE-Detector"

ra='39.99708'
dec='-34.55083'
radiuss='2.0'
pmcuts='3.'  # std of pm, 0: no cut
plxcuts='3.'  # std of parallax, 0: no cut

kernels='gaussian'
# kernels='poisson'
s_ones='0.005  0.01'
s_twos='0.05  0.1'
pixels='0.002'
r_areas='1.0  2.0'


for kernel in $kernels;  do
for s_one in $s_ones;  do
for s_two in $s_twos;  do
for pixel in $pixels;  do
for r_area in $r_areas;  do
for radius in $radiuss;  do
for pmcut in $pmcuts;  do
for plxcut in $plxcuts;  do

  if [ $pmcut != "0" ] && [ $plxcut != "0" ]; then
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_r"$r_area"_pmcut_"$pmcut"_plxcut_"$plxcut
    else
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_pmcut_"$pmcut"_plxcut_"$plxcut
    fi
  elif [ $pmcut = "0" ] && [ $plxcut != "0" ]; then
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_r"$r_area"_plxcut_"$plxcut
    else
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_plxcut_"$plxcut
    fi
  elif [ $pmcut != "0" ] && [ $plxcut = "0" ]; then
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_r"$r_area"_pmcut_"$pmcut
    else
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_pmcut_"$pmcut
    fi
  else
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two"_r"$r_area
    else
      run_name="R"$radius"_"$kernel"_p"$pixel"_s"$s_one"_s"$s_two
    fi
  fi

  cd  $job_path
  cp  -r  $source_path  $run_name
  cd  $run_name

  sed  "s/PIXEL_SIZE = .*/PIXEL_SIZE = $pixel    # pixel size/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  sed  "s/SIGMA1 = .*/SIGMA1 = $s_one    # searching scale in deg/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  sed  "s/SIGMA2 = .*/SIGMA2 = $s_two    # background scale in deg/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  sed  "s/KERNEL_BG = .*/KERNEL_BG = '$kernel'    # kernel/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  sed  "s/DR_FROM_S2 = .*/DR_FROM_S2 = $r_area    # delta distance in deg/g"  param.py > tmp.py  &&  mv  tmp.py  param.py

  sed  "s/RA = .*/RA = $ra    # ra of target deg/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  sed  "s/DEC = .*/DEC = $dec    # dec of target in deg/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  sed  "s/RADIUS = .*/RADIUS = $radius    # querying radius in deg/g"  param.py > tmp.py  &&  mv  tmp.py  param.py

  if [ $pmcut != "0" ]; then
    sed  "s/PM_CUT = .*/PM_CUT = True    # True:on, Flase:off/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  else
    sed  "s/PM_CUT = .*/PM_CUT = False    # True:on, Flase:off/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
    sed  "s/    PM_CUT_STD = .*/    PM_CUT_STD = $pmcut/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  fi

  if [ $plxcut != "0" ]; then
    sed  "s/PARALLAX_CUT = .*/PARALLAX_CUT = True    # True:on, Flase:off/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  else
    sed  "s/PARALLAX_CUT = .*/PARALLAX_CUT = False    # True:on, Flase:off/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
    sed  "s/    PARALLAX_CUT_STD = .*/    PARALLAX_CUT_STD = $plxcut/g"  param.py > tmp.py  &&  mv  tmp.py  param.py
  fi

  # qsub  submit-coma.job  -q  bigmem
  qsub  submit-coma.job


done  #  for pmcut
done  #  for pkxcut
done  #  for radius
done  #  for r_area
done  #  for pixel
done  #  for s_two
done  #  for s_one
done  #  for kernel
