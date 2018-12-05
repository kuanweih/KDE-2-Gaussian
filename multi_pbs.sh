#!/bin/bash


job_path="/Users/kwhuang/bash-multi-pbs/Fornax"
source_path="/Users/kwhuang/bash-multi-pbs/KDE-Detecor"

ra='39.99708'
dec='-34.55083'
radiuss='2.0'
pmnans='True  False'
pmcuts='0  1'  # std of pm

kernels='gaussian  poisson'
s_ones='0.005  0.01'
s_twos='0.05  0.1'
ngs='400  600'
r_areas='1.0  2.0'
# r_areas='1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0  10.0'


for kernel in $kernels;  do
for s_one in $s_ones;  do
for s_two in $s_twos;  do
for ng in $ngs;  do
for r_area in $r_areas;  do
for radius in $radiuss;  do
for pmnan in $pmnans;  do
for pmcut in $pmcuts;  do

  if [ $pmnan = "True" ] && [ $pmcut != "0" ]; then
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_N"$ng"_s"$s_one"_s"$s_two"_r"$r_area"_pmcut_"$pmcut
    else
      run_name="R"$radius"_"$kernel"_N"$ng"_s"$s_one"_s"$s_two"_pmcut_"$pmcut
    fi
  elif [ $pmnan = "True" ] && [ $pmcut = "0" ]; then
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_N"$ng"_s"$s_one"_s"$s_two"_r"$r_area"_pmnan"
    else
      run_name="R"$radius"_"$kernel"_N"$ng"_s"$s_one"_s"$s_two"_pmnan"
    fi
  elif [ $pmnan = "False" ] && [ $pmcut = "0" ]; then
    if [ $kernel = "poisson" ]; then
      run_name="R"$radius"_"$kernel"_N"$ng"_s"$s_one"_s"$s_two"_r"$r_area
    else
      run_name="R"$radius"_"$kernel"_N"$ng"_s"$s_one"_s"$s_two
    fi
  else
    break;
  fi

  cd  $job_path
  cp  -r  $source_path  $run_name
  cd  $run_name

  sed  "s/NUM_GRID = .*/NUM_GRID = $ng    # number of meshgrids/g"  param_den.py > tmp.py  &&  mv  tmp.py  param_den.py
  sed  "s/SIGMA1 = .*/SIGMA1 = $s_one    # searching scale in deg/g"  param_den.py > tmp.py  &&  mv  tmp.py  param_den.py
  sed  "s/SIGMA2 = .*/SIGMA2 = $s_two    # background scale in deg/g"  param_den.py > tmp.py  &&  mv  tmp.py  param_den.py
  sed  "s/KERNEL_BG = .*/KERNEL_BG = '$kernel'    # kernel/g"  param_den.py > tmp.py  &&  mv  tmp.py  param_den.py
  sed  "s/RATIO_AREA_TG_BG = .*/RATIO_AREA_TG_BG = $r_area    # ratio of area/g"  param_den.py > tmp.py  &&  mv  tmp.py  param_den.py

  sed  "s/RA = .*/RA = $ra    # ra of target deg/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py
  sed  "s/DEC = .*/DEC = $dec    # dec of target in deg/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py
  sed  "s/RADIUS = .*/RADIUS = $radius    # querying radius in deg/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py
  sed  "s/REMOVE_PM_NAN = .*/REMOVE_PM_NAN = $pmnan    # True:on, Flase:off/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py


  if [ $pmcut != "0" ]; then
    sed  "s/PM_CUT = .*/PM_CUT = True    # True:on, Flase:off/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py
  else
    sed  "s/PM_CUT = .*/PM_CUT = False    # True:on, Flase:off/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py
    sed  "s/    PM_CUT_STD = .*/    PM_CUT_STD = $pmcut/g"  param_get.py > tmp.py  &&  mv  tmp.py  param_get.py
  fi

  # qsub  submit-coma.job  -q  bigmem
  qsub  submit-coma.job


done  #  for pmcut
done  #  for pmnan
done  #  for radius
done  #  for r_area
done  #  for ng
done  #  for s_two
done  #  for s_one
done  #  for kernel
