#!/bin/bash

export LD_LIBRARY_PATH=/msu/data/t3work3/ivanp/gcc-7.2.0/hep/root-6.10.02/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/hep/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/gcc/lib64:/msu/data/t3work3/ivanp/gcc-7.2.0/gcc/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/lib64:/msu/data/t3work3/ivanp/gcc-7.2.0/lib

phi=${1:?}

# data=/home/ivanp/work/bh_analysis2/H1j_angles.root
data=/msu/data/t3work2/ivanp/H1j_cos_theta.root
# data=/home/ivanp/work/angles_hj/data/H1j_mtop_unweighted.root

/home/ivanp/work/angles_hj/bin/fit $data \
  -o ${phi}.root \
  -M 12:250:550 -p 0:0:0:${phi} \
  -n 3 -r 0.8 --use-chi2-pars --nbins=50

/home/ivanp/work/angles_hj/bin/draw \
  ${phi}.root -o ${phi}.pdf -y 0.25:1.25

