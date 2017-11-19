#!/bin/bash

dir=1
mkdir -p $dir

for m in 150_200 200_250 250_350 350_450
do
  ./bin/select1 data/H1j_mtop_unweighted.dat -o $dir/${m}.root -m `sed 's/_/:/' <<< $m`
  ./bin/fit $dir/${m}.root -o $dir/${m}_fit.root --nbins 50 -n 3 -l 2:2:2
  ./bin/draw $dir/${m}_fit.root -o $dir/${m}_fit.pdf #--logy
done
