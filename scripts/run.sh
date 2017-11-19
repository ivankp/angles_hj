#!/bin/bash

for m in 150_200 200_250 250_350 350_450
do
  ./bin/select1 H1j_mtop_unweighted.dat -o 1/${m}.root -m `sed 's/_/:/' <<< $m`
  ./bin/fit 1/${m}.root -o 1/${m}_fit.root --nbins 50 -n 3 -l 2:2:2
  ./bin/draw 1/${m}_fit.root -o 1/${m}_fit.pdf #--logy
done
