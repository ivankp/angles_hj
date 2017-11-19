#!/bin/bash

# for f in 150_200 200_250 250_350 350_450
# do scripts/hist.py 1/${f}.root 1/${f}_h.root
# done

rxplot $1/*_fit.root -o $1/abs_cos_theta.pdf -ltl \
  -r 'fl|(.*/)?(.+)_fit\.root|M = \2 GeV' 't/^.*/ ' 'ng///norm=1' 'x/^.*/|cos #theta|'

