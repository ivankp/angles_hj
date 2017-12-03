#!/bin/bash

do_fit=false
use_unweighted=false

while getopts 'fu' flag; do
  case "${flag}" in
    f) do_fit=true ;;
    u) use_unweighted=true ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

for r in 8
do

all_fits=()

for p in 3 4
do

if $use_unweighted; then
  suf=unw
  y="0:2"
else
  suf=full
  case "$r" in
    5) y="0.35:0.8" ;;
    8) y="0.25:1.25" ;;
    *) y="0:2" ;;
  esac
fi

fits=fits_${suf}_${p}_${r}.root
all_fits+=($fits)

if $do_fit; then
  if $use_unweighted; then
    ./bin/fit data/H1j_mtop_unweighted.root -o $fits \
      -M 12:250:550 -n${p} -r 0.${r} --use-chi2-pars -l 2:2:2
  else
    ./bin/fit ../bh_analysis2/H1j_angles.root -o $fits \
      -M 12:250:550 -n${p} -r 0.${r} --use-chi2-pars --nbins=50
      # -M 30:250:550 -n${p} -r 0.${r} --use-chi2-pars
  fi
fi

./bin/draw $fits -y$y

pars=pars_${suf}_${p}_${r}.root

./bin/pars $fits -o $pars

rxplot $pars -ltr \
  -r 'gg/^.*' nl \
     "tt/^.*/${p} par fit,  cos**theta **in [-0.${r},0.${r})" \
     'xx/^.*/Hj mass [GeV]'

done

./bin/llr ${all_fits[*]} -o llr_${suf}_${r}.root

rxplot llr_${suf}_${r}.root \
  -r 'sn/^P.*' 'xx/^.*/Hj mass [GeV]'

done
