#!/usr/bin/env python3

import sys

if (len(sys.argv)!=3): sys.exit(1)

from ROOT import TFile, TTree

f = TFile(sys.argv[1])
t = f.Get("angles")
print("{:,}".format( t.GetEntries() ))

f2 = TFile(sys.argv[2],"recreate")
t.Draw("abs(cos_theta)>>h_cos_theta(50,0,1)")

f2.Write()
