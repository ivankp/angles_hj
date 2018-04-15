#!/usr/bin/env python3

import sys, re
from collections import defaultdict
from ROOT import gROOT, TFile, TCanvas, TH1, TF1

_d = lambda: defaultdict(_d)
d = defaultdict(_d)

for arg in sys.argv[1:]:
    print(arg)
    f = TFile(arg)
    phi = float(re.sub(r"\.root$","",arg))
    for h in f.GetListOfKeys():
        h = h.ReadObj()
        name = h.GetName()
        print(name)
        M = tuple([ float(x) for x in
            re.search(r"hj_mass\[([0-9,\.]*)\)",name).group(1).split(',') ])

        fs = h.GetListOfFunctions()
        d[M][phi] = float(
            next( f for f in h.GetListOfFunctions() if "-logl" in f.GetName() )
            .GetTitle().split('=')[1])

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pages = PdfPages('scan.pdf')

for M, d_phi in d.items():
    print("hj_mass:",M)
    logl0 = None
    x = []
    y = []
    for phi, logl in sorted(d_phi.items()):
        if logl0 is None: logl0 = logl
        logl -= logl0
        print("{}: {}".format(phi,logl))
        # if phi!=1.1:
        # if logl < 1e3:
        if logl < 1e1:
            x.append(phi)
            y.append(logl)

    fig = plt.figure()
    plt.plot(x,y,'ro-')
    plt.title(r"hj_mass $\in [{:.0f},{:.0f})$".format(*M))
    plt.xlabel(r"$\phi_2$", horizontalalignment='right', x=1.0)
    plt.ylabel(r"$\Delta L = L_0 - L(\phi_2)$", horizontalalignment='right', y=1.0)

    plt.figtext(0.75,0.95,r"$L = -2 log \sum w_i f(x_i)$")
    plt.figtext(0.75,0.90,r"$L_0 = {:.5E}$".format(logl0))

    pages.savefig(fig)

pages.close()
