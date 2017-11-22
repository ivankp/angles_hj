#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

#include "program_options.hh"
#include "tc_msg.hh"

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;

int main(int argc, char* argv[]) {
  const char *ifname, *ofname;
  bool logy = false;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        (ofname,'o',"output file",req())
        (logy,"--logy")
        .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TFile fin(ifname);
  info("Input file",fin.GetName());
  if (fin.IsZombie()) return 1;

  TH1D *h;
  TF1  *f;
  fin.GetObject("abs_cos_theta",h);
  fin.GetObject("fit",f);

  h->Scale(1./h->Integral("width"));

  TCanvas canv;
  if (logy) canv.SetLogy();
  gStyle->SetOptStat(0);
  h->SetLineWidth(2);
  h->SetTitle("");
  h->SetXTitle("|cos #theta|");
  h->Draw();
  f->Draw("same");

  TLatex latex;
  latex.SetTextSize(0.025);

  auto l = [&](double x, double y, const char* name){
    TNamed *named;
    fin.GetObject(name,named);
    latex.DrawLatexNDC(x,y,cat(name," = ",named->GetTitle()).c_str());
  };
  l(0.15,0.85,"c2");
  l(0.15,0.80,"c4");
  l(0.15,0.75,"c6");
  l(0.15,0.70,"phi2");
  l(0.15,0.65,"-2LogL");

  canv.SaveAs(ofname);
}
