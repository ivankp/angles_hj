#include <iostream>
#include <vector>

#include <TF1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

#include "program_options.hh"
#include "tc_msg.hh"
#include "Legendre.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;

#define NPAR 4
const char* pars_names[NPAR] = {"c2","c4","c6","#phi2"};

int main(int argc, char* argv[]) {
  std::array<double,NPAR> pars{};
  std::vector<double> phi2;
  bool logy = false;
  const char* ofname;

  try {
    using namespace ivanp::po;
    if (program_options()
        (pars,{},"function parameters",req(),pos(1))
        (ofname,'o',"output file",req())
        (phi2,{},"φ2",pos())
        (logy,"--logy")
        .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TF1 f("Legendre",Legendre,-1.,1.,NPAR);
  for (int i=0; i<NPAR; ++i)
    f.SetParName(i,pars_names[i]);
  f.SetParameters(pars.data());

  TCanvas canv;
  if (logy) canv.SetLogy();
  gStyle->SetOptStat(0);

  // f.GetYaxis()->SetRangeUser(0.,2.);
  f.Draw();

  for (double x : phi2) {
    f.SetParameter(3,x);
    f.DrawCopy("SAME")->GetYaxis()->SetRangeUser(0.,2.);
  }

  canv.SaveAs(ofname);
}
