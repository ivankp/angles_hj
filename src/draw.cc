#include <iostream>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

#include "program_options.hh"
#include "tc_msg.hh"
#include "tkey.hh"
#include "ordered_map.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;

ordered_map<std::vector<TObject*>> hj_mass_bins;

void loop(TDirectory* dir) { // LOOP
  for (TKey& key : get_keys(dir)) {
    const TClass* key_class = get_class(key);
    if (
      key_class->InheritsFrom(TH1::Class()) ||
      key_class->InheritsFrom(TF1::Class())
    ) {
      hj_mass_bins[strchr(key.GetName(),'[')].push_back(key.ReadObj());
    } else if (key_class->InheritsFrom(TDirectory::Class())) { // DIR
      loop(read_key<TDirectory>(key));
    }
  }
}

int main(int argc, char* argv[]) {
  std::string ifname, ofname;
  bool logy = false;
  boost::optional<std::array<double,2>> y_range;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        (ofname,'o',"output file")
        (y_range,'y',"y-axis range")
        (logy,"--logy")
        .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TFile fin(ifname.c_str());
  info("Input file",fin.GetName());
  if (fin.IsZombie()) return 1;

  loop(&fin);
  hj_mass_bins.sort();

  TCanvas canv;
  if (logy) canv.SetLogy();
  gStyle->SetOptStat(0);

  TLatex latex;
  latex.SetTextSize(0.025);

  if (ofname.empty()) {
    ofname = ifname.substr(ifname.rfind('/')+1);
    ofname = ofname.substr(0,ofname.rfind(".root"))+".pdf";
  }
  info("Output file",ofname);
  ofname += '(';
  unsigned page_cnt = hj_mass_bins.size();
  bool first_page = true;
  for (const auto& bin : hj_mass_bins) {
    --page_cnt;
    const auto& name = bin.first;

    int fi = 0;
    double scale = 1.;

    for (auto* p : bin.second) {
      if (p->InheritsFrom(TH1::Class())) { // HIST
        TH1 *h = static_cast<TH1*>(p);

        // h->Scale(1./h->Integral("width")); // normalize
        h->SetLineWidth(2);
        h->SetTitle(cat("hj_mass #in ",strchr(h->GetName(),'[')).c_str());

        scale = 1./h->Integral("width");
        h->Scale(scale);

        if (y_range) h->GetYaxis()->SetRangeUser((*y_range)[0],(*y_range)[1]);
        h->Draw();

      } else if (p->InheritsFrom(TF1::Class())) { // TF1
        TF1 *f = static_cast<TF1*>(p);
        TH1 *h = new TH1D("","",100,0,1);
        h->Add(f);
        h->SetLineWidth(2);
        h->SetLineColor(f->GetLineColor());

        if (strstr(f->GetName(),"-logl-")) {
          h->Scale(1./h->Integral("width"));
        } else {
          h->Scale(scale);
        }

        h->Draw("C SAME");

        auto l = [&](double x, double y, int i){
          latex.DrawLatexNDC(x,y,cat(
              f->GetParName(i)," = ",f->GetParameter(i)
            ).c_str());
        };
        const int npar = f->GetNpar();
        for (int i=0; i<npar; ++i)
          l(0.15+0.2*fi,0.85-0.05*i,i);
        latex.DrawLatexNDC(0.15+0.2*fi,0.85-0.05*npar,f->GetTitle());
        ++fi;
      }
    }

    if (!page_cnt) ofname += ')';
    canv.Print(ofname.c_str(),("Title:"+name).c_str());
    if (first_page) ofname.pop_back(), first_page = false;
  }

}
