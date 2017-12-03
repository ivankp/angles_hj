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
#include <TLine.h>

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

std::vector<std::string> split(const char* str, char d) {
  std::vector<std::string> tok;
  unsigned i = 0;
  for (char c; (c=str[i])!='\0'; ++i) {
    if (c==d) {
      tok.emplace_back(str,i);
      str += i+1;
      i = 0;
    }
  }
  tok.emplace_back(str,i);
  return tok;
}

int main(int argc, char* argv[]) {
  std::string ifname, ofname;
  bool logy = false, more_logy = false;
  boost::optional<std::array<double,2>> y_range;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifname,'i',"input file",req(),pos())
      (ofname,'o',"output file")
      (y_range,'y',"y-axis range")
      (more_logy,"--more-logy","more y-axis log labels")
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
  if (hj_mass_bins.size()==0) {
    error("nothing to draw");
    return 1;
  } else info("hj_mass bins",hj_mass_bins.size());
  hj_mass_bins.sort();

  TCanvas canv("","",700,600);
  if (logy) canv.SetLogy();
  gStyle->SetOptStat(0);

  TPad *pad1 = new TPad("","",0,0.25,1,1);
  pad1->SetMargin(0.05,0.05,0,0.1);
  TPad *pad2 = new TPad("","",0,0,1,0.25);
  pad2->SetMargin(0.05,0.05,0.25,0);

  TLatex latex;
  latex.SetTextSize(0.025);

  if (ofname.empty()) {
    ofname = ifname.substr(ifname.rfind('/')+1);
    ofname = ofname.substr(0,ofname.rfind(".root"))+".pdf";
  }
  info("Output file",ofname);
  if (hj_mass_bins.size()>1) ofname += '(';
  unsigned page_cnt = hj_mass_bins.size();
  bool first_page = true;
  TH1 *h_mc = nullptr, *h_fit = nullptr;
  for (const auto& bin : hj_mass_bins) {
    --page_cnt;
    const std::string& name = bin.first;

    int fi = 0;
    double scale = 1.;

    pad1->cd();
    for (auto* p : bin.second) {
      if (p->InheritsFrom(TH1::Class())) { // HIST
        TH1 *h = static_cast<TH1*>(p);

        // h->Scale(1./h->Integral("width")); // normalize
        h->SetLineWidth(2);
        h->SetTitle(cat("hj_mass #in ",name).c_str());

        scale = 1./h->Integral("width");
        h->Scale(scale);
        h_mc = static_cast<TH1*>(h->Clone());
        h->SetXTitle("");

        auto* ya = h->GetYaxis();
        if (y_range) ya->SetRangeUser((*y_range)[0],(*y_range)[1]);
        if (more_logy) ya->SetMoreLogLabels();
        h->Draw();
        latex.DrawLatexNDC(0.70,0.85,cat("Events: ",h->GetEntries()).c_str());

      } else if (p->InheritsFrom(TF1::Class())) { // TF1
        TF1 *f = static_cast<TF1*>(p);
        TH1 *h = new TH1D("","",f->GetNpx(),f->GetXmin(),f->GetXmax());
        h->Add(f);
        h->SetLineWidth(2);
        h->SetLineColor(f->GetLineColor());

        if (strstr(f->GetName(),"-logl-")) {
          h->Scale(1./h->Integral("width"));
          h_fit = static_cast<TH1*>(h->Clone());
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
        int line = 0;
        for (; line<npar; ++line)
          l(0.15+0.2*fi,0.85-0.04*line,line);
        bool first = true;
        for (const auto& str : split(f->GetTitle(),',')) {
          if (first)
            latex.SetTextColor(f->GetLineColor());
          latex.DrawLatexNDC(0.15+0.2*fi,0.85-0.04*line,str.c_str());
          if (first) {
            latex.SetTextColor(1);
            first = false;
          }
          ++line;
        }
        ++fi;
      }
    }

    pad2->cd();
    h_mc->SetTitle("");
    h_fit->SetTitle("");
    TAxis *ax = h_mc->GetXaxis();
    TAxis *ay = h_mc->GetYaxis();
    ax->SetLabelSize(0.1);
    ay->SetLabelSize(0.1);
    ax->SetTitleSize(0.1);
    ax->SetTickLength(0.09);
    ay->SetTickLength(0.02);
    const int factor = h_fit->GetNbinsX()/h_mc->GetNbinsX();
    h_fit->Rebin(factor);
    h_fit->Scale(1./factor);
    h_mc->Divide(h_fit);
    ay->SetRangeUser(0.95,1.05);
    h_mc->Draw();
    static TLine line1(-1,1,1,1);
    line1.Draw();

    canv.cd();
    pad1->Draw();
    pad2->Draw();

    if (!page_cnt && !first_page) ofname += ')';
    canv.Print(ofname.c_str(),("Title:"+std::string(name,1,name.size()-2)).c_str());
    if (first_page) ofname.pop_back(), first_page = false;
  }

}
