#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>

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

ordered_map<std::vector<double>> hj_mass_bins;

void loop(TDirectory* dir) { // LOOP
  for (TKey& key : get_keys(dir)) {
    const TClass* key_class = get_class(key);
    const char* key_name = key.GetName();
    if (inherits_from<TF1>(key_class)) { // TF1
      if (!strstr(key_name,"-logl-")) continue;
      const char* bin1 = strchr(key_name,'[');
      const char* bin2 = strchr(bin1+1,')');

      hj_mass_bins[std::string(bin1,bin2-bin1+1)].push_back(
        atof(strchr(read_key<TF1>(key)->GetTitle(),'=')+1));

    } else if (inherits_from<TDirectory>(key_class)) { // DIR
      loop(read_key<TDirectory>(key));
    }
  }
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char *ofname;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input file\nfirst is fixed",req(),pos())
      (ofname,'o',"output file",req())
      .parse(argc,argv,true)) return 0;
    if (ifnames.size()<2) throw std::runtime_error(
      "need at least 2 input files");
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  for (const char* ifname : ifnames) {
    TFile f(ifname);
    info("Input file",f.GetName());
    if (f.IsZombie()) return 1;
    loop(&f);
  }
  if (hj_mass_bins.size()==0) {
    error("no function found");
    return 1;
  } else info("hj_mass bins",hj_mass_bins.size());
  hj_mass_bins.sort();

  TFile fout(ofname,"recreate");
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  std::vector<double> bins;
  std::set<unsigned> skip;
  unsigned i = 1;
  for (const auto& b : hj_mass_bins) {
    const auto& str = b.first;
    const auto d = str.find(',');
    std::array<double,2> bin {
      stod(str.substr(1,d-1)),
      stod(str.substr(d+1,str.size()-d-2))
    };

    if (bins.empty() || bins.back()!=bin[0]) {
      if (!bins.empty()) skip.insert(i);
      bins.push_back(bin[0]);
    }
    bins.push_back(bin[1]);
    ++i;
  }

  for (unsigned f=1, nf=ifnames.size(); f<nf; ++f) {
    TH1D* h_llr = new TH1D(cat("LLR",f).c_str(),"LLR",bins.size()-1,bins.data());
    TH1D* h_p = new TH1D(cat("P",f).c_str(),"P-value",bins.size()-1,bins.data());
    i = 1;
    for (const auto& b : hj_mass_bins) {
      if (skip.count(i)) ++i;
      const double llr = b.second.front()-b.second[f];
      h_llr->SetBinContent(i,llr);
      h_p->SetBinContent(i,TMath::Prob(llr,1));
      ++i;
    }
  }

  fout.Write();
}
