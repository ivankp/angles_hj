#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>

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

ordered_map<TF1*> hj_mass_bins;

void loop(TDirectory* dir) { // LOOP
  for (TKey& key : get_keys(dir)) {
    const TClass* key_class = get_class(key);
    const char* key_name = key.GetName();
    if (inherits_from<TH1>(key_class)) { // TH1
      const char* bin1 = strchr(key_name,'[');
      const char* bin2 = strchr(bin1+1,')');
      for (TObject* obj : *key_cast<TH1>(key)->GetListOfFunctions()) {
        if (!obj->InheritsFrom(TF1::Class())) continue;
        if (!strstr(obj->GetName(),"-logl")) continue;

        hj_mass_bins[std::string(bin1,bin2-bin1+1)] =
          static_cast<TF1*>(obj);
      }
    } else if (inherits_from<TDirectory>(key_class)) { // DIR
      loop(key_cast<TDirectory>(key));
    }
  }
}

int main(int argc, char* argv[]) {
  const char *ifname, *ofname;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifname,'i',"input file",req(),pos())
      (ofname,'o',"output file",req())
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TFile fin(ifname);
  info("Input file",fin.GetName());
  if (fin.IsZombie()) return 1;

  loop(&fin);
  if (hj_mass_bins.size()==0) {
    error("nothing to draw");
    return 1;
  } else info("hj_mass bins",hj_mass_bins.size());
  hj_mass_bins.sort();

  TFile fout(ofname,"recreate");
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TF1* f1 = hj_mass_bins.front().second;

  std::vector<double> bins;
  std::set<unsigned> skip;
  for (const auto& b : hj_mass_bins) {
    static unsigned i = 1;
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

  for (unsigned p=0, np=f1->GetNpar(); p<np; ++p) {
    const char* name = f1->GetParName(p);
    TH1D* h = new TH1D(name,name,bins.size()-1,bins.data());
    unsigned i = 1;
    for (const auto& b : hj_mass_bins) {
      if (skip.count(i)) ++i;
      h->SetBinContent(i,b.second->GetParameter(p));
      // h.SetBinError(i,b.second->GetParError(p));
      ++i;
    }
  }

  fout.Write();
}
