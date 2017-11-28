#include <iostream>
#include <array>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "math.hh"
#include "minuit.hh"
#include "Legendre.hh"

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;
using namespace ivanp::math;

template <typename... Args>
inline void write(const char* name, Args&&... args) {
  TNamed(name,cat(args...).c_str()).Write();
}

#define NPAR 4

int main(int argc, char* argv[]) {
  const char *ifname, *ofname;
  // std::array<double,2> range {0,1};
  unsigned npar = NPAR, nbins = 100;
  std::array<double,NPAR> pars_init {0,0,0,0}, pars_lim {1,1,1,M_PI};

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        (ofname,'o',"output file",req())
        // (range,'r',cat("|cos_theta Θ| range [",range[0],',',range[1],']'))
        (npar,'n',cat("number of fit parameters [",npar,']'))
        (pars_init,'p',"parameters' initial values")
        (pars_lim,'l',"parameters' limits")
        (nbins,"--nbins",cat('[',nbins,']'))
        .parse(argc,argv,true)) return 0;
    if (npar>NPAR) throw std::runtime_error("npar > " STR(NPAR));
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TFile fin(ifname);
  info("Input file",fin.GetName());
  if (fin.IsZombie()) return 1;

  TTree *tin;
  fin.GetObject("angles",tin);
  if (!tin) return 1;

  double cos_theta;
  tin->SetBranchAddress("cos_theta",&cos_theta);

  TFile fout(ofname,"recreate");
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TH1D *h = new TH1D("abs_cos_theta","|cos_theta #theta|",nbins,0,1);

  std::vector<double> vals;
  const Long64_t nent = tin->GetEntries();
  vals.reserve(nent);

  for (timed_counter<Long64_t> ent(nent); !!ent; ++ent) {
    tin->GetEntry(ent);
    const double x = std::abs(cos_theta);
    // if (x>range[1]) continue;
    vals.push_back(x);
    h->Fill(x);
  }
  info("Number of points",vals.size());

  auto LogL = [&](double* c){
    double logl = 0.;
    for (double x : vals) logl += std::log(Legendre(&x,c));
    return -2.*logl;
  };

  minuit<decltype(LogL)> m(NPAR,LogL);

  // parameter number
  // parameter name
  // start value
  // step size
  // mininum
  // maximum
  m.DefineParameter(0,"c2",  pars_init[0],0.1,-pars_lim[0],pars_lim[0]);
  m.DefineParameter(1,"c4",  pars_init[1],0.1,-pars_lim[1],pars_lim[1]);
  m.DefineParameter(2,"c6",  pars_init[2],0.1,-pars_lim[2],pars_lim[2]);
  m.DefineParameter(3,"phi2",pars_init[3],0.1,-pars_lim[3],pars_lim[3]);

  switch (npar) {
    case 0: m.FixParameter(0);
    case 1: m.FixParameter(1);
    case 2: m.FixParameter(2);
    case 3: m.FixParameter(3);
    default: break;
  }

  m.Migrad();
  double pars[NPAR], errs[NPAR];
  for (unsigned i=0; i<NPAR; ++i)
    m.GetParameter(i,pars[i],errs[i]);

  TF1* tf = new TF1("fit", Legendre, 0, 1, npar);
  tf->SetParameters(pars);
  tf->Write();
  // h->GetListOfFunctions()->Add(tf);

  write("-2LogL",LogL(pars));
  write("c2",  pars[0]," #pm ",errs[0]);
  write("c4",  pars[1]," #pm ",errs[1]);
  write("c6",  pars[2]," #pm ",errs[2]);
  write("phi2",pars[3]," #pm ",errs[3]);

  fout.Write(0,TObject::kOverwrite);
}
