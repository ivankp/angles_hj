#include <iostream>
#include <array>
#include <vector>
#include <tuple>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "minuit.hh"
#include "binner.hh"
#include "math.hh"

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;
using namespace ivanp::math;

double weight, hj_mass, cos_theta, abs_cos_theta;

struct wx { double w, x; };
struct mass_bin {
  std::vector<wx> v;
  TH1D * const h;
  mass_bin(): v(), h(new TH1D("","",100,0,1)) { v.reserve(1024); }
  inline void operator()() {
    v.push_back({weight,abs_cos_theta});
    h->Fill(abs_cos_theta,weight);
  }
  inline auto begin() const noexcept { return v.begin(); }
  inline auto   end() const noexcept { return v.  end(); }
};

#define NPAR 4

double fitf(const double* x, const double* c) {
  const double x2 = sq(*x), x4 = x2*x2, x6 = x4*x2;

  const double p2 = 1.5*x2 - 0.5;
  const double p4 = 4.375*x4 - 3.75*x2 + 0.375;
  const double p6 = 14.4375*x6 - 19.6875*x4 + 6.5625*x2 - 0.3125;

  const double c0 = std::sqrt(
    1. - (0.2*sq(c[0]) + (1./9.)*sq(c[1]) + (1./13.)*sq(c[2])) );

  const auto phase = exp(std::polar<double>(1.,c[3]));

  return norm( c0 + c[0]*phase*p2 + c[1]*p4 + c[2]*p6 );
}

int main(int argc, char* argv[]) {
  const char *ifname, *ofname;
  unsigned npar = NPAR, nbins = 100;
  std::array<double,NPAR> pars_init {0,0,0,0}, pars_lim {1,1,1,M_PI};
  std::tuple<unsigned,double,double> hj_mass_binning;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        (ofname,'o',"output file",req())
        (hj_mass_binning,'M',"Higgs+jet mass binning",req())
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

  tin->SetBranchAddress("weight",&weight);
  tin->SetBranchAddress("hj_mass",&hj_mass);
  tin->SetBranchAddress("cos_theta",&cos_theta);

  TH1::AddDirectory(false);

  TFile fout(ofname,"recreate");
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  binner<mass_bin, std::tuple<
    axis_spec<uniform_axis<double>, false, false> >
  > hj_mass_bins(hj_mass_binning);

  // tree LOOP ======================================================
  for (timed_counter<Long64_t> ent(tin->GetEntries()); !!ent; ++ent) {
    tin->GetEntry(ent);
    abs_cos_theta = std::abs(cos_theta);
    hj_mass_bins(hj_mass);
  }

  // Fit in mass bins ===============================================
  auto range = [i=0u,&ax=hj_mass_bins.axis()]() mutable {
    return ++i, cat('[',ax.lower(i),',',ax.upper(i),')');
  };

  for (const auto& bin : hj_mass_bins) {
    const auto hj_mass_range = range();
    info("Fitting hj_mass",hj_mass_range);
    info("Events",bin.v.size());

    bin.h->SetName(("abs_cos_theta-hj_mass"+hj_mass_range).c_str());
    bin.h->SetTitle(("hj_mass "+hj_mass_range).c_str());
    bin.h->SetXTitle("|cos #theta|");
    bin.h->Write();

    auto LogL = [&v = bin.v](const double* c){
      double logl = 0.;
      for (const auto& e : v) logl += e.w*log(fitf(&e.x,c));
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

    TF1* fit = new TF1(("fit-hj_mass"+hj_mass_range).c_str(), fitf, 0, 1, npar);
    fit->SetParameters(pars);
    fit->SetParErrors(errs);
    fit->Write();
  }

  fout.Write(0,TObject::kOverwrite);
}
