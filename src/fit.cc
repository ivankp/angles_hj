#include <iostream>
#include <array>
#include <vector>
#include <tuple>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "minuit.hh"
#include "binner.hh"
#include "math.hh"
#include "Legendre.hh"

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;
using namespace ivanp::math;

unsigned nbins = 100;
double weight = 1., hj_mass, cos_theta;
double fit_scale = 1.;

struct wx { double w, x; };
struct mass_bin {
  std::vector<wx> v;
  TH1D * const h;
  mass_bin(): v(), h(new TH1D("","",nbins,-1.,1.)) { v.reserve(1024); }
  inline void operator()() {
    cos_theta *= fit_scale;
    if (std::abs(cos_theta)>1.) return;
    v.push_back({weight,cos_theta});
    h->Fill(cos_theta,weight);
  }
  inline auto begin() const noexcept { return v.begin(); }
  inline auto   end() const noexcept { return v.  end(); }
};

#define NPAR 4
const char* pars_names[NPAR] = {"c2","c4","c6","#phi2"};

int main(int argc, char* argv[]) {
  const char *ifname, *ofname;
  unsigned npar = NPAR;
  std::array<double,NPAR> pars_init{0,0,0,0}, pars_lim{1,1,1,M_PI};
  std::tuple<unsigned,double,double> hj_mass_binning;
  double fit_range = 1.;
  int print_level = 0;
  bool use_chi2_pars = false;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        (ofname,'o',"output file",req())
        (hj_mass_binning,'M',"Higgs+jet mass binning",req())
        (npar,'n',cat("number of fit parameters [",npar,']'))
        (fit_range,'r',cat("max |cos θ| fit range [",fit_range,']'))
        (pars_init,'p',"parameters' initial values")
        (pars_lim,'l',"parameters' limits")
        (use_chi2_pars,"--use-chi2-pars")
        (nbins,"--nbins",cat('[',nbins,']'))
        (print_level,"--print-level",
         "-1 - quiet (also suppress all warnings)\n"
         " 0 - normal (default)\n"
         " 1 - verbose")
        .parse(argc,argv,true)) return 0;
    if (npar>NPAR) throw std::runtime_error("npar > " STR(NPAR));
    if (fit_range > 1.) throw std::runtime_error("fit range > 1");
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  fit_scale = 1./fit_range;

  TFile fin(ifname);
  info("Input file",fin.GetName());
  if (fin.IsZombie()) return 1;

  TTree *tin;
  fin.GetObject("angles",tin);
  if (!tin) return 1;

  tin->SetBranchAddress("hj_mass",&hj_mass);
  tin->SetBranchAddress("cos_theta",&cos_theta);
  for (auto b : *tin->GetListOfBranches()) {
    if (!strcmp(b->GetName(),"weight")) {
      tin->SetBranchAddress("weight",&weight);
      break;
    }
  }

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
    hj_mass_bins(hj_mass);
    // if (ent > 1e6) break;
  }

  // Fit in mass bins ===============================================
  TF1 *fit = new TF1("fit",Legendre,-1.,1.,npar);
  // fit->SetNpx(nbins);
  for (unsigned i=0; i<NPAR; ++i)
    fit->SetParName(i,pars_names[i]);

  TF1 *fit2 = new TF1("fit2",
    [](const double* x, const double* c){ return c[0]*Legendre(x,c+1); },
    -1.,1.,npar+1);
  // fit2->SetNpx(nbins);
  fit2->SetParName(0,"A");
  for (unsigned i=0; i<NPAR; ++i)
    fit2->SetParName(i+1,pars_names[i]);

  for (const auto& bin : hj_mass_bins) {
    static unsigned bin_i = 0;
    const std::string hj_mass_bin = hj_mass_bins.bin_str(++bin_i);
    info("Fitting hj_mass",hj_mass_bin);
    info("Events",bin.v.size());

    bin.h->SetName(("cos_theta-hj_mass"+hj_mass_bin).c_str());
    bin.h->SetTitle(("hj_mass "+hj_mass_bin).c_str());
    bin.h->SetXTitle("|cos #theta|");
    bin.h->Write();

    info("χ² fit");
    fit2->SetParameter(0,bin.h->Integral(1,bin.h->GetNbinsX()+1));
    for (unsigned i=0; i<NPAR; ++i)
      fit2->SetParameter(i+1,pars_init[i]);
    auto result = bin.h->Fit(fit2,"S","",-1.,1.);
    fit2->SetLineColor(8);
    fit2->SetTitle(cat("#chi^{2} = ",result->Chi2()).c_str());
    fit2->SetName(("fit-chi2-hj_mass"+hj_mass_bin).c_str());
    fit2->Write();
    fit2->SetName("fit2");

    if (use_chi2_pars) {
      for (unsigned i=0; i<npar; ++i)
        pars_init[i] = fit2->GetParameter(i+1);
      mod_phi_ref(pars_init[3]);
    }

    info("LogL fit");
    auto LogL = [=,&v=bin.v](const double* c){
      double logl = 0.;
      // for (const auto& e : v)
      const unsigned n = v.size();
      #pragma omp parallel for reduction(+:logl)
      for (unsigned i=0; i<n; ++i)
        logl += v[i].w*log(Legendre(&v[i].x,c));
      return -2.*logl;
    };

    minuit<decltype(LogL)> m(NPAR,LogL);
    m.SetPrintLevel(print_level);

    for (unsigned i=0; i<NPAR; ++i)
      m.DefineParameter(
        i,             // parameter number
        pars_names[i], // parameter name
        pars_init[i],  // start value
        0.1,           // step size
        -pars_lim[i],  // mininum
        pars_lim[i]    // maximum
      );

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

    fit->SetParameters(pars);
    fit->SetParErrors(errs);
    fit->SetLineColor(2);
    fit->SetTitle(cat("-2LogL = ",LogL(pars)).c_str());
    fit->SetName(("fit-logl-hj_mass"+hj_mass_bin).c_str());
    fit->Write();
    fit->SetName("fit");
  }

  info("Saving",fout.GetName());
  fout.Write(0,TObject::kOverwrite);
}
