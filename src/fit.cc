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
  std::array<double,NPAR> pars_init{0,0,0,0};
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
        (fit_range,'r',cat("max cosθ fit range [",fit_range,']'))
        (pars_init,'p',"parameters' initial values")
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

  TF1 *fit = new TF1("fit-logl",Legendre,-1.,1.,NPAR);
  fit->SetLineColor(2);
  for (unsigned i=0; i<NPAR; ++i) {
    fit->SetParName(i,pars_names[i]);
    // fit->SetParLimits(i,-pars_lim[i],pars_lim[i]);
  }
  for (unsigned i=npar; i<NPAR; ++i)
    fit->FixParameter(i,pars_init[i]);

  TF1 *fit2 = new TF1("fit-chi2",
    [](const double* x, const double* c){ return c[0]*Legendre(x,c+1); },
    -1.,1.,NPAR+1);
  fit2->SetLineColor(418);
  fit2->SetParName(0,"A");
  for (unsigned i=0; i<NPAR; ++i) {
    fit2->SetParName(i+1,pars_names[i]);
    fit2->SetParameter(i+1,pars_init[i]);
  }
  fit2->SetParLimits(4,-1e-4,M_PI);
  for (unsigned i=npar; i<NPAR; ++i)
    fit2->FixParameter(i+1,pars_init[i]);

  if (nbins<100) {
    const int n = ceil(91./nbins)*nbins;
    fit ->SetNpx(n);
    fit2->SetNpx(n);
  }

  // Fit in mass bins ===============================================
  double pars[NPAR], errs[NPAR];

  for (const auto& bin : hj_mass_bins) {
    static unsigned bin_i = 0;
    const std::string hj_mass_bin = hj_mass_bins.bin_str(++bin_i);
    info("Fitting hj_mass",hj_mass_bin);
    info("Events",bin.v.size());

    bin.h->SetName(("cos_theta-hj_mass"+hj_mass_bin).c_str());
    bin.h->SetTitle(("hj_mass "+hj_mass_bin).c_str());
    bin.h->SetXTitle(cat("cos #theta / ",fit_range).c_str());

    auto LogL = [&v=bin.v](const double* c) -> double {
      long double logl = 0.;
      const unsigned n = v.size();
      #pragma omp parallel for reduction(+:logl)
      for (unsigned i=0; i<n; ++i)
        logl += v[i].w*log(Legendre(&v[i].x,c));
      return -2.*logl;
    };

    info("χ² fit");
    fit2->SetParameter(0,bin.h->Integral(1,bin.h->GetNbinsX()+1));
    auto result = bin.h->Fit(fit2,"SR0");
    TF1 *f = static_cast<TF1*>(bin.h->GetListOfFunctions()->At(0));
    f->SetTitle(cat(
        std::setprecision(15),std::scientific,
        "#chi^{2} = ",result->Chi2(),","
        "-2LogL = ",LogL(fit2->GetParameters()+1)
      ).c_str());
    // f->SetParameter(4, npar>3 ? mod_phi(f->GetParameter(4)) : 0.);

    info("LogL fit");
    minuit<decltype(LogL)> m(NPAR,LogL);
    m.SetPrintLevel(print_level);

    for (unsigned i=0; i<NPAR; ++i)
      m.DefineParameter(
        i,             // parameter number
        pars_names[i], // parameter name
        use_chi2_pars ? f->GetParameter(i+1) : pars_init[i],  // start value
        0.01,          // step size
        i==3 ? -1e-4 : 0, // mininum
        i==3 ?  M_PI : 0  // maximum
      );

    switch (npar) {
      case 0: m.FixParameter(0);
      case 1: m.FixParameter(1);
      case 2: m.FixParameter(2);
      case 3: m.FixParameter(3);
      default: break;
    }

    m.Migrad();
    for (unsigned i=0; i<NPAR; ++i)
      m.GetParameter(i,pars[i],errs[i]);

    fit->SetParameters(pars);
    fit->SetParErrors(errs);
    fit->SetTitle(cat(
      std::setprecision(17),std::scientific,
      "-2LogL = ",LogL(pars)).c_str());
    bin.h->GetListOfFunctions()->Add(fit);

    bin.h->Write();
  }

  info("Saving",fout.GetName());
  fout.Write(0,TObject::kOverwrite);
}
