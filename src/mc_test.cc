#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <random>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "minuit.hh"
#include "math.hh"
#include "random.hh"
#include "Legendre.hh"

#define TEST(VAR) \
  std::cout << tc::cyan << #VAR << tc::reset << " = " << VAR << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;
using namespace ivanp::math;

#define NPAR 4
const char* pars_names[NPAR] = {"c2","c4","c6","#phi2"};

/*
double testf(const double* x, const double* c) {
  const double _x = *x;
  return (1. + c[0]*_x + c[1]*_x*_x)/(2. + (2./3.)*c[1]);
}
double testf2(const double* x, const double* c) { return c[0]*testf(x,c+1); }
*/

int main(int argc, char* argv[]) {
  const char* ofname;
  long unsigned nevents = 100000;
  unsigned npar = NPAR, nbins = 100;
  std::array<double,NPAR> coeffs{}, pars_init{}, pars_lim{1,1,1,M_PI};
  double fit_range = 1.;
  int print_level = 0;
  auto seed = std::mt19937::default_seed;
  bool use_chi2_pars = false;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ofname,'o',"output file",req())
        (coeffs,'c',"coefficients",pos(1),req())
        (nevents,'N',cat("number of MC events [",nevents,']'),pos(1))
        (npar,'n',cat("number of fit parameters [",npar,']'))
        (fit_range,'r',cat("max |cos θ| fit range [",fit_range,']'))
        (pars_init,'p',"parameters' initial values")
        (pars_lim,'l',"parameters' limits")
        (nbins,"--nbins",cat('[',nbins,']'))
        (print_level,"--print-level",
         "-1 - quiet (also suppress all warnings)\n"
         " 0 - normal (default)\n"
         " 1 - verbose")
        (seed,"--seed",cat('[',seed,']'))
        (use_chi2_pars,"--use-chi2-pars")
        .parse(argc,argv,true)) return 0;
    if (npar>NPAR) throw std::runtime_error("npar > " STR(NPAR));
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  // Output file ====================================================
  TFile fout(ofname,"recreate");
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;

  TH1::AddDirectory(false);

  TH1D *h = new TH1D("cos_theta-[test)","test",nbins,-1,1);
  h->SetXTitle("|cos #theta|");

  // generate =======================================================
  info("SEED",seed);
  std::mt19937 gen(seed); // mersenne twister random number generator
  auto dist = sample(-fit_range,fit_range,10.,
    [=](double x){ return Legendre(&x,&coeffs.front()); });

  std::vector<double> v;
  v.reserve(nevents);
  for (timed_counter<decltype(nevents)> i(nevents); !!i; ++i) {
    const auto x = dist(gen);
    v.push_back(x);
    h->Fill(x);
  }

  h->Write();

  // Fit in mass bins ===============================================
  TF1 *fit = new TF1("fit-logl",Legendre,-1.,1.,npar);
  // fit->SetNpx(nbins);
  for (unsigned i=0; i<NPAR; ++i)
    fit->SetParName(i,pars_names[i]);

  TF1 *fit2 = new TF1("fit2-chi2",
    [](const double* x, const double* c){ return c[0]*Legendre(x,c+1); },
    -1.,1.,npar+1);
  // fit2->SetNpx(nbins);
  fit2->SetParName(0,"A");
  for (unsigned i=0; i<NPAR; ++i)
    fit2->SetParName(i+1,pars_names[i]);

  info("χ² fit");
  fit2->SetParameter(0,h->Integral(1,h->GetNbinsX()+1));
  for (unsigned i=0; i<npar; ++i)
    fit2->SetParameter(i+1,pars_init[i]);
  auto result = h->Fit(fit2,"S","",-1.,1.);
  fit2->SetLineColor(3);
  fit2->SetTitle(cat("#chi^{2} = ",result->Chi2()).c_str());
  fit2->SetName("fit-chi2-[test)");
  fit2->Write();
  fit2->SetName("fit2");

  if (use_chi2_pars) {
    for (unsigned i=0; i<npar; ++i)
      pars_init[i] = fit2->GetParameter(i+1);
    mod_phi_ref(pars_init[3]);
  }

  info("LogL fit");
  auto LogL = [=,&v](const double* c){
    double logl = 0.;
    for (double x : v)
      logl += log(Legendre(&x,c));
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
  fit->SetName("fit-logl-[test)");
  fit->Write();
  fit->SetName("fit");

  info("Saving",fout.GetName());
  fout.Write(0,TObject::kOverwrite);
}
