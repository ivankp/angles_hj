#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <tuple>
#include <memory>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TLorentzVector.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "minuit.hh"
#include "binner.hh"
#include "category_bin.hh"
#include "math.hh"
#include "Legendre.hh"
#include "float_or_double_reader.hh"

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
using std::get;
using boost::optional;
using namespace ivanp;
using namespace ivanp::math;

double weight = 1., cos_theta;

MAKE_ENUM(isp,(all)(gg)(gq)(qq))
isp get_isp(Int_t id1, Int_t id2) noexcept {
  const bool g1 = (id1 == 21), g2 = (id2 == 21);
  if (g1 == g2) return g1 ? isp::gg : isp::qq;
  else return isp::gq;
}

struct xw { double x, w; };
struct mass_bin {
  std::vector<xw> v;
  mass_bin(): v() { v.reserve(1<<10); }
  inline void operator()() {
    if (std::abs(cos_theta)>1.) return;
    v.push_back({cos_theta,weight});
  }
  inline auto begin() const noexcept { return v.begin(); }
  inline auto   end() const noexcept { return v.  end(); }
};

struct lo_bin {
  double w = 0, w2 = 0;
  long unsigned n = 0;
  inline void operator()(double weight) noexcept {
    w  += weight;
    w2 += weight*weight;
    ++n;
  }
};

template <typename T>
struct total { T all=0, selected=0; };

#define NPAR 4

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char *ofname, *cfname;
  const char* tree_name = "t3";
  int print_level = 0;
  unsigned prec = 10;

  struct {
    using _v = std::tuple<unsigned,double,double>;
    struct _p { std::string name; double init, step, a, b; };
    std::map<std::string,_v> v;
    std::array<_p,NPAR> p;
  } cfg;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input file",req(),pos())
      (ofname,'o',"output file",req())
      (cfname,'c',"config file",req())
      (tree_name,{"-t","--tree"},cat("input TTree name [",tree_name,']'))
      (prec,"--prec",cat("float precision [",prec,']'))
      (print_level,{"-v","--print-level"},
       "-1 - quiet (also suppress all warnings)\n"
       " 0 - normal (default)\n"
       " 1 - verbose",
       switch_init(1))
      .parse(argc,argv,true)) return 0;

    try {
      std::ifstream f(cfname);
      for (std::string str; std::getline(f,str); ) {
        if (str.empty() || str[0]=='#') continue;
        std::stringstream ss(str);
        ss >> str;
        if (str=="p") {
          ss >> str;
          auto& p = [&]() -> auto& {
            for (auto& p : cfg.p)
              if (p.name.empty()) { p.name = str; return p; }
            throw error("extra parameter ",str);
          }();
          ss >> p.init >> p.step >> p.a >> p.b;
        } else {
          auto& v = cfg.v[str];
          ss >> get<0>(v) >> get<1>(v) >> get<2>(v);
        }
      }
      if (cfg.p.back().name.empty()) throw error("missing parameters");
    } catch (const std::exception& e) {
      throw error("in ",cfname,": ",e.what());
    }
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  info("Input ntuples");
  for (const char* name : ifnames) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  // Set up branches for reading
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  float_or_double_value_reader _weight(reader,"weight2");

  optional<TTreeReaderValue<Int_t>> _ncount;

  TTreeReaderValue<Int_t> _id1(reader,"id1"), _id2(reader,"id2");

#define OPT_BRANCH(NAME) \
  if (!strcmp(bo->GetName(),#NAME)) _##NAME.emplace(reader,#NAME);

  for (auto bo : *reader.GetTree()->GetListOfBranches()) {
    OPT_BRANCH(ncount)
    if (_ncount) break;
  }

#undef OPT_BRANCH

  binner<category_bin<mass_bin,isp>, std::tuple<
    axis_spec<uniform_axis<double>, false, false> >
  > hj_mass_bins(cfg.v.at("M"));

  TLorentzVector Higgs;
  std::vector<TLorentzVector> jets;

  total<double> total_weight;
  total<long unsigned> total_entries, total_ncount;
  unsigned ncount;

  // LOOP ===========================================================
  using cnt = ivanp::timed_counter<Long64_t>;
  for (cnt ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    total_weight.all += (weight = *_weight);
    total_ncount.all += _ncount ? (ncount = **_ncount) : 1;
    ++total_entries.all;

    // Read particles -----------------------------------------------
    const unsigned np = *_nparticle;
    jets.clear();
    if (jets.capacity()==0) jets.reserve(np-1);
    for (unsigned i=0; i<np; ++i) {
      if (_kf[i]==25) {
        Higgs.SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
      } else {
        jets.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }

    decltype(hj_mass_bins)::bin_type::id<isp>()
      = (unsigned)get_isp(*_id1,*_id2);
    // --------------------------------------------------------------

    // Jets ---------------------------------------------------------
    std::sort( jets.begin(), jets.end(), [](const auto& a, const auto& b){
      return a.Pt() > b.Pt();
    });
    const auto& jet1 = jets.front();

    // Angles -------------------------------------------------------
    const auto Q = Higgs + jet1;
    const double Q2 = Q*Q;
    const double hj_mass = std::sqrt(Q2);

    const TLorentzVector Z(0,0,Q.E(),Q.Pz());
    const auto ell = ((Q*jet1)/Q2)*Higgs - ((Q*Higgs)/Q2)*jet1;

    cos_theta = (ell*Z) / std::sqrt(sq(ell)*sq(Z));

    // Fill ---------------------------------------------------------
    hj_mass_bins(hj_mass);
  } // end event loop

  // OUTPUT FILE ####################################################
  std::ofstream out(ofname);
  out.precision(prec);
  out << "[{"
    "\"weight\":["<<total_weight.all<<"],"
    "\"entries\":["<<total_entries.all<<"],"
    "\"ncount\":["<<total_ncount.all<<"]"
    "},{";

  { bool first = true;
  for (const auto& var : cfg.v) {
    if (!first) out << ',';
    else first = false;
    out << '\"' << var.first << "\":["
        << get<0>(var.second) << ','
        << get<1>(var.second) << ','
        << get<2>(var.second) << ']';
  }}

  out << "},";

  decltype(hj_mass_bins)::bin_type::id<isp>() = 0;

  // FITTING ########################################################
  double pars[NPAR+1], errs[NPAR+1];

  for (const auto& bin : hj_mass_bins) { // loop over bins
    static unsigned bin_i = 0;
    if (bin_i) out << ',';
    const std::string hj_mass_bin = hj_mass_bins.bin_str(++bin_i);
    info("Fitting hj_mass",hj_mass_bin);
    const auto& v = bin->v;
    info("Events",v.size());

    out << "[[" << hj_mass_bins.axis().lower(bin_i)
        << ',' << hj_mass_bins.axis().upper(bin_i) << "],";

    binner<lo_bin, std::tuple<
      axis_spec<uniform_axis<double>, false, false> >
    > h(cfg.v.at("cos"));
    const auto& axis = h.axis();

    const unsigned nbins = axis.nbins();

    std::vector<double> mid(nbins);
    for (unsigned i=0; i<nbins; ++i) {
      const double a = axis.lower(i+1).get(),
                   b = axis.upper(i+1).get();
      mid[i] = (a+b)*0.5;
    }

    for (const auto& e : v) h(e.x,e.w);

    double total_w = 0;
    for (const auto& b : h) total_w += b.w;

    out << '{';

    info("χ² fit"); // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    auto fChi2 = [&b=h.bins(),&mid](const double* c) -> double {
      double chi2 = 0.;
      const unsigned n = b.size();
      for (unsigned i=0; i<n; ++i)
        chi2 += sq(b[i].w - c[NPAR]*Legendre(&mid[i],c))/b[i].w2;
      return chi2;
    };

    minuit<decltype(fChi2)> mChi2(NPAR,fChi2);
    mChi2.SetPrintLevel(print_level);

    for (unsigned i=0; i<NPAR; ++i) {
      const auto& p = cfg.p[i];
      mChi2.DefineParameter(i, p.name.c_str(), p.init, p.step, p.a, p.b);
    }
    mChi2.DefineParameter(
      NPAR, "A", total_w, total_w*1e-2, total_w*0.1, total_w*10);

    mChi2.Migrad();

    out << "\"chi2\":{";
    { bool first = true;
    for (unsigned i=0; i<=NPAR; ++i) {
      if (!first) out << ',';
      else first = false;
      mChi2.GetParameter(i,pars[i],errs[i]);
      out <<'\"'<< mChi2.fCpnam[i] << "\":["
          << pars[i] <<','<< errs[i] << "]";
    }}
    out << "},";

    info("LogL fit"); // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    auto fLogL = [&v](const double* c) -> double {
      long double logl = 0.;
      const unsigned n = v.size();
      #pragma omp parallel for reduction(+:logl)
      for (unsigned i=0; i<n; ++i)
        logl += v[i].w*log(Legendre(&v[i].x,c));
      return -2.*logl;
    };

    minuit<decltype(fLogL)> mLogL(NPAR,fLogL);
    mLogL.SetPrintLevel(print_level);

    for (unsigned i=0; i<NPAR; ++i) {
      const auto& p = cfg.p[i];
      mLogL.DefineParameter(i, p.name.c_str(), pars[i], p.step, p.a, p.b);
    }

    mLogL.Migrad();

    out << "\"logl\":{";
    { bool first = true;
    for (unsigned i=0; i<NPAR; ++i) {
      if (!first) out << ',';
      else first = false;
      mChi2.GetParameter(i,pars[i],errs[i]);
      out <<'\"'<< mChi2.fCpnam[i] << "\":["
          << pars[i] <<','<< errs[i] << "]";
    }}
    out << "}";

    out << "},[";
    { bool first = true;
    for (const auto& b : h) {
      if (!first) out << ',';
      else first = false;
      out << '[' << b.w <<','<< b.w2 <<','<< b.n << ']';
    }}
    out << "]]";
  }
  out << ']';
}
