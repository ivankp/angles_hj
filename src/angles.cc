#include <iostream>
#include <array>

#include <TFile.h>
#include <TTree.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "math.hh"
#include "lorentz_vector.hh"

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;
using namespace ivanp::math;

template <typename... Args>
inline void write(const char* name, Args&&... args) {
  TNamed(name,cat(args...).c_str()).Write();
}

int main(int argc, char* argv[]) {
  const char *ifname, *ofname;
  std::array<double,2> Hj_mass_range;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        (ofname,'o',"output file",req())
        (Hj_mass_range,{"-m","--mass"},"Hj mass range",req())
        .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TFile fin(ifname);
  info("Input  file",fin.GetName());
  if (fin.IsZombie()) return 1;

  TTree *tin;
  fin.GetObject("events",tin);
  if (!tin) return 1;

  double px[2], py[2], pz[2], E[2];
  tin->SetBranchAddress("px",px);
  tin->SetBranchAddress("py",py);
  tin->SetBranchAddress("pz",pz);
  tin->SetBranchAddress("E",E);

  TFile fout(ofname,"recreate","",109);
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TTree *tout = new TTree("angles","");
  double cos_theta;
  tout->Branch("cos_theta",&cos_theta);

  for (timed_counter<Long64_t> ent(tin->GetEntries()); !!ent; ++ent) {
    tin->GetEntry(ent);

    const lorentz_vector pH {E[0],px[0],py[0],pz[0]},
                         pj {E[1],px[1],py[1],pz[1]};
    const auto Q = pH + pj;
    const double Q2 = Q*Q;
    const double M = sqrt(Q2);

    if (M < Hj_mass_range[0]) continue;
    if (Hj_mass_range[1] <= M) continue;

    const lorentz_vector Z {Q.z,0,0,Q.t};
    const auto ell = ((Q*pj)/Q2)*pH - ((Q*pH)/Q2)*pj;

    cos_theta = (ell*Z) / sqrt((ell*ell)*(Z*Z));
    tout->Fill();
  }

  write("M range",'[',Hj_mass_range[0],',',Hj_mass_range[1],')');

  fout.Write(0,TObject::kOverwrite);
}
