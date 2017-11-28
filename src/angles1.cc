#include <iostream>
#include <fstream>
#include <array>

#include <TFile.h>
#include <TTree.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "math.hh"
#include "lorentz_vector.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

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

  int pid;
  lorentz_vector pH, pj;
  // constexpr lorentz_vector p1{1,0,0,1}, p2{1,0,0,-1};

  TFile fout(ofname,"recreate","",109);
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TTree *tout = new TTree("angles","");
  double hj_mass, cos_theta;
  tout->Branch("hj_mass",&hj_mass);
  tout->Branch("cos_theta",&cos_theta);

  info("Input file",ifname);
  write("Input file",ifname);
  std::ifstream fin(ifname);
  for (
    timed_counter<unsigned> cnt;
    fin >> pid >> pH.x >> pH.y >> pH.z >> pH.t
        >> pid >> pj.x >> pj.y >> pj.z >> pj.t;
    ++cnt
  ) {
    const auto Q = pH + pj;
    const double Q2 = Q*Q;
    hj_mass = sqrt(Q2);

    // const auto Z = ((Q*p2)/Q2)*p1 - ((Q*p1)/Q2)*p2;
    const lorentz_vector Z {Q.z,0,0,Q.t};
    const auto ell = ((Q*pj)/Q2)*pH - ((Q*pH)/Q2)*pj;

    cos_theta = (ell*Z) / sqrt(sq(ell)*sq(Z));

    tout->Fill();
  }

  info("Saving",fout.GetName());
  fout.Write(0,TObject::kOverwrite);
}
