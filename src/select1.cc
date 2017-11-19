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

  int pid;
  lorentz_vector pH, pj;
  // constexpr lorentz_vector p1{1,0,0,1}, p2{1,0,0,-1};

  TFile fout(ofname,"recreate","",109);
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TTree *tout = new TTree("angles","");
  double cos_theta;
  tout->Branch("cos_theta",&cos_theta);

  std::ifstream fin(ifname);
  for (
    timed_counter<unsigned> cnt;
    fin >> pid >> pH.x >> pH.y >> pH.z >> pH.t
        >> pid >> pj.x >> pj.y >> pj.z >> pj.t;
    ++cnt
  ) {
    const auto Q = pH + pj;
    const double Q2 = Q*Q;
    const double M = sqrt(Q2);

    if (M < Hj_mass_range[0]) continue;
    if (Hj_mass_range[1] <= M) continue;

    // const auto Z   = ((Q*p2)/Q2)*p1 - ((Q*p1)/Q2)*p2;
    const lorentz_vector Z {Q.z,0,0,Q.t};
    const auto ell = ((Q*pj)/Q2)*pH - ((Q*pH)/Q2)*pj;

    cos_theta = (ell*Z) / sqrt((ell*ell)*(Z*Z));

    tout->Fill();
  }

  fout.Write(0,TObject::kOverwrite);
}
