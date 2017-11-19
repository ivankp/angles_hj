#include <iostream>
#include <array>

#include <TFile.h>
#include <TTree.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "tc_msg.hh"
#include "math.hh"

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

  TTree *tout = new TTree("cos","");
  double cos;
  tout->Branch("cos",&cos);

  using cnt = timed_counter<Long64_t>;
  for (cnt ent(tin->GetEntries()); !!ent; ++ent) {
  // for (Long64_t nentries=tin->GetEntries(), i=0; i<nentries; ++i) {
    tin->GetEntry(ent);

    const double M = std::sqrt(
      sq(E[0]+E[1]) - sq(px[0]+px[1],py[0]+py[1],pz[0]+pz[1]) );

    if (M < Hj_mass_range[0]) continue;
    if (Hj_mass_range[1] <= M) continue;

    cos = pz[0]/std::sqrt(sq(px[0],py[0],pz[0]));
    tout->Fill();
  }

  fout.Write(0,TObject::kOverwrite);
}
