// Author: Yipeng Sun
// Last Change: Thu Oct 08, 2020 at 10:28 PM +0800
// Description: FF reweighting for R(D(*)) run 1, step 1 ntuples.
// Based on:
//   https://github.com/ZishuoYang/my-hammer-reweighting/blob/master/Bc2JpsiMuNu.cc

#include <Hammer/Hammer.hh>
#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>
#include <Hammer/Process.hh>
#include <Hammer/Tools/HammerRoot.hh>

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <string>

using namespace std;

auto particle(Double_t pe, Double_t px, Double_t py, Double_t pz, Int_t pid) {
  auto four_mom = Hammer::FourMomentum(pe, px, py, pz);
  auto part_id  = static_cast<Hammer::PdgId>(pid);

  return Hammer::Particle(four_mom, part_id);
}

void reweight_dst(TFile* input_file, TFile* output_file,
                  const char* tree        = "mc_dst_tau_aux",
                  const char* tree_output = "mc_dst_tau_ff_w") {
  TTreeReader reader(tree, input_file);
  TTree       output(tree_output, tree_output);

  // Read input branches ///////////////////////////////////////////////////////
  // General
  TTreeReaderValue<ULong64_t> eventNumber(reader, "eventNumber");
  TTreeReaderValue<UInt_t>    runNumber(reader, "runNumber");
  // B
  TTreeReaderValue<Int_t>    b_id(reader, "b_id");
  TTreeReaderValue<Double_t> b_true_pe(reader, "b_true_pe");
  TTreeReaderValue<Double_t> b_true_px(reader, "b_true_px");
  TTreeReaderValue<Double_t> b_true_py(reader, "b_true_py");
  TTreeReaderValue<Double_t> b_true_pz(reader, "b_true_pz");
  // D*
  TTreeReaderValue<Int_t>    dst_id(reader, "dst_id");
  TTreeReaderValue<Double_t> dst_true_pe(reader, "dst_true_pe");
  TTreeReaderValue<Double_t> dst_true_px(reader, "dst_true_px");
  TTreeReaderValue<Double_t> dst_true_py(reader, "dst_true_py");
  TTreeReaderValue<Double_t> dst_true_pz(reader, "dst_true_pz");
  // D0
  TTreeReaderValue<Int_t>    d0_id(reader, "d0_id");
  TTreeReaderValue<Double_t> d0_true_pe(reader, "d0_true_pe");
  TTreeReaderValue<Double_t> d0_true_px(reader, "d0_true_px");
  TTreeReaderValue<Double_t> d0_true_py(reader, "d0_true_py");
  TTreeReaderValue<Double_t> d0_true_pz(reader, "d0_true_pz");
  // Mu
  TTreeReaderValue<Int_t>    mu_id(reader, "mu_id");
  TTreeReaderValue<Double_t> mu_true_pe(reader, "mu_true_pe");
  TTreeReaderValue<Double_t> mu_true_px(reader, "mu_true_px");
  TTreeReaderValue<Double_t> mu_true_py(reader, "mu_true_py");
  TTreeReaderValue<Double_t> mu_true_pz(reader, "mu_true_pz");
  // K
  TTreeReaderValue<Int_t>    k_id(reader, "k_id");
  TTreeReaderValue<Double_t> k_true_pe(reader, "k_true_pe");
  TTreeReaderValue<Double_t> k_true_px(reader, "k_true_px");
  TTreeReaderValue<Double_t> k_true_py(reader, "k_true_py");
  TTreeReaderValue<Double_t> k_true_pz(reader, "k_true_pz");
  // Pi
  TTreeReaderValue<Int_t>    pi_id(reader, "pi_id");
  TTreeReaderValue<Double_t> pi_true_pe(reader, "pi_true_pe");
  TTreeReaderValue<Double_t> pi_true_px(reader, "pi_true_px");
  TTreeReaderValue<Double_t> pi_true_py(reader, "pi_true_py");
  TTreeReaderValue<Double_t> pi_true_pz(reader, "pi_true_pz");
  // Slow Pi
  TTreeReaderValue<Int_t>    spi_id(reader, "spi_id");
  TTreeReaderValue<Double_t> spi_true_pe(reader, "spi_true_pe");
  TTreeReaderValue<Double_t> spi_true_px(reader, "spi_true_px");
  TTreeReaderValue<Double_t> spi_true_py(reader, "spi_true_py");
  TTreeReaderValue<Double_t> spi_true_pz(reader, "spi_true_pz");
  // Tau
  TTreeReaderValue<Int_t>    tau_id(reader, "tau_id");
  TTreeReaderValue<Double_t> tau_true_pe(reader, "tau_true_pe");
  TTreeReaderValue<Double_t> tau_true_px(reader, "tau_true_px");
  TTreeReaderValue<Double_t> tau_true_py(reader, "tau_true_py");
  TTreeReaderValue<Double_t> tau_true_pz(reader, "tau_true_pz");
  // Anti-Nu_Tau
  TTreeReaderValue<Int_t>    anu_tau_id(reader, "anu_tau_id");
  TTreeReaderValue<Double_t> anu_tau_true_pe(reader, "anu_tau_true_pe");
  TTreeReaderValue<Double_t> anu_tau_true_px(reader, "anu_tau_true_px");
  TTreeReaderValue<Double_t> anu_tau_true_py(reader, "anu_tau_true_py");
  TTreeReaderValue<Double_t> anu_tau_true_pz(reader, "anu_tau_true_pz");
  // Nu_Tau
  TTreeReaderValue<Int_t>    nu_tau_id(reader, "nu_tau_id");
  TTreeReaderValue<Double_t> nu_tau_true_pe(reader, "nu_tau_true_pe");
  TTreeReaderValue<Double_t> nu_tau_true_px(reader, "nu_tau_true_px");
  TTreeReaderValue<Double_t> nu_tau_true_py(reader, "nu_tau_true_py");
  TTreeReaderValue<Double_t> nu_tau_true_pz(reader, "nu_tau_true_pz");
  // Anti-Nu_Mu
  TTreeReaderValue<Int_t>    anu_mu_id(reader, "anu_mu_id");
  TTreeReaderValue<Double_t> anu_mu_true_pe(reader, "anu_mu_true_pe");
  TTreeReaderValue<Double_t> anu_mu_true_px(reader, "anu_mu_true_px");
  TTreeReaderValue<Double_t> anu_mu_true_py(reader, "anu_mu_true_py");
  TTreeReaderValue<Double_t> anu_mu_true_pz(reader, "anu_mu_true_pz");

  // Define output branches ////////////////////////////////////////////////////
  ULong64_t eventNumber_out;
  output.Branch("eventNumber", &eventNumber_out);
  UInt_t runNumber_out;
  output.Branch("runNumber", &runNumber_out);
  Int_t ham_proc_id_out;
  output.Branch("ham_proc_id", &ham_proc_id_out);

  // Setup HAMMER //////////////////////////////////////////////////////////////
  Hammer::Hammer   ham{};
  Hammer::IOBuffer ham_buf;

  auto semi_tau_decay = vector<string>{"BD*TauNu", "TauEllNuNu"};
  // auto semi_tau_decay = vector<string>{"BD*TauNu"};

  ham.includeDecay(semi_tau_decay);
  ham.addFFScheme("Scheme1", {{"BD*", "BGL"}});
  ham.setOptions("BctoJpsiBGL: {dvec: [0., 0., 0.] }");
  ham.setFFInputScheme({{"BD*", "CLN"}});

  ham.setUnits("GeV");  // This is the default, but let's be explicit.

  ham.initRun();

  while (reader.Next()) {
    eventNumber_out = *eventNumber;
    runNumber_out   = *runNumber;

    auto B   = particle(*b_true_pe, *b_true_px, *b_true_py, *b_true_pz, *b_id);
    auto Dst = particle(*dst_true_pe, *dst_true_px, *dst_true_py, *dst_true_pz,
                        *dst_id);
    auto SlowPi = particle(*spi_true_pe, *spi_true_px, *spi_true_py,
                           *spi_true_pz, *spi_id);
    auto D0 =
        particle(*d0_true_pe, *d0_true_px, *d0_true_py, *d0_true_pz, *d0_id);
    auto K = particle(*k_true_pe, *k_true_px, *k_true_py, *k_true_pz, *k_id);
    auto Pi =
        particle(*pi_true_pe, *pi_true_px, *pi_true_py, *pi_true_pz, *pi_id);
    auto Mu =
        particle(*mu_true_pe, *mu_true_px, *mu_true_py, *mu_true_pz, *mu_id);
    auto Tau = particle(*tau_true_pe, *tau_true_px, *tau_true_py, *tau_true_pz,
                        *tau_id);
    auto Anti_Nu_Mu = particle(*anu_mu_true_pe, *anu_mu_true_px,
                               *anu_mu_true_py, *anu_mu_true_pz, *anu_mu_id);
    auto Anti_Nu_Tau =
        particle(*anu_tau_true_pe, *anu_tau_true_px, *anu_tau_true_py,
                 *anu_tau_true_pz, *anu_tau_id);
    auto Nu_Tau = particle(*nu_tau_true_pe, *nu_tau_true_px, *nu_tau_true_py,
                           *nu_tau_true_pz, *nu_tau_id);

    Hammer::Process proc;

    auto B_idx           = proc.addParticle(B);
    auto Dst_idx         = proc.addParticle(Dst);
    auto SlowPi_idx      = proc.addParticle(SlowPi);
    auto D0_idx          = proc.addParticle(D0);
    auto K_idx           = proc.addParticle(K);
    auto Pi_idx          = proc.addParticle(Pi);
    auto Mu_idx          = proc.addParticle(Mu);
    auto Tau_idx         = proc.addParticle(Tau);
    auto Anti_Nu_Mu_idx  = proc.addParticle(Anti_Nu_Mu);
    auto Anti_Nu_Tau_idx = proc.addParticle(Anti_Nu_Tau);
    auto Nu_Tau_idx      = proc.addParticle(Nu_Tau);

    proc.addVertex(B_idx, {Dst_idx, Tau_idx, Anti_Nu_Tau_idx});
    proc.addVertex(Tau_idx, {Mu_idx, Nu_Tau_idx, Anti_Nu_Mu_idx});
    proc.addVertex(Dst_idx, {D0_idx, SlowPi_idx});
    proc.addVertex(D0_idx, {K_idx, Pi_idx});

    ham.initEvent();
    auto proc_id = ham.addProcess(proc);

    if (proc_id != 0) {
      ham_proc_id_out = proc_id;

      ham.processEvent();
      ham_buf = ham.saveEventWeights();

      output.Fill();
    }
  }

  output_file->Write("", TObject::kOverwrite);
}

int main(int, char** argv) {
  TFile* input_file  = new TFile(argv[1], "read");
  TFile* output_file = new TFile(argv[2], "recreate");

  reweight_dst(input_file, output_file);

  delete input_file;
  delete output_file;

  return 0;
}
