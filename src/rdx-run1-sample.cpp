// Author: Yipeng Sun
// Last Change: Mon Oct 05, 2020 at 09:24 PM +0800
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
                  const char* tree = "mc_dst_tau_aux") {
  TTreeReader reader(tree, input_file);
  TTree       output(tree, tree);

  // Read input branches ///////////////////////////////////////////////////////
  // General
  TTreeReaderValue<ULong64_t> eventNumber(reader, "eventNumber");
  TTreeReaderValue<UInt_t>    runNumber(reader, "runNumber");
  // B
  TTreeReaderValue<Double_t> b_true_pe(reader, "b_true_pe");
  TTreeReaderValue<Double_t> b_true_px(reader, "b_true_px");
  TTreeReaderValue<Double_t> b_true_py(reader, "b_true_py");
  TTreeReaderValue<Double_t> b_true_pz(reader, "b_true_pz");
  // D*
  TTreeReaderValue<Double_t> dst_true_pe(reader, "dst_true_pe");
  TTreeReaderValue<Double_t> dst_true_px(reader, "dst_true_px");
  TTreeReaderValue<Double_t> dst_true_py(reader, "dst_true_py");
  TTreeReaderValue<Double_t> dst_true_pz(reader, "dst_true_pz");
  // D0
  TTreeReaderValue<Double_t> d0_true_pe(reader, "d0_true_pe");
  TTreeReaderValue<Double_t> d0_true_px(reader, "d0_true_px");
  TTreeReaderValue<Double_t> d0_true_py(reader, "d0_true_py");
  TTreeReaderValue<Double_t> d0_true_pz(reader, "d0_true_pz");
  // Mu
  TTreeReaderValue<Double_t> mu_true_pe(reader, "mu_true_pe");
  TTreeReaderValue<Double_t> mu_true_px(reader, "mu_true_px");
  TTreeReaderValue<Double_t> mu_true_py(reader, "mu_true_py");
  TTreeReaderValue<Double_t> mu_true_pz(reader, "mu_true_pz");
  // K
  TTreeReaderValue<Double_t> k_true_pe(reader, "k_true_pe");
  TTreeReaderValue<Double_t> k_true_px(reader, "k_true_px");
  TTreeReaderValue<Double_t> k_true_py(reader, "k_true_py");
  TTreeReaderValue<Double_t> k_true_pz(reader, "k_true_pz");
  // Pi
  TTreeReaderValue<Double_t> pi_true_pe(reader, "pi_true_pe");
  TTreeReaderValue<Double_t> pi_true_px(reader, "pi_true_px");
  TTreeReaderValue<Double_t> pi_true_py(reader, "pi_true_py");
  TTreeReaderValue<Double_t> pi_true_pz(reader, "pi_true_pz");
  // Slow Pi
  TTreeReaderValue<Double_t> spi_true_pe(reader, "spi_true_pe");
  TTreeReaderValue<Double_t> spi_true_px(reader, "spi_true_px");
  TTreeReaderValue<Double_t> spi_true_py(reader, "spi_true_py");
  TTreeReaderValue<Double_t> spi_true_pz(reader, "spi_true_pz");
  // Anti-Nu_Tau
  TTreeReaderValue<Double_t> anu_tau_true_pe(reader, "anu_tau_true_pe");
  TTreeReaderValue<Double_t> anu_tau_true_px(reader, "anu_tau_true_px");
  TTreeReaderValue<Double_t> anu_tau_true_py(reader, "anu_tau_true_py");
  TTreeReaderValue<Double_t> anu_tau_true_pz(reader, "anu_tau_true_pz");
  // Nu_Tau
  TTreeReaderValue<Double_t> nu_tau_true_pe(reader, "nu_tau_true_pe");
  TTreeReaderValue<Double_t> nu_tau_true_px(reader, "nu_tau_true_px");
  TTreeReaderValue<Double_t> nu_tau_true_py(reader, "nu_tau_true_py");
  TTreeReaderValue<Double_t> nu_tau_true_pz(reader, "nu_tau_true_pz");
  // Anti-Nu_Mu
  TTreeReaderValue<Double_t> anu_mu_true_pe(reader, "anu_mu_true_pe");
  TTreeReaderValue<Double_t> anu_mu_true_px(reader, "anu_mu_true_px");
  TTreeReaderValue<Double_t> anu_mu_true_py(reader, "anu_mu_true_py");
  TTreeReaderValue<Double_t> anu_mu_true_pz(reader, "anu_mu_true_pz");

  // Define output branches ////////////////////////////////////////////////////
}

int main(int, char** argv) {
  TFile* input_file  = new TFile(argv[1], "read");
  TFile* output_file = new TFile(argv[2], "recreate");

  delete input_file;
  delete output_file;

  return 0;
}
