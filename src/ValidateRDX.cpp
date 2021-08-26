// Author: Yipeng Sun
// Last Change: Thu Aug 26, 2021 at 07:36 PM +0200

#include <iostream>
#include <map>
#include <vector>

#include <math.h>

#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>
#include <TTree.h>

#include <Hammer/Hammer.hh>
#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>
#include <Hammer/Process.hh>

#include <ff_dstaunu.hpp>

#include "utils.h"

using namespace std;

///////////////////
// Configurables //
///////////////////

const Double_t B_MASS  = 5.27932;
const Double_t B0_MASS = 5.27963;

const Double_t Dst_MASS = 2.01026;
const Double_t D0_MASS  = 1.86483;

const Double_t TAU_MASS = 1.77682;

// clang-format off
const auto LEGAL_B_MESON_IDS = vector<Int_t>{511, 521};

void set_input_ff(Hammer::Hammer& ham) {
  ham.setFFInputScheme({
    {"BD", "ISGW2"},
    {"BD*", "ISGW2"},
  });
}

void set_output_ff(Hammer::Hammer& ham) {
  ham.addFFScheme("OutputFF", {
    {"BD", "CLN_1"},
    {"BD*", "CLN_2"},
  });

  // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
  ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET3
  // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
  ham.setOptions("BtoD*CLN_2: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2
}
// clang-format on

void set_decays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BD*TauNu");
}

/////////////////////
// General helpers //
/////////////////////

typedef map<Int_t, Hammer::FourMomentum> PartEmu;

// Everything's in GeV!
auto gen_B_decay(Int_t B_id, Double_t B_mass, Int_t D_id, Double_t D_mass,
                 Int_t l_id, Double_t l_mass, Int_t nu_id, Double_t q2,
                 TRandom& rng) {
  // Neutrino mass is set to 0
  PartEmu result{};

  auto B_pz = rng.Uniform(20, 200);
  auto B_p = Hammer::FourMomentum(TMath::Sqrt(B_pz * B_pz + B_mass * B_mass), 0,
                                  0, B_pz);

  auto fact_neg = B_mass * B_mass - D_mass * D_mass;
  auto fact_pos = B_mass * B_mass + D_mass * D_mass;
  auto D_p_mag =
      TMath::Sqrt(1 / (4 * B_mass * B_mass) *
                  (TMath::Power(fact_neg, 2) + q2 * q2 - 2 * q2 * fact_pos));
  auto D_theta = rng.Uniform(0.3, 6.2);
  auto D_phi   = rng.Uniform(0.3, 3.1);

  auto D_p_rest =
      Hammer::FourMomentum(TMath::Sqrt(D_p_mag * D_p_mag + D_mass * D_mass),
                           D_p_mag * TMath::Sin(D_phi) * TMath::Cos(D_theta),
                           D_p_mag * TMath::Sin(D_phi) * TMath::Sin(D_theta),
                           D_p_mag * TMath::Cos(D_phi));
  auto D_p = D_p_rest.boostFromRestFrameOf(B_p);

  auto l_sys_p = B_p - D_p;
  auto l_p_mag = (q2 - l_mass * l_mass) / (2 * TMath::Sqrt(q2));
  auto l_theta = rng.Uniform(0.3, 6.2);
  auto l_phi   = rng.Uniform(0.3, 3.1);

  auto l_p_rest =
      Hammer::FourMomentum(TMath::Sqrt(l_mass * l_mass + l_p_mag * l_p_mag),
                           l_p_mag * TMath::Sin(l_phi) * TMath::Cos(l_theta),
                           l_p_mag * TMath::Sin(l_phi) * TMath::Sin(l_theta),
                           l_p_mag * TMath::Cos(l_phi));
  auto nu_p_rest = Hammer::FourMomentum(l_p_mag, -l_p_rest.px(), -l_p_rest.py(),
                                        -l_p_rest.pz());

  auto l_p  = l_p_rest.boostFromRestFrameOf(l_sys_p);
  auto nu_p = nu_p_rest.boostFromRestFrameOf(l_sys_p);

  result[B_id]  = B_p;
  result[D_id]  = D_p;
  result[l_id]  = l_p;
  result[nu_id] = nu_p;
  result[-1]    = Hammer::FourMomentum(q2, 0, 0, 0);  // store q2 as a vector

  return result;
}

auto gen_BDstTau_decay(Double_t q2, TRandom& rng) {
  return gen_B_decay(511, B0_MASS, -413, Dst_MASS, -15, TAU_MASS, 16, q2, rng);
}

auto gen_BDTau_decay(Double_t q2, TRandom& rng) {
  return gen_B_decay(521, B_MASS, -421, D0_MASS, -15, TAU_MASS, 16, q2, rng);
}

/////////////////
// Reweighting //
/////////////////

void weight_gen(vector<PartEmu> cands, Int_t B_key, Int_t D_key,
                TFile* output_ntp, TString tree_name, Hammer::Hammer& ham) {
  auto output_tree = new TTree(tree_name, tree_name);
  auto calc_BDst   = BToDstaunu{};

  if (TMath::Abs(B_key) == 511) calc_BDst.SetMasses(0);  // neutral B

  Bool_t ham_ok;
  output_tree->Branch("ham_ok", &ham_ok);

  Double_t q2_out;
  output_tree->Branch("q2_true", &q2_out);

  Double_t ff_out;
  output_tree->Branch("w_ff", &ff_out);
  Double_t ff_calc_out;
  output_tree->Branch("w_ff_calc", &ff_calc_out);

  Int_t b_id_out;
  output_tree->Branch("b_id", &b_id_out);
  Int_t d_id_out;
  output_tree->Branch("d_id", &d_id_out);

  Double_t b_pe_out;
  output_tree->Branch("b_pe", &b_pe_out);
  Double_t b_px_out;
  output_tree->Branch("b_px", &b_px_out);
  Double_t b_py_out;
  output_tree->Branch("b_py", &b_py_out);
  Double_t b_pz_out;
  output_tree->Branch("b_pz", &b_pz_out);

  Double_t d_pe_out;
  output_tree->Branch("d_pe", &d_pe_out);
  Double_t d_px_out;
  output_tree->Branch("d_px", &d_px_out);
  Double_t d_py_out;
  output_tree->Branch("d_py", &d_py_out);
  Double_t d_pz_out;
  output_tree->Branch("d_pz", &d_pz_out);

  for (auto& cand : cands) {
    Hammer::Process proc;
    ham_ok = false;

    auto q2 = cand[-1].E();
    q2_out  = q2;

    b_id_out = B_key;
    d_id_out = D_key;

    b_pe_out = cand[B_key].E();
    b_px_out = cand[B_key].px();
    b_py_out = cand[B_key].py();
    b_pz_out = cand[B_key].pz();

    d_pe_out = cand[D_key].E();
    d_px_out = cand[D_key].px();
    d_py_out = cand[D_key].py();
    d_pz_out = cand[D_key].pz();

    auto part_B     = particle(cand[B_key], B_key);
    auto part_B_idx = proc.addParticle(part_B);

    auto part_D     = particle(cand[D_key], D_key);
    auto part_D_idx = proc.addParticle(part_D);

    auto part_L     = particle(cand[-15], -15);
    auto part_L_idx = proc.addParticle(part_L);

    auto part_NuL     = particle(cand[16], 16);
    auto part_NuL_idx = proc.addParticle(part_NuL);

    proc.addVertex(part_B_idx, {part_D_idx, part_L_idx, part_NuL_idx});

    ham.initEvent();
    auto proc_id = ham.addProcess(proc);

    ff_out      = 1.0;
    ff_calc_out = 1.0;

    if (proc_id != 0) {
      ham.processEvent();
      ff_out = ham.getWeight("OutputFF");

      if (!isnan(ff_out) && !isinf(ff_out)) {
        ham_ok = true;

        Double_t calc_isgw2 = 1.;
        Double_t calc_cln   = 1.;

        if (TMath::Abs(D_key) == 413) {
          calc_isgw2 = calc_BDst.Compute(q2_out, false, TAU_MASS);
          calc_cln   = calc_BDst.Compute(q2_out, true, TAU_MASS);
        }

        ff_calc_out = calc_cln / calc_isgw2;
      } else
        ff_out = 1.0;
    }

    output_tree->Fill();
  }

  output_ntp->Write();
  delete output_tree;
}

//////////
// Main //
//////////

int main(int, char** argv) {
  auto output_ntp = new TFile(argv[1], "recreate");

  Hammer::Hammer ham{};
  auto           rng = TRandom(42);

  set_decays(ham);
  set_input_ff(ham);
  set_output_ff(ham);

  ham.setUnits("GeV");
  ham.initRun();

  auto q2s = vector<Double_t>{};
  for (auto i = 3.2; i <= 12.2; i += 0.001) {
    q2s.push_back(i);
  }

  auto cands_BDst = vector<PartEmu>{};
  auto cands_BD   = vector<PartEmu>{};
  for (auto& q2 : q2s) {
    cands_BDst.push_back(gen_BDstTau_decay(q2, rng));
    cands_BD.push_back(gen_BDTau_decay(q2, rng));
  }

  weight_gen(cands_BDst, 511, -413, output_ntp, "tree_BDst", ham);

  delete output_ntp;
}
