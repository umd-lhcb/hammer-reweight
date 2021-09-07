// Author: Yipeng Sun
// Last Change: Tue Sep 07, 2021 at 04:14 PM +0200

#include <any>
#include <iostream>
#include <map>
#include <string>
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

typedef Hammer::FourMomentum hp4;
typedef map<string, any>     PartEmu;

Double_t compute_p(Double_t m2_mom, Double_t m2_dau1, Double_t m2_dau2) {
  auto denom = 2 * TMath::Sqrt(m2_mom);
  auto nom   = TMath::Sqrt(
      m2_mom * m2_mom + m2_dau1 * m2_dau1 + m2_dau2 * m2_dau2 -
      2 * (m2_mom * m2_dau1 + m2_mom * m2_dau2 + m2_dau1 * m2_dau2));
  return nom / denom;
}

// Everything's in GeV!
PartEmu gen_B_decay(Int_t B_id, Double_t B_mass, Int_t D_id, Double_t D_mass,
                    Int_t l_id, Double_t l_mass, Int_t nu_id, Double_t q2,
                    TRandom& rng) {
  // Neutrino mass is set to 0
  PartEmu result{};

  // Remember that we are in the B rest frame
  // No need to boost back from B rest frame
  auto B_p = Hammer::FourMomentum(B_mass, 0, 0, 0);

  auto fac_neg = B_mass * B_mass - D_mass * D_mass;
  auto fac_pos = B_mass * B_mass + D_mass * D_mass;
  auto D_p_mag =
      TMath::Sqrt(1 / (4 * B_mass * B_mass) *
                  (TMath::Power(fac_neg, 2) + q2 * q2 - 2 * q2 * fac_pos));

  // Say D is flying in the z direction
  auto D_p = Hammer::FourMomentum(
      TMath::Sqrt(D_p_mag * D_p_mag + D_mass * D_mass), 0, 0, D_p_mag);

  auto l_sys_p = B_p - D_p;
  auto l_p_mag = compute_p(q2, l_mass * l_mass, 0);

  // Only θ_l is physical
  auto theta_l = rng.Uniform(0.1, 3.1);

  // Leptons are in the x-z plane
  // Angles are defined in the rest frame of the lepton pair, so rotate first
  // before boosting back to the B rest frame
  auto l_p_rest = Hammer::FourMomentum(
      TMath::Sqrt(l_mass * l_mass + l_p_mag * l_p_mag),
      l_p_mag * TMath::Sin(theta_l), 0, l_p_mag * TMath::Cos(theta_l));
  auto nu_p_rest = Hammer::FourMomentum(l_p_mag, -l_p_rest.px(), -l_p_rest.py(),
                                        -l_p_rest.pz());

  auto l_p  = l_p_rest.boostFromRestFrameOf(l_sys_p);
  auto nu_p = nu_p_rest.boostFromRestFrameOf(l_sys_p);

  result["B_id"] = B_id;
  result["B_p"]  = B_p;

  result["D_id"] = D_id;
  result["D_p"]  = D_p;

  result["theta_l"] = theta_l;  // the only physical angle
  result["l_id"]    = l_id;
  result["l_p"]     = l_p;
  result["nu_id"]   = nu_id;
  result["nu_p"]    = nu_p;

  result["q2"] = q2;

  return result;
}

PartEmu gen_B_decay(Int_t B_id, Double_t B_mass, Int_t D_id, Double_t D_mass,
                    Int_t l_id, Double_t l_mass, Int_t nu_id, Int_t D_dau_id,
                    Double_t D_dau_mass, Int_t pi_id, Double_t pi_mass,
                    Double_t q2, TRandom& rng) {
  auto result = gen_B_decay(B_id, B_mass, D_id, D_mass, l_id, l_mass, nu_id,
                            D_dau_id, D_dau_mass, pi_id, pi_mass, q2, rng);

  auto D_dau_p_rest =
      compute_p(D_mass * D_mass, D_dau_mass * D_dau_mass, pi_mass * pi_mass);

  // Now we have 2 additional physical angles: θ_v and χ
  auto theta_v = rng.Uniform(0.1, 3.1);
  auto chi     = -rng.Uniform(
      0.1, 6.2);  // From Bernlochner's definition, this is always negative

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
