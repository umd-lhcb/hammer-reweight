// Author: Yipeng Sun

#include <any>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
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
#include <ff_dtaunu.hpp>

#include "utils_general.h"
#include "utils_ham.h"

using namespace std;
using TMath::Abs;
using TMath::Cos;
using TMath::Power;
using TMath::Sin;
using TMath::Sqrt;

///////////////////
// Configurables //
///////////////////
// Masses are defined in GeV!

const Double_t B_MASS  = 5.27932;
const Double_t B0_MASS = 5.27963;

const Double_t DST_MASS = 2.01026;
const Double_t D0_MASS  = 1.86483;

const Double_t TAU_MASS = 1.77682;
const Double_t PI_MASS  = 0.13957;

const bool LEPTON_POSITIVE = true;

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
  ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET2
  // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
  ham.setOptions("BtoD*CLN_2: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2
}
// clang-format on

void set_decays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BD*TauNu");
}

/////////////////////////////////
// Event generation: Interface //
/////////////////////////////////

typedef Hammer::FourMomentum hp4;
typedef map<string, any>     PartEmu;

class IRandGenerator {
 public:
  virtual vector<Double_t> get() = 0;
  virtual PartEmu          gen() = 0;

  virtual ~IRandGenerator() = 0;  // Just in case to avoid potential memory leak
};

IRandGenerator::~IRandGenerator() {
  cout << "Generator itself doesn't provide a destructor or it was deleted via "
          "a base pointer."
       << endl;
}

//////////////////////////
// Event generation: D0 //
//////////////////////////

class BToDUniformGenerator : public IRandGenerator {
  Double_t _q2_min  = 0.;
  Double_t _q2_step = 0.01;
  Double_t _q2, _q2_max, _theta_l_min, _theta_l_max;

 public:
  BToDUniformGenerator(Double_t q2_min, Double_t q2_max, Double_t theta_l_min,
                       Double_t theta_l_max, TRandom* rng);

  vector<Double_t> get() override;
  PartEmu          gen() override;

  void setStepInQ(Double_t step) { _q2_step = step; };
  void reset() { _q2 = _q2_min; };

 protected:
  TRandom* _rng;

  Double_t computeP(Double_t m2_mom, Double_t m2_dau1, Double_t m2_dau2);
  PartEmu  genBD(Int_t B_id, Double_t B_mass, Int_t D_id, Double_t D_mass,
                 Int_t l_id, Double_t l_mass, Int_t nu_id, Double_t q2,
                 Double_t theta_l);
};

BToDUniformGenerator::BToDUniformGenerator(Double_t q2_min, Double_t q2_max,
                                           Double_t theta_l_min,
                                           Double_t theta_l_max, TRandom* rng)
    : _q2_min(q2_min),
      _q2_max(q2_max),
      _theta_l_min(theta_l_min),
      _theta_l_max(theta_l_max),
      _rng(rng) {
  reset();
}

vector<Double_t> BToDUniformGenerator::get() {
  vector<Double_t> result{};
  if (_q2 <= _q2_max) {
    result.push_back(_q2);
    result.push_back(_rng->Uniform(_theta_l_min, _theta_l_max));
    _q2 += _q2_step;
  } else
    throw(
        domain_error("The q2 is out of its upper limit. Reset the generator to "
                     "start over!"));

  return result;
}

PartEmu BToDUniformGenerator::gen() {
  auto rand_num = get();
  auto q2       = rand_num[0];
  auto theta_l  = rand_num[1];

  return genBD(521, B_MASS, -421, D0_MASS, -15, TAU_MASS, 16, q2, theta_l);
}

// Protected methods ///////////////////////////////////////////////////////////
Double_t BToDUniformGenerator::computeP(Double_t m2_mom, Double_t m2_dau1,
                                        Double_t m2_dau2) {
  auto denom = 2 * Sqrt(m2_mom);
  auto nom =
      Sqrt(m2_mom * m2_mom + m2_dau1 * m2_dau1 + m2_dau2 * m2_dau2 -
           2 * (m2_mom * m2_dau1 + m2_mom * m2_dau2 + m2_dau1 * m2_dau2));
  return nom / denom;
}

// Everything's in GeV!
PartEmu BToDUniformGenerator::genBD(Int_t B_id, Double_t B_mass, Int_t D_id,
                                    Double_t D_mass, Int_t l_id,
                                    Double_t l_mass, Int_t nu_id, Double_t q2,
                                    Double_t theta_l) {
  PartEmu result{};

  // Remember that we are in the B rest frame
  // No need to boost back from B rest frame
  auto B_p = hp4{B_mass, 0, 0, 0};

  auto fac_neg = B_mass * B_mass - D_mass * D_mass;
  auto fac_pos = B_mass * B_mass + D_mass * D_mass;
  auto D_p_mag = Sqrt(1 / (4 * B_mass * B_mass) *
                      (Power(fac_neg, 2) + q2 * q2 - 2 * q2 * fac_pos));

  // Say D is flying in the z direction
  auto D_p = hp4(Sqrt(D_p_mag * D_p_mag + D_mass * D_mass), 0, 0, D_p_mag);

  auto l_sys_p = B_p - D_p;
  auto l_p_mag = computeP(q2, l_mass * l_mass, 0);

  // Leptons are in the x-z plane
  // Angles are defined in the rest frame of the lepton pair, so rotate first
  // before boosting back to the B rest frame
  auto l_p_rest  = hp4(Sqrt(l_mass * l_mass + l_p_mag * l_p_mag),
                      l_p_mag * Sin(theta_l), 0, l_p_mag * Cos(theta_l));
  auto nu_p_rest = hp4(l_p_mag, -l_p_rest.px(), -l_p_rest.py(), -l_p_rest.pz());

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

//////////////////////////
// Event generation: D* //
//////////////////////////

class BToDstUniformGenerator : BToDUniformGenerator {
  Double_t _theta_v_min, _theta_v_max, _chi_min, _chi_max;

 public:
  BToDstUniformGenerator(Double_t q2_min, Double_t q2_max, Double_t theta_l_min,
                         Double_t theta_l_max, Double_t theta_v_min,
                         Double_t theta_v_max, Double_t chi_min,
                         Double_t chi_max, TRandom* rng);

  vector<Double_t> get() override;
  PartEmu          gen() override;

 protected:
  PartEmu genBDst(Int_t B_id, Double_t B_mass, Int_t D_id, Double_t D_mass,
                  Int_t l_id, Double_t l_mass, Int_t nu_id, Int_t D_dau_id,
                  Double_t D_dau_mass, Int_t pi_id, Double_t pi_mass,
                  Double_t q2, Double_t theta_l, Double_t theta_v,
                  Double_t chi);
};

BToDstUniformGenerator::BToDstUniformGenerator(
    Double_t q2_min, Double_t q2_max, Double_t theta_l_min,
    Double_t theta_l_max, Double_t theta_v_min, Double_t theta_v_max,
    Double_t chi_min, Double_t chi_max, TRandom* rng)
    : BToDUniformGenerator(q2_min, q2_max, theta_l_min, theta_l_max, rng),
      _theta_v_min(theta_v_min),
      _theta_v_max(theta_v_max),
      _chi_min(chi_min),
      _chi_max(chi_max) {}

vector<Double_t> BToDstUniformGenerator::get() {
  auto result = BToDUniformGenerator::get();
  result.push_back(_rng->Uniform(_theta_v_min, _theta_v_max));
  result.push_back(_rng->Uniform(_chi_min, _chi_max));

  return result;
}

PartEmu BToDstUniformGenerator::gen() {
  auto rand_num = get();
  auto q2       = rand_num[0];
  auto theta_l  = rand_num[1];
  auto theta_v  = rand_num[2];
  auto chi      = rand_num[3];

  return genBDst(511, B0_MASS, -413, DST_MASS, -15, TAU_MASS, 16, -421, D0_MASS,
                 -211, PI_MASS, q2, theta_l, theta_v, chi);
}

// Protected methods ///////////////////////////////////////////////////////////
PartEmu BToDstUniformGenerator::genBDst(Int_t B_id, Double_t B_mass, Int_t D_id,
                                        Double_t D_mass, Int_t l_id,
                                        Double_t l_mass, Int_t nu_id,
                                        Int_t D_dau_id, Double_t D_dau_mass,
                                        Int_t pi_id, Double_t pi_mass,
                                        Double_t q2, Double_t theta_l,
                                        Double_t theta_v, Double_t chi) {
  auto result =
      genBD(B_id, B_mass, D_id, D_mass, l_id, l_mass, nu_id, q2, theta_l);

  auto D_dau_p_mag =
      computeP(D_mass * D_mass, D_dau_mass * D_dau_mass, pi_mass * pi_mass);

  // These are defined in the D* rest frame
  auto D_dau_p_rest =
      hp4(Sqrt(D_dau_mass * D_dau_mass + D_dau_p_mag * D_dau_p_mag),
          D_dau_p_mag * Sin(theta_v) * Cos(chi),
          D_dau_p_mag * Sin(theta_v) * Sin(chi), D_dau_p_mag * Cos(theta_v));
  auto pi_p_rest = hp4(D_mass, 0, 0, 0) - D_dau_p_rest;

  // Boost back to B rest frame from D* rest frame
  auto D_p     = any_cast<hp4>(result["D_p"]);
  auto D_dau_p = D_dau_p_rest.boostFromRestFrameOf(D_p);
  auto pi_p    = pi_p_rest.boostFromRestFrameOf(D_p);

  result["theta_v"] = theta_v;
  result["chi"]     = chi;

  result["D_dau_id"] = D_dau_id;
  result["D_dau_p"]  = D_dau_p;
  result["pi_id"]    = pi_id;
  result["pi_p"]     = pi_p;

  return result;
}

/////////////////
// Reweighting //
/////////////////

void weight_gen(IRandGenerator* rng, TFile* output_ntp, TString tree_name,
                Hammer::Hammer& ham, Int_t max_entries = 10000) {
  auto output_tree = new TTree(tree_name, tree_name);
  auto calc_BDst   = BToDstaunu{};
  auto calc_BD     = BToDtaunu{};

  vector<PartEmu> cands{};

  for (auto i = 0; i < max_entries; i++) {
    try {
      cands.push_back(rng->gen());
    } catch (const domain_error& e) {
      break;
    }
  }

  auto   B_key  = any_cast<Int_t>(cands[0]["B_id"]);
  auto   D_key  = any_cast<Int_t>(cands[0]["D_id"]);
  Bool_t is_Dst = (Abs(D_key) == 413);

  if (Abs(B_key) == 511) calc_BDst.SetMasses(0);  // neutral B

  Bool_t ham_ok;
  output_tree->Branch("ham_ok", &ham_ok);
  Bool_t ff_calc_ok;
  output_tree->Branch("ff_calc_ok", &ff_calc_ok);

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
  Int_t d_dau_id_out;
  output_tree->Branch("d_dau_id", &d_dau_id_out);
  Int_t pi_id_out;
  output_tree->Branch("pi_id", &pi_id_out);

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

  Double_t d_dau_pe_out;
  output_tree->Branch("d_dau_pe", &d_dau_pe_out);
  Double_t d_dau_px_out;
  output_tree->Branch("d_dau_px", &d_dau_px_out);
  Double_t d_dau_py_out;
  output_tree->Branch("d_dau_py", &d_dau_py_out);
  Double_t d_dau_pz_out;
  output_tree->Branch("d_dau_pz", &d_dau_pz_out);

  Double_t l_pe_out;
  output_tree->Branch("l_pe", &l_pe_out);
  Double_t l_px_out;
  output_tree->Branch("l_px", &l_px_out);
  Double_t l_py_out;
  output_tree->Branch("l_py", &l_py_out);
  Double_t l_pz_out;
  output_tree->Branch("l_pz", &l_pz_out);

  Double_t nu_pe_out;
  output_tree->Branch("nu_pe", &nu_pe_out);
  Double_t nu_px_out;
  output_tree->Branch("nu_px", &nu_px_out);
  Double_t nu_py_out;
  output_tree->Branch("nu_py", &nu_py_out);
  Double_t nu_pz_out;
  output_tree->Branch("nu_pz", &nu_pz_out);

  Double_t pi_pe_out;
  output_tree->Branch("pi_pe", &pi_pe_out);
  Double_t pi_px_out;
  //      $Id: RateCalc.hh,v 1.2 2021/09/09 00:43:50 yipengsun Exp $
  output_tree->Branch("pi_px", &pi_px_out);
  Double_t pi_py_out;
  output_tree->Branch("pi_py", &pi_py_out);
  Double_t pi_pz_out;
  output_tree->Branch("pi_pz", &pi_pz_out);

  Double_t theta_l_out;
  output_tree->Branch("theta_l", &theta_l_out);

  Double_t theta_v_out;
  output_tree->Branch("theta_v", &theta_v_out);
  Double_t chi_out;
  output_tree->Branch("chi", &chi_out);

  for (auto& cand : cands) {
    Hammer::Process proc;
    ham_ok     = false;
    ff_calc_ok = false;

    q2_out = any_cast<Double_t>(cand["q2"]);

    b_id_out = any_cast<Int_t>(cand["B_id"]);
    auto b_p = any_cast<hp4>(cand["B_p"]);
    b_pe_out = b_p.E();
    b_px_out = b_p.px();
    b_py_out = b_p.py();
    b_pz_out = b_p.pz();

    d_id_out = any_cast<Int_t>(cand["D_id"]);
    auto d_p = any_cast<hp4>(cand["D_p"]);
    d_pe_out = d_p.E();
    d_px_out = d_p.px();
    d_py_out = d_p.py();
    d_pz_out = d_p.pz();

    auto l_p = any_cast<hp4>(cand["l_p"]);
    l_pe_out = l_p.E();
    l_px_out = l_p.px();
    l_py_out = l_p.py();
    l_pz_out = l_p.pz();

    auto nu_p = any_cast<hp4>(cand["nu_p"]);
    nu_pe_out = nu_p.E();
    nu_px_out = nu_p.px();
    nu_py_out = nu_p.py();
    nu_pz_out = nu_p.pz();

    auto part_B     = particle(b_p, b_id_out);
    auto part_B_idx = proc.addParticle(part_B);

    auto part_D     = particle(d_p, d_id_out);
    auto part_D_idx = proc.addParticle(part_D);

    auto part_L     = particle(l_p, -15);
    auto part_L_idx = proc.addParticle(part_L);

    auto part_NuL     = particle(nu_p, 16);
    auto part_NuL_idx = proc.addParticle(part_NuL);

    proc.addVertex(part_B_idx, {part_D_idx, part_L_idx, part_NuL_idx});

    // Angles
    theta_l_out = any_cast<Double_t>(cand["theta_l"]);

    if (is_Dst) {
      d_dau_id_out = any_cast<Int_t>(cand["D_dau_id"]);
      auto d_dau_p = any_cast<hp4>(cand["D_dau_p"]);
      d_dau_pe_out = d_dau_p.E();
      d_dau_px_out = d_dau_p.px();
      d_dau_py_out = d_dau_p.py();
      d_dau_pz_out = d_dau_p.pz();

      pi_id_out = any_cast<Int_t>(cand["pi_id"]);
      auto pi_p = any_cast<hp4>(cand["pi_p"]);
      pi_pe_out = pi_p.E();
      pi_px_out = pi_p.px();
      pi_py_out = pi_p.py();
      pi_pz_out = pi_p.pz();

      auto part_D_dau     = particle(d_dau_p, d_dau_id_out);
      auto part_D_dau_idx = proc.addParticle(part_D_dau);

      auto part_Pi     = particle(pi_p, pi_id_out);
      auto part_Pi_idx = proc.addParticle(part_Pi);

      proc.addVertex(part_D_idx, {part_D_dau_idx, part_Pi_idx});

      // Additional angles
      theta_v_out = any_cast<Double_t>(cand["theta_v"]);
      chi_out     = any_cast<Double_t>(cand["chi"]);
    } else {
      d_dau_id_out = pi_id_out = 0;
      d_dau_pe_out = d_dau_px_out = d_dau_py_out = d_dau_pz_out = 0.;
      pi_pe_out = pi_px_out = pi_py_out = pi_pz_out = 0.;

      theta_v_out = chi_out = 0.;
    }

    ham.initEvent();
    auto proc_id = ham.addProcess(proc);

    if (proc_id != 0) {
      ham.processEvent();
      ff_out = ham.getWeight("OutputFF");

      if (!isnan(ff_out) && !isinf(ff_out)) {
        ham_ok = true;
        // Compute FF weights w/ Manuel's calculator
        Double_t calc_isgw2, calc_cln, a1, v, a2, a0, fplus, fminus;
        auto     cos_theta_l = Cos(theta_l_out);
        auto     cos_theta_v = Cos(theta_v_out);
        if (is_Dst) {
          calc_BDst.ComputeISGW2(q2_out, a1, v, a2, a0);
          calc_isgw2 = calc_BDst.Gamma_q2Angular(
              q2_out, cos_theta_l, cos_theta_v, chi_out, false, LEPTON_POSITIVE,
              a1, v, a2, a0, TAU_MASS);

          calc_BDst.ComputeCLN(q2_out, a1, v, a2, a0);
          calc_cln = calc_BDst.Gamma_q2Angular(q2_out, cos_theta_l, cos_theta_v,
                                               chi_out, false, LEPTON_POSITIVE,
                                               a1, v, a2, a0, TAU_MASS);

        } else {
          calc_BD.ComputeISGW2(q2_out, fplus, fminus);
          calc_isgw2 =
              calc_BD.Gamma_q2tL(q2_out, theta_l_out, fplus, fminus, TAU_MASS);

          calc_BD.ComputeCLN(q2_out, fplus, fminus);
          calc_cln =
              calc_BD.Gamma_q2tL(q2_out, theta_l_out, fplus, fminus, TAU_MASS);
        }
        ff_calc_out = calc_cln / calc_isgw2;

        if (!isnan(ff_calc_out) && !isinf(ff_calc_out)) ff_calc_ok = true;
        // DEBUG
        // cout << "CLN: " << calc_cln << "; ISGW2: " << calc_isgw2 << endl;

      } else {
        ff_out      = 1.0;
        ff_calc_out = 1.0;
      }
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
  auto rng        = new TRandom(42);

  Hammer::Hammer ham{};

  set_decays(ham);
  set_input_ff(ham);
  set_output_ff(ham);

  ham.setUnits("GeV");
  ham.initRun();

  auto q2_min = TAU_MASS * TAU_MASS;

  auto gen_B_to_D = new BToDUniformGenerator(
      q2_min, Power(B0_MASS - DST_MASS, 2), 0, 3.14159, rng);

  // weight_gen(cands_BDst, output_ntp, "tree_BDst", ham);
  weight_gen(gen_B_to_D, output_ntp, "tree_BD", ham);

  delete gen_B_to_D;
  delete rng;
  delete output_ntp;
}
