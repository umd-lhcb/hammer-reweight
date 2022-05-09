// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon May 09, 2022 at 05:19 PM -0400

#include <any>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <math.h>

#include <TFile.h>
#include <TH2D.h>
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

#include "const.h"
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

const double PI = 3.141592;

const double B_MASS  = 5.27932;
const double B0_MASS = 5.27963;

const double DST_MASS = 2.01026;
const double D0_MASS  = 1.86483;

const double TAU_MASS = 1.77682;
const double PI_MASS  = 0.13957;

const bool LEPTON_POSITIVE = true;

// clang-format off

void setInputFF(Hammer::Hammer& ham) {
  ham.setFFInputScheme({
    {"BD", "ISGW2"},
    {"BD*", "ISGW2"},
  });
}

void setOutputFF(Hammer::Hammer& ham) {
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

void setDecays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BD*TauNu");
}

/////////////////////////////////
// Event generation: Interface //
/////////////////////////////////

typedef Hammer::FourMomentum HFM;
typedef map<string, any>     PartEmu;

class IRandGenerator {
 public:
  virtual vector<double> get() = 0;
  virtual PartEmu        gen() = 0;

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
 public:
  BToDUniformGenerator(double q2Min, double q2Max, double thetaLMin,
                       double thetaLMax, TRandom* rng);

  vector<double> get() override;
  PartEmu        gen() override;

  void setStepInQ(double step) { q2Step = step; };
  void reset() { q2 = q2Min; };

 protected:
  double   q2Min  = 0.;
  double   q2Step = 0.01;
  double   q2, q2Max, thetaLMin, thetaLMax;
  TRandom* rng;

  double  computeP(double m2_mom, double m2_dau1, double m2_dau2);
  PartEmu genBD(int B_id, double B_mass, int D_id, double D_mass, int l_id,
                double l_mass, int nu_id, double q2, double theta_l);
};

BToDUniformGenerator::BToDUniformGenerator(double q2Min, double q2Max,
                                           double thetaLMin, double thetaLMax,
                                           TRandom* rng)
    : q2Min(q2Min),
      q2Max(q2Max),
      thetaLMin(thetaLMin),
      thetaLMax(thetaLMax),
      rng(rng) {
  reset();
}

vector<double> BToDUniformGenerator::get() {
  vector<double> result{};
  if (q2 <= q2Max) {
    result.push_back(q2);
    result.push_back(rng->Uniform(thetaLMin, thetaLMax));
    q2 += q2Step;
  } else
    throw(
        domain_error("The q2 is out of its upper limit. Reset the generator to "
                     "start over!"));

  return result;
}

PartEmu BToDUniformGenerator::gen() {
  auto inputs = get();
  auto q2     = inputs[0];
  auto thetaL = inputs[1];
  return genBD(521, B_MASS, -421, D0_MASS, -15, TAU_MASS, 16, q2, thetaL);
}

// Protected methods ///////////////////////////////////////////////////////////
double BToDUniformGenerator::computeP(double m2Mom, double m2Dau1,
                                      double m2Dau2) {
  auto denom = 2 * Sqrt(m2Mom);
  auto nom   = Sqrt(m2Mom * m2Mom + m2Dau1 * m2Dau1 + m2Dau2 * m2Dau2 -
                  2 * (m2Mom * m2Dau1 + m2Mom * m2Dau2 + m2Dau1 * m2Dau2));
  return nom / denom;
}

// Everything's in GeV!
PartEmu BToDUniformGenerator::genBD(int bId, double mB, int dId, double mD,
                                    int lId, double mL, int nuId, double q2,
                                    double thetaL) {
  PartEmu result{};

  // Remember that we are in the B rest frame
  // No need to boost back from B rest frame
  auto pB = HFM{mB, 0, 0, 0};

  auto facNeg = mB * mB - mD * mD;
  auto facPos = mB * mB + mD * mD;
  auto pDMag =
      Sqrt(1 / (4 * mB * mB) * (Power(facNeg, 2) + q2 * q2 - 2 * q2 * facPos));

  // Say D is flying in the z direction
  auto pD = HFM(Sqrt(pDMag * pDMag + mD * mD), 0, 0, pDMag);

  auto pLSys = pB - pD;
  auto pLMag = computeP(q2, mL * mL, 0);

  // Leptons are in the x-z plane
  // Angles are defined in the rest frame of the lepton pair, so rotate first
  // before boosting back to the B rest frame
  auto pLRest  = HFM(Sqrt(mL * mL + pLMag * pLMag), pLMag * Sin(thetaL), 0,
                    pLMag * Cos(thetaL));
  auto pNuRest = HFM(pLMag, -pLRest.px(), -pLRest.py(), -pLRest.pz());

  auto pL  = pLRest.boostFromRestFrameOf(pLSys);
  auto pNu = pNuRest.boostFromRestFrameOf(pLSys);

  result["B_id"] = bId;
  result["B_p"]  = pB;

  result["D_id"] = dId;
  result["D_p"]  = pD;

  result["theta_l"] = thetaL;  // the only physical angle
  result["l_id"]    = lId;
  result["l_p"]     = pL;
  result["nu_id"]   = nuId;
  result["nu_p"]    = pNu;

  result["q2"] = q2;

  return result;
}

/////////////////////////////////////////////
// Event generation: D0, real distribution //
/////////////////////////////////////////////

class BToDRealGenerator : public BToDUniformGenerator {
 public:
  BToDRealGenerator(double q2_min, double q2_max, double theta_l_min,
                    double theta_l_max, TRandom* rng, string ff_mode = "ISGW2",
                    int xbins = 200, int ybins = 200);
  ~BToDRealGenerator();

  vector<double> get() override;
  PartEmu        gen() override;

  void setFF(string ff_mode);

 protected:
  double _theta_l_step;
  string _ff_mode;
  TH2D*  _histo;

  void buildHisto();

 private:
};

BToDRealGenerator::BToDRealGenerator(double q2_min, double q2_max,
                                     double theta_l_min, double theta_l_max,
                                     TRandom* rng, string ff_mode, int xbins,
                                     int ybins)
    : BToDUniformGenerator(q2_min, q2_max, theta_l_min, theta_l_max, rng),
      _ff_mode(ff_mode) {
  _histo = new TH2D("histo_BD", "histo_BD", xbins, q2_min, q2_max, ybins,
                    theta_l_min, theta_l_max);

  _q2_step      = (q2_max - q2_min) / xbins;
  _theta_l_step = (theta_l_max - theta_l_min) / ybins;

  buildHisto();
}

BToDRealGenerator::~BToDRealGenerator() { delete _histo; }

void BToDRealGenerator::setFF(string ff_mode) { _ff_mode = ff_mode; }

vector<double> BToDRealGenerator::get() {
  double q2, theta_l;
  _histo->GetRandom2(q2, theta_l, _rng);
  return vector<double>{q2, theta_l};
}

PartEmu BToDRealGenerator::gen() { return BToDUniformGenerator::gen(); }

// Protected methods ///////////////////////////////////////////////////////////

// NOTE: We hard-code to use B0 and Tau
void BToDRealGenerator::buildHisto() {
  auto   ff_model = BToDtaunu{};
  double fplus, fminus;

  if (_ff_mode == "ISGW2") {
    for (auto q2 = _q2_min + _q2_step / 2; q2 <= _q2_max - _q2_step / 2;
         q2 += _q2_step) {
      ff_model.ComputeISGW2(q2, fplus, fminus);
      for (auto theta_l = _theta_l_min + _theta_l_step / 2;
           theta_l <= _theta_l_max - _theta_l_step / 2;
           theta_l += _theta_l_step) {
        auto ff_val = ff_model.Gamma_q2tL(q2, theta_l, fplus, fminus, TAU_MASS);
        _histo->Fill(q2, theta_l, ff_val);

        // DEBUG
        // cout << "q2: " << q2 << " theta_l: " << theta_l << " ff val: " <<
        // ff_val
        //<< endl;
      }
    }
  } else if (_ff_mode == "CLN") {
    for (auto q2 = _q2_min; q2 <= _q2_max; q2 += _q2_step) {
      ff_model.ComputeCLN(q2, fplus, fminus);
      for (auto theta_l = _theta_l_min; theta_l <= _theta_l_max;
           theta_l += _theta_l_step) {
        auto ff_val = ff_model.Gamma_q2tL(q2, theta_l, fplus, fminus, TAU_MASS);
        _histo->Fill(q2, theta_l, ff_val);
      }
    }
  } else
    throw(domain_error("Unknown FF parameterization: " + _ff_mode));
}

//////////////////////////
// Event generation: D* //
//////////////////////////

class BToDstUniformGenerator : public BToDUniformGenerator {
  double _theta_v_min, _theta_v_max, _chi_min, _chi_max;

 public:
  BToDstUniformGenerator(double q2_min, double q2_max, double theta_l_min,
                         double theta_l_max, double theta_v_min,
                         double theta_v_max, double chiE_min, double chi_max,
                         TRandom* rng);

  vector<double> get() override;
  PartEmu        gen() override;

 protected:
  PartEmu genBDst(int B_id, double B_mass, int D_id, double D_mass, int l_id,
                  double l_mass, int nu_id, int D_dau_id, double D_dau_mass,
                  int pi_id, double pi_mass, double q2, double theta_l,
                  double theta_v, double chi);
};

BToDstUniformGenerator::BToDstUniformGenerator(
    double q2_min, double q2_max, double theta_l_min, double theta_l_max,
    double theta_v_min, double theta_v_max, double chi_min, double chi_max,
    TRandom* rng)
    : BToDUniformGenerator(q2_min, q2_max, theta_l_min, theta_l_max, rng),
      _theta_v_min(theta_v_min),
      _theta_v_max(theta_v_max),
      _chi_min(chi_min),
      _chi_max(chi_max) {}

vector<double> BToDstUniformGenerator::get() {
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
PartEmu BToDstUniformGenerator::genBDst(
    int B_id, double B_mass, int D_id, double D_mass, int l_id, double l_mass,
    int nu_id, int D_dau_id, double D_dau_mass, int pi_id, double pi_mass,
    double q2, double theta_l, double theta_v, double chi) {
  auto result =
      genBD(B_id, B_mass, D_id, D_mass, l_id, l_mass, nu_id, q2, theta_l);

  auto D_dau_p_mag =
      computeP(D_mass * D_mass, D_dau_mass * D_dau_mass, pi_mass * pi_mass);

  // These are defined in the D* rest frame
  auto D_dau_p_rest =
      HFM(Sqrt(D_dau_mass * D_dau_mass + D_dau_p_mag * D_dau_p_mag),
          D_dau_p_mag * Sin(theta_v) * Cos(chi),
          D_dau_p_mag * Sin(theta_v) * Sin(chi), D_dau_p_mag * Cos(theta_v));
  auto pi_p_rest = HFM(D_mass, 0, 0, 0) - D_dau_p_rest;

  // Boost back to B rest frame from D* rest frame
  auto D_p     = any_cast<HFM>(result["D_p"]);
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
                Hammer::Hammer& ham, int max_entries = 100000) {
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

  auto   B_key  = any_cast<int>(cands[0]["B_id"]);
  auto   D_key  = any_cast<int>(cands[0]["D_id"]);
  Bool_t is_Dst = (Abs(D_key) == 413);

  if (Abs(B_key) == 511) {
    calc_BDst.SetMasses(0);  // neutral B
  }

  Bool_t ham_ok;
  output_tree->Branch("ham_ok", &ham_ok);
  Bool_t ff_calc_ok;
  output_tree->Branch("ff_calc_ok", &ff_calc_ok);

  double q2_out;
  output_tree->Branch("q2_true", &q2_out);

  double ff_out;
  output_tree->Branch("wff", &ff_out);
  double ff_calc_out;
  output_tree->Branch("wff_calc", &ff_calc_out);

  int b_id_out;
  output_tree->Branch("b_id", &b_id_out);
  int d_id_out;
  output_tree->Branch("d_id", &d_id_out);
  int d_dau_id_out;
  output_tree->Branch("d_dau_id", &d_dau_id_out);
  int pi_id_out;
  output_tree->Branch("pi_id", &pi_id_out);

  double b_pe_out;
  output_tree->Branch("b_pe", &b_pe_out);
  double b_px_out;
  output_tree->Branch("b_px", &b_px_out);
  double b_py_out;
  output_tree->Branch("b_py", &b_py_out);
  double b_pz_out;
  output_tree->Branch("b_pz", &b_pz_out);

  double d_pe_out;
  output_tree->Branch("d_pe", &d_pe_out);
  double d_px_out;
  output_tree->Branch("d_px", &d_px_out);
  double d_py_out;
  output_tree->Branch("d_py", &d_py_out);
  double d_pz_out;
  output_tree->Branch("d_pz", &d_pz_out);

  double d_dau_pe_out;
  output_tree->Branch("d_dau_pe", &d_dau_pe_out);
  double d_dau_px_out;
  output_tree->Branch("d_dau_px", &d_dau_px_out);
  double d_dau_py_out;
  output_tree->Branch("d_dau_py", &d_dau_py_out);
  double d_dau_pz_out;
  output_tree->Branch("d_dau_pz", &d_dau_pz_out);

  double l_pe_out;
  output_tree->Branch("l_pe", &l_pe_out);
  double l_px_out;
  output_tree->Branch("l_px", &l_px_out);
  double l_py_out;
  output_tree->Branch("l_py", &l_py_out);
  double l_pz_out;
  output_tree->Branch("l_pz", &l_pz_out);

  double nu_pe_out;
  output_tree->Branch("nu_pe", &nu_pe_out);
  double nu_px_out;
  output_tree->Branch("nu_px", &nu_px_out);
  double nu_py_out;
  output_tree->Branch("nu_py", &nu_py_out);
  double nu_pz_out;
  output_tree->Branch("nu_pz", &nu_pz_out);

  double pi_pe_out;
  output_tree->Branch("pi_pe", &pi_pe_out);
  double pi_px_out;
  //      $Id: RateCalc.hh,v 1.2 2021/09/09 00:43:50 yipengsun Exp $
  output_tree->Branch("pi_px", &pi_px_out);
  double pi_py_out;
  output_tree->Branch("pi_py", &pi_py_out);
  double pi_pz_out;
  output_tree->Branch("pi_pz", &pi_pz_out);

  double theta_l_out;
  output_tree->Branch("theta_l", &theta_l_out);

  double theta_v_out;
  output_tree->Branch("theta_v", &theta_v_out);
  double chi_out;
  output_tree->Branch("chi", &chi_out);

  for (auto& cand : cands) {
    Hammer::Process proc;
    ham_ok     = false;
    ff_calc_ok = false;

    q2_out = any_cast<double>(cand["q2"]);

    b_id_out = any_cast<int>(cand["B_id"]);
    auto b_p = any_cast<HFM>(cand["B_p"]);
    b_pe_out = b_p.E();
    b_px_out = b_p.px();
    b_py_out = b_p.py();
    b_pz_out = b_p.pz();

    d_id_out = any_cast<int>(cand["D_id"]);
    auto d_p = any_cast<HFM>(cand["D_p"]);
    d_pe_out = d_p.E();
    d_px_out = d_p.px();
    d_py_out = d_p.py();
    d_pz_out = d_p.pz();

    auto l_p = any_cast<HFM>(cand["l_p"]);
    l_pe_out = l_p.E();
    l_px_out = l_p.px();
    l_py_out = l_p.py();
    l_pz_out = l_p.pz();

    auto nu_p = any_cast<HFM>(cand["nu_p"]);
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
    theta_l_out = any_cast<double>(cand["theta_l"]);

    if (is_Dst) {
      d_dau_id_out = any_cast<int>(cand["D_dau_id"]);
      auto d_dau_p = any_cast<HFM>(cand["D_dau_p"]);
      d_dau_pe_out = d_dau_p.E();
      d_dau_px_out = d_dau_p.px();
      d_dau_py_out = d_dau_p.py();
      d_dau_pz_out = d_dau_p.pz();

      pi_id_out = any_cast<int>(cand["pi_id"]);
      auto pi_p = any_cast<HFM>(cand["pi_p"]);
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
      theta_v_out = any_cast<double>(cand["theta_v"]);
      chi_out     = any_cast<double>(cand["chi"]);
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
        double calc_isgw2, calc_cln, a1, v, a2, a0, fplus, fminus;
        auto   cos_theta_l = Cos(theta_l_out);
        auto   cos_theta_v = Cos(theta_v_out);
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

  setDecays(ham);
  setInputFF(ham);
  setOutputFF(ham);

  ham.setUnits("GeV");
  ham.initRun();

  auto q2_min = TAU_MASS * TAU_MASS;

  auto gen_B_to_D =
      new BToDRealGenerator(q2_min, Power(B0_MASS - D0_MASS, 2), 0, PI, rng);
  auto gen_B_to_Dst = new BToDstUniformGenerator(
      q2_min, Power(B0_MASS - DST_MASS, 2), 0, PI, 0, PI, 0, 2 * PI, rng);

  weight_gen(gen_B_to_Dst, output_ntp, "tree_BDst", ham);
  weight_gen(gen_B_to_D, output_ntp, "tree_BD", ham);

  delete gen_B_to_D;
  delete rng;
  delete output_ntp;
}
