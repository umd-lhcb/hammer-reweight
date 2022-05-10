// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon May 09, 2022 at 08:31 PM -0400

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

  ham.addFFScheme("OutputFFBGL", {
    {"BD", "BGL_1"},
    {"BD*", "BGL_2"},
  });

  ham.addFFScheme("OutputFFBGLVarRef", {
    {"BD", "BGLVar_Ref1"},
    {"BD*", "BGLVar_Ref2"},
  });

  ham.addFFScheme("OutputFFBGLVar", {
    {"BD", "BGLVar_1"},
    {"BD*", "BGLVar_2"},
  });

  // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
  ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET2
  // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
  ham.setOptions("BtoD*CLN_2: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2

  ham.setFFEigenvectors("BtoD", "BGLVar_1", {{"delta_a1", 1.0}});  // HQET2
  ham.setFFEigenvectors("BtoD*", "BGLVar_1", {{"delta_a1", 1.0}});  // HQET2
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
  BToDRealGenerator(double q2Min, double q2Max, double thetaLMin,
                    double thetaLMax, TRandom* rng, string ffMode = "ISGW2",
                    int xBins = 200, int yBins = 200);
  ~BToDRealGenerator();

  vector<double> get() override;
  PartEmu        gen() override;

  void setFF(string ffMode);

 protected:
  double thetaLStep;
  string ffMode;
  TH2D*  histo;

  void buildHisto();

 private:
};

BToDRealGenerator::BToDRealGenerator(double q2Min, double q2Max,
                                     double thetaLMin, double thetaLMax,
                                     TRandom* rng, string ffMode, int xBins,
                                     int yBins)
    : BToDUniformGenerator(q2Min, q2Max, thetaLMin, thetaLMax, rng),
      ffMode(ffMode) {
  histo = new TH2D("histo_BD", "histo_BD", xBins, q2Min, q2Max, yBins,
                   thetaLMin, thetaLMax);

  q2Step     = (q2Max - q2Min) / xBins;
  thetaLStep = (thetaLMax - thetaLMin) / yBins;
  buildHisto();
}

BToDRealGenerator::~BToDRealGenerator() { delete histo; }

void BToDRealGenerator::setFF(string ffMode) { this->ffMode = ffMode; }

vector<double> BToDRealGenerator::get() {
  double q2, thetaL;
  histo->GetRandom2(q2, thetaL, rng);
  return vector<double>{q2, thetaL};
}

PartEmu BToDRealGenerator::gen() { return BToDUniformGenerator::gen(); }

// Protected methods ///////////////////////////////////////////////////////////

// NOTE: We hard-code to use B0 and Tau
void BToDRealGenerator::buildHisto() {
  auto   ffModel = BToDtaunu{};
  double fPlus, fMinus;

  if (ffMode == "ISGW2") {
    // using the center of each bin
    for (auto q2 = q2Min + q2Step / 2; q2 <= q2Max - q2Step / 2; q2 += q2Step) {
      ffModel.ComputeISGW2(q2, fPlus, fMinus);
      for (auto thetaL = thetaLMin + thetaLStep / 2;
           thetaL <= thetaLMax - thetaLStep / 2; thetaL += thetaLStep) {
        auto ffVal = ffModel.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);
        histo->Fill(q2, thetaL, ffVal);

        // DEBUG
        // cout << "q2: " << q2 << " theta_l: " << thetaL << " ff val: " <<
        // ffVal
        //<< endl;
      }
    }
  } else if (ffMode == "CLN") {
    for (auto q2 = q2Min; q2 <= q2Max; q2 += q2Step) {
      ffModel.ComputeCLN(q2, fPlus, fMinus);
      for (auto thetaL = thetaLMin; thetaL <= thetaLMax; thetaL += thetaLStep) {
        auto ffVal = ffModel.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);
        histo->Fill(q2, thetaL, ffVal);
      }
    }
  } else
    throw(domain_error("Unknown FF parameterization: " + ffMode));
}

//////////////////////////
// Event generation: D* //
//////////////////////////

class BToDstUniformGenerator : public BToDUniformGenerator {
  double thetaVMin, thetaVMax, chiMin, chiMax;

 public:
  BToDstUniformGenerator(double q2Min, double q2Max, double thetaLMin,
                         double thetaLMax, double thetaVMin, double thetaVMax,
                         double chiMin, double chiMax, TRandom* rng);

  vector<double> get() override;
  PartEmu        gen() override;

 protected:
  PartEmu genBDst(int bId, double mB, int dId, double mD, int lId, double mL,
                  int nuId, int dDauId, double mDDau, int piId, double mPi,
                  double q2, double thetaL, double thetaV, double chi);
};

BToDstUniformGenerator::BToDstUniformGenerator(double q2Min, double q2Max,
                                               double thetaLMin,
                                               double thetaLMax,
                                               double thetaVMin,
                                               double thetaVMax, double chiMin,
                                               double chiMax, TRandom* rng)
    : BToDUniformGenerator(q2Min, q2Max, thetaLMin, thetaLMax, rng),
      thetaVMin(thetaVMin),
      thetaVMax(thetaVMax),
      chiMin(chiMin),
      chiMax(chiMax) {}

vector<double> BToDstUniformGenerator::get() {
  auto result = BToDUniformGenerator::get();
  result.push_back(rng->Uniform(thetaVMin, thetaVMax));
  result.push_back(rng->Uniform(chiMin, chiMax));
  return result;
}

PartEmu BToDstUniformGenerator::gen() {
  auto inputs = get();
  auto q2     = inputs[0];
  auto thetaL = inputs[1];
  auto thetaV = inputs[2];
  auto chi    = inputs[3];

  return genBDst(511, B0_MASS, -413, DST_MASS, -15, TAU_MASS, 16, -421, D0_MASS,
                 -211, PI_MASS, q2, thetaL, thetaV, chi);
}

// Protected methods ///////////////////////////////////////////////////////////
PartEmu BToDstUniformGenerator::genBDst(int bId, double mB, int dId, double mD,
                                        int lId, double mL, int nuId,
                                        int dDauId, double mDDau, int piId,
                                        double mPi, double q2, double thetaL,
                                        double thetaV, double chi) {
  auto result = genBD(bId, mB, dId, mD, lId, mL, nuId, q2, thetaL);

  auto pDDauMag = computeP(mD * mD, mDDau * mDDau, mPi * mPi);

  // These are defined in the D* rest frame
  auto pDDauRest =
      HFM(Sqrt(mDDau * mDDau + pDDauMag * pDDauMag),
          pDDauMag * Sin(thetaV) * Cos(chi), pDDauMag * Sin(thetaV) * Sin(chi),
          pDDauMag * Cos(thetaV));
  auto pPiRest = HFM(mD, 0, 0, 0) - pDDauRest;

  // Boost back to B rest frame from D* rest frame
  auto pD    = any_cast<HFM>(result["D_p"]);
  auto pDDau = pDDauRest.boostFromRestFrameOf(pD);
  auto pPi   = pPiRest.boostFromRestFrameOf(pD);

  result["theta_v"] = thetaV;
  result["chi"]     = chi;

  result["D_dau_id"] = dDauId;
  result["D_dau_p"]  = pDDau;
  result["pi_id"]    = piId;
  result["pi_p"]     = pPi;

  return result;
}

/////////////////
// Reweighting //
/////////////////

void weightGen(IRandGenerator* rng, TFile* outputNtp, TString treeName,
               Hammer::Hammer& ham, int maxEntries = 100000) {
  auto outputTree = new TTree(treeName, treeName);
  auto calcBDst   = BToDstaunu{};
  auto calcBD     = BToDtaunu{};

  vector<PartEmu> cands{};

  for (auto i = 0; i < maxEntries; i++) {
    try {
      cands.push_back(rng->gen());
    } catch (const domain_error& e) {
      break;
    }
  }

  auto   bKey  = any_cast<int>(cands[0]["B_id"]);
  auto   dKey  = any_cast<int>(cands[0]["D_id"]);
  Bool_t isDst = (Abs(dKey) == 413);

  if (Abs(bKey) == 511) {
    calcBDst.SetMasses(0);  // neutral B
  }

  Bool_t ham_ok;
  outputTree->Branch("ham_ok", &ham_ok);
  Bool_t ff_calc_ok;
  outputTree->Branch("ff_calc_ok", &ff_calc_ok);

  double q2_out;
  outputTree->Branch("q2_true", &q2_out);

  double ff_out;
  outputTree->Branch("wff", &ff_out);
  double ff_bgl_out;
  outputTree->Branch("wff_bgl", &ff_bgl_out);
  double ff_bgl_var_ref_out;
  outputTree->Branch("wff_bgl_var_ref", &ff_bgl_var_ref_out);
  double ff_bgl_var_out;
  outputTree->Branch("wff_bgl_var", &ff_bgl_var_out);
  double ff_calc_out;
  outputTree->Branch("wff_calc", &ff_calc_out);

  int b_id_out;
  outputTree->Branch("b_id", &b_id_out);
  int d_id_out;
  outputTree->Branch("d_id", &d_id_out);
  int d_dau_id_out;
  outputTree->Branch("d_dau_id", &d_dau_id_out);
  int pi_id_out;
  outputTree->Branch("pi_id", &pi_id_out);

  double b_pe_out;
  outputTree->Branch("b_pe", &b_pe_out);
  double b_px_out;
  outputTree->Branch("b_px", &b_px_out);
  double b_py_out;
  outputTree->Branch("b_py", &b_py_out);
  double b_pz_out;
  outputTree->Branch("b_pz", &b_pz_out);

  double d_pe_out;
  outputTree->Branch("d_pe", &d_pe_out);
  double d_px_out;
  outputTree->Branch("d_px", &d_px_out);
  double d_py_out;
  outputTree->Branch("d_py", &d_py_out);
  double d_pz_out;
  outputTree->Branch("d_pz", &d_pz_out);

  double d_dau_pe_out;
  outputTree->Branch("d_dau_pe", &d_dau_pe_out);
  double d_dau_px_out;
  outputTree->Branch("d_dau_px", &d_dau_px_out);
  double d_dau_py_out;
  outputTree->Branch("d_dau_py", &d_dau_py_out);
  double d_dau_pz_out;
  outputTree->Branch("d_dau_pz", &d_dau_pz_out);

  double l_pe_out;
  outputTree->Branch("l_pe", &l_pe_out);
  double l_px_out;
  outputTree->Branch("l_px", &l_px_out);
  double l_py_out;
  outputTree->Branch("l_py", &l_py_out);
  double l_pz_out;
  outputTree->Branch("l_pz", &l_pz_out);

  double nu_pe_out;
  outputTree->Branch("nu_pe", &nu_pe_out);
  double nu_px_out;
  outputTree->Branch("nu_px", &nu_px_out);
  double nu_py_out;
  outputTree->Branch("nu_py", &nu_py_out);
  double nu_pz_out;
  outputTree->Branch("nu_pz", &nu_pz_out);

  double pi_pe_out;
  outputTree->Branch("pi_pe", &pi_pe_out);
  double pi_px_out;
  //      $Id: RateCalc.hh,v 1.2 2021/09/09 00:43:50 yipengsun Exp $
  outputTree->Branch("pi_px", &pi_px_out);
  double pi_py_out;
  outputTree->Branch("pi_py", &pi_py_out);
  double pi_pz_out;
  outputTree->Branch("pi_pz", &pi_pz_out);

  double theta_l_out;
  outputTree->Branch("theta_l", &theta_l_out);

  double theta_v_out;
  outputTree->Branch("theta_v", &theta_v_out);
  double chi_out;
  outputTree->Branch("chi", &chi_out);

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

    auto part_B     = buildHamPart(b_p, b_id_out);
    auto part_B_idx = proc.addParticle(part_B);

    auto part_D     = buildHamPart(d_p, d_id_out);
    auto part_D_idx = proc.addParticle(part_D);

    auto part_L     = buildHamPart(l_p, -15);
    auto part_L_idx = proc.addParticle(part_L);

    auto part_NuL     = buildHamPart(nu_p, 16);
    auto part_NuL_idx = proc.addParticle(part_NuL);

    proc.addVertex(part_B_idx, {part_D_idx, part_L_idx, part_NuL_idx});

    // Angles
    theta_l_out = any_cast<double>(cand["theta_l"]);

    if (isDst) {
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

      auto part_D_dau     = buildHamPart(d_dau_p, d_dau_id_out);
      auto part_D_dau_idx = proc.addParticle(part_D_dau);

      auto part_Pi     = buildHamPart(pi_p, pi_id_out);
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
      ff_out             = ham.getWeight("OutputFF");
      ff_bgl_out         = ham.getWeight("OutputFFBGL");
      ff_bgl_var_ref_out = ham.getWeight("OutputFFBGLVarRef");
      ff_bgl_var_out     = ham.getWeight("OutputFFBGLVar");

      if (!isnan(ff_out) && !isinf(ff_out)) {
        ham_ok = true;
        // Compute FF weights w/ Manuel's calculator
        double calc_isgw2, calc_cln, a1, v, a2, a0, fPlus, fMinus;
        auto   cos_theta_l = Cos(theta_l_out);
        auto   cos_theta_v = Cos(theta_v_out);
        if (isDst) {
          calcBDst.ComputeISGW2(q2_out, a1, v, a2, a0);
          calc_isgw2 = calcBDst.Gamma_q2Angular(
              q2_out, cos_theta_l, cos_theta_v, chi_out, false, LEPTON_POSITIVE,
              a1, v, a2, a0, TAU_MASS);

          calcBDst.ComputeCLN(q2_out, a1, v, a2, a0);
          calc_cln = calcBDst.Gamma_q2Angular(q2_out, cos_theta_l, cos_theta_v,
                                              chi_out, false, LEPTON_POSITIVE,
                                              a1, v, a2, a0, TAU_MASS);

        } else {
          calcBD.ComputeISGW2(q2_out, fPlus, fMinus);
          calc_isgw2 =
              calcBD.Gamma_q2tL(q2_out, theta_l_out, fPlus, fMinus, TAU_MASS);

          calcBD.ComputeCLN(q2_out, fPlus, fMinus);
          calc_cln =
              calcBD.Gamma_q2tL(q2_out, theta_l_out, fPlus, fMinus, TAU_MASS);
        }
        ff_calc_out = calc_cln / calc_isgw2;

        if (!isnan(ff_calc_out) && !isinf(ff_calc_out)) ff_calc_ok = true;
        // DEBUG
        // cout << "CLN: " << calc_cln << "; ISGW2: " << calc_isgw2 << endl;

      } else {
        ff_out             = 1.0;
        ff_bgl_out         = 1.0;
        ff_bgl_var_ref_out = 1.0;
        ff_bgl_var_out     = 1.0;
        ff_calc_out        = 1.0;
      }
    }

    outputTree->Fill();
  }

  outputNtp->Write();
  delete outputTree;
}

//////////
// Main //
//////////

int main(int, char** argv) {
  auto outputNtp = new TFile(argv[1], "recreate");
  auto rng       = new TRandom(42);

  Hammer::Hammer ham{};

  setDecays(ham);
  setInputFF(ham);
  setOutputFF(ham);

  ham.setUnits("GeV");
  ham.initRun();

  auto q2Min = TAU_MASS * TAU_MASS;

  auto genBToD =
      new BToDRealGenerator(q2Min, Power(B0_MASS - D0_MASS, 2), 0, PI, rng);
  auto genBToDst = new BToDstUniformGenerator(
      q2Min, Power(B0_MASS - DST_MASS, 2), 0, PI, 0, PI, 0, 2 * PI, rng);

  weightGen(genBToDst, outputNtp, "tree_BDst", ham);
  weightGen(genBToD, outputNtp, "tree_BD", ham);

  delete genBToD;
  delete genBToDst;
  delete rng;
  delete outputNtp;
}
