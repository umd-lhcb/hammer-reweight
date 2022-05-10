// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Tue May 10, 2022 at 05:05 AM -0400

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
    // {"BD*", "BGLVar_2"},
  });

  // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
  ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET2
  // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
  ham.setOptions("BtoD*CLN_2: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2

  ham.setFFEigenvectors("BtoD", "BGLVar_1", {{"delta_ap2", 0.1}});
  // ham.setFFEigenvectors("BtoD*", "BGLVar_2", {{"delta_a2", 10.0}});
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

  result["bId"] = bId;
  result["pB"]  = pB;

  result["dId"] = dId;
  result["pD"]  = pD;

  result["thetaL"] = thetaL;  // the only physical angle
  result["lId"]    = lId;
  result["pL"]     = pL;
  result["nuId"]   = nuId;
  result["pNu"]    = pNu;

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
  auto pD    = any_cast<HFM>(result["pD"]);
  auto pDDau = pDDauRest.boostFromRestFrameOf(pD);
  auto pPi   = pPiRest.boostFromRestFrameOf(pD);

  result["thetaV"] = thetaV;
  result["chi"]    = chi;

  result["dDauId"] = dDauId;
  result["pDDau"]  = pDDau;
  result["piId"]   = piId;
  result["pPi"]    = pPi;

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
  for (auto i = 0; i < maxEntries; i++) try {
      cands.push_back(rng->gen());
    } catch (const domain_error& e) {
      break;
    };

  auto bKey  = any_cast<int>(cands[0]["bId"]);
  auto dKey  = any_cast<int>(cands[0]["dId"]);
  bool isDst = (Abs(dKey) == 413);

  if (Abs(bKey) == 511) calcBDst.SetMasses(0);  // neutral B

  bool hamOk;
  outputTree->Branch("ham_ok", &hamOk);
  bool ffCalcOk;
  outputTree->Branch("ff_calc_ok", &ffCalcOk);

  double q2;
  outputTree->Branch("q2_true", &q2);

  double ff;
  outputTree->Branch("wff", &ff);
  double ffBgl;
  outputTree->Branch("wff_bgl", &ffBgl);
  double ffBglVarRef;
  outputTree->Branch("wff_bgl_var_ref", &ffBglVarRef);
  double ffBglVar;
  outputTree->Branch("wff_bgl_var", &ffBglVar);
  double ffCalc;
  outputTree->Branch("wff_calc", &ffCalc);

  int bId;
  outputTree->Branch("b_id", &bId);
  int dId;
  outputTree->Branch("d_id", &dId);
  int dDauId;
  outputTree->Branch("d_dau_id", &dDauId);
  int piId;
  outputTree->Branch("pi_id", &piId);

  double bPE;
  outputTree->Branch("b_pe", &bPE);
  double bPX;
  outputTree->Branch("b_px", &bPX);
  double bPY;
  outputTree->Branch("b_py", &bPY);
  double bPZ;
  outputTree->Branch("b_pz", &bPZ);

  double dPE;
  outputTree->Branch("d_pe", &dPE);
  double dPX;
  outputTree->Branch("d_px", &dPX);
  double dPY;
  outputTree->Branch("d_py", &dPY);
  double dPZ;
  outputTree->Branch("d_pz", &dPZ);

  double dDauPE;
  outputTree->Branch("d_dau_pe", &dDauPE);
  double dDauPX;
  outputTree->Branch("d_dau_px", &dDauPX);
  double dDauPY;
  outputTree->Branch("d_dau_py", &dDauPY);
  double dDauPZ;
  outputTree->Branch("d_dau_pz", &dDauPZ);

  double lPE;
  outputTree->Branch("l_pe", &lPE);
  double lPX;
  outputTree->Branch("l_px", &lPX);
  double lPY;
  outputTree->Branch("l_py", &lPY);
  double lPZ;
  outputTree->Branch("l_pz", &lPZ);

  double nuPE;
  outputTree->Branch("nu_pe", &nuPE);
  double nuPX;
  outputTree->Branch("nu_px", &nuPX);
  double nuPY;
  outputTree->Branch("nu_py", &nuPY);
  double nuPZ;
  outputTree->Branch("nu_pz", &nuPZ);

  double piPE;
  outputTree->Branch("pi_pe", &piPE);
  double piPX;
  //      $Id: RateCalc.hh,v 1.2 2021/09/09 00:43:50 yipengsun Exp $
  outputTree->Branch("pi_px", &piPX);
  double piPY;
  outputTree->Branch("pi_py", &piPY);
  double piPZ;
  outputTree->Branch("pi_pz", &piPZ);

  double thetaL;
  outputTree->Branch("theta_l", &thetaL);

  double thetaV;
  outputTree->Branch("theta_v", &thetaV);
  double chi;
  outputTree->Branch("chi", &chi);

  for (auto& cand : cands) {
    Hammer::Process proc;
    hamOk    = true;
    ffCalcOk = true;

    q2 = any_cast<double>(cand["q2"]);

    bId     = bKey;
    auto pB = any_cast<HFM>(cand["pB"]);
    bPE     = pB.E();
    bPX     = pB.px();
    bPY     = pB.py();
    bPZ     = pB.pz();

    dId     = dKey;
    auto pD = any_cast<HFM>(cand["pD"]);
    dPE     = pD.E();
    dPX     = pD.px();
    dPY     = pD.py();
    dPZ     = pD.pz();

    auto pL = any_cast<HFM>(cand["pL"]);
    lPE     = pL.E();
    lPX     = pL.px();
    lPY     = pL.py();
    lPZ     = pL.pz();

    auto pNu = any_cast<HFM>(cand["pNu"]);
    nuPE     = pNu.E();
    nuPX     = pNu.px();
    nuPY     = pNu.py();
    nuPZ     = pNu.pz();

    auto partB    = buildHamPart(pB, bId);
    auto partBIdx = proc.addParticle(partB);

    auto partD    = buildHamPart(pD, dId);
    auto partDIdx = proc.addParticle(partD);

    auto partL    = buildHamPart(pL, -15);
    auto partLIdx = proc.addParticle(partL);

    auto partNuL    = buildHamPart(pNu, 16);
    auto partNuLIdx = proc.addParticle(partNuL);

    proc.addVertex(partBIdx, {partDIdx, partLIdx, partNuLIdx});

    // Angles
    thetaL = any_cast<double>(cand["thetaL"]);

    if (isDst) {
      dDauId     = any_cast<int>(cand["dDauId"]);
      auto pDDau = any_cast<HFM>(cand["pDDau"]);
      dDauPE     = pDDau.E();
      dDauPX     = pDDau.px();
      dDauPY     = pDDau.py();
      dDauPZ     = pDDau.pz();

      piId     = any_cast<int>(cand["piId"]);
      auto pPi = any_cast<HFM>(cand["pPi"]);
      piPE     = pPi.E();
      piPX     = pPi.px();
      piPY     = pPi.py();
      piPZ     = pPi.pz();

      auto partDDau    = buildHamPart(pDDau, dDauId);
      auto partDDauIdx = proc.addParticle(partDDau);

      auto partPi    = buildHamPart(pPi, piId);
      auto partPiIdx = proc.addParticle(partPi);

      proc.addVertex(partDIdx, {partDDauIdx, partPiIdx});

      // Additional angles
      thetaV = any_cast<double>(cand["thetaV"]);
      chi    = any_cast<double>(cand["chi"]);
    } else {
      dDauId = piId = 0;
      dDauPE = dDauPX = dDauPY = dDauPZ = 0.;
      piPE = piPX = piPY = piPZ = 0.;
      thetaV = chi = 0.;
    }

    ham.initEvent();
    auto procId = ham.addProcess(proc);
    if (procId != 0) {
      try {
        ham.processEvent();
        ff    = ham.getWeight("OutputFF");
        ffBgl = ham.getWeight("OutputFFBGL");
      } catch (const std::exception& e) {
        hamOk = false;
      }

      if (hamOk)
        if (isnan(ff) || isinf(ff) || isnan(ffBgl) || isinf(ffBgl))
          hamOk = false;

      if (hamOk) {
        // Compute FF variations
        auto mapFFOut = map<string, double*>{
            {"OutputFFBGLVarRef", &ffBglVarRef},
            {"OutputFFBGLVar", &ffBglVar},
        };

        for (auto [ffScheme, outBrPtr] : mapFFOut) {
          try {
            *outBrPtr = ham.getWeight(ffScheme);
          } catch (const std::exception& e) {
            *outBrPtr = 1.0;
          }
        }
        // Compute FF weights w/ Manuel's calculator
        double calcIsgw2, calcCln, a1, v, a2, a0, fPlus, fMinus;
        auto   cosThetaL = Cos(thetaL);
        auto   cosThetaV = Cos(thetaV);
        if (isDst) {
          calcBDst.ComputeISGW2(q2, a1, v, a2, a0);
          calcIsgw2 = calcBDst.Gamma_q2Angular(q2, cosThetaL, cosThetaV, chi,
                                               false, LEPTON_POSITIVE, a1, v,
                                               a2, a0, TAU_MASS);

          calcBDst.ComputeCLN(q2, a1, v, a2, a0);
          calcCln = calcBDst.Gamma_q2Angular(q2, cosThetaL, cosThetaV, chi,
                                             false, LEPTON_POSITIVE, a1, v, a2,
                                             a0, TAU_MASS);
        } else {
          calcBD.ComputeISGW2(q2, fPlus, fMinus);
          calcIsgw2 = calcBD.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);

          calcBD.ComputeCLN(q2, fPlus, fMinus);
          calcCln = calcBD.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);
        }
        ffCalc = calcCln / calcIsgw2;
        if (isnan(ffCalc) || isinf(ffCalc)) ffCalcOk = false;
      } else {
        ff          = 1.0;
        ffBgl       = 1.0;
        ffBglVarRef = 1.0;
        ffBglVar    = 1.0;
        ffCalc      = 1.0;
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
