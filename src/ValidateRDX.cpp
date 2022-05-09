// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon May 09, 2022 at 02:30 PM -0400

#include <any>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <math.h>

#include <TH2D.h>
#include <TMath.h>
#include <TRandom.h>
#include <ROOT/RDataFrame.hxx>

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
using ROOT::RDataFrame;
using ROOT::RDF::RNode;
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
const double mPi      = 0.13957;

const bool LEPTON_POSITIVE = true;

const string D_TREENAME   = "tree_BD";
const string DST_TREENAME = "tree_BDst";

// clang-format off

void setInputFf(Hammer::Hammer& ham) {
  ham.setFFInputScheme({
    {"BD", "ISGW2"},
    {"BD*", "ISGW2"},
  });
}

void setOutputFf(Hammer::Hammer& ham) {
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
  // virtual RNode          gen(RNode df)      = 0;
  virtual vector<string> getOutputBrsName() = 0;

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
  // RNode          gen(RNode df) override;
  vector<string> getOutputBrsName() override;

 protected:
  double   q2Min, q2Max, thetaLMin, thetaLMax;
  TRandom* rng;
  // variable name, output branch name
  map<string, string> outputBrsInt = {
      {"bId", "B_id"},
      {"dId", "D_id"},
      {"lId", "l_id"},
      {"nuId", "nu_id"},
  };
  map<string, string> outputBrsDouble = {{"pB", "B_p"}, {"pD", "D_p"},
                                         {"pL", "l_p"}, {"pNu", "nu_p"},
                                         {"q2", "q2"},  {"thetaL", "theta_l"}};

  double  computeP(double m2Mom, double m2Dau1, double m2Dau2);
  PartEmu genBD(int bId, double mB, int dId, double mD, int lId, double mL,
                int nuId, double q2, double thetaL);
};

BToDUniformGenerator::BToDUniformGenerator(double q2Min, double q2Max,
                                           double thetaLMin, double thetaLMax,
                                           TRandom* rng)
    : q2Min(q2Min),
      q2Max(q2Max),
      thetaLMin(thetaLMin),
      thetaLMax(thetaLMax),
      rng(rng) {}

vector<double> BToDUniformGenerator::get() {
  return {rng->Uniform(q2Min, q2Max), rng->Uniform(thetaLMin, thetaLMax)};
}

PartEmu BToDUniformGenerator::gen() {
  auto inputs = get();
  auto q2     = inputs[0];
  auto thetaL = inputs[1];
  return genBD(521, B_MASS, -421, D0_MASS, -15, TAU_MASS, 16, q2, thetaL);
}

// RNode BToDUniformGenerator::gen(RNode df) {
//   auto result = gen();

//   for (auto [var, br] : outputBrsInt)
//     df = df.Define(br, any_cast<int>(result[var]));
//   for (auto [var, br] : outputBrsDouble)
//     df = df.Define(br, any_cast<double>(result[var]));

//   return df;
// }

vector<string> BToDUniformGenerator::getOutputBrsName() {
  vector<string> result{};

  for (auto [var, br] : outputBrsInt) result.emplace_back(br);
  for (auto [var, br] : outputBrsDouble) result.emplace_back(br);

  return result;
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

  auto lSysP = pB - pD;
  auto pLMag = computeP(q2, mL * mL, 0);

  // Leptons are in the x-z plane
  // Angles are defined in the rest frame of the lepton pair, so rotate first
  // before boosting back to the B rest frame
  auto pLRest  = HFM(Sqrt(mL * mL + pLMag * pLMag), pLMag * Sin(thetaL), 0,
                    pLMag * Cos(thetaL));
  auto pNuRest = HFM(pLMag, -pLRest.px(), -pLRest.py(), -pLRest.pz());

  auto pL  = pLRest.boostFromRestFrameOf(lSysP);
  auto pNu = pNuRest.boostFromRestFrameOf(lSysP);

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
                    double thetaLMax, TRandom* rng, string ff_mode = "ISGW2",
                    int xBins = 200, int yBins = 200);
  ~BToDRealGenerator() override;

  vector<double> get() override;
  PartEmu        gen() override;
  // RNode          gen(RNode df) override;

  void setFF(string ff_mode);

 protected:
  string ffMode;
  TH2D*  histo;
  double q2Min, q2Max, thetaLMin, thetaLMax;
  double q2Step, thetaLStep;

  void buildHisto();
};

BToDRealGenerator::BToDRealGenerator(double q2Min, double q2Max,
                                     double thetaLMin, double thetaLMax,
                                     TRandom* rng, string ffMode, int xBins,
                                     int yBins)
    : BToDUniformGenerator(q2Min, q2Max, thetaLMin, thetaLMax, rng),
      q2Min(q2Min),
      q2Max(q2Max),
      thetaLMin(thetaLMin),
      thetaLMax(thetaLMax),
      ffMode(ffMode) {
  histo      = new TH2D("histo_BD", "histo_BD", xBins, q2Min, q2Max, yBins,
                   thetaLMin, thetaLMax);
  q2Step     = (q2Max - q2Min) / xBins;
  thetaLStep = (thetaLMax - thetaLMin) / yBins;
  buildHisto();
}

BToDRealGenerator::~BToDRealGenerator() { delete histo; }

void BToDRealGenerator::setFF(string ffMode) { ffMode = ffMode; }

vector<double> BToDRealGenerator::get() {
  double q2, thetaL;
  histo->GetRandom2(q2, thetaL, rng);
  return {q2, thetaL};
}

PartEmu BToDRealGenerator::gen() { return BToDUniformGenerator::gen(); }

// RNode BToDRealGenerator::gen(RNode df) { return
// BToDUniformGenerator::gen(df); }

// Protected methods ///////////////////////////////////////////////////////////

// NOTE: We hard-code to use B0 and Tau
void BToDRealGenerator::buildHisto() {
  auto   ffModel = BToDtaunu{};
  double fPlus, fMinus;

  if (ffMode == "ISGW2") {
    for (auto q2 = q2Min + q2Step / 2; q2 <= q2Max - q2Step / 2; q2 += q2Step) {
      ffModel.ComputeISGW2(q2, fPlus, fMinus);
      for (auto thetaL = thetaLMin + thetaLStep / 2;
           thetaL <= thetaLMax - thetaLStep / 2; thetaL += thetaLStep) {
        auto ffVal = ffModel.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);
        histo->Fill(q2, thetaL, ffVal);

        // DEBUG
        // cout << "q2: " << q2 << " thetaL: " << thetaL << " ff val: " <<
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
  // RNode          gen(RNode df) override;

 protected:
  // variable name, output branch name
  map<string, string> outputBrsInt = {
      {"bId", "B_id"},   {"dId", "D_id"},        {"lId", "l_id"},
      {"nuId", "nu_id"}, {"dDauId", "D_dau_id"}, {"piId", "pi_id"},
  };
  map<string, string> outputBrsDouble = {
      {"pB", "B_p"},   {"pD", "D_p"},         {"pL", "l_p"},
      {"pNu", "nu_p"}, {"pDDau", "D_dau_p"},  {"pPi", "pi_p"},
      {"q2", "q2"},    {"thetaL", "theta_l"}, {"thetaV", "theta_v"},
      {"chi", "chi"}};

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
  result.emplace_back(rng->Uniform(thetaVMin, thetaVMax));
  result.emplace_back(rng->Uniform(chiMin, chiMax));

  return result;
}

PartEmu BToDstUniformGenerator::gen() {
  auto inputs = get();
  auto q2     = inputs[0];
  auto thetaL = inputs[1];
  auto thetaV = inputs[2];
  auto chi    = inputs[3];

  return genBDst(511, B0_MASS, -413, DST_MASS, -15, TAU_MASS, 16, -421, D0_MASS,
                 -211, mPi, q2, thetaL, thetaV, chi);
}

// RNode BToDstUniformGenerator::gen(RNode df) {
//   auto result = gen();

//   for (auto [var, br] : outputBrsInt)
//     df = df.Define(br, any_cast<int>(result[var]));
//   for (auto [var, br] : outputBrsDouble)
//     df = df.Define(br, any_cast<double>(result[var]));

//   return df;
// }

// Protected methods
// ///////////////////////////////////////////////////////////
PartEmu BToDstUniformGenerator::genBDst(int bId, double mB, int dId, double mD,
                                        int lId, double mL, int nuId,
                                        int dDauId, double mDDau, int piId,
                                        double mPi, double q2, double thetaL,
                                        double thetaV, double chi) {
  auto result   = genBD(bId, mB, dId, mD, lId, mL, nuId, q2, thetaL);
  auto pDDauMag = computeP(mD * mD, mDDau * mDDau, mPi * mPi);

  // These are defined in the D* rest frame
  auto pDDau_rest =
      HFM(Sqrt(mDDau * mDDau + pDDauMag * pDDauMag),
          pDDauMag * Sin(thetaV) * Cos(chi), pDDauMag * Sin(thetaV) * Sin(chi),
          pDDauMag * Cos(thetaV));
  auto pPiRest = HFM(mD, 0, 0, 0) - pDDau_rest;

  // Boost back to B rest frame from D* rest frame
  auto pD    = any_cast<HFM>(result["pD"]);
  auto pDDau = pDDau_rest.boostFromRestFrameOf(pD);
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

double calcD0FFLib(int bId, int dId, int lId, int nuId, double pB, double pD,
                   double pL, double pNu, double q2, double thetaL) {
  double a1, v, a2, a0, fPlus, fMinus;
  auto   calc = BToDtaunu();
  if (Abs(bId) == 511) calc.SetMasses(0);  // neutral B

  calc.ComputeISGW2(q2, fPlus, fMinus);
  auto ffIsgw2 = calc.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);

  calc.ComputeCLN(q2, fPlus, fMinus);
  auto ffCln = calc.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);

  return ffCln / ffIsgw2;
}

double calcDstFFLib(int bId, int dId, int lId, int nuId, double pB, double pD,
                    double pL, double pNu, double q2, double thetaL,
                    double thetaV, double chi) {
  double a1, v, a2, a0;
  auto   calc = BToDstaunu();
  if (Abs(bId) == 511) calc.SetMasses(0);  // neutral B

  calc.ComputeISGW2(q2, a1, v, a2, a0);
  auto ffIsgw2 = calc.Gamma_q2Angular(q2, Cos(thetaL), Cos(thetaV), chi, false,
                                      LEPTON_POSITIVE, a1, v, a2, a0, TAU_MASS);

  calc.ComputeCLN(q2, a1, v, a2, a0);
  auto ffCln = calc.Gamma_q2Angular(q2, Cos(thetaL), Cos(thetaV), chi, false,
                                    LEPTON_POSITIVE, a1, v, a2, a0, TAU_MASS);

  return ffCln / ffIsgw2;
}

vector<double> ffGen(Hammer::Hammer& ham, int bId, int dId, int lId, int nuId) {
  vector<double> result{};

  return result;
}

template <typename F>
RNode treeGen(IRandGenerator* rng, Hammer::Hammer& ham, F calcAux,
              vector<string> auxInputBrs, int maxEntries = 1e6) {
  auto df = static_cast<RNode>(RDataFrame(maxEntries));
  // df      = rng.gen(df);

  df = df.Define("test", []() { return 1; }, {});
  // df = df.Define("result", [=]() { return rng->gen(); }, {});
  // df      = df.Define("wff_calc", calcAux, auxInputBrs);

  return df;
}

/*
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

  auto   B_key  = any_cast<int>(cands[0]["bId"]);
  auto   D_key  = any_cast<int>(cands[0]["dId"]);
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
  int piId_out;
  output_tree->Branch("piId", &piId_out);

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

  double pLe_out;
  output_tree->Branch("pLe", &pLe_out);
  double pLx_out;
  output_tree->Branch("pLx", &pLx_out);
  double pLy_out;
  output_tree->Branch("pLy", &pLy_out);
  double pLz_out;
  output_tree->Branch("pLz", &pLz_out);

  double pNue_out;
  output_tree->Branch("pNue", &pNue_out);
  double pNux_out;
  output_tree->Branch("pNux", &pNux_out);
  double pNuy_out;
  output_tree->Branch("pNuy", &pNuy_out);
  double pNuz_out;
  output_tree->Branch("pNuz", &pNuz_out);

  double pPie_out;
  output_tree->Branch("pPie", &pPie_out);
  double pPix_out;
  //      $Id: RateCalc.hh,v 1.2 2021/09/09 00:43:50 yipengsun Exp $
  output_tree->Branch("pPix", &pPix_out);
  double pPiy_out;
  output_tree->Branch("pPiy", &pPiy_out);
  double pPiz_out;
  output_tree->Branch("pPiz", &pPiz_out);

  double thetaL_out;
  output_tree->Branch("thetaL", &thetaL_out);

  double thetaV_out;
  output_tree->Branch("thetaV", &thetaV_out);
  double chi_out;
  output_tree->Branch("chi", &chi_out);

  for (auto& cand : cands) {
    Hammer::Process proc;
    ham_ok     = false;
    ff_calc_ok = false;

    q2_out = any_cast<double>(cand["q2"]);

    b_id_out = any_cast<int>(cand["bId"]);
    auto b_p = any_cast<HFM>(cand["pB"]);
    b_pe_out = b_p.E();
    b_px_out = b_p.px();
    b_py_out = b_p.py();
    b_pz_out = b_p.pz();

    d_id_out = any_cast<int>(cand["dId"]);
    auto d_p = any_cast<HFM>(cand["pD"]);
    d_pe_out = d_p.E();
    d_px_out = d_p.px();
    d_py_out = d_p.py();
    d_pz_out = d_p.pz();

    auto pL = any_cast<HFM>(cand["pL"]);
    pLe_out = pL.E();
    pLx_out = pL.px();
    pLy_out = pL.py();
    pLz_out = pL.pz();

    auto pNu = any_cast<HFM>(cand["pNu"]);
    pNue_out = pNu.E();
    pNux_out = pNu.px();
    pNuy_out = pNu.py();
    pNuz_out = pNu.pz();

    auto part_B    = particle(b_p, b_id_out);
    auto part_bIdx = proc.addParticle(part_B);

    auto part_D    = particle(d_p, d_id_out);
    auto part_dIdx = proc.addParticle(part_D);

    auto part_L    = particle(pL, -15);
    auto part_lIdx = proc.addParticle(part_L);

    auto part_NuL    = particle(pNu, 16);
    auto part_NulIdx = proc.addParticle(part_NuL);

    proc.addVertex(part_bIdx, {part_dIdx, part_lIdx, part_NulIdx});

    // Angles
    thetaL_out = any_cast<double>(cand["thetaL"]);

    if (is_Dst) {
      d_dau_id_out = any_cast<int>(cand["dDauId"]);
      auto d_dau_p = any_cast<HFM>(cand["pDDau"]);
      d_dau_pe_out = d_dau_p.E();
      d_dau_px_out = d_dau_p.px();
      d_dau_py_out = d_dau_p.py();
      d_dau_pz_out = d_dau_p.pz();

      piId_out = any_cast<int>(cand["piId"]);
      auto pPi = any_cast<HFM>(cand["pPi"]);
      pPie_out = pPi.E();
      pPix_out = pPi.px();
      pPiy_out = pPi.py();
      pPiz_out = pPi.pz();

      auto part_D_dau   = particle(d_dau_p, d_dau_id_out);
      auto part_dDauIdx = proc.addParticle(part_D_dau);

      auto part_Pi    = particle(pPi, piId_out);
      auto part_piIdx = proc.addParticle(part_Pi);

      proc.addVertex(part_dIdx, {part_dDauIdx, part_piIdx});

      // Additional angles
      thetaV_out = any_cast<double>(cand["thetaV"]);
      chi_out    = any_cast<double>(cand["chi"]);
    } else {
      d_dau_id_out = piId_out = 0;
      d_dau_pe_out = d_dau_px_out = d_dau_py_out = d_dau_pz_out = 0.;
      pPie_out = pPix_out = pPiy_out = pPiz_out = 0.;

      thetaV_out = chi_out = 0.;
    }

    ham.initEvent();
    auto proc_id = ham.addProcess(proc);

    if (proc_id != 0) {
      ham.processEvent();
      ff_out = ham.getWeight("OutputFF");

      if (!isnan(ff_out) && !isinf(ff_out)) {
        ham_ok = true;
        // Compute FF weights w/ Manuel's calculator
        double calc_isgw2, calc_cln, a1, v, a2, a0, fPlus, fMinus;
        auto   cos_thetaL = Cos(thetaL_out);
        auto   cos_thetaV = Cos(thetaV_out);
        if (is_Dst) {
          calc_BDst.ComputeISGW2(q2_out, a1, v, a2, a0);
          calc_isgw2 = calc_BDst.Gamma_q2Angular(
              q2_out, cos_thetaL, cos_thetaV, chi_out, false, LEPTON_POSITIVE,
              a1, v, a2, a0, TAU_MASS);

          calc_BDst.ComputeCLN(q2_out, a1, v, a2, a0);
          calc_cln = calc_BDst.Gamma_q2Angular(q2_out, cos_thetaL, cos_thetaV,
                                               chi_out, false, LEPTON_POSITIVE,
                                               a1, v, a2, a0, TAU_MASS);

        } else {
          calc_BD.ComputeISGW2(q2_out, fPlus, fMinus);
          calc_isgw2 =
              calc_BD.Gamma_q2tL(q2_out, thetaL_out, fPlus, fMinus, TAU_MASS);

          calc_BD.ComputeCLN(q2_out, fPlus, fMinus);
          calc_cln =
              calc_BD.Gamma_q2tL(q2_out, thetaL_out, fPlus, fMinus, TAU_MASS);
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
*/

//////////
// Main //
//////////

int main(int argc, char** argv) {
  // auto output_ntp = new TFile(argv[1], "recreate");

  auto           rng = new TRandom(42);
  Hammer::Hammer ham{};

  setDecays(ham);
  setInputFf(ham);
  setOutputFf(ham);

  ham.setUnits("GeV");
  ham.initRun();

  auto q2Min    = TAU_MASS * TAU_MASS;
  auto mB       = B0_MASS;
  auto q2MaxD0  = Power(mB - D0_MASS, 2);
  auto q2MaxDst = Power(mB - DST_MASS, 2);

  auto genBToD = BToDRealGenerator(q2Min, q2MaxD0, 0, PI, rng);
  auto genBToDst =
      BToDstUniformGenerator(q2Min, q2MaxDst, 0, PI, 0, PI, 0, 2 * PI, rng);

  auto brsD   = genBToD.getOutputBrsName();
  auto brsDst = genBToDst.getOutputBrsName();

  auto dfD = treeGen(&genBToD, ham, calcD0FFLib, brsD);

  dfD.Snapshot(argv[1], D_TREENAME);

  delete rng;
}
