// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sun Jun 12, 2022 at 03:11 AM -0400

#include <any>
#include <chrono>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <boost/histogram.hpp>
#include <boost/histogram/fwd.hpp>

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
using namespace std::chrono;
using namespace boost::histogram;
using namespace TMath;

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

void setBtoDBGLDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {ChiT: 0.0005131}");
  ham.setOptions(scheme + ": {ChiL: 0.006332}");
  ham.setOptions(scheme + ": {BcStatesp: [6.329, 6.92, 7.02]}");
  ham.setOptions(scheme + ": {BcStates0: [6.716, 7.121]}");
  ham.setOptions(scheme + ": {ap: [0.01566, -0.0342, -0.09, 0.0]}");
  ham.setOptions(scheme + ": {a0: [0.07935, -0.205, -0.23, 0.0]}");
}

void setBtoDstarBGLDefault(Hammer::Hammer& ham, const string scheme) {
  ham.setOptions(scheme + ": {Vcb: 0.0415}");
  ham.setOptions(scheme + ": {Chim: 0.0003068}");
  ham.setOptions(scheme + ": {Chip: 0.000528}");
  ham.setOptions(scheme + ": {ChimL: 0.002466}");
  ham.setOptions(scheme + ": {avec: [0.00133258741, -0.0060989894, -0.02506434]}");
  ham.setOptions(scheme + ": {bvec: [0.0005188318380000001, 0.00015456343000000002, 0.0008354780000000001]}");
  ham.setOptions(scheme + ": {cvec: [6.266085e-06, 0.0032583642]}");
  ham.setOptions(scheme + ": {dvec: [0.00220148453, -0.0081459105]}");
}

void setOutputFF(Hammer::Hammer& ham) {
  ham.addFFScheme("OutputFF", {{"BD", "CLN_1"}, {"BD*", "CLN_2"}});
  ham.addFFScheme("OutputFFBGL", {{"BD", "BGL_1"}, {"BD*", "BGL_1"}});
  ham.addFFScheme("OutputFFBGLVarShift", {{"BD", "BGL_2"}});
  ham.addFFScheme("OutputFFBGLVar", {{"BD", "BGLVar_1"}});
  ham.addFFScheme("OutputFFBGLN3", {{"BD", "BGL_3"}, {"BD*", "BGL_3"}});

  // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
  ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET2
  // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
  ham.setOptions("BtoD*CLN_2: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2

  // D0, BGL settings
  setBtoDBGLDefault(ham, "BtoDBGL_1");
  // +1 variation by shifting nominal
  setBtoDBGLDefault(ham, "BtoDBGL_2");
  ham.setOptions("BtoDBGL_2: {ap: [0.015642660612597052, -0.034035685234997386, -0.0898070271651518, 0.0]}");
  ham.setOptions("BtoDBGL_2: {a0: [0.0858147593040351, -0.3050702334404309, -0.22752112980378628, 0.0]}");
  // +1 variation w/ BGLVar
  setBtoDBGLDefault(ham, "BtoDBGLVar_1");

  // D*, BGL settings
  setBtoDstarBGLDefault(ham, "BtoD*BGL_1");
}

// FF variations for B -> D, in u1p, u1m, u2p, u2m, ... order
auto varParamsB2DBGL = vector<map<string, double>> {
  {{"delta_ap0", -1.7339387402949836e-05}, {"delta_ap1", 0.0001643147650026172}, {"delta_ap2", 0.0001929728348481839}, {"delta_a01", -0.10007023344043088}, {"delta_a02", 0.0024788701962137416}, {"delta_a00", 0.006464759304035082}},
  {{"delta_ap0", 1.7339387402949836e-05}, {"delta_ap1", -0.0001643147650026172}, {"delta_ap2", -0.0001929728348481839}, {"delta_a01", 0.10007023344043088}, {"delta_a02", -0.0024788701962137416}, {"delta_a00", -0.006464759304035082}},
  {{"delta_ap0", -0.0001238406964858629}, {"delta_ap1", 0.0014520126238037846}, {"delta_ap2", 0.003897645034731049}, {"delta_a01", -0.0005198293137311275}, {"delta_a02", -0.02138567598323755}, {"delta_a00", 4.7861119207580826e-05}},
  {{"delta_ap0", 0.0001238406964858629}, {"delta_ap1", -0.0014520126238037846}, {"delta_ap2", -0.003897645034731049}, {"delta_a01", 0.0005198293137311275}, {"delta_a02", 0.02138567598323755}, {"delta_a00", -4.7861119207580826e-05}},
  {{"delta_ap0", -0.0004745338176480161}, {"delta_ap1", -0.013866861740418688}, {"delta_ap2", 0.0003080116100249672}, {"delta_a01", -4.393049104058412e-05}, {"delta_a02", -0.0008815590134471655}, {"delta_a00", -0.006792469545215629}},
  {{"delta_ap0", 0.0004745338176480161}, {"delta_ap1", 0.013866861740418688}, {"delta_ap2", -0.0003080116100249672}, {"delta_a01", 4.393049104058412e-05}, {"delta_a02", 0.0008815590134471655}, {"delta_a00", 0.006792469545215629}},
  {{"delta_ap0", 2.865149761619912e-06}, {"delta_ap1", 8.988500416478102e-07}, {"delta_ap2", 9.41795313532076e-05}, {"delta_a01", 6.085188006960159e-07}, {"delta_a02", 1.7194330288888333e-05}, {"delta_a00", 1.6454168439027227e-05}},
  {{"delta_ap0", -2.865149761619912e-06}, {"delta_ap1", -8.988500416478102e-07}, {"delta_ap2", -9.41795313532076e-05}, {"delta_a01", -6.085188006960159e-07}, {"delta_a02", -1.7194330288888333e-05}, {"delta_a00", -1.6454168439027227e-05}},
  {{"delta_ap0", -0.0010139800624278266}, {"delta_ap1", 3.447442492334894e-05}, {"delta_ap2", 2.808344227692626e-05}, {"delta_a01", 6.163074193685879e-07}, {"delta_a02", 1.3315838706715554e-05}, {"delta_a00", -0.005048232266588653}},
  {{"delta_ap0", 0.0010139800624278266}, {"delta_ap1", -3.447442492334894e-05}, {"delta_ap2", -2.808344227692626e-05}, {"delta_a01", -6.163074193685879e-07}, {"delta_a02", -1.3315838706715554e-05}, {"delta_a00", 0.005048232266588653}},
};
// clang-format on

map<string, complex<double>> specializedWC = {
    {"SM", 1},     {"S_qLlL", 0}, {"S_qRlL", 0}, {"V_qLlL", 0},
    {"V_qRlL", 0}, {"T_qLlL", 0}, {"S_qLlR", 0}, {"S_qRlR", 0},
    {"V_qLlR", 0}, {"V_qRlR", 0}, {"T_qRlR", 0}};

void setDecays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BD*TauNu");
}

//////////////////
// Stat helpers //
//////////////////

// Get histogram bin edges
template <typename Axes>
auto getBinEdges(const histogram<Axes>& histo) {
  vector<pair<double, double>> binEdges;
  const size_t                 ndim = histo.rank();
  for (int axisIdx = 0; axisIdx < ndim; axisIdx++) {
    auto   axis   = histo.axis(axisIdx);
    double binLow = axis.bin(0).lower();
    double binHi  = axis.bin(axis.size() - 1).upper();
    binEdges.emplace_back(pair<double, double>{binLow, binHi});
  }
  return binEdges;
}

template <typename Axes>
double getMaxBinCount(const histogram<Axes>& histo) {
  double max = 0;
  for (auto x : indexed(histo)) {
    const double binCount = *x;
    if (binCount >= max) max = binCount;
  }
  return max;
}

template <typename Axes, typename T>
auto getBinIdx(const histogram<Axes>& histo, vector<T> vals) {
  vector<unsigned int> binIdx{};
  for (int axisIdx = 0; axisIdx < vals.size(); axisIdx++)
    binIdx.emplace_back(histo.axis(axisIdx).index(vals[axisIdx]));
  return binIdx;
}

// Naive MC sampling
template <typename Axes>
auto getRand(const histogram<Axes>& histo, TRandom* rng, int maxTry = 1000) {
  auto         binEdges    = getBinEdges(histo);
  auto         maxBinCount = getMaxBinCount(histo);
  const size_t ndim        = histo.rank();

  int            numTry = 0;
  vector<double> randVals;
  while (numTry <= maxTry) {
    randVals = {};
    // generate some random numbers
    for (const auto& [min, max] : binEdges)
      randVals.emplace_back(rng->Uniform(min, max));

    auto   binIdx = getBinIdx(histo, randVals);
    double freq   = rng->Uniform(0, maxBinCount);
    if (freq <= histo.at(binIdx)) {
      // for (const auto& elem : randVals) cout << elem << " ";
      // cout << endl;
      return randVals;
    }

    numTry += 1;
  }
  return randVals;
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
};

/////////////////////////////////////////////
// Event generation: D0, real distribution //
/////////////////////////////////////////////

typedef histogram<vector<axis::regular<double>>> HistoND;

class BToDRealGenerator : public IRandGenerator {
 public:
  BToDRealGenerator(){};
  BToDRealGenerator(double q2Min, double q2Max, double thetaLMin,
                    double thetaLMax, TRandom* rng, string ffMode = "ISGW2",
                    int xBins = 300, int yBins = 300);

  vector<double> get() override;
  PartEmu        gen() override;

  void setFF(string ffMode);

 protected:
  TRandom* rng;
  double   q2Min, q2Max, q2Step;
  double   thetaLMin, thetaLMax, thetaLStep;
  int      xBins, yBins;
  string   ffMode;
  HistoND  histo;

  void    buildHisto();
  double  computeP(double m2Mom, double m2Dau1, double m2Dau2);
  PartEmu genBD(int bID, double mB, int dID, double mD, int lId, double mL,
                int nuId, double q2, double thetaL);
};

BToDRealGenerator::BToDRealGenerator(double q2Min, double q2Max,
                                     double thetaLMin, double thetaLMax,
                                     TRandom* rng, string ffMode, int xBins,
                                     int yBins)
    : q2Min(q2Min),
      q2Max(q2Max),
      thetaLMin(thetaLMin),
      thetaLMax(thetaLMax),
      xBins(yBins),
      yBins(yBins),
      rng(rng),
      ffMode(ffMode) {
  q2Step     = (q2Max - q2Min) / xBins;
  thetaLStep = (thetaLMax - thetaLMin) / yBins;
  buildHisto();
}

void BToDRealGenerator::setFF(string ffMode) { this->ffMode = ffMode; }

vector<double> BToDRealGenerator::get() { return getRand(histo, rng); }

PartEmu BToDRealGenerator::gen() {
  auto inputs = get();
  auto q2     = inputs[0];
  auto thetaL = inputs[1];
  return genBD(521, B_MASS, -421, D0_MASS, -15, TAU_MASS, 16, q2, thetaL);
}

// Protected methods ///////////////////////////////////////////////////////////

double BToDRealGenerator::computeP(double m2Mom, double m2Dau1, double m2Dau2) {
  auto denom = 2 * Sqrt(m2Mom);
  auto nom   = Sqrt(m2Mom * m2Mom + m2Dau1 * m2Dau1 + m2Dau2 * m2Dau2 -
                  2 * (m2Mom * m2Dau1 + m2Mom * m2Dau2 + m2Dau1 * m2Dau2));
  return nom / denom;
}

// Everything's in GeV!
PartEmu BToDRealGenerator::genBD(int bId, double mB, int dId, double mD,
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

// NOTE: We hard-code to use B0 and Tau
void BToDRealGenerator::buildHisto() {
  auto axes = vector<axis::regular<double>>{
      axis::regular<double>(xBins, q2Min, q2Max),
      axis::regular<double>(yBins, thetaLMin, thetaLMax),
  };
  histo          = make_histogram(std::move(axes));
  auto   ffModel = BToDtaunu{};
  double fPlus, fMinus;

  // using the center of each bin
  for (auto q2 = q2Min + q2Step / 2; q2 <= q2Max - q2Step / 2; q2 += q2Step) {
    if (ffMode == "ISGW2")
      ffModel.ComputeISGW2(q2, fPlus, fMinus);
    else if (ffMode == "CLN")
      ffModel.ComputeCLN(q2, fPlus, fMinus);

    for (auto thetaL = thetaLMin + thetaLStep / 2;
         thetaL <= thetaLMax - thetaLStep / 2; thetaL += thetaLStep) {
      auto ffVal = ffModel.Gamma_q2tL(q2, thetaL, fPlus, fMinus, TAU_MASS);
      histo(q2, thetaL, weight(ffVal));
    }
  }
}

//////////////////////////
// Event generation: D* //
//////////////////////////

class BToDstRealGenerator : public BToDRealGenerator {
 public:
  BToDstRealGenerator(double q2Min, double q2Max, double thetaLMin,
                      double thetaLMax, double thetaVMin, double thetaVMax,
                      double chiMin, double chiMax, TRandom* rng, string ffMode,
                      int xBins, int yBins, int zBins, int wBins);

  vector<double> get() override;
  PartEmu        gen() override;

 protected:
  TRandom* rng;
  double   q2Min, q2Max, q2Step;
  double   thetaLMin, thetaLMax, thetaLStep;
  double   thetaVMin, thetaVMax, thetaVStep;
  double   chiMin, chiMax, chiStep;
  int      xBins, yBins, zBins, wBins;
  string   ffMode;
  HistoND  histo;

  void    buildHisto();
  PartEmu genBDst(int bId, double mB, int dId, double mD, int lId, double mL,
                  int nuId, int dDauId, double mDDau, int piId, double mPi,
                  double q2, double thetaL, double thetaV, double chi);

 private:
  friend class BToDRealGenerator;
};

BToDstRealGenerator::BToDstRealGenerator(double q2Min, double q2Max,
                                         double thetaLMin, double thetaLMax,
                                         double thetaVMin, double thetaVMax,
                                         double chiMin, double chiMax,
                                         TRandom* rng, string ffMode = "ISGW2",
                                         int xBins = 100, int yBins = 50,
                                         int zBins = 50, int wBins = 50)
    : q2Min(q2Min),
      q2Max(q2Max),
      thetaLMin(thetaLMin),
      thetaLMax(thetaLMax),
      thetaVMin(thetaVMin),
      thetaVMax(thetaVMax),
      chiMin(chiMin),
      chiMax(chiMax),
      xBins(xBins),
      yBins(yBins),
      zBins(zBins),
      wBins(wBins),
      rng(rng),
      ffMode(ffMode) {
  q2Step = (q2Max - q2Min) / xBins;

  thetaLStep = (thetaLMax - thetaLMin) / yBins;
  thetaVStep = (thetaVMax - thetaVMin) / zBins;
  chiStep    = (chiMax - chiMin) / wBins;
  buildHisto();
}

vector<double> BToDstRealGenerator::get() { return getRand(histo, rng); }

PartEmu BToDstRealGenerator::gen() {
  auto inputs = get();
  auto q2     = inputs[0];
  auto thetaL = inputs[1];
  auto thetaV = inputs[2];
  auto chi    = inputs[3];
  return genBDst(511, B0_MASS, -413, DST_MASS, -15, TAU_MASS, 16, -421, D0_MASS,
                 -211, PI_MASS, q2, thetaL, thetaV, chi);
}

// Protected methods ///////////////////////////////////////////////////////////
PartEmu BToDstRealGenerator::genBDst(int bId, double mB, int dId, double mD,
                                     int lId, double mL, int nuId, int dDauId,
                                     double mDDau, int piId, double mPi,
                                     double q2, double thetaL, double thetaV,
                                     double chi) {
  auto result = genBD(bId, mB, dId, mD, lId, mL, nuId, q2, thetaL);
  auto pDDauMag =
      BToDRealGenerator::computeP(mD * mD, mDDau * mDDau, mPi * mPi);

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

void BToDstRealGenerator::buildHisto() {
  auto axes = vector<axis::regular<double>>{
      axis::regular<double>(xBins, q2Min, q2Max),
      axis::regular<double>(yBins, thetaLMin, thetaLMax),
      axis::regular<double>(zBins, thetaVMin, thetaVMax),
      axis::regular<double>(wBins, chiMin, chiMax)};
  histo          = make_histogram(std::move(axes));
  auto   ffModel = BToDstaunu{};
  double A1, V, A2, A0;

  // using the center of each bin
  for (auto q2 = q2Min + q2Step / 2; q2 <= q2Max - q2Step / 2; q2 += q2Step) {
    if (ffMode == "ISGW2")
      ffModel.ComputeISGW2(q2, A1, V, A2, A0);
    else if (ffMode == "CLN")
      ffModel.ComputeCLN(q2, A1, V, A2, A0);

    for (auto thetaL = thetaLMin + thetaLStep / 2;
         thetaL <= thetaLMax - thetaLStep / 2; thetaL += thetaLStep) {
      for (auto thetaV = thetaVMin + thetaVStep / 2;
           thetaV <= thetaVMax - thetaLStep / 2; thetaV += thetaVStep) {
        for (auto chi = chiMin + chiStep / 2; chi <= chiMax - chiStep / 2;
             chi += chiStep) {
          auto ffVal =
              ffModel.Gamma_q2Angular(q2, Cos(thetaL), Cos(thetaV), chi, 0,
                                      false, A1, V, A2, A0, TAU_MASS);
          histo(q2, thetaL, thetaV, chi, weight(ffVal));
        }
      }
    }
  }
}

/////////////////
// Reweighting //
/////////////////

void weightGen(IRandGenerator* rng, TFile* outputNtp, TString treeName,
               Hammer::Hammer& ham, int maxEntries = 5e4) {
  cout << "Start generating tree: " << treeName << endl;
  auto outputTree = new TTree(treeName, treeName);
  auto calcBDst   = BToDstaunu{};
  auto calcBD     = BToDtaunu{};

  auto timeNoVar     = microseconds(0);
  auto timeVarP      = microseconds(0);
  auto timeVarPShift = microseconds(0);

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
  double ffBglN3;
  outputTree->Branch("wff_bgl_n3", &ffBglN3);
  double ffBglVarP;
  outputTree->Branch("wff_bgl_var_p", &ffBglVarP);
  double ffBglVarPShift;
  outputTree->Branch("wff_bgl_var_p_shift", &ffBglVarPShift);
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
        ff = ham.getWeight("OutputFF");
      } catch (const std::exception& e) {
        hamOk = false;
        ff    = 1.0;
      }

      try {
        // reference weight at N=3
        ffBglN3 = ham.getWeight("OutputFFBGLN3");

        auto startNoVar = high_resolution_clock::now();
        ffBgl           = ham.getWeight("OutputFFBGL");
        auto stopNovar  = high_resolution_clock::now();
        timeNoVar += duration_cast<microseconds>(stopNovar - startNoVar);
      } catch (const std::exception& e) {
        ffBgl = ffBglN3 = 1.0;
      }

      if (hamOk && (isnan(ff) || isinf(ff) || isnan(ffBgl) || isinf(ffBgl)))
        hamOk = false;

      if (hamOk) {
        // Compute FF variations: u1p
        try {
          auto varParams = varParamsB2DBGL[0];
          auto startVarP = high_resolution_clock::now();

          ham.setFFEigenvectors("BtoD", "BGLVar_1", varParams);
          ffBglVarP = ham.getWeight("OutputFFBGLVar");
          ham.resetFFEigenvectors("BtoD", "BGLVar_1");

          auto stopVarP = high_resolution_clock::now();
          timeVarP += duration_cast<microseconds>(stopVarP - startVarP);
        } catch (const std::exception& e) {
          ffBglVarP = 1.0;
        }
        // Compute FF variations u1p by shifting the central values
        try {
          auto startVarPShift = high_resolution_clock::now();
          ffBglVarPShift      = ham.getWeight("OutputFFBGLVarShift");
          auto stopVarPShift  = high_resolution_clock::now();
          timeVarPShift +=
              duration_cast<microseconds>(stopVarPShift - startVarPShift);
        } catch (const std::exception& e) {
          ffBglVarPShift = 1.0;
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
        ff             = 1.0;
        ffBgl          = 1.0;
        ffBglVarP      = 1.0;
        ffBglVarPShift = 1.0;
        ffCalc         = 1.0;
      }
    }
    outputTree->Fill();
  }

  outputNtp->Write();
  delete outputTree;

  cout << "The no variation BGL took " << timeNoVar.count() << " us to execute."
       << endl;
  cout << "The +1 variation BGL took " << timeVarP.count() << " us to execute."
       << endl;
  cout << "The +1 by shift nominal BGL took " << timeVarPShift.count()
       << " us to execute." << endl;
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

  // only use SM Wilson coefficients
  ham.specializeWCInWeights("BtoCTauNu", specializedWC);
  ham.specializeWCInWeights("BtoCMuNu", specializedWC);

  auto q2Min = TAU_MASS * TAU_MASS;

  auto genBToD =
      new BToDRealGenerator(q2Min, Power(B0_MASS - D0_MASS, 2), 0, PI, rng);
  auto genBToDst = new BToDstRealGenerator(q2Min, Power(B0_MASS - DST_MASS, 2),
                                           0, PI, 0, PI, 0, 2 * PI, rng);

  weightGen(genBToD, outputNtp, "tree_BD", ham, 1e5);
  weightGen(genBToDst, outputNtp, "tree_BDst", ham, 3e4);

  delete genBToD;
  delete genBToDst;
  delete rng;
  delete outputNtp;
}
