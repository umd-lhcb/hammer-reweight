// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Fri Jun 17, 2022 at 11:08 PM -0400

#include <algorithm>
#include <array>
#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <math.h>

#include <TMath.h>
#include <TString.h>
#include <ROOT/RDataFrame.hxx>

#include <Hammer/Hammer.hh>
#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>
#include <Hammer/Process.hh>

#include <cxxopts.hpp>

#include "const.h"
#include "utils_general.h"
#include "utils_ham.h"

using namespace std;
using ROOT::RDataFrame;
using ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

#define DEBUG_OUT cout << __LINE__ << endl;

///////////////////
// Configurables //
///////////////////

//#define DEBUG_CLI
//#define FORCE_MOMENTUM_CONSERVATION_LEPTONIC
#define RADIATIVE_CORRECTION
#define SOFT_PHOTON_THRESH 0.1

typedef map<vector<Int_t>, unsigned long> DecayFreq;

void setInputFF(Hammer::Hammer& ham, TString run) {
  if (run == "run1") {
    ham.setFFInputScheme({
        {"BD", "ISGW2"},  // 12573010
        {"BD*", "ISGW2"}  // 11574010
        // for run 1, B -> D*MuNu is modelled w/ HQET, which is not implemented
        // in HAMMER
    });
  } else if (run == "run2") {
    ham.setFFInputScheme({{"BD", "CLN_1"},
                          {"BD*", "CLN_1"},
                          {"BD**0*", "ISGW2"},
                          {"BD**1", "ISGW2"},
                          {"BD**1*", "ISGW2"},
                          {"BD**2*", "ISGW2"}});

    // 12573001, 12573012
    // clang-format off
    // HQET2(hqetrho2, hqetv1_1, indelta): 1.131 1.035 0.38
    ham.setOptions("BtoDCLN_1: {RhoSq: 1.131, Delta: 0.38, G1: 1.035}");  // HQET3
    // 11574011, 11574021
    // HQET2(hqetrho2, hqetha1_1, hqetr1_1, hqetr2_1, hqetr0_1): 1.122 0.908 1.270 0.852 1.15
    ham.setOptions("BtoD*CLN_1: {RhoSq: 1.122, F1: 0.908, R1: 1.270, R2: 0.852, R0: 1.15}");  // HQET2 11874440, rest of D**
    // ISGW2, which has no configurable parameter
    // clang-format on
  }
}

map<string, complex<double>> specializedWC = {
    {"SM", 1},     {"S_qLlL", 0}, {"S_qRlL", 0}, {"V_qLlL", 0},
    {"V_qRlL", 0}, {"T_qLlL", 0}, {"S_qLlR", 0}, {"S_qRlR", 0},
    {"V_qLlR", 0}, {"V_qRlR", 0}, {"T_qRlR", 0}};

// clang-format off
// +, -, +, -, ...
vector<vector<string>> BtoDVars = {
  // shifting 1-th param in + direction....
  {
    "{ap: [0.015642660612597052, -0.034035685234997386, -0.0898070271651518, 0.0]}",
    "{a0: [0.0858147593040351, -0.3050702334404309, -0.22752112980378628, 0.0]}"
  },
  //shifting 1-th param in - direction....
  {
    "{ap: [0.01567733938740295, -0.03436431476500262, -0.09019297283484819, 0.0]}",
    "{a0: [0.07288524069596491, -0.10492976655956911, -0.23247887019621374, 0.0]}"
  },
  // shifting 2-th param in + direction....
  {
    "{ap: [0.015536159303514138, -0.03274798737619622, -0.08610235496526895, 0.0]}",
    "{a0: [0.07939786111920759, -0.2055198293137311, -0.2513856759832376, 0.0]}"
  },
  // shifting 2-th param in - direction....
  {
    "{ap: [0.015783840696485862, -0.03565201262380378, -0.09389764503473104, 0.0]}",
    "{a0: [0.07930213888079242, -0.20448017068626886, -0.20861432401676247, 0.0]}"
  },
  // shifting 3-th param in + direction....
  {
    "{ap: [0.015185466182351984, -0.04806686174041869, -0.08969198838997502, 0.0]}",
    "{a0: [0.07255753045478437, -0.20504393049104058, -0.23088155901344717, 0.0]}"
  },
  // shifting 3-th param in - direction....
  {
    "{ap: [0.016134533817648015, -0.020333138259581315, -0.09030801161002497, 0.0]}",
    "{a0: [0.08614246954521564, -0.2049560695089594, -0.22911844098655285, 0.0]}"
  },
  // shifting 4-th param in + direction....
  {
    "{ap: [0.01566286514976162, -0.034199101149958354, -0.08990582046864679, 0.0]}",
    "{a0: [0.07936645416843903, -0.20499939148119928, -0.22998280566971113, 0.0]}"
  },
  // shifting 4-th param in - direction....
  {
    "{ap: [0.01565713485023838, -0.03420089885004165, -0.09009417953135321, 0.0]}",
    "{a0: [0.07933354583156098, -0.2050006085188007, -0.2300171943302889, 0.0]}"
  },
  //shifting 5-th param in + direction....
  {
    "{ap: [0.014646019937572174, -0.034165525575076655, -0.08997191655772308, 0.0]}",
    "{a0: [0.07430176773341135, -0.2049993836925806, -0.22998668416129328, 0.0]}"
  },
  // shifting 5-th param in - direction....
  {
    "{ap: [0.016673980062427828, -0.03423447442492335, -0.09002808344227692, 0.0]}",
    "{a0: [0.08439823226658866, -0.20500061630741936, -0.23001331583870674, 0.0]}"
  }
};

// Various FF config
map<string, vector<vector<string>>> ffVarSpecs = {
  {"BD", BtoDVars}
};

map<string, string> ffSchemeByDecay = {
  {"BD", "BGL"},
  {"BD*", "BGL"},
  {"BD**0", "BLR"},
  {"BD**1", "BLR"},
  {"BD**1*", "BLR"},
  {"BD**2", "BLR"},
};

void setBtoDBGLDefault(Hammer::Hammer& ham, const string scheme) {
  cout << "Configuring default BtoDBGL for scheme " << scheme << endl;
  ham.setOptions(scheme + ": {ChiT: 0.0005131}");
  ham.setOptions(scheme + ": {ChiL: 0.006332}");
  ham.setOptions(scheme + ": {BcStatesp: [6.329, 6.92, 7.02]}");
  ham.setOptions(scheme + ": {BcStates0: [6.716, 7.121]}");
  ham.setOptions(scheme + ": {ap: [0.01566, -0.0342, -0.09, 0.0]}");
  ham.setOptions(scheme + ": {a0: [0.07935, -0.205, -0.23, 0.0]}");
}

void setBtoDstarBGLDefault(Hammer::Hammer& ham, const string scheme) {
  // placeholder
}

void setGenericFFDefault(Hammer::Hammer& ham, const string scheme) {};

map<string, function<void(Hammer::Hammer&, const string)>>
ffSchemeDefaultsByDecay = {
  {"BD", setBtoDBGLDefault},
  {"BD*", setBtoDstarBGLDefault},
  {"BD**0", setGenericFFDefault},
  {"BD**1", setGenericFFDefault},
  {"BD**1*", setGenericFFDefault},
  {"BD**2", setGenericFFDefault}
};
// clang-format on

const int numOfFFVar = 20;

string decayDescr(const string decay) {
  auto daughter = decay.substr(1);
  return "Bto" + daughter;
}

vector<string> setOutputFF(Hammer::Hammer& ham) {
  vector<string> hamFFSchemes{"OutputFF"};

  ham.addFFScheme("OutputFF", ffSchemeByDecay);
  // Set defaults
  for (auto [decay, ffName] : ffSchemeByDecay) {
    auto fullDescr = decayDescr(decay) + ffName;
    cout << "Decay: " << decay << "; default FF: " << fullDescr << endl;
    ffSchemeDefaultsByDecay[decay](ham, fullDescr);
  }

  // Set variations
  for (int i = 1; i <= numOfFFVar; i++) {
    cout << "Configuring FF scheme: OutputFFVar" + to_string(i) << endl;
    map<string, string> schemes{};
    for (auto const& [decay, vars] : ffVarSpecs) {
      if (i <= vars.size()) {
        // Need to reweight the decay in this HAMMER scheme
        auto ffName    = ffSchemeByDecay[decay] + "_" + to_string(i);
        auto descr     = decayDescr(decay);
        schemes[decay] = ffName;
        // Configure the FF scheme defaults for this decay
        cout << "  Variation for decay: " << decay
             << "; with FF: " << descr + ffName << endl;
        ffSchemeDefaultsByDecay[decay](ham, descr + ffName);
        for (auto const& shift : vars[i - 1])
          ham.setOptions(descr + ffName + ": " +
                         shift);  // configure FF variations
      }
    }
    ham.addFFScheme("OutputFFVar" + to_string(i), schemes);
  }

  return hamFFSchemes;
}

void setDecays(Hammer::Hammer& ham) {
  ham.includeDecay("BDTauNu");
  ham.includeDecay("BDMuNu");

  ham.includeDecay("BD*TauNu");
  ham.includeDecay("BD*MuNu");

  ham.includeDecay("BD**0*TauNu");
  ham.includeDecay("BD**1TauNu");
  ham.includeDecay("BD**1*TauNu");
  ham.includeDecay("BD**2*TauNu");

  ham.includeDecay("BD**0*MuNu");
  ham.includeDecay("BD**1MuNu");
  ham.includeDecay("BD**1*MuNu");
  ham.includeDecay("BD**2*MuNu");
}

/////////////
// Filters //
/////////////

bool truthMatchOk(double q2True, bool isTauDecay, int bMesonId, int dau1Id,
                  int dau2Id, int muID) {
  double q2Min = 100 * 100;
  if (isTauDecay) q2Min = 1700 * 1700;

  // we require that there's ONE and ONLY ONE D meson
  return findIn(LEGAL_B_MESON_IDS, TMath::Abs(bMesonId)) && q2True > q2Min &&
         isDMeson(TMath::Abs(dau1Id)) && !isDMeson(TMath::Abs(dau2Id)) &&
         TMath::Abs(muID) == 13;
}

/////////////////////////
// Branch name helpers //
/////////////////////////

vector<string> getBTrueP(string particle, string bId) {
  return setBrPrefix(particle, {"TRUEP_E", "TRUEP_X", "TRUEP_Y", "TRUEP_Z"},
                     {bId});
}

vector<string> getDauTrueP(string particle, string dau) {
  return setBrPrefix(particle + "_TrueHadron_" + dau,
                     {"PE", "PX", "PY", "PZ", "ID"});
}

////////////////////
// FSR correction //
////////////////////

typedef tuple<double, double, double, double, int> HamPartCtn;

bool isSoftPhoton(Hammer::Particle part) {
  if (part.p().E() < SOFT_PHOTON_THRESH) return true;
  return false;
}

bool isSoftPhoton(HamPartCtn part) {
  auto pe = get<0>(part);
  if (pe < SOFT_PHOTON_THRESH) return true;
  return false;
}

vector<HamPartCtn> buildPhotonVec(RVec<float>& arrPe, RVec<float>& arrPx,
                                  RVec<float>& arrPy, RVec<float>& arrPz,
                                  RVec<float>& arrId, int size) {
  vector<HamPartCtn> result{};
  for (auto idx = 0; idx != size; idx++) {
    auto part = buildPartVec(arrPe[idx], arrPx[idx], arrPy[idx], arrPz[idx],
                             static_cast<int>(arrId[idx]));
    if (!isSoftPhoton(part)) result.emplace_back(part);
  }
  return result;
}

string printP(double pe, double px, double py, double pz) {
  stringstream buffer;
  buffer << pe << "," << px << "," << py << "," << pz;
  return buffer.str();
}

string printP(Hammer::Particle part) {
  stringstream buffer;
  auto         p = part.p();
  buffer << p.E() << "," << p.px() << "," << p.py() << "," << p.pz();
  return buffer.str();
}

string addRadiativePhotons(Hammer::Process& proc, Hammer::ParticleIndices& idx,
                           int refMomId, vector<HamPartCtn>& photons) {
  stringstream buffer;

  for (const auto& p : photons) {
    auto [pe, px, py, pz, photonMomId] = p;
    if (photonMomId == refMomId || photonMomId == -refMomId) {
      auto partPhoton = buildHamPart(pe, px, py, pz, 22);
      idx.emplace_back(proc.addParticle(partPhoton));
      buffer << "  Adding photon: " << printP(pe, px, py, pz) << " to "
             << refMomId << endl;
    }
  }

  return buffer.str();
}

/////////////////
// Reweighting //
/////////////////

pair<RNode, vector<string>> prepAuxOutput(RNode df, string bMesonName) {
  auto outputBrs = vector<string>{};

  // truth variables
  df = df.Define("q2_true", bMesonName + "_True_Q2 / 1000 / 1000");
  outputBrs.emplace_back("q2_true");

  df = df.Define("is_tau", bMesonName + "_True_IsTauDecay");
  outputBrs.emplace_back("is_tau");

  auto kinematicSuffix = {"_PE", "_PX", "_PY", "_PZ"};
  for (int i = 0; i < 2; i++) {
    auto partName = "d_meson" + to_string(i + 1);

    df = df.Define(partName + "_true_id",
                   bMesonName + "_TrueHadron_D" + to_string(i) + "_ID");
    outputBrs.emplace_back(partName + "_true_id");

    auto sIdx         = to_string(i);
    auto kinematicBrs = vector<string>{};
    for (const auto& suf : kinematicSuffix) {
      kinematicBrs.emplace_back("TrueHadron_D" + sIdx + suf);
    }

    df = df.Define(partName + "_true_m", invM,
                   setBrPrefix(bMesonName, kinematicBrs));
    outputBrs.emplace_back(partName + "_true_m");
  }

  // truth-matching ok
  df = df.Define("ham_tm_ok", truthMatchOk,
                 setBrPrefix(bMesonName,
                             {"True_Q2", "True_IsTauDecay", "TRUEID",
                              "TrueHadron_D0_ID", "TrueHadron_D1_ID"},
                             {"mu_TRUEID"}));
  outputBrs.emplace_back("ham_tm_ok");

  return {df, outputBrs};
}

// NOTE: Refine branch names if input ntuples have different tree structrue!
pair<RNode, vector<string>> prepHamInput(RNode df, string bMesonName) {
  auto outputBrs = vector<string>{};

  // B meson
  df = df.Define("b_id_fixed", bIdFix,
                 setBrPrefix(bMesonName, {"TRUEID", "TrueHadron_D0_ID"}));
  df = df.Define("part_B", buildPartVec, getBTrueP(bMesonName, "b_id_fixed"));

  // D meson (can be a D**, D*, or D)
  df = df.Define("part_D", buildPartVec, getDauTrueP(bMesonName, "D0"));

  // D meson daughters
  df =
      df.Define("part_D_dau0", buildPartVec, getDauTrueP(bMesonName, "D0_GD0"));
  df =
      df.Define("part_D_dau1", buildPartVec, getDauTrueP(bMesonName, "D0_GD1"));
  df =
      df.Define("part_D_dau2", buildPartVec, getDauTrueP(bMesonName, "D0_GD2"));

#ifdef RADIATIVE_CORRECTION
  // the radiative photons container
  df = df.Define("part_photon_arr", buildPhotonVec,
                 setBrPrefix(bMesonName, {"MCTrue_gamma_E", "MCTrue_gamma_PX",
                                          "MCTrue_gamma_PY", "MCTrue_gamma_PZ",
                                          "MCTrue_gamma_mother_ID",
                                          "MCTrue_gamma_ArrayLength"}));

  // Tau/Mu, Nu_Tau/Nu_Mu associated w/ B -> D decay
  df = df.Define("part_Tau_id", tauIdFix, {"mu_TRUEID"});
  df = df.Define(
      "part_Tau", buildPartVec,
      setBrPrefix(bMesonName,
                  {"TrueTau_PE", "TrueTau_PX", "TrueTau_PY", "TrueTau_PZ"},
                  {"part_Tau_id"}));
  df = df.Define("part_Mu_id", muIdFix, {"mu_TRUEID"});
  df = df.Define(
      "part_Mu", buildPartVec,
      setBrPrefix(bMesonName,
                  {"TrueMu_PE", "TrueMu_PX", "TrueMu_PY", "TrueMu_PZ"},
                  {"part_Mu_id"}));
  df = df.Define("part_L",
                 [](HamPartCtn pTau, HamPartCtn pMu, bool isTau) {
                   if (isTau) return pTau;
                   return pMu;
                 },
                 {"part_Tau", "part_Mu", bMesonName + "_True_IsTauDecay"});

  df = df.Define("part_NuL_id", nuIdFix,
                 {"mu_TRUEID", bMesonName + "_True_IsTauDecay"});
  df = df.Define("part_NuL", buildPartVec,
                 setBrPrefix(bMesonName,
                             {"TrueNeutrino_PE", "TrueNeutrino_PX",
                              "TrueNeutrino_PY", "TrueNeutrino_PZ"},
                             {"part_NuL_id"}));
  // for leptons produced by primary B -> D decays, ALWAYS use part_L and
  // part_NuL

  // TauNu, MuNu associated w/ Tau -> Mu decay
  df = df.Define("part_NuTau_id", tauNuTauIdFix, {"mu_TRUEID"});
  df = df.Define("part_NuTau", buildPartVec,
                 setBrPrefix(bMesonName,
                             {"TrueTauNuTau_PE", "TrueTauNuTau_PX",
                              "TrueTauNuTau_PY", "TrueTauNuTau_PZ"},
                             {"part_NuTau_id"}));
  df = df.Define("part_NuMu_id", tauNuMuIdFix, {"mu_TRUEID"});
  df = df.Define("part_NuMu", buildPartVec,
                 setBrPrefix(bMesonName,
                             {"TrueTauNuMu_PE", "TrueTauNuMu_PX",
                              "TrueTauNuMu_PY", "TrueTauNuMu_PZ"},
                             {"part_NuMu_id"}));
#else
  // this is just a place holder in case we don't add radiative photons back!
  df = df.Define("part_photon_arr", []() { return vector<HamPartCtn>{}; }, {});
#endif

  return {df, outputBrs};
}

auto reweightWrapper(Hammer::Hammer& ham, unsigned long& numOfEvt,
                     unsigned long& numOfEvtOk, vector<string>& schemes) {
  return [&](bool truthMatchOk, bool isTau, HamPartCtn pB, HamPartCtn pD,
             HamPartCtn pDDau0, HamPartCtn pDDau1, HamPartCtn pDDau2,
             HamPartCtn pL, HamPartCtn pNuL, HamPartCtn pMu, HamPartCtn pNuMu,
             HamPartCtn pNuTau, vector<HamPartCtn> pPhotons) {
    bool                      hamOk = true;
    double                    wtFF  = 1.0;
    array<double, numOfFFVar> wtFFVars;
    fill_n(wtFFVars.begin(), numOfFFVar, 1.0);  // all default to 1.0

    string debugMsg = "====\n";
    numOfEvt += 1;
    if (!truthMatchOk)
      // FIXME: it's ugly! but not sure RDataFrame frame can deduce type if we
      // use anything fancy.
      return tuple<bool, double, double, double, double, double, double, double,
                   double, double, double, double, double, double, double,
                   double, double, double, double, double, double, double>{
          false,        wtFF,         wtFFVars[0],  wtFFVars[1],  wtFFVars[2],
          wtFFVars[3],  wtFFVars[4],  wtFFVars[5],  wtFFVars[6],  wtFFVars[7],
          wtFFVars[8],  wtFFVars[9],  wtFFVars[10], wtFFVars[11], wtFFVars[12],
          wtFFVars[13], wtFFVars[14], wtFFVars[15], wtFFVars[16], wtFFVars[17],
          wtFFVars[18], wtFFVars[19]};

    Hammer::Process proc;
    auto            partB   = buildHamPart(pB);
    auto            partD   = buildHamPart(pD);
    auto            partL   = buildHamPart(pL);
    auto            partNuL = buildHamPart(pNuL);
    auto particles = vector<Hammer::Particle>{partB, partD, partL, partNuL};

    // add B meson
    auto partBIdx = proc.addParticle(partB);

    debugMsg += "  B meson 4-mom: " + printP(partB) + '\n';
    debugMsg += "  D meson 4-mom: " + printP(partD) + '\n';
    debugMsg += "  primary charged lepton 4-mom: " + printP(partL) + '\n';
    debugMsg += "  primary neutrino 4-mom: " + printP(partNuL) + '\n';

    // add direct B daughters
    auto                    partDIdx   = proc.addParticle(partD);
    auto                    partLIdx   = proc.addParticle(partL);
    auto                    partNuLIdx = proc.addParticle(partNuL);
    Hammer::ParticleIndices partBDauIdx{partDIdx, partLIdx, partNuLIdx};
#ifdef RADIATIVE_CORRECTION
    debugMsg += addRadiativePhotons(proc, partBDauIdx, partB.pdgId(), pPhotons);
#endif
    proc.addVertex(partBIdx, partBDauIdx);

    // in case of a D*, add its daughters as well
    Hammer::ParticleIndices partDDauIdx{};
#ifdef RADIATIVE_CORRECTION
    debugMsg += addRadiativePhotons(proc, partDDauIdx, partD.pdgId(), pPhotons);
#endif

    bool isDst = isDstMeson(get<4>(pD));
    if (isDst) {
      auto partDDaus = {buildHamPart(pDDau0), buildHamPart(pDDau1),
                        buildHamPart(pDDau2)};
      for (const auto& p : partDDaus) {
        // don't add placeholder particle
        // don't add photons as it will be a duplicate
        if (p.pdgId() == 0 || p.pdgId() == 22) continue;
        partDDauIdx.emplace_back(proc.addParticle(p));
        particles.emplace_back(p);
        debugMsg += "  D daughters: " + printP(p) + "\n";
      }
    }
    if (partDDauIdx.size()) proc.addVertex(partDIdx, partDDauIdx);

    // in case of a Tau, add its daughters
    Hammer::ParticleIndices partLDauIdx{};
#ifdef RADIATIVE_CORRECTION
    debugMsg += addRadiativePhotons(proc, partLDauIdx, partL.pdgId(), pPhotons);
#endif

    if (isTau) {
      auto partLDaus = {buildHamPart(pMu), buildHamPart(pNuMu),
                        buildHamPart(pNuTau)};
      for (const auto& p : partLDaus) {
        partLDauIdx.emplace_back(proc.addParticle(p));
        particles.emplace_back(p);
        debugMsg += "  secondary leptons: " + printP(p) + "\n";
      }
    }
    if (partLDauIdx.size()) proc.addVertex(partLIdx, partLDauIdx);

    // make sure invariant mass is not negative
    auto partInvMOk = vector<bool>{};
    for (const auto& p : particles) {
      if (p.p().mass() >= 0)
        partInvMOk.emplace_back(true);
      else
        partInvMOk.emplace_back(false);
    }
    bool allPartsOk =
        find(partInvMOk.begin(), partInvMOk.end(), false) == partInvMOk.end();
    if (!allPartsOk) {
      cout << "  WARN: Bad kinematics for candidate: " << numOfEvt << endl;
      hamOk = false;
    }

#ifdef DEBUG_CLI
    cout << debugMsg;
#endif

    // add the whole decay chain to hammer and see if it likes it
    if (hamOk) {
      ham.initEvent();
      int procId = 0;
      try {
        procId = ham.addProcess(proc);
      } catch (const exception& e) {
        cout << "  WARN: HAMMER doesn't add process properly: " << numOfEvt
             << endl;
        cout << e.what() << endl;
        hamOk = false;
      }
      if (procId == 0) hamOk = false;
    }

    // compute FF weight
    auto nominalFFScheme = schemes[0];
    if (hamOk) {
      try {
        ham.processEvent();
        wtFF = ham.getWeight(nominalFFScheme);
      } catch (const exception& e) {
        cout << "  WARN: HAMMER doesn't like candidate for reweighting: "
             << numOfEvt << endl;
        cout << e.what() << endl;
        hamOk = false;
      }

      if (isnan(wtFF) || isinf(wtFF)) hamOk = false;

      if (hamOk) {
        numOfEvtOk += 1;
        // compute various FF variation weights here
        // failure's allowed because not all FF schemes require 20 variations!
        try {
          for (int i = 0; i < numOfFFVar; i++) {
            auto ffName = "OutputFFVar" + to_string(i + 1);
            wtFFVars[i] = ham.getWeight(ffName);
          }
        } catch (const exception& e) {
        }
      }
    }

#ifdef DEBUG_CLI
    if (hamOk) cout << "  FF weight: " << wtFF << endl;
#endif

    // FIXME: it's ugly! but not sure RDataFrame frame can deduce type if we
    // use anything fancy.
    return tuple<bool, double, double, double, double, double, double, double,
                 double, double, double, double, double, double, double, double,
                 double, double, double, double, double, double>{
        hamOk,        wtFF,         wtFFVars[0],  wtFFVars[1],  wtFFVars[2],
        wtFFVars[3],  wtFFVars[4],  wtFFVars[5],  wtFFVars[6],  wtFFVars[7],
        wtFFVars[8],  wtFFVars[9],  wtFFVars[10], wtFFVars[11], wtFFVars[12],
        wtFFVars[13], wtFFVars[14], wtFFVars[15], wtFFVars[16], wtFFVars[17],
        wtFFVars[18], wtFFVars[19]};
  };
}

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("ReweightRDX", "Reweight RDX FF w/ HAMMER.");

  // clang-format off
  argOpts.add_options()
    // positional
    ("ntpIn", "specify input ntuple.", cxxopts::value<string>())
    ("ntpOut", "specify output ntuple.", cxxopts::value<string>())
    ("extra", "unused.", cxxopts::value<vector<string>>())
    // keyword
    ("h,help", "print help.")
    ("t,trees", "specify tree name.",
     cxxopts::value<vector<string>>()
     ->default_value("TupleBminus/DecayTree,TupleB0/DecayTree"))
    ("b,bMesons", "specify B meson name.",
     cxxopts::value<vector<string>>()->default_value("b,b0"))
    ("r,run", "specify run.", cxxopts::value<string>()->default_value("run2"))
  ;
  // setup positional argument
  argOpts.parse_positional({"ntpIn", "ntpOut", "extra"});
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  auto ntpIn   = parsedArgs["ntpIn"].as<string>();
  auto ntpOut  = parsedArgs["ntpOut"].as<string>();
  auto trees   = parsedArgs["trees"].as<vector<string>>();
  auto bMesons = parsedArgs["bMesons"].as<vector<string>>();
  auto run     = parsedArgs["run"].as<string>();

  Hammer::Hammer ham{};

  setDecays(ham);
  setInputFF(ham, run);
  auto ffSchemes = setOutputFF(ham);

  ham.setUnits("MeV");
  ham.setOptions("ProcessCalc: {CheckForNaNs: true}");
  ham.initRun();

  // only use SM Wilson coefficients
  ham.specializeWCInWeights("BtoCTauNu", specializedWC);
  ham.specializeWCInWeights("BtoCMuNu", specializedWC);

  // output option
  auto writeOpts  = ROOT::RDF::RSnapshotOptions{};
  writeOpts.fMode = "UPDATE";

  for (int idx = 0; idx != trees.size(); idx++) {
    RNode          df     = static_cast<RNode>(RDataFrame(trees[idx], ntpIn));
    auto           bMeson = bMesons[idx];
    vector<string> outputBrs{"runNumber", "eventNumber"};
    unsigned long  numOfEvt   = 0;
    unsigned long  numOfEvtOk = 0;

    cout << "Handling " << trees[idx] << " with B meson name " << bMeson
         << endl;

    // prepare aux output branches like q2_true
    auto [dfAux, outputBrsAux] = prepAuxOutput(df, bMeson);
    df                         = dfAux;
    for (const auto& br : outputBrsAux) outputBrs.emplace_back(br);

    // prepare HAMMER particles
    tie(df, ignore) = prepHamInput(df, bMeson);

    // reweight FF
    auto reweight = reweightWrapper(ham, numOfEvt, numOfEvtOk, ffSchemes);
    df            = df.Define("ff_result", reweight,
                   {"ham_tm_ok", "is_tau", "part_B", "part_D", "part_D_dau0",
                    "part_D_dau1", "part_D_dau2", "part_L", "part_NuL",
                    "part_Mu", "part_NuMu", "part_NuTau", "part_photon_arr"});
    df            = df.Define("ham_ok", "get<0>(ff_result)");
    df            = df.Define("wff", "get<1>(ff_result)");
    for (int i = 0; i < numOfFFVar; i++) {
      auto outputBrName = "wff_var" + to_string(i + 1);
      df = df.Define(outputBrName, "get<" + to_string(i + 2) + ">(ff_result)");
      outputBrs.emplace_back(outputBrName);
    }
    outputBrs.emplace_back("ham_ok");
    outputBrs.emplace_back("wff");

    df.Snapshot(trees[idx], ntpOut, outputBrs, writeOpts);

    cout << "Total number of candidates: " << numOfEvt << endl;
    cout << "Hammer reweighted candidates: " << numOfEvtOk << endl;
    cout << "Reweighted fraction: "
         << static_cast<float>(numOfEvtOk) / static_cast<float>(numOfEvt)
         << endl;
  }
}
