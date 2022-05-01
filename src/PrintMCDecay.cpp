// Author: Yipeng Sun
// Last Change: Sat Apr 30, 2022 at 09:05 PM -0400

#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <ROOT/RDataFrame.hxx>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <cxxopts.hpp>

#include "const.h"
#include "utils_general.h"

using namespace std;

using ROOT::RDataFrame;

///////////////////
// Configurables //
///////////////////

// clang-format off
const auto DECAY_NAMES = vector<string_view>{
  "B meson ID: ",
  "First daughter ID: ",
  "  First G-daughter ID: ",
  "  Second G-daughter ID: ",
  "  Third G-daughter ID: ",
  "Second daughter meson ID: ",
  "  First G-daughter ID: ",
  "  Second G-daughter ID: ",
  "  Third G-daughter ID: ",
  "Third daughter meson ID: ",
  "  First G-daughter ID: ",
  "  Second G-daughter ID: ",
  "  Third G-daughter ID: "
};
// clang-format on

const auto BRANCH_ALIAES = vector<pair<string, string>>{
    {"q2_true", "True_Q2"},
    {"is_tau", "True_IsTauDecay"},
    {"b_id", "TRUEID"},
    {"dau0_id", "TrueHadron_D0_ID"},
    {"dau1_id", "TrueHadron_D1_ID"},
    {"dau2_id", "TrueHadron_D2_ID"},
    {"dau0_gd0_id", "TrueHadron_D0_GD0_ID"},
    {"dau0_gd1_id", "TrueHadron_D0_GD1_ID"},
    {"dau0_gd2_id", "TrueHadron_D0_GD2_ID"},
    {"dau1_gd0_id", "TrueHadron_D1_GD0_ID"},
    {"dau1_gd1_id", "TrueHadron_D1_GD1_ID"},
    {"dau1_gd2_id", "TrueHadron_D1_GD2_ID"},
    {"dau2_gd0_id", "TrueHadron_D2_GD0_ID"},
    {"dau2_gd1_id", "TrueHadron_D2_GD1_ID"},
    {"dau2_gd2_id", "TrueHadron_D2_GD2_ID"},
};

const auto DECAY_SIGNATURE = vector<string>{
    "is_tau",      "b_id",        "dau0_id",     "dau0_gd0_id", "dau0_gd1_id",
    "dau0_gd2_id", "dau1_id",     "dau1_gd0_id", "dau1_gd1_id", "dau1_gd2_id",
    "dau2_id",     "dau2_gd0_id", "dau2_gd1_id", "dau2_gd2_id"};

/////////////
// Helpers //
/////////////

typedef map<vector<int>, unsigned long> DecayFreq;

template <typename A, typename B>
pair<B, A> flipPair(const pair<A, B>& p) {
  return pair<B, A>(p.second, p.first);
}

template <typename A, typename B>
multimap<B, A> flipMap(const map<A, B>& src) {
  multimap<B, A> dst;
  transform(src.begin(), src.end(), inserter(dst, dst.begin()), flipPair<A, B>);
  return dst;
}

/////////////
// Filters //
/////////////

bool truthMatchOk(double q2True, bool isTauDecay, int bMesonId, int dMesonId) {
  double q2Min = 100 * 100;
  if (isTauDecay) q2Min = 1700 * 1700;

  return findIn(LEGAL_B_MESON_IDS, TMath::Abs(bMesonId)) && q2True > q2Min &&
         isDMeson(TMath::Abs(dMesonId));  // NOTE: Taking abs is crucial!
}

//////////////
// Printers //
//////////////

void countDecayFreq(DecayFreq& freq, unsigned long& numOfEvt,
                    unsigned long& numOfEvtWithB, bool truthMatch,
                    vector<int> truthSignature) {
  numOfEvt += 1;

  if (truthMatch) {
    numOfEvtWithB += 1;
    vector<int> key = {};
    for (auto v : truthSignature) key.emplace_back(TMath::Abs(v));

    if (freq.find(key) == freq.end())
      freq[key] = 1l;
    else
      freq[key] += 1;
  }
}

void printDecayFreq(DecayFreq freq, TDatabasePDG* db) {
  auto sorted = flipMap(freq);
  for (auto const& [val, key] : boost::adaptors::reverse(sorted)) {
    cout << "======" << endl;
    cout << "The following decay has " << val << " candidates." << endl;
    cout << "Is Tau decay: " << key[0] << endl;
    for (auto idx = 1; idx < key.size(); idx++) {
      if (key[idx]) {
        cout << DECAY_NAMES[idx - 1] << getParticleName(key[idx], db, true)
             << endl;
      }
    }
  }
}

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("PrintMCDecay", "print decays of valid B mesons.");

  // clang-format off
  argOpts.add_options()
    // positional
    ("ntp", "specify input ntuple.", cxxopts::value<string>())
    ("extra", "unused.", cxxopts::value<vector<string>>())
    // keyword
    ("h,help", "print help")
    ("t,tree", "specify tree name.",
     cxxopts::value<string>()->default_value("TupleBminus/DecayTree"))
    ("p,particle", "specify B meson name.",
     cxxopts::value<string>()->default_value("b"))
  ;
  // setup positional argument
  argOpts.parse_positional({"ntp", "extra"});
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  auto tree   = parsedArgs["tree"].as<string>();
  auto ntp    = parsedArgs["ntp"].as<string>();
  auto bMeson = parsedArgs["particle"].as<string>();

  unsigned long numOfEvt      = 0;
  unsigned long numOfEvtWithB = 0;
  auto          freq          = DecayFreq{};
  auto          db            = make_unique<TDatabasePDG>();
  auto          dfInit        = RDataFrame(tree, ntp);

  // functions to be JIT'ed
  gInterpreter->Declare(
      "auto makeVecInt = [](auto...args) { return vector<int>{args...}; };");

  auto df = defineBranch(dfInit, BRANCH_ALIAES, bMeson);
  df      = df.Define("truthmatch", truthMatchOk,
                 {"q2_true", "is_tau", "b_id", "dau0_id"})
           .Define("signature",
                   "makeVecInt(" +
                       boost::algorithm::join(DECAY_SIGNATURE, ",") + ")");
  df.Foreach(
      [&](bool truthMatch, vector<int> truthSignature) {
        countDecayFreq(freq, numOfEvt, numOfEvtWithB, truthMatch,
                       truthSignature);
      },
      {"truthmatch", "signature"});

  printDecayFreq(freq, db.get());
  cout << endl;
  cout << "Total number of candidates: " << numOfEvt << endl;
  cout << "Truth-matched candidates: " << numOfEvtWithB << endl;
  cout << "Truth-matched fraction: "
       << static_cast<float>(numOfEvtWithB) / static_cast<float>(numOfEvt)
       << endl;
}
