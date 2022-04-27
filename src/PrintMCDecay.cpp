// Author: Yipeng Sun
// Last Change: Wed Apr 27, 2022 at 01:19 AM -0400

#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <boost/range/adaptor/reversed.hpp>
#include <cxxopts.hpp>

#include "utils_general.h"

using namespace std;

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

const auto LEGAL_B_MESON_IDS = vector<int>{511, 521};

const auto BRANCH_ALIAES = vector<pair<string, string>>{
    {"q2_true", "_TRUE_Q2"},
    {"is_tau", "_True_IsTauDecay"},
    {"b_id", "_TRUEID"},
    {"dau0_id", "_TrueHadron_D0_ID"},
    {"dau1_id", "_TrueHadron_D1_ID"},
    {"dau2_id", "_TrueHadron_D2_ID"},
    {"dau0_gd0_id", "_TrueHadron_D0_GD0_ID"},
    {"dau0_gd1_id", "_TrueHadron_D0_GD1_ID"},
    {"dau0_gd2_id", "_TrueHadron_D0_GD2_ID"},
    {"dau1_gd0_id", "_TrueHadron_D1_GD0_ID"},
    {"dau1_gd1_id", "_TrueHadron_D1_GD1_ID"},
    {"dau1_gd2_id", "_TrueHadron_D1_GD2_ID"},
    {"dau2_gd0_id", "_TrueHadron_D2_GD0_ID"},
    {"dau2_gd1_id", "_TrueHadron_D2_GD1_ID"},
    {"dau2_gd2_id", "_TrueHadron_D2_GD2_ID"},
};

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
         isDMeson(dMesonId);
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
    if (freq.find(truthSignature) == freq.end())
      freq[truthSignature] = 1l;
    else
      freq[truthSignature] += 1;
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

  cout << parsedArgs["ntp"].as<string>();

  unsigned long numOfEvt      = 0;
  unsigned long numOfEvtWithB = 0;

  // auto freq = print_id(ntp, tree_name);
  // print_decay_freq(freq, db);

  // cout << "Total number of candidates: " << num_of_evt << endl;
  // cout << "Truth-matched candidates: " << num_of_evt_w_b_meson << endl;
  // cout << "Truth-matched fraction: "
  //      << static_cast<float>(num_of_evt_w_b_meson) /
  //             static_cast<float>(num_of_evt)
  //      << endl;

  // delete ntp;
  // delete db;
}
