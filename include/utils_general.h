// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sat Jun 11, 2022 at 03:11 PM -0400

#pragma once

#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <TDataType.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>

using ROOT::RDF::RNode;
using std::map;
using std::pair;
using std::string;
using std::vector;

/////////////////////
// General helpers //
/////////////////////

template <template <typename, typename> class Iterable, typename T,
          typename Allocator>
bool findIn(Iterable<T, Allocator> iter, T elem) {
  if (find(iter.begin(), iter.end(), elem) != iter.end()) return true;
  return false;
}

vector<string> split(const string& s, char delim) {
  std::stringstream ss(s);
  string            item;
  vector<string>    elems;

  while (getline(ss, item, delim)) {
    elems.push_back(move(item));
  }

  return elems;
}

TString dirname(string s) {
  auto splitted = split(s, '/');
  if (splitted.size() == 1) return TString("");

  TString dir = "";
  for (int idx = 0; idx < splitted.size() - 1; idx++) {
    dir += splitted[idx];
    if (idx < splitted.size() - 2) dir += "/";
  }

  return dir;
}

TString basename(string s) { return TString(split(s, '/').back()); }

int digitIs(Int_t num, int digit, int base = 10) {
  int fac = num / TMath::Power(base, digit - 1);
  return fac % 10;
}

string getParticleName(int id, TDatabasePDG* db, bool useAbsId = false) {
  if (!id) return string("None");

  auto absId = id;
  if (useAbsId) absId = TMath::Abs(id);

  char buf[50];
  sprintf(buf, " (%d)", absId);
  auto strId    = string(buf);
  auto particle = db->GetParticle(absId);

  if (particle != nullptr) return string(particle->GetName()) + buf;
  return string("Unknown") + buf;
}

template <typename T, size_t... Indices>
auto createTupleHelper(const vector<T>& vec, std::index_sequence<Indices...>) {
  return std::make_tuple(vec[Indices]...);
}

template <size_t N, typename T>
auto createTuple(const vector<T>& vec) {
  assert(vec.size() >= N);
  return createTupleHelper(vec, std::make_index_sequence<N>{});
}

////////////////////////
// RDataframe helpers //
////////////////////////

RNode defineBranch(RNode df, const vector<pair<string, string>>& rules,
                   const string particle, int idx = 0) {
  if (rules.size() == idx) return df;

  auto inputBrName = rules[idx].second;
  if (particle != "") inputBrName = particle + "_" + inputBrName;

  return defineBranch(df.Define(rules[idx].first, inputBrName), rules, particle,
                      idx + 1);
}

vector<string> setBrPrefix(const string prefix, const vector<string>& vars,
                           const vector<string>& varsAppend = {}) {
  vector<string> result{};
  for (const auto& v : vars) result.emplace_back(prefix + "_" + v);
  for (const auto& v : varsAppend) result.emplace_back(v);
  return result;
}

////////////////////////
// Kinematics helpers //
////////////////////////

double invM(double pe, double px, double py, double pz) {
  return TMath::Sqrt(pe * pe - px * px - py * py - pz * pz);
}

/////////////////////////
// Particle ID helpers //
/////////////////////////

typedef pair<bool, TString>    DMesonPack;
typedef map<TString, Int_t>    PartIdMap;
typedef map<TString, Double_t> PartMomMap;

DMesonPack isDMeson(const PartIdMap parts) {
  for (const auto [key, val] : parts) {
    auto id = TMath::Abs(val);
    if (digitIs(id, 3) == 4) return DMesonPack{true, key};
  }

  return DMesonPack{false, "none"};
}

bool isDMeson(const int id) {
  if (digitIs(TMath::Abs(id), 3) == 4) return true;
  return false;
}

bool isHadron(const int id) {
  if (TMath::Abs(id) > 100) return true;
  return false;
}

bool isDstMeson(const int id) {
  return TMath::Abs(id) == 413 || TMath::Abs(id) == 423;
}

// We fix particle IDs based on Muon's true ID
int bIdFix(const int bId, const int dId) {
  if (bId * dId > 0) return -bId;
  return bId;
}

int idFix(const int muId, const int trueId) {
  Int_t sign = muId / TMath::Abs(muId);
  return sign * trueId;
}

int muIdFix(const int muId) { return idFix(muId, 13); }

int tauIdFix(const int muId) { return idFix(muId, 15); }

int nuIdFix(const int muId, bool isTau) {
  if (isTau) return idFix(muId, -16);
  return idFix(muId, -14);
}

int tauNuMuIdFix(const int muId) { return idFix(muId, -14); }

int tauNuTauIdFix(const int muId) { return idFix(muId, 16); }
