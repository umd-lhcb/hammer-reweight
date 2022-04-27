// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Tue Apr 26, 2022 at 11:44 PM -0400

#pragma once

#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <TDataType.h>
#include <TDatabasePDG.h>
#include <TMath.h>

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

////////////////////////
// RDataframe helpers //
////////////////////////

////////////////////////
// Kinematics helpers //
////////////////////////

double invM(double pe, double px, double py, double pz) {
  return TMath::Sqrt(pe * pe - px * px - py * py - pz * pz);
}

/////////////////////////
// Particle ID helpers //
/////////////////////////

typedef pair<bool, TString>    dMesonPack;
typedef map<TString, Int_t>    partIdMap;
typedef map<TString, Double_t> partMomMap;

dMesonPack isDMeson(const partIdMap parts) {
  for (const auto [key, val] : parts) {
    auto id = TMath::Abs(val);
    if (digitIs(id, 3) == 4) return dMesonPack{true, key};
  }

  return dMesonPack{false, "none"};
}

bool isDMeson(const Int_t id) {
  if (digitIs(id, 3) == 4) return true;
  return false;
}

bool isHadron(const Int_t id) {
  if (TMath::Abs(id) > 100) return true;
  return false;
}

// We fix particle IDs based on Muon's true ID

int bIdFix(const int bId, const int dId) {
  if (bId * dId > 0) return -bId;
  return bId;
}

int muIdFix(int muId, int trueId = 13) {
  Int_t sign = muId / TMath::Abs(muId);
  return sign * trueId;
}

int tauIdFix(int muId) { return muIdFix(muId, 15); }

int nuIdFix(int muId, bool isTau) {
  if (isTau) return muIdFix(muId, -16);
  return muIdFix(muId, -14);
}

int tauNuMuIdFix(int muId) { return muIdFix(muId, -14); }

int tauNuTauIdFix(int muId) { return muIdFix(muId, 16); }
