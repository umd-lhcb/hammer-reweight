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
Bool_t findIn(Iterable<T, Allocator> iter, T elem) {
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

Int_t digitIs(Int_t num, Int_t digit, Int_t base = 10) {
  Int_t fac = num / TMath::Power(base, digit - 1);
  return fac % 10;
}

string getParticleName(Int_t id, TDatabasePDG* db, Bool_t useAbsId = false) {
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

////////////////
// Kinematics //
////////////////

auto invM(Double_t pe, Double_t px, Double_t py, Double_t pz) {
  return TMath::Sqrt(pe * pe - px * px - py * py - pz * pz);
}

/////////////////////////
// Particle ID helpers //
/////////////////////////

typedef pair<Bool_t, TString>  dMesonPack;
typedef map<TString, Int_t>    partIdMap;
typedef map<TString, Double_t> partMomMap;

dMesonPack isDMeson(const partIdMap parts) {
  for (const auto [key, val] : parts) {
    auto id = TMath::Abs(val);
    if (digitIs(id, 3) == 4) return dMesonPack{true, key};
  }

  return dMesonPack{false, "none"};
}

Bool_t isDMeson(const Int_t id) {
  if (digitIs(id, 3) == 4) return true;
  return false;
}

Bool_t isHadron(const Int_t id) {
  if (TMath::Abs(id) > 100) return true;
  return false;
}

// We fix particle IDs based on Muon's true ID

Int_t bIdFix(const Int_t bId, const Int_t dId) {
  if (bId * dId > 0) return -bId;
  return bId;
}

int muIdFix(Int_t muId, Int_t trueId = 13) {
  Int_t sign = muId / TMath::Abs(muId);
  return sign * trueId;
}

int tauIdFix(Int_t muId) { return muIdFix(muId, 15); }

int nuIdFix(Int_t muId, Bool_t isTau) {
  if (isTau) return muIdFix(muId, -16);
  return muIdFix(muId, -14);
}

int tauNuMuIdFix(Int_t muId) { return muIdFix(muId, -14); }

int tauNuTauIdFix(Int_t muId) { return muIdFix(muId, 16); }
