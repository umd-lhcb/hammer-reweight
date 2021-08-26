#ifndef _HAM_RWT_UTILS_
#define _HAM_RWT_UTILS_

#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <TDataType.h>
#include <TMath.h>

#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>

using std::map;
using std::pair;
using std::string;
using std::vector;

/////////////////////
// General helpers //
/////////////////////

template <template <typename, typename> class Iterable, typename T,
          typename Allocator>
Bool_t find_in(Iterable<T, Allocator> iter, T elem) {
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

Int_t digit_is(Int_t num, Int_t digit, Int_t base = 10) {
  Int_t fac = num / TMath::Power(base, digit - 1);
  return fac % 10;
}

TString print_p(const Hammer::FourMomentum& p) {
  char tmp[80];
  sprintf(tmp, "%.2f, %.2f, %.2f, %.2f", p.E(), p.px(), p.py(), p.pz());
  return TString(tmp);
}

/////////////////////////
// Particle ID helpers //
/////////////////////////

typedef pair<Bool_t, TString>  DMesonPack;
typedef map<TString, Int_t>    PartIdMap;
typedef map<TString, Double_t> PartMomMap;

DMesonPack is_D_meson(const PartIdMap parts) {
  for (const auto [key, val] : parts) {
    auto id = TMath::Abs(val);
    if (digit_is(id, 3) == 4) return DMesonPack{true, key};
  }

  return DMesonPack{false, "none"};
}

Bool_t is_hadron(const Int_t id) {
  if (TMath::Abs(id) > 100) return true;
  return false;
}

Int_t B_id_fix(const Int_t B_id, const Int_t D_id) {
  if (B_id * D_id > 0) return -B_id;
  return B_id;
}

int Mu_id(Int_t mu_id, Int_t true_id = 13) {
  Int_t sign = mu_id / TMath::Abs(mu_id);
  return sign * true_id;
}

int Tau_id(Int_t mu_id) { return Mu_id(mu_id, 15); }

int Nu_id(Int_t mu_id, Bool_t is_tau) {
  if (is_tau) return Mu_id(mu_id, -16);
  return Mu_id(mu_id, -14);
}

int Tau_NuMu_id(Int_t mu_id) { return Mu_id(mu_id, -14); }

int Tau_NuTau_id(Int_t mu_id) { return Mu_id(mu_id, 16); }

////////////////////////////
// HAMMER-related helpers //
////////////////////////////

auto particle(Double_t pe, Double_t px, Double_t py, Double_t pz, Int_t pid) {
  auto four_mom = Hammer::FourMomentum(pe, px, py, pz);
  auto part_id  = static_cast<Hammer::PdgId>(pid);

  return Hammer::Particle(four_mom, part_id);
}

auto particle(Hammer::FourMomentum four_mom, Int_t pid) {
  auto part_id = static_cast<Hammer::PdgId>(pid);
  return Hammer::Particle(four_mom, part_id);
}

#endif
