#ifndef _HAM_RWT_UTILS_HAM_
#define _HAM_RWT_UTILS_HAM_

#include <TString.h>

#include <Hammer/Math/FourMomentum.hh>
#include <Hammer/Particle.hh>

////////////////////////////
// HAMMER-related helpers //
////////////////////////////

TString print_p(const Hammer::FourMomentum& p) {
  char tmp[80];
  sprintf(tmp, "%.2f, %.2f, %.2f, %.2f", p.E(), p.px(), p.py(), p.pz());
  return TString(tmp);
}

auto inv_m(Double_t pe, Double_t px, Double_t py, Double_t pz) {
  return pe * pe - px * px - py * py - pz * pz;
}

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
